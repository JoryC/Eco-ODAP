library(tximport)
library(dplyr)
library(readr)

# Function to write table with specific formatting
write_custom_table <- function(data_matrix, file_path) {
  # Get column and row names
  col_names <- colnames(data_matrix)
  row_names <- rownames(data_matrix)
  
  # Open file connection
  con <- file(file_path, "w")
  
  # Write header row - tab followed by quoted column names
  header <- paste0("\t", paste0("\"", col_names, "\"", collapse = "\t"))
  cat(header, "\n", file = con, sep = "")
  
  # Write each data row - quoted gene ID followed by unquoted counts
  for (i in 1:nrow(data_matrix)) {
    row_data <- paste0("\"",
                       row_names[i],
                       "\"",
                       "\t",
                       paste0(data_matrix[i, ], collapse = "\t"))
    cat(row_data, "\n", file = con, sep = "")
  }
  
  # Close connection
  close(con)
  cat("Table written to", file_path, "\n")
}

# Define argument inputs
args <- commandArgs(trailingOnly = TRUE)
salmon_dir <- args[1]
tx2gene_file <- args[2]
output_dir <- args[3]
Sequencing_type <- args[4]

#==============================================================================
# PART 1: Process transcript-level files (quant.sf)
#==============================================================================

files <- list.files(
  path = salmon_dir,
  pattern = "quant.sf$",
  full.names = TRUE,
  recursive = TRUE
)
# Exclude quant.genes.sf files from this list
files <- files[!grepl("quant.genes.sf", files)]

if (length(files) > 0) {
  sample_names <- gsub(".*/([^/]+)_quant/quant.sf$", "\\1", files)
  names(files) <- sample_names
  
  # Print sample names for verification
  cat("Found", length(files), "transcript-evel quant samples:\n")
  cat(paste(sample_names, collapse = "\n"), "\n")
  
  # Read in tx2gene file and clean up colnames
  tx2gene <- read_tsv(tx2gene_file, col_names = FALSE)
  colnames(tx2gene) <- c("TXNAME", "GENEID", "GENENAME", "BIOTYPE")
  tx2gene <- tx2gene[which(duplicated(tx2gene) == FALSE), ]
  tx2gene <- tx2gene %>%
    select(TXNAME, GENEID)
  
  #Import the tx2gene file
  if (Sequencing_type == "3PrimeRNASeq") {
    txi <- tximport(
      files,
      type = "salmon",
      tx2gene = tx2gene,
      countsFromAbundance = "no"
    )
  } else {
    txi <- tximport(
      files,
      type = "salmon",
      tx2gene = tx2gene,
      countsFromAbundance = "lengthScaledTPM"
    )
  }
  
  # Remove NAs
  txi$counts[is.na(txi$counts)] <- 0
  # Round the counts to nearest integer value
  counts_matrix <- round(txi$counts)
  
  # Write the count matrix to a file
  output_file <- file.path(output_dir, "genes.data.tsv")
  write_custom_table(counts_matrix, output_file)  # This function copies the format of rsem-generate-data-matrix
  
  cat("Transcript-to-gene-level count matrix created at",
      output_file,
      "\n")
}

#==============================================================================
# PART 2: Process direct gene-level files (quant.genes.sf)
#==============================================================================

# Find all quant.genes.sf files
gene_files <- list.files(
  path = salmon_dir,
  pattern = "quant.genes.sf$",
  full.names = TRUE,
  recursive = TRUE
)

if (length(gene_files) > 0) {
  # Extract sample names
  gene_sample_names <- gsub(".*/([^/]+)_quant/quant.genes.sf$", "\\1", gene_files)
  
  # Print sample names for verification
  cat("Found", length(gene_files), "gene-level samples:\n")
  cat(paste(gene_sample_names, collapse = "\n"), "\n")
  
  # Initialize a list to store count data
  gene_counts_list <- list()
  
  # Process each gene-level file
  for (i in 1:length(gene_files)) {
    file <- gene_files[i]
    sample <- gene_sample_names[i]
    
    # Read the quant.genes.sf file
    gene_data <- read_tsv(file)
    
    # Extract gene IDs and counts (NumReads column, usually column 5)
    counts <- gene_data %>%
      select(Name, NumReads) %>%
      rename(gene_id = Name)
    
    # Store in list
    gene_counts_list[[sample]] <- counts
  }
  
  # Merge all samples into a single data frame
  gene_counts_all <- gene_counts_list[[1]]
  colnames(gene_counts_all)[2] <- gene_sample_names[1]
  
  if (length(gene_files) > 1) {
    for (i in 2:length(gene_files)) {
      sample <- gene_sample_names[i]
      counts <- gene_counts_list[[sample]]
      gene_counts_all <- left_join(gene_counts_all, counts, by = "gene_id")
      colnames(gene_counts_all)[ncol(gene_counts_all)] <- sample
    }
  }
  
  # Replace NAs with zeros
  gene_counts_all[is.na(gene_counts_all)] <- 0
  
  # Round to integers
  gene_counts_all[, gene_sample_names] <- round(gene_counts_all[, gene_sample_names])
  
  # Create the exact format matrix
  direct_gene_matrix <- as.matrix(gene_counts_all[, gene_sample_names])
  rownames(direct_gene_matrix) <- gene_counts_all$gene_id
  
  write_custom_table(direct_gene_matrix,
                     file.path(output_dir, "alternative.genes.data.tsv"))
  
  cat("Direct gene-level expression matrices created in",
      output_dir,
      "\n")
}

# #Checking validity of each method:
# # Get all unique gene IDs from both matrices and sort them
# all_genes <- sort(unique(c(rownames(counts_matrix), rownames(direct_gene_matrix))))
# 
# # Subset and reorder both matrices to include all genes
# new_counts_matrix <- matrix(0, nrow = length(all_genes), ncol = ncol(counts_matrix))
# rownames(new_counts_matrix) <- all_genes
# colnames(new_counts_matrix) <- colnames(counts_matrix)
# # Fill in values for genes that exist in counts_matrix
# common_genes_tx <- intersect(rownames(counts_matrix), all_genes)
# new_counts_matrix[common_genes_tx, ] <- counts_matrix[common_genes_tx, ]
# 
# new_direct_matrix <- matrix(0, nrow = length(all_genes), ncol = ncol(direct_gene_matrix))
# rownames(new_direct_matrix) <- all_genes
# colnames(new_direct_matrix) <- colnames(direct_gene_matrix)
# # Fill in values for genes that exist in direct_gene_matrix
# common_genes_direct <- intersect(rownames(direct_gene_matrix), all_genes)
# new_direct_matrix[common_genes_direct, ] <- direct_gene_matrix[common_genes_direct, ]
# 
# # Replace original matrices with new sorted ones
# counts_matrix <- new_counts_matrix
# direct_gene_matrix <- new_direct_matrix
# 
# # Basic check: Are the rownames identical and in the same order?
# if (identical(rownames(counts_matrix), rownames(direct_gene_matrix))) {
#   cat("✓ Matrices have identical gene order with", nrow(counts_matrix), "genes\n")
# } else {
#   cat("✗ Matrices have different gene orders!\n")
#   
#   # Check if they have the same number of rows
#   if (nrow(counts_matrix) != nrow(direct_gene_matrix)) {
#     cat("  → Different number of rows:", nrow(counts_matrix), "vs", nrow(direct_gene_matrix), "\n")
#   }
#   
#   # Check if they contain the same genes
#   genes_in_tx_not_direct <- setdiff(rownames(counts_matrix), rownames(direct_gene_matrix))
#   genes_in_direct_not_tx <- setdiff(rownames(direct_gene_matrix), rownames(counts_matrix))
#   
#   if (length(genes_in_tx_not_direct) > 0) {
#     cat("  → Genes in transcript-derived matrix but not in direct-gene matrix:", 
#         length(genes_in_tx_not_direct), "\n")
#   }
#   
#   if (length(genes_in_direct_not_tx) > 0) {
#     cat("  → Genes in direct-gene matrix but not in transcript-derived matrix:", 
#         length(genes_in_direct_not_tx), "\n")
#   }
#   
#   # Check first few positions to see where they differ
#   if (length(rownames(counts_matrix)) > 0 && length(rownames(direct_gene_matrix)) > 0) {
#     cat("  → First 5 genes in transcript-derived matrix:\n    ", 
#         paste(head(rownames(counts_matrix), 5), collapse=", "), "\n")
#     cat("  → First 5 genes in direct-gene matrix:\n    ", 
#         paste(head(rownames(direct_gene_matrix), 5), collapse=", "), "\n")
#   }
# }
# 
# # Additional check: Is the row ordering strictly alphabetical?
# if (all(rownames(counts_matrix) == sort(rownames(counts_matrix)))) {
#   cat("✓ Gene order is alphabetically sorted\n")
# } else {
#   cat("✗ Gene order is not alphabetically sorted\n")
# }