library('AnnotationHub')
library('biomaRt')
library(GenomicFeatures)
library(rtracklayer)
library(rentrez)
library(dplyr)
library(readr)
library(httr)
library(purrr)
library(furrr)
library(jsonlite)
setAnnotationHubOption("ASK", FALSE)
# source(file = "C:/Users/curryj/OneDrive - EC-EC/Desktop/tPOD_Curation_Project/tPOD_curation/Handy_scripts/OrthoDBFunctions.R")

#############################################

# Define the species you would like to annotate
# Choose from: mouse, zebrafish, rainbowtrout, chinesemedaka, fatheadminnow
species <- "fatheadminnow"
#Define the name of the GTF file which should be in this /genome directory.
mygtffile <- "./FHM_genomic.gtf"
mysummarizedblastnresults <- "summarized_annotated_results_FHM_with_ZF_homologs.txt" #For Fathead Minnow only, since reference and annotations are not available on Ensembl!!!

#############
# Mouse #
#############
#mouse
if (species == "mouse") {
  ensembl_species <- "mmusculus_gene_ensembl"
  species_gene_symbol <- "mgi_symbol"
}
#############
# Zebrafish # tested
#############
#zebrafish
if (species == "zebrafish") {
  ensembl_species <- "drerio_gene_ensembl"
  species_gene_symbol <- "external_gene_name"
}
#################
# Rainbow trout # tested
#################
#rainbowtrout
if (species == "rainbowtrout") {
  ensembl_species <- "omykiss_gene_ensembl"
  species_gene_symbol <- "external_gene_name"
  version <- "112" # Version 112 is required!!!
  ensembl <- useEnsembl(biomart = "genes", version = version)
  datasets <- listDatasets(ensembl)
  searchDatasets(mart = ensembl, pattern = "omykiss")
  ensembl <- useEnsembl(biomart = "genes", dataset = ensembl_species, version = "112")
  attributes <- listAttributes(ensembl)
}
##################
# Chinese Medaka # tested
##################
#chinesemedaka
if (species == "chinesemedaka") {
  ensembl_species <- "osinensis_gene_ensembl"
  species_gene_symbol <- "external_gene_name"
  version <- "112" # Version 112 is required!!!
  ensembl <- useEnsembl(biomart = "genes", version = version)
  datasets <- listDatasets(ensembl)
  searchDatasets(mart = ensembl, pattern = "osinensis")
  ensembl <- useEnsembl(biomart = "genes", dataset = ensembl_species, version = "112")
  attributes <- listAttributes(ensembl)
}
##################
# Fathead Minnow # tested
##################
#fatheadminnow
if (species == "fatheadminnow") {
  ensembl_species <- "pimephales_promelas_gene_ensembl"
  species_gene_symbol <- "external_gene_name"
}
##################
# Daphnia Magna   # tested
##################
#daphnia
if (species == "daphnia") {
  # Connect to Ensembl Metazoa and list marts
  marts <- listMarts(host = "https://metazoa.ensembl.org")
  metazoa_mart <- as.character(marts[[1]][1])
  dataset_list <- useMart(metazoa_mart, host = "https://metazoa.ensembl.org")
  dataset_list <- listDatasets(dataset_list)
  ensembl_species <- "dmgca020631705v2_eg_gene" # Daphnia Magna Ensembl dataset
  ensembl <- useMart(metazoa_mart, dataset = ensembl_species, host = "https://metazoa.ensembl.org")
  attributes <- listAttributes(ensembl)
  species_gene_symbol <- "external_gene_name" # Daphnia Magna uses the same attribute for gene symbols
}
#####################################
# Algae (Chlamydomonas_reinhardtii) # tested
#####################################
#algae
if (species == "algae") {
  # Connect to Ensembl plants and list marts
  marts <- listMarts(host = "https://plants.ensembl.org")
  plant_mart <- as.character(marts[[1]][1])
  dataset_list <- useMart(plant_mart, host = "https://plants.ensembl.org")
  dataset_list <- listDatasets(dataset_list)
  ensembl_species <- "creinhardtii_eg_gene" # Algae Ensembl dataset
  ensembl <- useMart(plant_mart, dataset = ensembl_species, host = "https://plants.ensembl.org")
  attributes <- listAttributes(ensembl)
  species_gene_symbol <- "external_gene_name" # Algae uses the same attribute for gene symbols
}
#############################################

# species <- if_else(
#   ensembl_species == "drerio_gene_ensembl",
#   "zebrafish",
#   if_else(
#     ensembl_species == "omykiss_gene_ensembl",
#     "rainbowtrout",
#     if_else(
#       ensembl_species == "osinensis_gene_ensembl",
#       "chinesemedaka",
#       if_else(
#         ensembl_species == "pimephales_promelas_gene_ensembl",
#         "fatheadminnow",
#         if_else(
#           ensembl_species == "mmusculus_gene_ensembl",
#           "mouse",
#           NA
#         )
#       )
#     )
#   )
# )

biomart_filter <- "ensembl_gene_id"
#Note - useMart() will not work for Fathead Minnow
#Note - Chinese Medaka is not included int the latest version (version 113) of Ensembl for some reason...
#Note - Rainbow Trout is not included int the latest version (version 113) of Ensembl for some reason...

if (!species %in% c("rainbowtrout", "chinesemedaka", "fatheadminnow")) {
  #all_ensembl <- useMart("ensembl")
  #all_datasets <- listDatasets(all_ensembl)
  ensembl <- useMart("ensembl", dataset = ensembl_species, host = "https://useast.ensembl.org") #West is no longer working
  attributes <- listAttributes(ensembl)
}

# Example generic gene IDs (for testing purposes)
#genes <- c("ENSDARG00000000001", "ENSDARG00000000002", "ENSDARG00000000003", "ENSDARG00000000018", "ENSDARG00000000019", "ENSDARG00000000068", "ENSDARG00000000069", "ENSOMYG00000042616", "ENSOMYG00000042616", "ENSOMYG00000042617", "ENSOMYG00000042618", "ENSOSIG00000001234", "ENSOSIG00000005678", "ENSOSIG00000009123", "ENSPPRG00000012345", "ENSPPRG00000067890", "ENSPPRG00000011223")
#genes <- lapply(resListAll,
#                function(x) row.names(as.data.frame(x)))
#genes <- unlist(genes) %>% unique()

if (species == "daphnia") {
  gtf_file <- mygtffile
  # Import the GFF/GTF file
  annotations <- import(gtf_file)
  
  # Convert annotations to a data frame
  annotations_df <- as.data.frame(annotations) %>%
    mutate(prefix = sub("_.*", "", transcript_id))
  
  #Count occurences of each prefix
  prefix_counts <- annotations_df %>%
    group_by(prefix) %>%
    summarize(count = n())
  
  genes <- annotations_df %>%
    dplyr::pull(gene_id) %>%
    unique()
  #SOME OF THE DAPHNIA GENE NAMES DON'T MATCH EXCACTLY WIOTH WHAT IS IN BIOMART SO I AM MANIPULATING THE GENE ID STRINGS TO MATCH THE BIOMART ENTRIES SO WE CAN GET THE FULL LIST OF ANNOTATED GENES.
  # Function to change only the specific pattern
  fix_trna_names <- function(gene_names) {
    # Only modify names that match the Trna pattern with underscore
    trna_pattern <- "^Trna.*-.*_[0-9]+$"
    
    # Find which genes match the pattern
    matches <- grepl(trna_pattern, gene_names)
    
    # For matching genes, replace underscore with hyphen
    gene_names[matches] <- gsub("_([0-9]+)$", "-\\1", gene_names[matches])
    
    return(gene_names)
  }
  
  # Usage
  corrected_gene_names <- fix_trna_names(genes)
  
  id_table_entrez <- getBM(
    #filters = biomart_filter,
    attributes = c(
      biomart_filter,
      species_gene_symbol,
      "description",
      "entrezgene_id",
      "entrezgene_accession"
    ),
    values = corrected_gene_names, #genes/all the genes in the genome
    mart = ensembl
  )
  
  fix_trna_names <- function(gene_names) {
    # Only modify names that match the Trna pattern with underscore
    trna_pattern <- "^Trna.*-.*-[0-9]+$"
    
    # Find which genes match the pattern
    matches <- grepl(trna_pattern, gene_names)
    
    # For matching genes, replace underscore with hyphen
    gene_names[matches] <- gsub("-([0-9]+)$", "_\\1", gene_names[matches])
    
    return(gene_names)
  }
  
  id_table_entrez <- id_table_entrez %>%
    mutate(ensembl_gene_id = if_else(stringr::str_starts(ensembl_gene_id, "GeneID"), external_gene_name, ensembl_gene_id)) %>%
    mutate(ensembl_gene_id = fix_trna_names(ensembl_gene_id))
  
  #############################################################################
  #Drosophilla homologues
  id_table_entrez <- getBM(
    #filters = biomart_filter,
    attributes = c(
      biomart_filter,
      species_gene_symbol,
      "description",
      "dmelanogaster_eg_homolog_ensembl_gene",
      "dmelanogaster_eg_homolog_associated_gene_name"
    ),
    values = corrected_gene_names, #genes/all the genes in the genome
    mart = ensembl
  )
  
  id_table_entrez <- id_table_entrez %>%
    unique()
  
  readr::write_csv(id_table_entrez, "./daphnia_annotated_id_table_entrez_wholegenome.csv")
}

if (species == "algae") {
  gtf_file <- mygtffile
  # Import the GFF/GTF file
  annotations <- import(gtf_file)
  
  # Convert annotations to a data frame
  annotations_df <- as.data.frame(annotations) %>%
    mutate(prefix = sub("_.*", "", transcript_id))
  
  #Count occurences of each prefix
  prefix_counts <- annotations_df %>%
    group_by(prefix) %>%
    summarize(count = n())
  
  genes <- annotations_df %>%
    dplyr::pull(gene_id) %>%
    unique()
  
  id_table_entrez <- getBM(
    #filters = biomart_filter,
    attributes = c(
      biomart_filter,
      species_gene_symbol,
      "description",
      "entrezgene_id",
      "entrezgene_accession"
    ),
    values = genes, #genes/all the genes in the genome
    mart = ensembl
  )
  
  id_table_entrez <- id_table_entrez %>%
    unique()
  
  readr::write_csv(id_table_entrez, "./algae_annotated_id_table_entrez_wholegenome.csv")
}

if (species %in% c("zebrafish", "hamster", "rat", "mouse", "human")) {
  ############################################################
  ### Zebrafish, Rainbow trout, Hamster, Rat, Mouse, Human ###
  ############################################################
  
  if (species == "zebrafish") {
    # Path to your local GFF/GTF file
    gtf_file <- mygtffile
  } else if (species == "mouse") {
    gtf_file <- mygtffile
  } else {
    stop("No gtf annotation files in directory... please download")
  }
  
  # Import the GFF/GTF file
  annotations <- import(gtf_file)
  
  # Convert annotations to a data frame
  annotations_df <- as.data.frame(annotations) %>%
    mutate(prefix = sub("_.*", "", transcript_id))
  
  #Count occurences of each prefix
  prefix_counts <- annotations_df %>%
    group_by(prefix) %>%
    summarize(count = n())
  
  genes <- annotations_df %>%
    dplyr::pull(gene_id) %>%
    unique()
  
  id_table_entrez <- getBM(
    filters = biomart_filter,
    attributes = c(
      biomart_filter,
      species_gene_symbol,
      "description",
      "entrezgene_id",
      "entrezgene_accession"
    ),
    values = genes, #genes/all the genes in the genome
    mart = ensembl
  )
  
  if (species == "zebrafish") {
    readr::write_csv(id_table_entrez, "./zebrafish_annotated_id_table_entrez_wholegenome.csv")    
  } else if (species == "mouse") {
    readr::write_csv(id_table_entrez, "./mouse_annotated_id_table_entrez_wholegenome.csv") 
  } else {
    print("No files written")
  }

} else if (species %in% c("chinesemedaka")) {
  ######################
  ### Chinese Medaka ###
  ######################
  
  ## Get zebrafish homolog genes ##
  
  if (species == "chinesemedaka") {
    # Path to your local GFF/GTF file
    gtf_file <- mygtffile
  }
  
  # Import the GFF/GTF file
  annotations <- import(gtf_file)
  
  # Convert annotations to a data frame
  annotations_df <- as.data.frame(annotations)
  
  genes <- annotations_df %>%
    dplyr::pull(gene_id) %>%
    unique()
  
  id_table_entrez <- getBM(filters = biomart_filter,
                           attributes = c(biomart_filter,
                                          species_gene_symbol,
                                          "description",
                                          #"zfin_id_id",
                                          #"zfin_id_symbol",
                                          "drerio_homolog_ensembl_gene",
                                          "drerio_homolog_associated_gene_name"),
                           values = genes,
                           mart = ensembl)
  
  if (species == "chinesemedaka") {
    readr::write_csv(id_table_entrez, "./chinesemedaka_annotated_id_table_entrez_wholegenome.csv")
  } 
  
} else if (species %in% c("fatheadminnow", "rainbowtrout")) {

  if (species == "fatheadminnow") {
    ######################
    ### Fathead Minnow ###
    ######################
    # Path to your local GFF/GTF file
    gtf_file <-
      mygtffile
  } else if (species == "rainbowtrout") {
    #####################
    ### Rainbow trout ###
    #####################
    # Path to your local GFF/GTF file
    gtf_file <-
      mygtffile
  }
  
  if (species == "fatheadminnow") {
    ######################
    ### Fathead Minnow ###
    ######################
    
    #Note: No BioMart is available so all information is strictly from the .gtf file and BLAST
    
    # Import the GFF/GTF file
    annotations <- import(gtf_file)
    
    # Convert annotations to a data frame
    annotations_df <- as.data.frame(annotations)
    
    # Extract product information for transcripts and exons
    transcript_exon_products <- annotations_df %>%
      dplyr::filter(type %in% c("transcript", "exon", "CDS", "start_codon", "stop_codon")) %>%
      dplyr::select(gene_id, product, transcript_id) %>%
      dplyr::filter(!is.na(product)) %>%
      dplyr::distinct()
    
    # Aggregate product information at the gene level
    gene_products <- transcript_exon_products %>%
      dplyr::group_by(gene_id, transcript_id) %>%
      dplyr::summarize(product = paste(unique(product), collapse = "; "), .groups = 'drop')
    
    # Add product information to gene annotations
    df <- annotations_df %>%
      dplyr::left_join(gene_products, by = "gene_id") %>%
      dplyr::mutate(product = dplyr::coalesce(product.y, product.x)) %>%
      dplyr::mutate(transcript_id = dplyr::coalesce(transcript_id.y, transcript_id.x)) %>%
      dplyr::select(-product.x, -product.y, -transcript_id.y, -transcript_id.x)
    
    # Filter for genes and select relevant columns
    gene_annotations_df <- df %>%
      dplyr::filter(type == "gene") %>%
      dplyr::select(seqnames, start, end, strand, gene_id, transcript_id, db_xref, product)
    
    annotations_df <- gene_annotations_df %>%
      dplyr::mutate(ensembl_gene_id = gsub("GeneID:", "", db_xref)) %>%
      dplyr::mutate(entrezgene_id = ensembl_gene_id) %>%
      dplyr::mutate(external_gene_name = gene_id) %>%
      dplyr::mutate(entrezgene_accession = gene_id) %>%
      dplyr::mutate(prefix = sub("_.*", "", transcript_id)) %>%
      dplyr::rename(query_id = transcript_id) %>%
      dplyr::rename(description = product) %>%
      dplyr::select(ensembl_gene_id, external_gene_name, description, entrezgene_id, entrezgene_accession, query_id, prefix)
    
    ##########################
    #     BLASTN method      #
    ##########################
    
    # In linux terminal:
    
    # Download and install blast+
    # sudo apt-get install ncbi-blast+
    
    # Download the Zebrafish transcript FASTA and FHM Transcript FASTA
    
    # Create a blast database for zebrafish
    # makeblastdb -in danRer11.rna.fna -dbtype nucl -out zebrafish_nucletoide_db
    
    # Run Blastn
    # blastn -query fathead_minnow_nucleotides.fna -db zebrafish_nucleotide_db -outfmt 6 -out blastn_results.txt
    # -oufmt 6 =  tab delimited
    
    # blast_results <-
    #   read.table("FHM_ZF_blastn_results.txt",
    #              header = FALSE,
    #              sep = "\t")
    # colnames(blast_results) <-
    #   c(
    #     "query_id",
    #     "subject_id",
    #     "percent_identity",
    #     "alignment_length",
    #     "mismatches",
    #     "gap_opens",
    #     "q_start",
    #     "q_end",
    #     "s_start",
    #     "s_end",
    #     "evalue",
    #     "bit_score"
    #   )
    # 
    # # Remove version numbers
    # blast_results$Subject_ID <-
    #   gsub("\\..*", "", blast_results$subject_id)
    # 
    # # Extract the prefixes (first 2 characters before '_')
    # blast_results <- blast_results %>%
    #   mutate(prefix = sub("_.*", "", Subject_ID))
    # blast_results <- blast_results %>%
    #   mutate(prefix_query = sub("_.*", "", query_id))
    # 
    # # Count the occurrences of each prefix
    # prefix_counts <- blast_results %>%
    #   group_by(prefix) %>%
    #   summarize(count = n())
    # prefix_counts_query <- blast_results %>%
    #   group_by(prefix_query) %>%
    #   summarize(count = n())
    # 
    # # Display the counts
    # print(prefix_counts)
    # print(prefix_counts_query)
    # 
    # # Split by RefSeq type
    # nm_ids <-
    #   blast_results[grep("^NM_", blast_results$Subject_ID), "Subject_ID"] #Manually curated mRNA sequences for protein-coding genes.
    # nr_ids <-
    #   blast_results[grep("^NR_", blast_results$Subject_ID), "Subject_ID"] #Manually curated non-coding RNA sequences, such as rRNA, tRNA, or other regulatory RNAs
    # xm_ids <-
    #   blast_results[grep("^XM_", blast_results$Subject_ID), "Subject_ID"] #Computationally predicted mRNA sequences for protein-coding genes
    # xr_ids <-
    #   blast_results[grep("^XR_", blast_results$Subject_ID), "Subject_ID"] #Computationally predicted non-coding RNA sequences
    # 
    # #Connect to Ensembl
    # ensembl <- useMart("ensembl", dataset = "drerio_gene_ensembl")
    # 
    # attributes <- listAttributes(ensembl)
    # 
    # # Retrieve annotations for each type
    # nm_annotations <-
    #   getBM(
    #     attributes = c('refseq_mrna', 'ensembl_gene_id', 'external_gene_name'),
    #     filters = 'refseq_mrna',
    #     values = nm_ids,
    #     mart = ensembl
    #   )
    # nm_annotations <- nm_annotations %>%
    #   rename(subject_id = 'refseq_mrna')
    # 
    # nr_annotations <-
    #   getBM(
    #     attributes = c('refseq_ncrna', 'ensembl_gene_id', 'external_gene_name'),
    #     filters = 'refseq_ncrna',
    #     values = nm_ids,
    #     mart = ensembl
    #   )
    # nr_annotations <- nr_annotations %>%
    #   rename(subject_id = 'refseq_ncrna')
    # 
    # xm_annotations <-
    #   getBM(
    #     attributes = c(
    #       'refseq_mrna_predicted',
    #       'ensembl_gene_id',
    #       'external_gene_name'
    #     ),
    #     filters = 'refseq_mrna_predicted',
    #     values = xm_ids,
    #     mart = ensembl
    #   )
    # xm_annotations <- xm_annotations %>%
    #   rename(subject_id = 'refseq_mrna_predicted')
    # 
    # xr_annotations <-
    #   getBM(
    #     attributes = c(
    #       'refseq_ncrna_predicted',
    #       'ensembl_gene_id',
    #       'external_gene_name'
    #     ),
    #     filters = 'refseq_ncrna_predicted',
    #     values = xr_ids,
    #     mart = ensembl
    #   )
    # xr_annotations <- xr_annotations %>%
    #   rename(subject_id = 'refseq_ncrna_predicted')
    # 
    # # Combine all annotations
    # annotations <-
    #   rbind(nm_annotations,
    #         nr_annotations,
    #         xm_annotations,
    #         xr_annotations)
    # 
    # # Merge with original BLAST results
    # annotated_results <-
    #   merge(
    #     blast_results,
    #     annotations,
    #     by.x = "Subject_ID",
    #     by.y = "subject_id",
    #     all.x = TRUE
    #   ) %>% distinct()
    # 
    # write_tsv(annotated_results, file = "Annotated_FHM_genes_with_ZF_homologs.txt")
    # 
    # #annotated_results <- readr::read_tsv("Annotated_FHM_genes_with_ZF_homologs.txt")
    # 
    # ranked_annotated_results <- annotated_results %>%
    #   dplyr::group_by(query_id) %>%
    #   dplyr::mutate(
    #     evalue_rank = rank(evalue, ties.method = "min"),
    #     identity_rank = rank(-percent_identity, ties.method = "min"),
    #     bit_score_rank = rank(-bit_score, ties.method = "min"),
    #     combined_rank = evalue_rank + identity_rank + bit_score_rank
    #   ) %>%
    #   dplyr::ungroup()
    # 
    # # Function to get the best non-NA or the best-ranked NA entry
    # get_best_annotation <- function(data) {
    #   # Check if there is any non-NA ensembl_gene_id
    #   if (any(!is.na(data$ensembl_gene_id))) {
    #     # Select the best non-NA ensembl_gene_id
    #     data %>%
    #       filter(!is.na(ensembl_gene_id)) %>%
    #       arrange(combined_rank) %>%
    #       slice(1)
    #   } else {
    #     # If all are NA, select the best-ranked entry
    #     data %>%
    #       arrange(combined_rank) %>%
    #       slice(1)
    #   }
    # }
    # 
    # # Apply the function across all query_ids
    # best_annotated_results <- ranked_annotated_results %>%
    #   dplyr::group_by(query_id) %>%
    #   group_modify( ~ get_best_annotation(.x)) %>%
    #   ungroup()
    # 
    # summarized_annotated_results <-
    #   best_annotated_results %>% dplyr::select(Subject_ID,
    #                                            subject_id,
    #                                            query_id,
    #                                            ensembl_gene_id,
    #                                            external_gene_name,
    #                                            prefix) %>% distinct()
    # 
    # write_tsv(summarized_annotated_results, file = "summarized_annotated_results_FHM_with_ZF_homologs.txt")
    # 
    
    summarized_annotated_results <- readr::read_tsv(mysummarizedblastnresults)
    
    zf_homolog_annotated_results <- summarized_annotated_results %>%
      dplyr::mutate(prefix_query = sub("_.*", "", query_id)) %>%
      dplyr::filter(prefix_query %in% c("XM", "XR")) %>%
      dplyr::left_join(annotations_df, by = "query_id") %>%
      dplyr::select(ensembl_gene_id.x, external_gene_name.x, description, entrezgene_id, entrezgene_accession, query_id) %>%
      dplyr::rename(ensembl_gene_id = ensembl_gene_id.x, external_gene_name = external_gene_name.x)
    
    # Print the updated dataframe
    print(zf_homolog_annotated_results)
  }
  else if (species == "rainbowtrout") {
    #####################
    ### Rainbow trout ###
    #####################
    
    #Note: BioMart is now available!
    
    # Import the GFF/GTF file
    annotations <- import(gtf_file)
    
    # Convert annotations to a data frame
    annotations_df <- as.data.frame(annotations)
    
    genes <- annotations_df %>%
      dplyr::pull(gene_id) %>%
      unique()
    
    id_table_entrez <- getBM(filters = biomart_filter,
                             attributes = c(biomart_filter,
                                            species_gene_symbol,
                                            "description",
                                            #"zfin_id_id",
                                            #"zfin_id_symbol",
                                            "drerio_homolog_ensembl_gene",
                                            "drerio_homolog_associated_gene_name"),
                             values = genes,
                             mart = ensembl)
    
    # # Extract product information for transcripts and exons
    # transcript_exon_products <- annotations_df %>%
    #   dplyr::filter(type %in% c("transcript", "exon", "CDS", "start_codon", "stop_codon")) %>%
    #   dplyr::select(gene_id, product, transcript_id) %>%
    #   dplyr::filter(!is.na(product)) %>%
    #   dplyr::distinct()
    # 
    # # Aggregate product information at the gene level
    # gene_products <- transcript_exon_products %>%
    #   dplyr::group_by(gene_id, transcript_id) %>%
    #   dplyr::summarize(product = paste(unique(product), collapse = "; "), .groups = 'drop')
    # 
    # # Add product information to gene annotations
    # df <- annotations_df %>%
    #   dplyr::left_join(gene_products, by = "gene_id") %>%
    #   dplyr::mutate(product = dplyr::coalesce(product.y, product.x)) %>%
    #   dplyr::mutate(transcript_id = dplyr::coalesce(transcript_id.y, transcript_id.x)) %>%
    #   dplyr::select(-product.x, -product.y, -transcript_id.y, -transcript_id.x)
    # 
    # # Filter for genes and select relevant columns
    # gene_annotations_df <- df %>%
    #   dplyr::filter(type == "gene") %>%
    #   dplyr::select(seqnames, start, end, strand, gene_id, transcript_id, db_xref, product)
    # 
    # annotations_df <- gene_annotations_df %>%
    #   dplyr::mutate(ensembl_gene_id = gsub("GeneID:", "", db_xref)) %>%
    #   dplyr::mutate(entrezgene_id = ensembl_gene_id) %>%
    #   dplyr::mutate(external_gene_name = gene_id) %>%
    #   dplyr::mutate(entrezgene_accession = gene_id) %>%
    #   dplyr::mutate(prefix = sub("_.*", "", transcript_id)) %>%
    #   dplyr::rename(query_id = transcript_id) %>%
    #   dplyr::rename(description = product) %>%
    #   dplyr::select(ensembl_gene_id, external_gene_name, description, entrezgene_id, entrezgene_accession, query_id, prefix)
    
    ################
    # OrgDB method #
    ################
    
    # Add a new column to store RefSeq zebrafish orthologs
    #annotations_df$zebrafish_refseq <- NA
    
    # Loop over each row to get RefSeq IDs
    #for (i in 1:nrow(annotations_df)) {
    #  entrez_id <- annotations_df$entrezgene_id[i]
    #  refseq_ids <- tryCatch({
    #    get_all_refseq_ids(entrez_id)
    #  }, error = function(e) {
    #    message("Error with Entrez ID: ", entrez_id)
    #    return(NA)
    #  })
    
    # Add the RefSeq ID to the dataframe, if any are found
    #  if (length(refseq_ids) > 0) {
    #    annotations_df$zebrafish_refseq[i] <- paste(refseq_ids, collapse = "; ")
    #  }
    #}
    
    ##########################
    # Preferred BLAST method #
    ##########################
    
    # In linux terminal:
    
    # Download and install blast+
    # sudo apt-get install ncbi-blast+
    
    # Download the Zebrafish transcript FASTA and FHM Transcript FASTA
    
    # Create a blast database for zebrafish
    # makeblastdb -in danRer11.rna.fna -dbtype nucl -out zebrafish_nucletoide_db
    
    # Run Blastn
    # blastn -query fathead_minnow_nucleotides.fna -db zebrafish_nucleotide_db -outfmt 6 -out blastn_results.txt
    # -oufmt 6 =  tab delimited
    
    # # Install and load biomaRt
    # if (!requireNamespace("biomaRt", quietly = TRUE)) {
    #   BiocManager::install("biomaRt")
    # }
    # library(biomaRt)
    # library(tidyverse)
    # blast_results <-
    #   read.table("RT_ZF_blastn_results.txt",
    #              header = FALSE,
    #              sep = "\t")
    # colnames(blast_results) <-
    #   c(
    #     "query_id",
    #     "subject_id",
    #     "percent_identity",
    #     "alignment_length",
    #     "mismatches",
    #     "gap_opens",
    #     "q_start",
    #     "q_end",
    #     "s_start",
    #     "s_end",
    #     "evalue",
    #     "bit_score"
    #   )
    # 
    # # Remove version numbers
    # blast_results$Subject_ID <-
    #   gsub("\\..*", "", blast_results$subject_id)
    # 
    # # Extract the prefixes (first 2 characters before '_')
    # blast_results <- blast_results %>%
    #   mutate(prefix = sub("_.*", "", Subject_ID))
    # blast_results <- blast_results %>%
    #   mutate(prefix_query = sub("_.*", "", query_id))
    # 
    # # Count the occurrences of each prefix
    # prefix_counts <- blast_results %>%
    #   group_by(prefix) %>%
    #   summarize(count = n())
    # prefix_counts_query <- blast_results %>%
    #   group_by(prefix_query) %>%
    #   summarize(count = n())
    # 
    # # Display the counts
    # print(prefix_counts)
    # print(prefix_counts_query)
    # 
    # # Split by RefSeq type
    # nm_ids <-
    #   blast_results[grep("^NM_", blast_results$Subject_ID), "Subject_ID"] #Manually curated mRNA sequences for protein-coding genes.
    # nr_ids <-
    #   blast_results[grep("^NR_", blast_results$Subject_ID), "Subject_ID"] #Manually curated non-coding RNA sequences, such as rRNA, tRNA, or other regulatory RNAs
    # xm_ids <-
    #   blast_results[grep("^XM_", blast_results$Subject_ID), "Subject_ID"] #Computationally predicted mRNA sequences for protein-coding genes
    # xr_ids <-
    #   blast_results[grep("^XR_", blast_results$Subject_ID), "Subject_ID"] #Computationally predicted non-coding RNA sequences
    # 
    # #Connect to Ensembl
    # ensembl <- useMart("ensembl", dataset = "drerio_gene_ensembl")
    # 
    # #attributes <- listAttributes(ensembl)
    # 
    # # Retrieve annotations for each type
    # nm_annotations <-
    #   getBM(
    #     attributes = c('refseq_mrna', 'ensembl_gene_id', 'external_gene_name'),
    #     filters = 'refseq_mrna',
    #     values = nm_ids,
    #     mart = ensembl
    #   )
    # nm_annotations <- nm_annotations %>%
    #   rename(subject_id = 'refseq_mrna')
    # 
    # nr_annotations <-
    #   getBM(
    #     attributes = c('refseq_ncrna', 'ensembl_gene_id', 'external_gene_name'),
    #     filters = 'refseq_ncrna',
    #     values = nm_ids,
    #     mart = ensembl
    #   )
    # nr_annotations <- nr_annotations %>%
    #   rename(subject_id = 'refseq_ncrna')
    # 
    # xm_annotations <-
    #   getBM(
    #     attributes = c(
    #       'refseq_mrna_predicted',
    #       'ensembl_gene_id',
    #       'external_gene_name'
    #     ),
    #     filters = 'refseq_mrna_predicted',
    #     values = xm_ids,
    #     mart = ensembl
    #   )
    # xm_annotations <- xm_annotations %>%
    #   rename(subject_id = 'refseq_mrna_predicted')
    # 
    # xr_annotations <-
    #   getBM(
    #     attributes = c(
    #       'refseq_ncrna_predicted',
    #       'ensembl_gene_id',
    #       'external_gene_name'
    #     ),
    #     filters = 'refseq_ncrna_predicted',
    #     values = xr_ids,
    #     mart = ensembl
    #   )
    # xr_annotations <- xr_annotations %>%
    #   rename(subject_id = 'refseq_ncrna_predicted')
    # 
    # # Combine all annotations
    # annotations <-
    #   rbind(nm_annotations,
    #         nr_annotations,
    #         xm_annotations,
    #         xr_annotations)
    # 
    # # Merge with original BLAST results
    # annotated_results <-
    #   merge(
    #     blast_results,
    #     annotations,
    #     by.x = "Subject_ID",
    #     by.y = "subject_id",
    #     all.x = TRUE
    #   ) %>% distinct()
    # 
    # write_tsv(annotated_results, file = "Annotated_RT_genes_with_ZF_homologs.txt")
    # 
    # annotated_results <-
    #   readr::read_tsv("Annotated_RT_genes_with_ZF_homologs.txt")
    # 
    # ranked_annotated_results <- annotated_results %>%
    #   dplyr::group_by(query_id) %>%
    #   dplyr::mutate(
    #     evalue_rank = rank(evalue, ties.method = "min"),
    #     identity_rank = rank(-percent_identity, ties.method = "min"),
    #     bit_score_rank = rank(-bit_score, ties.method = "min"),
    #     combined_rank = evalue_rank + identity_rank + bit_score_rank
    #   ) %>%
    #   dplyr::ungroup()
    # 
    # # Function to get the best non-NA or the best-ranked NA entry
    # get_best_annotation <- function(data) {
    #   # Check if there is any non-NA ensembl_gene_id
    #   if (any(!is.na(data$ensembl_gene_id))) {
    #     # Select the best non-NA ensembl_gene_id
    #     data %>%
    #       filter(!is.na(ensembl_gene_id)) %>%
    #       arrange(combined_rank) %>%
    #       slice(1)
    #   } else {
    #     # If all are NA, select the best-ranked entry
    #     data %>%
    #       arrange(combined_rank) %>%
    #       slice(1)
    #   }
    # }
    # 
    # # Apply the function across all query_ids
    # best_annotated_results <- ranked_annotated_results %>%
    #   dplyr::group_by(query_id) %>%
    #   group_modify( ~ get_best_annotation(.x)) %>%
    #   ungroup()
    # 
    # summarized_annotated_results <-
    #   best_annotated_results %>% dplyr::select(Subject_ID,
    #                                            subject_id,
    #                                            query_id,
    #                                            ensembl_gene_id,
    #                                            external_gene_name,
    #                                            prefix) %>% distinct()
    # 
    # write_tsv(summarized_annotated_results, file = "summarized_annotated_results_RT_with_ZF_homologs.txt")
    
    # summarized_annotated_results <- readr::read_tsv("summarized_annotated_results_RT_with_ZF_homologs.txt")
    # 
    # zf_homolog_annotated_results <- summarized_annotated_results %>%
    #   dplyr::mutate(prefix_query = sub("_.*", "", query_id)) %>%
    #   dplyr::filter(prefix_query %in% c("XM", "XR")) %>%
    #   dplyr::left_join(id_table_entrez, by = "query_id") %>%
    #   dplyr::select(ensembl_gene_id.x, external_gene_name.x, description, entrezgene_id, entrezgene_accession, query_id) %>%
    #   dplyr::rename(ensembl_gene_id = ensembl_gene_id.x, external_gene_name = external_gene_name.x)
    print(id_table_entrez)
  }
  
  if (species == "fatheadminnow") {
    readr::write_csv(
      zf_homolog_annotated_results,
      "./fatheadminnow_annotated_id_table_entrez_wholegenome.csv"
    )
  } else if (species == "rainbowtrout") {
    readr::write_csv(
      id_table_entrez,
      "./rainbowtrout_annotated_id_table_entrez_wholegenome.csv"
    )
  }
}

###############################################################################

########################
# For Testing purposes #
########################

#Change paths if necessary
paths <- list()
paths$RData <- getwd() #file.path("../Study_dir/DEG_output/RData/")
paths$metadata <- file.path(getwd(), "..")
Platform <- "RNA-Seq"
#resListAll <- load("resListAll.RData")
params <- list()
params$use_cached_RData = FALSE

# Compile list of all potentially relevant genes
if (Platform == "RNA-Seq") {
  genes <- lapply(resListAll,
                  function(x) row.names(as.data.frame(x)))
  genes <- unlist(genes) %>% unique()
} else if (Platform == "TempO-Seq") {
  # NOTE: The original R-ODAF code for annotating TEMPO-SEQ data is old, outdated and does not work with manifests availbale on biospyder's website
  # It is not advisable to filter the probes by gene ID and gene name/ensembl id because probes come from transcripts. multiple probes can come from a single gene, so you lose that detail if you do what I've done here... I don't really care though because a DEG analysis is not my main priority here... it is BMDs and tPODs which will still have unique probe IDs.
  # Extract the result probes
  genes <- lapply(resListAll,
                  function(x) row.names(as.data.frame(x)))
  genes <- unlist(genes) %>% unique()
  # Extract gene symbols from the 'genes' list
  gene_ids <- sapply(genes, function(x) strsplit(x, "_")[[1]][1])
  gene_ids <- tibble(GENE_SYMBOL = gene_ids, PROBE = names(gene_ids))
  # load the manifest
  temposeq_genes <- readr::read_csv(normalizePath(file.path(projectdir, "genome", temposeq_manifest))) #Note - You will have to include this file in the genome directory. Original file copy can be found in project dir.
  # Filter the manifest by probe names that match the results in resList
  #gene_ids <- temposeq_genes %>%
  #	dplyr::filter(PROBE_NAME %in% genes) %>%
  #	pull(GENE_SYMBOL)
  # Now load the annotation file with the ensembl gene ids
  id_table_entrez_all_genes <- readr::read_csv(normalizePath(file.path(projectdir, "genome", paste0(species, "_annotated_id_table_entrez_wholegenome.csv")))) #Note - You will have to include this file in the genome directory. Original file copy can be found in project dir.
  # Filter the gene id column by those gene ids from the filtered manifest
  id_table_entrez <- gene_ids %>%
    dplyr::left_join(id_table_entrez_all_genes, by = c("GENE_SYMBOL" = species_gene_symbol), keep = TRUE) %>%
    dplyr::mutate(
      !!sym(biomart_filter) := if_else(is.na(!!sym(biomart_filter)), NA_character_, !!sym(biomart_filter)),
      !!sym(species_gene_symbol) := if_else(is.na(!!sym(species_gene_symbol)), NA_character_, !!sym(species_gene_symbol)),
      description = if_else(is.na(description), NA_character_, description),
      entrezgene_id = if_else(is.na(entrezgene_id), NA_integer_, entrezgene_id),
      entrezgene_accession = if_else(is.na(entrezgene_accession), NA_character_, entrezgene_accession),
      PROBE = if_else(is.na(PROBE), GENE_SYMBOL, PROBE)
    ) %>%
    dplyr::filter(PROBE %in% genes) %>%
    dplyr::select(!!sym(biomart_filter), !!sym(species_gene_symbol), description, entrezgene_id, entrezgene_accession, PROBE)
  save(id_table_entrez, file = normalizePath(file.path(paths$RData, "id_table.RData")))
  id_table <- distinct(id_table_entrez[, c(1,2,3,6)])
  missing_genes <- setdiff(genes, id_table$PROBE)
}


if (any(species %in% c("zebrafish", "rainbowtrout", "hamster", "rat", "mouse", "human", "chinesemedaka", "daphnia", "algae")) & Platform == "RNA-Seq") {
  
  if (file.exists(normalizePath(file.path(paths$RData, "id_table.RData"))) & params$use_cached_RData == TRUE) {
    
    print(paste("Already found ID table; loading from disk."))
    load(normalizePath(file.path(paths$RData, "id_table.RData")))
    id_table <- distinct(id_table_entrez[, 1:3])
    
  } else {
    
    id_table_entrez_all_genes <- readr::read_csv(normalizePath(file.path(paths$metadata, "genome", paste0(species, "_annotated_id_table_entrez_wholegenome.csv")))) #Note - You will have to include this file in the genome directory. Original file copy can be found in project dir.
    id_table_entrez <- id_table_entrez_all_genes %>%
      dplyr::filter(ensembl_gene_id %in% genes) %>%
      distinct()
    save(id_table_entrez, file = normalizePath(file.path(paths$RData, "id_table.RData")))
    id_table <- distinct(id_table_entrez[, 1:3])
    
  }
  
} else if (any(species %in% c("fatheadminnow")) & Platform == "RNA-Seq") {
  
  print("Fathead Minnow is not available in BiomaRt and Ensemble and Entrez IDs are not available... Only ZFIN IDs are available... Local annotation will begin.")
  
  if (file.exists(normalizePath(file.path(paths$RData, "id_table.RData"))) & params$use_cached_RData == TRUE) {
    
    print(paste("Already found ID table; loading from disk."))
    load(file.path(paths$RData, "id_table.RData"))
    id_table <- distinct(id_table_entrez[, 1:3])
    
  } else {
    
    id_table_entrez_all_genes <- readr::read_csv(normalizePath(file.path(paths$metadata, "genome", paste0(species, "_annotated_id_table_entrez_wholegenome.csv")))) #Note - You will have to include this file in the genome directory. Original file copy can be found in project dir.
    id_table_entrez <- id_table_entrez_all_genes %>%
      dplyr::filter(external_gene_name %in% genes) %>%
      distinct()
    save(id_table_entrez, file = file.path(paths$RData, "id_table.RData"))
    id_table <- distinct(id_table_entrez[, 1:3])
    id_table$ensembl_gene_id <- as.character(id_table$ensembl_gene_id)
    
  }
} else if (Platform == "RNA-Seq") {
  
  print("Species not recognized. Annotations will not be loaded")
  
} else {
  
  print("TempO-Seq data!")
  
}

# Check if id_table is empty
if (nrow(id_table) == 0) {
  print("Warning: id_table is empty.")
} 
if (nrow(id_table) != length(genes)) {
  print("Warning: nrow(id_table) is not equal to length(genes)")
}
if (Platform == "TempO-Seq") {
  print(paste("Warning: missing", length(missing_genes), "genes from id_table"))
}

print("id_table and id_table_entrez loaded successfully.")
