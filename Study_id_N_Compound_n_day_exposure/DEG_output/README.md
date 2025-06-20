# Study_dir - DEG_output

After executing Part 3 - Differential gene expression analysis, several files and folder will be saved here, including excel spreadhseets, comma-separated sheets, bmdexpress input files, etc.

The folder sctructure should appear as so:

<pre>
Study_dir
|____DEG_output
      |____1.Study_dir_DESeq2_RNA-Seq_DESeq_by_gene.xlsx    #wide DESeq2 results table for all genes
      |____2.Study_dir_DESeq2_RNA-Seq_DESeq_all.xlsx  #long DESeq2 results table for all genes
      |____3.Study_dir_DESeq2_RNA-Seq_DESeq_by_contrast.xlsx      #a table of the number DEGs passing all filters that were found per contrast
      |____4.Study_dir_DESeq2_RNA-Seq-CPM.xlsx  #A matrix of counts per million reads of each gene.
      |____5.Study_dir_DESeq2_RNA-Seq-IPA.xlsx  #Input data for Ingenuity Pathway Analysis
      |____Study_dir_DESeq2_RNA-Seq-DESeq_output_ALL.txt    #long DESeq2 results table for all genes
      |____Study_dir_DESeq2_RNA-Seq-DESeq_output_all_genes.txt    #wide DESeq2 results table for all genes 
      |____Study_dir_DESeq2_RNA-Seq-DESeq_output_significant.txt  #DESeq2 results table for features that passed all filters.
      |____Study_dir_DESeq2_RNA-Seq-Per_sample_CPM.txt      #A matrix of counts per million reads of each gene.
      |____Study_dir_DESeq2_RNA-Seq-Per_sample_normalized_counts.txt    #a matrix of counts per gene, normalized by DESeq2
      |____<b>bmdexpress_input_log2_transformed.txt</b>     #Key input for bmdexpress. Copy to `/BMD_and_tPOD_calculation/RNASeqData/DESeq2_and_log2_norm_counts`
      |____bmdexpress_relevance_filtered_input_log2_transformed.txt
      | 
      |____pathway_analysis 
      |    |____#If enabled, pathway analysis files will appear here... disabled by default
      | 
      |____plots 
      |    |____Various_plot.pdfs 
      | 
      |____RData
      |    |____various_key_R_objects_including_metadata_and_DEG_results
      |  
      |____R-ODAF
           |____text_files_containing_genes_and_expression_data_from_each_R-ODAF_filter_for_each_contrast
			 |____...FDR_0.05_DEG_table.txt #0.05 FDR by default
			 |____...FDR_0.05_failed_DEspikes_table.txt
			 |____...FDR_0.05_failed_quantile_table.txt
			 |____...FDR_0.05_Normalized_and_relvance_filtered_data.txt
			 |____...FDR_0.05_Normalized_data.txt
</pre>