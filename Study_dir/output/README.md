# Study_dir - output

Output files will appear here.

When Part 1 (data pre-processing) of the R-ODAF is completed the directory should appear as so:

<pre>
Study_dir
|____output 
      |____QC_per_sample.txt
      |____metadata.QC_applied.txt
      |____contrasts.txt
      |____samples_removed.txt 
      |____Plots_output_from_Sample_QC.pdf/.png 
      |____log.out
	  |____conda_environments.yml
      | 
      |____Trimmed_reads 
            |____SampleA_R1_trimmed.fq.gz 
            |____SampleA_R2_trimmed.fq.gz 
            | 
            |____fastQCoutput 
            |    |____Xfastp.html #used in MultiQC 
            |    |____X_fastp.json #used in MultiQC 
            | 
            |____RSEM 
            |    |____genes.data.tsv 
            |    |____isoforms.data.tsv 
            |    |____various X.genes.results and X.isoforms.results files #Each corresponding to one sample 
            |  
            |____MultiQC 
            |    |____MultiQC_Report.html #Review the report by opening in browser. Move to /Study_dir/reports if desired.
            |         | 
            |         |____MultiQC_Report_data 
            |              |____multiqc_general_stats.txt
            |              |____multiqc.log
            |              |____multiqc_data.json
            |              |____multiqc_sources.txt 
            |              |____multiqc_star.txt
            |              |____multiqc_rsem.txt
            |              |____multiqc_fastp.txt		
            | 
            |____STAR 
                 |____various_individual_STAR_files_used_by_MultiQC_and_RSEM_including_aligned_bam_files
</pre>