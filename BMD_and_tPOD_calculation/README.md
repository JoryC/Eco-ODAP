# BMD_and_tPOD_calculation

Once Part 3 - Differential gene expression analysis of the R-ODAF is complete, the final part of the pipeline is to calulate BMDs for dose-responsive genes and derive tPODs.

Dependencies:
- [BMDExpress3](https://github.com/auerbachs/BMDExpress-3/releases) - [command line tool](https://github.com/auerbachs/BMDExpress-3/wiki/Command-Line)
  - [Add bmdexpress3-cmd to your $PATH](https://www.youtube.com/watch?v=NZgHnV3ZXaw&list=PLX2Rd5DjtiTeR84Z4wRSUmKYMoAbilZEc&index=18)
- Rfast
- tidyverse
- purrr
- knitr
- kableExtra
- sessioninfo
- jsonlite
- vtree
- data.table

To calculate BMDs and derive tPODs:
- copy the `/Study_dir/DEG_output/bmdexpress_input_log2_transformed.txt` to the `/BMD_and_tPOD_calculation/RNASeqData/DESeq2_and_log2_norm_counts/` directory, and copy the `SraRunTable.csv` metadata file to the `/BMD_and_tPOD_calculation/RNASeqData/metadata/` directory.
- populate the paramaters in the YAML header of `/BMD_and_tPOD_calculation/BMD_and_tPOD_calculation.Rmd`.
- knit/render `/BMD_and_tPOD_calculation/BMD_and_tPOD_calculation.Rmd`.

The metadata file (either metadata.QC_applied.csv or SraRunTable.csv) must minimally have:
- A "Run" column (unique sample name).
- An "Organism" column.
- A "Days" column, corresponding to the number of exposure days.
- A compound/chemical column (flexible column name).
- A dose/concentration column (flexible column name).
- A "DoseUnits" column.

Once completed, the html report will appear here, and various outputs in the `/BMD_and_tPOD_calculation/Outputs/` directory. The current directory should appear as follows:
**Necessary inputs are bolded**.

<pre>
.
|____BMD_and_tPOD_calculation
	 |____<b>BMD_and_tPOD_calculation.Rmd</b> # R markdown script that exedcutes BMDExpress and calculates tPODs
	 |____BMD_and_tPOD_calculation.html # Output report
	 |
	 |____BMDExpressData
	 |    |____BMDExpressFile.bm2
	 |    |____BMDExpressFile_zf_homologs.bm2 (if applicable)
	 |    |____BMDExpressFile.json #Configuration file, created by `BMD_and_tPOD_calculation.Rmd`
	 |    |____BMDExpressFile_zf_homolgs.json (if applicable) #Configuration file, created by `BMD_and_tPOD_calculation.Rmd`
	 |    |
	 |    |____BMD
	 |    |    |____mybmddata.txt
	 |    |
	 |    |____GO_and_REACTOME
	 |    |____mycategoricaldata.txt
	 |    |____mycategoricaldata_zf_homolgs.txt (if applicable)
	 |	  |
	 |	  |____Williams
	 |    	   |____mywilliamsdata.txt
	 |
	 |____Functions
	 |	  |____<b>BMDExpressFunctions.R</b>
	 |	  |____<b>Bootstrapping_Functions.R</b>
	 |	  |____<b>mode_antimode.R</b>
	 |
	 |____Outputs
	 |	  |____all_BMD_list_logtransformed_R-ODAF_DEGs.RDS
	 |	  |____BMD_Data_for_Histogram_R-ODAF_DEGs.csv
	 |	  |____BMD_Histogram.pdf
	 |	  |____BMD_summary_table.RDS
	 |	  |____final_tPOD_table.csv
	 |	  |____GO_accumulation_plot.pdf
	 |	  |____go_BMD_list_logtransformed_R-ODAF_DEGsBPF.RDS
	 |	  |____GO_Data_R-ODAF_DEGs.csv
	 |	  |____go_raw_data_filtered_final_results_with_logBMDs.RDS
	 |	  |____go_raw_data_filtered_results.RDS
	 |	  |____GO_summary_table.RDS
	 |	  |____raw_data_filtered.RDS
	 |	  |____raw_data_filtered_results_.RDS
	 |	  |____REACTOME_accumulation_plot.pdf
	 |	  |____reactome_BMD_list_logtransformed_R-ODAF_DEGs.RDS
	 |	  |____Reactome_Data_R-ODAF_DEGs.csv
	 |	  |____reactome_raw_data_filtered_final_results_with_logBMDs.RDS
	 |	  |____reactome_raw_data_filtered_results.RDS
	 |	  |____Reactome_summary_table.RDS
	 |	  |____tPOD_Summary_Plot.pdf
	 |	  |____tpod_values.txt
	 |
	 |____RNASeqData
	 	  |____DESeq2_norm_and_log2_norm_counts
		  |	   |____<b>bmdexpress_input_log2_transformed.txt</b> # Input data
		  |
		  |____metadata
		  	   |____<b>metadata.csv</b> # Input metadata
</pre>

As a bonus, I have included `/BMD_and_tPOD_calculation/Render_multiple_reports.R` which is an RScript that can render multiple reports in a loop if you have many input files to analyze.

Note: BMDExpress is not a supported module on Digital Research Alliance of Canada (DRAC) servers. You will either need to contact someone and request a module build for a specific version of BMDExpress, or run it locally on a laptop/PC. For this reason I have not included a conda environment for this portion of the pipeline and none of the necessary program files and dependencies are contained in the apptainer image. 