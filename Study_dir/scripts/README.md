# Study_directory - scripts folder

This folder contains three executable script files:

- Part 1 of the R-ODAF - Preprocessing of the sequencing data (WARNING: not TempO-Seq capable) (`R-ODAF_1_sequencing_DataPreProcess.sh`)
  - Shell script (`.sh`)
  - executed by submitting `/submit_R_ODAF_job_1.sh` or `/submit_R_ODAF_job_1_array.sh` via the SLURM command `sbatch` or `srun`
- Part 2 of the R-ODAF - Quality control (`R-ODAF_2_render_studywide_QC_report.R`)
  - R script (`.R`)
  - renders an html report and outputs from `/Study_dir/Rmd/R-ODAF_2_Sample_QC.rmd`
  - executed by submitting `/submit_R_ODAF_job_2.sh` or `/submit_R_ODAF_job_2_array.sh` via the SLURM command `sbatch` or `srun`
  - Please note that it may be necessary for you to update the Shebang on line 1 which points to your `updated_my_r_pkgs` conda environment.
- Part 3 (final) of the R-ODAF - Differential gene expression analysis using DESeq2 (`R-ODAF_3_render_DESeq2_report.R`)
  - R script (`.R`)
  - renders an html report and outputs from `/Study_dir/Rmd/R-ODAF_3_DESeq2_report.rnaseq.Rmd`
  - executed by submitting `/submit_R_ODAF_job_3.sh` or `/submit_R_ODAF_job_3_array.sh` via the SLURM command `sbatch` or `srun`
  - Please note that it may be necessary for you to update the Shebang on line 1 which points to your `updated_my_r_pkgs` conda environment.