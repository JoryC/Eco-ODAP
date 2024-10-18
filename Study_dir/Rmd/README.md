# Study_dir - Rmd

This folder contains 5 files:

- Part 2 of the R-ODAF - Sample quality control (`R-ODAF_2_Sample_QC.rmd`)
  - R markdown file (`.rmd`)
  - rendered by `/Study_dir/scripts/R-ODAF_2_render_studywide_QC_report.R` which is executed by submitting the job `/submit_R_ODAF_job_2.sh` or `/submit_R_ODAF_job_2_array.sh` via the SLURM command `sbatch`.
  - customizable paramaters in the YAML header. Input params should match the configuration params in `/Study_dir/Rmd/config.qc.yml`
- The configuration file for sample quality control (`config.qc.yml`)
  - YAML file (`.yml`)
  - a crucial input for rendering the sample quality control report and outputs. Should match the params in `/Study_dir/Rmd/R-ODAF_2_Sample_QC.rmd`
- Part 3 of the R-ODAF - Differential gene expression analysis via DESeq2 (`R-ODAF_3_DESeq2_report.rnaseq.Rmd`)
 - R markdown file (`.rmd`)
 - rendered by `/Study_dir/scripts/R-ODAF_3_render_studywide_QC_report.R` which is executed by submitting the job `/submit_R_ODAF_job_2.sh` or `/submit_R_ODAF_job_2_array.sh` via the SLURM command `sbatch`.
- The configuration file for differential gene expression analysis (`config.yml`)
  - YAML file (`.yml`)
  - a crucial input for rendering the differential gene expression report and outputs. Should match the params in `/Study_dir/Rmd/R-ODAF_3_DESeq2_report.rnaseq.Rmd`
- The bibliography file for citing the original authors of the R-ODAF and critical package authors (`references.bib`)
  - citations are rendered at the end of both reports