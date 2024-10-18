#!/root/.conda/envs/updated_my_r_pkgs/bin/R
# Custom parameters for the report
library(tidyverse)
require(yaml)

#Make the RMD folder if not already present
#And populate the config.qc.yml file and the R-ODAF_2_Sample_QC.rmd file

base_dir <- paste0("/mnt/", Sys.getenv("STUDY_ID_DIR")) #Setting here to the /mnt directory

print(paste("Base Directory:", base_dir))

config <- yaml::read_yaml(file.path(base_dir,
                                    "Rmd/config.qc.yml"),
                          eval.expr = T)

print(paste("Config file contents:", config))

# Manually evaluate the R expressions in the config (When they're set with !r)
#config$params$projectdir <- eval(parse(text=config$params$projectdir))
#config$params$project_name <- eval(parse(text=config$params$project_name))
#config$params$Platform <- eval(parse(text=config$params$Platform))
# Use the string values directly from the config (When they're set with !expr)
projectdir <- config$params$projectdir
project_name <- config$params$project_name
print(paste("config params projectdir:", config$params$projectdir))
print(paste("config params project_name:", config$params$project_name))
#print(paste("config params Platform:", config$params$Platform))

title <- paste('Study-wide sample quality control:', gsub(pattern = '_', replacement = ' ', x = config$params$project_name))

print(paste("Title:", title))

bibliography <- file.path(config$params$projectdir, 'Rmd', 'references.bib')

print(paste("Bibliography file:", bibliography))

# Input file - Rmd
files <- list.files(path = file.path(config$params$projectdir, "Rmd"), pattern = "R-ODAF_2_Sample_QC.*\\.rmd$", full.names = TRUE)
print(paste("Input Rmd file(s):", files))
# Check if we have found the file
if (length(files) == 1) {
  inputFile <- files[1]
} else if (length(files) > 1) {
  stop("Found more than one R-ODAF_2_Sample_QC*\\.rmd$ file")
} else {
  stop("No matching Rmd file found in the directory")
}
# Define the directory path
reports_dir <- file.path(config$params$projectdir, "reports")
# Check if the directory exists, and create it if it doesn't
if (!dir.exists(reports_dir)) {
  dir.create(reports_dir, recursive = TRUE)
}
  message("Writing QC report for all samples in the experiment.")
  # Output file - HTML
  filename <- paste0("Study-wide_Sample_QC.",
                     config$params$Platform, "_",
                     config$params$project_name, "_",
                     format(Sys.time(),'%d-%m-%Y.%H.%M'),
                     ".html")
  outFile <- file.path(reports_dir,
                       filename)
  rmarkdown::render(input = inputFile,
                    encoding = "UTF-8",
                    output_file = outFile,
                    params = config$params,
                    envir = new.env())