# Load the necessary package
library(rmarkdown)
library(tidyverse)
library(vtree)
library(here) # For better path handling
library(grDevices)  # For controlling graphics devices

# Define the base filename format and chemical names
base_filename <- "BMD_and_tPOD_calculation_"
chemnames <- c("chem1", "chem2", "chem3", "chem4", "chem5", "chem6")

#Sys.sleep(2700) #delaying execution, this can be customized if need be. Not very complex but got the job done for me

# Set up a custom temporary directory for vtree
custom_tempdir <- tempdir()
options(vtree.tempdir = custom_tempdir)

# Function to safely close graphics devices
safe_dev_off <- function() {
  while (dev.cur() > 1) dev.off()
}

# Loop through each chemical name and render the corresponding .Rmd file
# Render each Rmd file in a loop
for (chem in chemnames) {
  filename <- paste0(base_filename, chem, ".Rmd")
  
  tryCatch({
    cat("Rendering:", filename, "\n")
    
    # Reset the graphics device before rendering
    safe_dev_off()  # Close any open graphics devices
    gc()       # Run garbage collection
    
    # Render the Rmd file
    render(
      input = filename,
      output_dir = here(),
      envir = new.env()  # Use a new environment for each render
    )
    
    Sys.sleep(10)  # Add delay to prevent file conflicts
    
  }, error = function(e) {
    cat("Error in rendering:", filename, "\n", conditionMessage(e), "\n")
    
    # Print diagnostic information
    cat("Using temp directory:", custom_tempdir, "\n")
    cat("Current working directory:", getwd(), "\n")
    
    # Capture and print the output of temp files
    temp_files <- list.files(custom_tempdir, full.names = TRUE)
    if (length(temp_files) > 0) {
      cat("Temp files:\n")
      print(temp_files)
    } else {
      cat("No temp files found.\n")
    }
    
    # Clean up graphics devices and temp files
    safe_dev_off()  # Ensure all graphics devices are closed
    unlink(file.path(custom_tempdir, "*"), recursive = TRUE)  # Clear temp files
  })
  
  # Clean up after each render to avoid conflicts
  safe_dev_off()
  gc()
  unlink(file.path(custom_tempdir, "*"), recursive = TRUE)
}
