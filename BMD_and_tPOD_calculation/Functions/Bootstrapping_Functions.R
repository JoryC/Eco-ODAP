####Libraries####
library(tidyverse)

#### nth Gene Bootstrap ####
nth_gene_bootstrap <- function(x, seed = 1, nth_gene = 20, repeats = 2000) {
  set.seed(seed)
  boot_nth_gene <- replicate(repeats, {
    sampleData <- sample(unlist(x), length(unlist(x)), replace = TRUE)
    sort(sampleData)[nth_gene]
  })
  return(quantile(boot_nth_gene, probs = c(0.025, 0.5, 0.975)))
}

#### nth Percentile Bootstrap ####
nth_percent_bootstrap <- function(x, seed = 1, nth_percent = 10, repeats = 2000) {
  set.seed(seed)
  boot_nth_percent <- replicate(repeats, {
    sampleData <- sample(unlist(x), length(unlist(x)), replace = TRUE)
    quantile(sort(sampleData), probs = (nth_percent / 100))
  })
  return(quantile(boot_nth_percent, probs = c(0.025, 0.5, 0.975)))
}

#### First Mode Bootstrap ####
mode_bootstrap <- function(x, seed = 1, repeats = 2000, min_size = 0.06, min_bw = 0.15) {
  #source("mode_antimode.R")
  set.seed(seed)
  boot_mode <- replicate(repeats, {
    sampleData <- sample(unlist(x), length(unlist(x)), replace = TRUE)
    dataMode <- mode.antimode(sampleData, min.size = min_size, bw = "SJ", min.bw = min_bw)
    dataMode$modes[[1]]
  })
  return(quantile(boot_mode, probs = c(0.025, 0.5, 0.975)))
}

# Define the function to calculate LCRD with the largest CRGB
LCRD <- function(bmc, probe, cut = 1.778, logbase = 10) {
  # Order BMC values and their corresponding probes
  ordered_indices <- order(bmc)
  sorted_bmc <- bmc[ordered_indices]
  sorted_probe <- probe[ordered_indices]
  
  # Calculate ratios
  ratios <- sorted_bmc[-1] / sorted_bmc[-length(sorted_bmc)]
  
  # Initialize variables to store CRGB information
  crgb_list <- list()
  current_start <- 1
  
  # Iterate through the ratios to find CRGBs
  for (i in seq_along(ratios)) {
    if (ratios[i] >= cut) {
      # Append each CRGB together into a list
      crgb_list <- append(crgb_list, list(list(start = current_start, end = i)))
      current_start <- i + 1
    }
  }
  # Append the last CRGB
  crgb_list <- append(crgb_list, list(list(start = current_start, end = length(sorted_bmc))))
  
  # Identify the largest CRGB
  crgb_lengths <- sapply(crgb_list, function(group) group$end - group$start + 1)
  largest_crgb_index <- which.max(crgb_lengths)
  largest_crgb <- crgb_list[[largest_crgb_index]]
  
  # Extract BMCs and probes of the largest CRGB
  largest_crgb_bmcs <- sorted_bmc[largest_crgb$start:largest_crgb$end]
  largest_crgb_probes <- sorted_probe[largest_crgb$start:largest_crgb$end]
  
  # Identify the lowest BMC in the largest CRGB
  lcrd_bmc <- min(largest_crgb_bmcs)
  lcrd_probe <- largest_crgb_probes[which.min(largest_crgb_bmcs)]
  
  # Extract the lowest BMC for each CRGB
  crgb_bmcs_df <- do.call(rbind, lapply(seq_along(crgb_list), function(i) {
    group <- crgb_list[[i]]
    crgb_bmcs <- sorted_bmc[group$start:group$end]
    crgb_probes <- sorted_probe[group$start:group$end]
    data.frame(
      Gene = crgb_probes[which.min(crgb_bmcs)],
      BMC = min(crgb_bmcs),
      logBMC = log(min(crgb_bmcs), base = logbase),
      CRGB_Number = i,
      CRGB_Size = length(crgb_bmcs)
    )
  }))
  
  # Return the result
  list(
    LCRD_Result = data.frame(
      Gene = lcrd_probe,
      BMC = lcrd_bmc,
      logBMC = log(lcrd_bmc, base = logbase),
      CRGB_Number = largest_crgb_index,
      CRGB_Size = length(largest_crgb_bmcs)
    ),
    CRGB_BMCs = crgb_bmcs_df
  )
}

### LCRD bootstrap ###
lcrd_bootstrap <- function(x, seed = 1, repeats = 2000, lcrdratiocut = 1.778, lcrdlogbase = 10) {
  set.seed(seed)
  bmds <- x[["BMD"]]
  names(bmds) <- x[["Probe ID"]]
  if (all(names(bmds) == x[["Probe ID"]], na.rm = TRUE) &&
      all(bmds == x[["BMD"]])) {
    boot_lcrd <- replicate(repeats, {
      sampleData <- sample(unlist(bmds), length(unlist(bmds)), replace = TRUE)
      result <-
        LCRD(bmc = sampleData,
             probe = names(sampleData),
             cut = lcrdratiocut,
             logbase = lcrdlogbase)
      result$LCRD_Result$logBMC
    })
    return(quantile(boot_lcrd, probs = c(0.025, 0.5, 0.975)))
  } else {
    print("Probe IDs do not match BMDs, please review Bootstrapping_Funtions.R")
  }
}

#### Pathway Bootstrap ####
pathway_bootstrap <- function(x, seed = 1, repeats = 2000) {
  set.seed(seed)
  boot_pathway <- replicate(repeats, {
    sampleData <- sample(unlist(x), length(unlist(x)), replace = TRUE)
    median(sampleData)
  })
  return(quantile(boot_pathway, probs = c(0.025, 0.5, 0.975)))
}

#### Average CI values ####
averageCI <- function(x) {
  mean_values <- colMeans(do.call(rbind, x), na.rm = TRUE)
  names(mean_values) <- c("2.5%", "50%", "97.5%")
  return(mean_values)
}