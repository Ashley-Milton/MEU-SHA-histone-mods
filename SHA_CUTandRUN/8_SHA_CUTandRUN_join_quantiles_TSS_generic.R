#Set working directory
setwd("/path/to/plottingdata")

#Load packages
library(data.table)
library(tidyverse)

#Load exp quantiles data
F_expquant <- fread("/path/to/female_exp_quantiles_geneid.tsv", header = F, col.names = c("gene_id", "quantile"))
M_expquant <- fread("/path/to/male_exp_quantiles_geneid.tsv", header = F, col.names = c("gene_id", "quantile"))

#Define histone variable
histones <- c("H3K27me3", "H3K9ac", "H3K9me3", "H4K20me3")

#Loop for females
for (histone in histones) {
  
  #Load CUT&RUN data - full dataset
  CUTRUN <- fread(paste0("F_TD_", histone, "_info_headers.tsv"))
  
  #Extract the gene_id values and create a new column
  CUTRUN[, gene_id := sub('.*gene_id "([^"]+)".*', '\\1', gene)]
  
  #Define expquant data to use
  expquant <- F_expquant
  
  #Left join
  joined_data <- merge(CUTRUN, expquant, by = "gene_id", all.x = TRUE)
  
  #Clear up memory
  rm(CUTRUN)
  gc()
  
  #Replace all NAs with Zero
  joined_data[is.na(quantile), quantile := "Zero"]
  
  #Update the "c" column based on the rules
  joined_data[, c := fifelse(grepl("NC_045432.1", position), "X",
                             fifelse(grepl("NC_045433.1", position), "Y",
                                     fifelse(grepl("NC_018788.1", position), "MT", "A")))]
  
  #Cut down columns
  reduced_joined_data <- joined_data[, .(relative_pos, height, c, quantile)]
  
  #Clear up memory
  rm(joined_data)
  gc()
  
  #Save cut down file
  fwrite(reduced_joined_data, file = paste0("F_TD_", histone, "_info_headers_quantiles.tsv"), quote = F)
  
  #Clear up memory
  rm(reduced_joined_data)
  gc()
}

#Loop for males
for (histone in histones) {
  
  #Load CUT&RUN data - full dataset
  CUTRUN <- fread(paste0("M_TD_", histone, "_info_headers.tsv"))
  
  #Extract the gene_id values and create a new column
  CUTRUN[, gene_id := sub('.*gene_id "([^"]+)".*', '\\1', gene)]
  
  #Define expquant data to use
  expquant <- M_expquant
  
  #Left join
  joined_data <- merge(CUTRUN, expquant, by = "gene_id", all.x = TRUE)
  
  #Clear up memory
  rm(CUTRUN)
  gc()
  
  #Replace all NAs with Zero
  joined_data[is.na(quantile), quantile := "Zero"]
  
  #Update the "c" column based on the rules
  joined_data[, c := fifelse(grepl("NC_045432.1", position), "X",
                             fifelse(grepl("NC_045433.1", position), "Y",
                                     fifelse(grepl("NC_018788.1", position), "MT", "A")))]
  
  #Cut down columns
  reduced_joined_data <- joined_data[, .(relative_pos, height, c, quantile)]
  
  #Clear up memory
  rm(joined_data)
  gc()
  
  #Save cut down file
  fwrite(reduced_joined_data, file = paste0("M_TD_", histone, "_info_headers_quantiles.tsv"), quote = F)
  
  #Clear up memory
  rm(reduced_joined_data)
  gc()
}