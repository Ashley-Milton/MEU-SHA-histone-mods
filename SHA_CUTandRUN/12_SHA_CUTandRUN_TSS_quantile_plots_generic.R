setwd("/path/to/normal/plottingdata")

library(tidyverse)
library(data.table)
library(scales)
library(readr)

#Define histone variable
histones <- c("H3K27me3", "H3K9ac", "H3K9me3", "H4K20me3")

#Read the file containing the number of mapped reads
mapped_reads <- read_delim("/path/to/subread_stats_final.txt", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE)



#Initialise an empty data frame to store the results
results <- data.frame(histone = character(), mapped_reads_F = numeric(), mapped_reads_M = numeric(), ratio = numeric())

#Define the window size and the range
window_size <- 50
range_start <- -10000
range_end <- 10000
windows <- seq(range_start, range_end, by = window_size)

#Function to calculate the ratios in windows
calculate_window_ratios <- function(data, data_bg, windows, window_size) {
  ratios <- data.frame() #Initialise an empty data frame to store the results
  
  for (w in windows) {
    #Filter data for the current window
    window_data <- data %>% filter(relative_pos >= w & relative_pos < w + window_size)
    window_data_bg <- data_bg %>% filter(relative_pos >= w & relative_pos < w + window_size)
    
    #Summarize data by sex_c and quantile
    window_sum <- window_data %>% group_by(sex_c, quantile) %>% summarise(sum_height = sum(mean_height), .groups = "drop")
    window_sum_bg <- window_data_bg %>% group_by(sex_c, quantile) %>% summarise(sum_height_bg = sum(mean_height), .groups = "drop")
    
    #Calculate the ratio
    window_ratio <- merge(window_sum, window_sum_bg, by = c("sex_c", "quantile")) %>% 
      mutate(ratio = sum_height / sum_height_bg) %>% 
      mutate(window_center = w + window_size / 2)
    
    #Append the results to the ratios data frame
    ratios <- rbind(ratios, window_ratio)
  }
  
  return(ratios)
}

for (histone in histones) {

  #Read data normal
  data_F <- fread(paste("F_TD_", histone, "_info_headers_quantiles.tsv", sep = ""), fill = T)
  data_M <- fread(paste("M_TD_", histone, "_info_headers_quantiles.tsv", sep = ""), fill = T)
  
  #Read data background
  setwd("/path/to/background/plottingdata")
  
  data_F_bg <- fread(paste("F_TD_", histone, "_info_headers_quantiles_lambda.tsv", sep = ""), fill = T)
  data_M_bg <- fread(paste("M_TD_", histone, "_info_headers_quantiles_lambda.tsv", sep = ""), fill = T)
  
  setwd("/path/to/normal/plottingdata")
  
  #Normalising
  
  #Get the number of mapped reads for the current histone and sex
  mapped_reads_F <- filter(mapped_reads, Histone_Mod == histone, Sex == "F")$Mapped
  mapped_reads_M <- filter(mapped_reads, Histone_Mod == histone, Sex == "M")$Mapped
  
  #Calculate the ratio of male to female reads
  ratio <- mapped_reads_M / mapped_reads_F
  
  #Add the results to the data frame
  results <- rbind(results, data.frame(histone = histone, mapped_reads_F = mapped_reads_F, mapped_reads_M = mapped_reads_M, ratio = ratio))
  
  #Calculate the scaling factor
  scaling_factor <- max(mapped_reads_F, mapped_reads_M) / min(mapped_reads_F, mapped_reads_M)
  
  #Scaling normal data
  #Determine which dataset has fewer reads
  if (mapped_reads_F < mapped_reads_M) {
    #If female data has fewer reads, scale up the female data
    data_F$height <- data_F$height * scaling_factor
  } else {
    #If male data has fewer reads, scale up the male data
    data_M$height <- data_M$height * scaling_factor
  }
  
  #Scaling background data
  if (mapped_reads_F < mapped_reads_M) {
    #If female data has fewer reads, scale up the female data
    data_F_bg$height <- data_F_bg$height * scaling_factor
  } else {
    #If male data has fewer reads, scale up the male data
    data_M_bg$height <- data_M_bg$height * scaling_factor
  }
  
  
  #TSS plots ---------------------------------------------------------------
  
  #Filtering for at least two datapoints per relative position, and
  #taking mean value at each position relative to TSS
  
  ##NORMAL DATA
  
  #Female
  meddata_F <- data_F %>% filter(c == "A" | c == "X") %>% group_by(relative_pos, c, quantile) %>% mutate(count=n()) %>% filter(count > 2) %>% summarise(mean_height=mean(height)) #Median height doesn't work here because there are too many low values
  meddata_F$mod <- histone
  meddata_F$sex <- "Female"
  meddata_F <- meddata_F %>% ungroup()
  
  #Male
  meddata_M <- data_M %>% filter(c == "A" | c == "X") %>% group_by(relative_pos, c, quantile) %>% mutate(count=n()) %>% filter(count > 2) %>% summarise(mean_height=mean(height)) #Median height doesn't work here because there are too many low values
  meddata_M$mod <- histone
  meddata_M$sex <- "Male"
  meddata_M <- meddata_M %>% ungroup()
  
  #Doubling the mean_height column for male X chromosomes to account for their hemizygosity
  meddata_M$mean_height <- ifelse(meddata_M$c == "X", meddata_M$mean_height*2, meddata_M$mean_height)
  
  #Putting female and male data together in one dataframe
  meddata <- rbind(meddata_F, meddata_M)
  meddata$sex_c <- paste(meddata$sex, meddata$c, sep = " ")
  
  meddata$relative_pos <- as.numeric(meddata$relative_pos)
  
  ##BACKGROUND DATA
  
  #Female
  meddata_F_bg <- data_F_bg %>% filter(c == "A" | c == "X") %>% group_by(relative_pos, c, quantile) %>% mutate(count=n()) %>% filter(count > 2) %>% summarise(mean_height=mean(height)) #Median height doesn't work here because there are too many low values
  meddata_F_bg$mod <- histone
  meddata_F_bg$sex <- "Female"
  meddata_F_bg <- meddata_F_bg %>% ungroup()
  
  #Male
  meddata_M_bg <- data_M_bg %>% filter(c == "A" | c == "X") %>% group_by(relative_pos, c, quantile) %>% mutate(count=n()) %>% filter(count > 2) %>% summarise(mean_height=mean(height)) #Median height doesn't work here because there are too many low values
  meddata_M_bg$mod <- histone
  meddata_M_bg$sex <- "Male"
  meddata_M_bg <- meddata_M_bg %>% ungroup()
  
  #Doubling the mean_height column for male X chromosomes to account for their hemizygosity
  meddata_M$mean_height <- ifelse(meddata_M$c == "X", meddata_M$mean_height*2, meddata_M$mean_height)
  
  #Putting female and male data together in one dataframe
  meddata_bg <- rbind(meddata_F_bg, meddata_M_bg)
  meddata_bg$sex_c <- paste(meddata_bg$sex, meddata_bg$c, sep = " ")
  
  meddata_bg$relative_pos <- as.numeric(meddata_bg$relative_pos)
  
  
  #Calculate ratios
  window_ratios <- calculate_window_ratios(meddata, meddata_bg, windows, window_size)
  
  #Modify the data to group autosomes by quantiles but not X chromosomes
  window_ratios <- window_ratios %>%
    mutate(group = ifelse(grepl("A", sex_c), paste(sex_c, quantile, sep = " "), sex_c))
  
  #Add a new column 'sex' with the first word of 'sex_c'
  window_ratios$sex <- sub(" .*", "", window_ratios$sex_c)
  
  #Set the desired order for the group column
  desired_group_order <- c(
    "Female A Zero", "Female A One", "Female A Two", "Female A Three", #Female autosomes by quantile
    "Female X",                                                       #Female X
    "Male A Zero", "Male A One", "Male A Two", "Male A Three",        #Male autosomes by quantile
    "Male X"                                                          #Male X
  )
  
  #Update the group column to be a factor with the desired order
  window_ratios <- window_ratios %>%
    mutate(group = factor(group, levels = desired_group_order))
  
  
  #Plot
  #Generate gradients for female (red-orange) and male (blue-teal)
  female_colors <- colorRampPalette(c("mediumorchid1", "deeppink","darkred", "orangered", "darkgoldenrod1"))(5) 
  male_colors <- colorRampPalette(c("gray15", "navy", "dodgerblue", "seagreen", "green3"))(5) 
  
  #Combine the gradients
  all_colors <- c(female_colors, male_colors)
  
  #Create separate plot for male
  window_ratios %>% filter(sex == "Male") %>% ggplot(aes(x = window_center, y = ratio, color = group)) +
    geom_smooth(aes(fill = group), se = T, method = "loess", span = 0.35) +
    labs(
      x = "Position relative to TSS (kb)",
      y = "Ratio (sample/background) of mean peak height"
    ) +
    coord_cartesian(x = c(-10000, 10000), y = c(0, 2)) +
    theme_classic() +
    scale_color_manual(values = male_colors, name = "") +
    scale_fill_manual(values = male_colors, name = "") +
    scale_x_continuous(labels = unit_format(unit = "", scale = 1e-3)) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 14),
      strip.text = element_text(size = 16),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 14),
      aspect.ratio = 1,
      legend.title = element_blank())
  
  ggsave(paste("/path/to/R_plots/TSS_", histone, "_backgroundratio_plot_mappedreads_SHA_maledoubled_quantiles_As_y_0_2_male.pdf", sep = ""), width = 7, height = 7, dpi = 800)
  
  #Create separate plot for female
  window_ratios %>% filter(sex == "Female") %>% ggplot(aes(x = window_center, y = ratio, color = group)) +
    geom_smooth(aes(fill = group), se = T, method = "loess", span = 0.35) +
    labs(
      x = "Position relative to TSS (kb)",
      y = "Ratio (sample/background) of mean peak height"
    ) +
    coord_cartesian(x = c(-10000, 10000), y = c(0, 2)) +
    theme_classic() +
    scale_color_manual(values = female_colors, name = "") +
    scale_fill_manual(values = female_colors, name = "") +
    scale_x_continuous(labels = unit_format(unit = "", scale = 1e-3)) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 14),
      strip.text = element_text(size = 16),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 14),
      aspect.ratio = 1,
      legend.title = element_blank())
  
  ggsave(paste("/path/to/R_plots/TSS_", histone, "_backgroundratio_plot_mappedreads_SHA_maledoubled_quantiles_As_y_0_2_female.pdf", sep = ""), width = 7, height = 7, dpi = 800)
  
  #Plot the data both sexes
  window_ratios %>% ggplot(aes(x = window_center, y = ratio, color = group)) +
    geom_smooth(aes(fill = group), se = TRUE, method = "loess", span = 0.35) +
    labs(x = "Position relative to TSS (kb)", y=paste("Ratio (sample/background) of mean peak height for", histone)) +
    coord_cartesian(x = c(-10000, 10000), y = c(0, 2)) +
    theme_classic() +
    scale_color_manual(values = all_colors, name = "") +
    scale_fill_manual(values = all_colors, name = "") +
    scale_x_continuous(labels = unit_format(unit = "", scale = 1e-3)) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 14),
      strip.text = element_text(size = 16),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 14),
      aspect.ratio = 1,
      legend.title = element_blank()
    )
  
  
  ggsave(paste("/path/to/R_plots/TSS_", histone, "_backgroundratio_plot_mappedreads_SHA_maledoubled_quantiles_As_y_0_2.pdf", sep = ""), width = 7, height = 7, dpi = 800)
}

write.csv(results, "/path/to/R_plots/mappedreads_scaling_backgroundratio_SHA_maledoubled_quantiles_As.csv", row.names = FALSE)
