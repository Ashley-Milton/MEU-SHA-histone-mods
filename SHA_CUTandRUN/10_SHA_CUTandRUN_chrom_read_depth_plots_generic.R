#Set working directory
setwd("/path/to/Chrom_wide_summary_stats_output")

#Load packages for plotting
library(tidyverse)

#Load data
data <- read.table("/path/to/chrom_wide_summary_stats.txt", header = T)

#Read the file containing the number of mapped reads
mapped_reads <- read_delim("/path/to/subread_stats_final.txt", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE)

#Update Sex column in mapped_reads
mapped_reads <- mapped_reads %>%
  mutate(Sex = case_when(
    Sex == "F" ~ "Female",
    Sex == "M" ~ "Male",
    TRUE ~ Sex  
  ))

#Merge the data and mapped_reads dataframes based on Histone_mod and Sex
data <- data %>%
  left_join(mapped_reads, by = c("Histone_mod" = "Histone_Mod", "Sex" = "Sex"))

#Normalise the Average_read_depth by the number of mapped reads
data <- data %>%
  mutate(Normalised_read_depth = Average_read_depth / Mapped)

#Order the histone_mod column manually
data$Histone_mod <- factor(data$Histone_mod, levels = c("H3K9ac", "H3K9me3", "H3K27me3", "H4K20me3"))

#Load devil lookup table of chromosome names
lookup <- read_tsv("/path/to/Devil_lookup_table.tsv")

#Replace chromosome names
#Create a named vector for mapping
chromosome_mapping <- setNames(lookup$`Chromosome name`, lookup$`RefSeq seq accession`)

#Replace chromosome names in the dataframe
data$Chromosome <- ifelse(data$Chromosome %in% names(chromosome_mapping),
                          chromosome_mapping[data$Chromosome],
                          data$Chromosome)

#Double male Xs to account for their hemizygosity
data$Normalised_read_depth <- ifelse(data$Chromosome == "X" & data$Sex == "Male", data$Normalised_read_depth*2, data$Normalised_read_depth)

#Plot
ggplot(data %>% filter(Chromosome != "Y"), 
       aes(x = Chromosome, y = Normalised_read_depth, group = Sex, color = Sex)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  facet_wrap(~ Histone_mod, nrow = 1) +
  theme_minimal() +
  labs(title = "Per chromosome mean read depth in Tasmanian devil",
       x = "Chromosome", y = "Normalised read depth") +
  scale_color_manual(values = c("mediumorchid4","darkorange"))+
  theme(aspect.ratio = 1, legend.position = "top")

ggsave("/path/to/outdir/SHA_mean_chrom_read_depth_dot_plot.pdf", dpi = 600, width = 30, height = 15, units = "cm")
