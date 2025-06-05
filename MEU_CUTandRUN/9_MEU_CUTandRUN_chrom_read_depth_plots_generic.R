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

#Filter for one tissue for now (Brain)
data_B <- data %>% filter(Tissue == "Brain") %>% select(!Tissue)

mapped_reads_B <- mapped_reads %>% filter(Tissue == "Br") %>% select(!Tissue)

#Merge the data_B and mapped_reads_B dataframes based on Histone_mod and Sex
data_B <- data_B %>%
  left_join(mapped_reads_B, by = c("Histone_mod" = "Histone_Mod", "Sex" = "Sex"))

#Normalise the Average_read_depth by the number of mapped reads
data_B <- data_B %>%
  mutate(Normalised_read_depth = Average_read_depth / Mapped)

#Order the histone_mod column manually
data_B$Histone_mod <- factor(data_B$Histone_mod, levels = c("H3K4me3", "H3K9me3", "H3K27me3", "H4K20me3"))

#Remove the 'chr' prefix from the Chromosome column
data_B$Chromosome <- sub("^chr", "", data_B$Chromosome)

#Double male Xs to account for their hemizygosity
data_B$Normalised_read_depth <- ifelse(data_B$Chromosome == "X" & data_B$Sex == "Male", data_B$Normalised_read_depth*2, data_B$Normalised_read_depth)
    
#Filter for one tissue for now (Kidney)
data_K <- data %>% filter(Tissue == "Kidney") %>% select(!Tissue)

mapped_reads_K <- mapped_reads %>% filter(Tissue == "Ki") %>% select(!Tissue)

#Merge the data_K and mapped_reads_K dataframes based on Histone_mod and Sex
data_K <- data_K %>%
  left_join(mapped_reads_K, by = c("Histone_mod" = "Histone_Mod", "Sex" = "Sex"))

#Normalise the Average_read_depth by the number of mapped reads
data_K <- data_K %>%
  mutate(Normalised_read_depth = Average_read_depth / Mapped)

#Order the histone_mod column manually
data_K$Histone_mod <- factor(data_K$Histone_mod, levels = c("H3K4me3", "H3K9me3", "H3K27me3", "H4K20me3"))

#Remove the 'chr' prefix from the Chromosome column
data_K$Chromosome <- sub("^chr", "", data_K$Chromosome)

#Double male Xs to account for their hemizygosity
data_K$Normalised_read_depth <- ifelse(data_K$Chromosome == "X" & data_K$Sex == "Male", data_K$Normalised_read_depth*2, data_K$Normalised_read_depth)

#Add tissue labels
data_B <- data_B %>% mutate(Tissue = "Brain")
data_K <- data_K %>% mutate(Tissue = "Kidney")

#Combine datasets
combined_data <- bind_rows(data_B, data_K)

#Filter combined data to remove the Y
combined_data_filtered <- combined_data %>% filter(Chromosome != "Y")

ggplot(combined_data_filtered, 
       aes(x = Chromosome, y = Normalised_read_depth, group = interaction(Sex, Tissue), shape = Tissue)) +
  geom_point(aes(color = Sex), size = 3) +
  geom_line(aes(color = Sex, linetype = Tissue), linewidth = 1) +
  facet_wrap(~ Histone_mod, nrow = 1) +
  theme_minimal() +
  labs(title = "Per chromosome mean read depth in tammar wallaby",
       x = "Chromosome",
       y = "Normalised read depth",
       color = "Sex",
       shape = "Tissue",
       linetype = "Tissue") +
  scale_color_manual(values = c("mediumorchid4","darkorange")) +
  theme(aspect.ratio = 1,
    legend.position = "top"
  )


ggsave("/path/to/outdir/MEU_mean_chrom_read_depth_dot_plot.pdf", dpi = 600, width = 30, height = 15, units = "cm")
