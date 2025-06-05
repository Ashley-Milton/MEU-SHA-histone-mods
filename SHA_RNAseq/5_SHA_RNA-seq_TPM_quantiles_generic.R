#Load required packages
library(tidyverse)

#Set working directory
setwd("/path/to/working/dir")

#Step 1: Read and process counts data --------------------------------------

#Read the counts data from featureCounts output
counts <- read.table("/path/to/featureCounts_s0.txt", header = TRUE) %>%
  select(-Start, -End)

#Step 2: Read and process GTF file for gene names --------------------------

#Read the GTF file containing gene annotations
GeneGFF <- read_delim("/path/to/genome/GCF_902635505.1_mSarHar1.11_genomic.gtf",
                      delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE, skip = 4) %>%
  filter(!grepl("##", X1))

#Extract relevant columns and separate gene information
GeneNames <- GeneGFF %>% filter(X3 == "gene") %>% select(X1, X4, X9)
colnames(GeneNames) <- c("Chromosome", "Start", "ID")

#Separate the ID column into key-value pairs and tidy up the data
GeneNames_separated <- GeneNames %>%
  separate_rows(ID, sep = "; ") %>%
  separate(ID, into = c("key", "value"), sep = " ", extra = "merge") %>%
  mutate(key = trimws(key), value = trimws(gsub('"', '', value))) %>%
  pivot_wider(names_from = key, values_from = value, values_fn = list(value = list)) %>%
  mutate(across(where(is.list), ~ sapply(., function(x) if (is.null(x)) NA else paste(unlist(x), collapse = "; ")))) %>%
  mutate(gene = ifelse(gene == "unknown" | gene == "Pc", gene_id, gene))

#Create a reference dataframe for gene IDs and names
GeneRef <- GeneNames_separated %>% select(gene_id, gene)

#Join counts data with gene names
SHACounts <- left_join(counts, GeneRef, by = c("Geneid" = "gene_id"))

#Set row names to the Geneid column and remove the Geneid column from data
rownames(SHACounts) <- make.names(SHACounts$Geneid, unique = TRUE)
SHACounts <- SHACounts %>% select(-Geneid)

#Rename specific columns for clarity
colnames(SHACounts)[c(4, 5, 6, 7)] <- c("Female_001", "Female_002", "Male_001", "Male_002")

#Separate counts data for females and males
female_counts <- SHACounts %>% select(Female_001, Female_002)
male_counts <- SHACounts %>% select(Male_001, Male_002)

#Extract gene lengths from the counts data
gene_lengths <- SHACounts$Length

#Step 3: Calculate TPM (Transcripts Per Million) ----------------------------

#Function to calculate TPM
calculate_tpm <- function(counts, gene_lengths) {
  #Calculate Reads Per Kilobase (RPK)
  rpk <- counts / gene_lengths
  #Calculate TPM by normalizing RPK values
  tpm <- t(t(rpk) / colSums(rpk) * 1e6)
  return(as.data.frame(tpm))
}

#Calculate TPM for females and males
female_tpm_df <- calculate_tpm(female_counts, gene_lengths)
male_tpm_df <- calculate_tpm(male_counts, gene_lengths)

#Step 4: Categorise genes into quantiles -----------------------------------

#Function to categorise genes
categorize_genes <- function(tpm_df) {
  tpm_df <- tpm_df %>%
    mutate(mean_tpm = rowMeans(tpm_df))
  
  quantiles <- quantile(tpm_df$mean_tpm[tpm_df$mean_tpm > 0], probs = seq(0, 1, length.out = 4))
  
  tpm_df <- tpm_df %>%
    mutate(category = case_when(
      mean_tpm == 0 ~ "Zero",
      mean_tpm > 0 & mean_tpm <= quantiles[2] ~ "One",
      mean_tpm > quantiles[2] & mean_tpm <= quantiles[3] ~ "Two",
      mean_tpm > quantiles[3] ~ "Three"
    ))
  
  return(tpm_df)
}

#Categorise genes for females and males
female_tpm_df <- categorize_genes(female_tpm_df)
male_tpm_df <- categorize_genes(male_tpm_df)

#Count the number of rows in each category
female_category_counts <- female_tpm_df %>%
  group_by(category) %>%
  summarize(count = n())

male_category_counts <- male_tpm_df %>%
  group_by(category) %>%
  summarize(count = n())

#Print the counts for each category
print(female_category_counts)
print(male_category_counts)

#Make minimal versions of the dataframes
female_exp_quantiles <- female_tpm_df %>% select(category)
male_exp_quantiles <- male_tpm_df %>% select(category)

#Save the minimal dataframes
write.table(female_exp_quantiles, file = "/path/to/outdir/female_exp_quantiles_geneid.tsv", row.names = TRUE, col.names = F, sep = "\t", quote = F)
write.table(male_exp_quantiles, file = "/path/to/outdir/male_exp_quantiles_geneid.tsv", row.names = TRUE, col.names = F, sep = "\t", quote = F)

#Print summaries of the TPM dataframes
print(summary(female_tpm_df))
print(summary(male_tpm_df))

