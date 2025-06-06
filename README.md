# MEU-SHA-histone-mods

## Overview

This repository contains scripts for processing and analysing CUT&RUN and RNA-seq data from tammar wallaby (*Macropus eugenii*, MEU) and Tasmanian devil (*Sarcophilus harrisii*, SHA). The pipelines use [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [MultiQC](https://multiqc.info/), [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), [Subread](http://subread.sourceforge.net/), [MACS3](https://github.com/macs3-project/MACS), [bedtools](https://bedtools.readthedocs.io/en/latest/), and custom R scripts for data processing and visualisation.

## CUT&RUN

### MEU (tammar wallaby)  
Scripts for MEU CUT&RUN data include:

- Quality control:  
  [`1_MEU_CUTandRUN_fastqc_multiqc_generic.pbs`](./MEU_CUTandRUN/1_MEU_CUTandRUN_fastqc_multiqc_generic.pbs)  
- Read trimming:  
  [`2_MEU_CUTandRUN_trimmomatic_generic.pbs`](./MEU_CUTandRUN/2_MEU_CUTandRUN_trimmomatic_generic.pbs)  
- Read alignment to reference genome:  
  [`3_MEU_CUTandRUN_parallel-subread-PE_generic.sh`](./MEU_CUTandRUN/3_MEU_CUTandRUN_parallel-subread-PE_generic.sh)  
- Peak calling:  
  [`4_MEU_CUTandRUN_macs3_generic.pbs`](./MEU_CUTandRUN/4_MEU_CUTandRUN_macs3_generic.pbs)  
- Setup for intersect directories:  
  [`5_MEU_CUTandRUN_setup_intersect_directories_generic.sh`](./MEU_CUTandRUN/5_MEU_CUTandRUN_setup_intersect_directories_generic.sh)  
- Bedtools intersection with MACS3 pileup files and GTF annotation:  
  [`6_MEU_CUTandRUN_bedtools_intersect_pileup_gtf_generic.sh`](./MEU_CUTandRUN/6_MEU_CUTandRUN_bedtools_intersect_pileup_gtf_generic.sh)  
- Combining TSV files (intersect output) and adding metadata:  
  [`7_MEU_CUTandRUN_combine_tsv_add_metadata_generic.sh`](./MEU_CUTandRUN/7_MEU_CUTandRUN_combine_tsv_add_metadata_generic.sh)  
- Chromosomal mean read depth calculation:  
  [`8_MEU_CUTandRUN_chrom_read_depth_calc_generic.pbs`](./MEU_CUTandRUN/8_MEU_CUTandRUN_chrom_read_depth_calc_generic.pbs)  
- Chromosomal mean read depth plotting:  
  [`9_MEU_CUTandRUN_chrom_read_depth_plots_generic.R`](./MEU_CUTandRUN/9_MEU_CUTandRUN_chrom_read_depth_plots_generic.R)  
- TSS plotting:  
  [`10_MEU_CUTandRUN_TSS_plots_generic.R`](./MEU_CUTandRUN/10_MEU_CUTandRUN_TSS_plots_generic.R)  

### SHA (Tasmanian devil)  
Scripts for SHA CUT&RUN data include:

- Quality control:  
  [`1_SHA_CUTandRUN_fastqc_multiqc_generic.pbs`](./SHA_CUTandRUN/1_SHA_CUTandRUN_fastqc_multiqc_generic.pbs)  
- Read trimming:  
  [`2_SHA_CUTandRUN_trimmomatic_generic.pbs`](./SHA_CUTandRUN/2_SHA_CUTandRUN_trimmomatic_generic.pbs)  
- Read alignment to reference genome:  
  [`3_SHA_CUTandRUN_parallel-subread-PE_generic.sh`](./SHA_CUTandRUN/3_SHA_CUTandRUN_parallel-subread-PE_generic.sh)  
- Peak calling:  
  [`4_SHA_CUTandRUN_macs3_generic.pbs`](./SHA_CUTandRUN/4_SHA_CUTandRUN_macs3_generic.pbs)  
- Setup for intersect directories:  
  [`5_SHA_CUTandRUN_setup_intersect_directories_generic.sh`](./SHA_CUTandRUN/5_SHA_CUTandRUN_setup_intersect_directories_generic.sh)  
- Bedtools intersection with MACS3 pileup files and GTF annotation:  
  [`6_SHA_CUTandRUN_bedtools_intersect_pileup_gtf_generic.sh`](./SHA_CUTandRUN/6_SHA_CUTandRUN_bedtools_intersect_pileup_gtf_generic.sh)  
- Combining TSV files (intersect output) and adding metadata:  
  [`7_SHA_CUTandRUN_combine_tsv_add_metadata_generic.sh`](./SHA_CUTandRUN/7_SHA_CUTandRUN_combine_tsv_add_metadata_generic.sh)  
- Join quantiles and intersect output:  
  [`8_SHA_CUTandRUN_join_quantiles_TSS_generic.R`](./SHA_CUTandRUN/8_SHA_CUTandRUN_join_quantiles_TSS_generic.R)  
- Chromosomal mean read depth calculation:  
  [`9_SHA_CUTandRUN_chrom_read_depth_calc_generic.pbs`](./SHA_CUTandRUN/9_SHA_CUTandRUN_chrom_read_depth_calc_generic.pbs)  
- Chromosomal mean read depth plotting:  
  [`10_SHA_CUTandRUN_chrom_read_depth_plots_generic.R`](./SHA_CUTandRUN/10_SHA_CUTandRUN_chrom_read_depth_plots_generic.R)  
- TSS plotting:  
  [`11_SHA_CUTandRUN_TSS_plots_generic.R`](./SHA_CUTandRUN/11_SHA_CUTandRUN_TSS_plots_generic.R)  
- TSS quantile plotting:  
  [`12_SHA_CUTandRUN_TSS_quantile_plots_generic.R`](./SHA_CUTandRUN/12_SHA_CUTandRUN_TSS_quantile_plots_generic.R)  

## RNA-seq (SHA Tasmanian devil)

SHA RNA-seq data processing and analysis scripts include:

- Quality control:  
  [`1_SHA_RNA-seq_fastqc_multiqc_generic.pbs`](./SHA_RNAseq/1_SHA_RNA-seq_fastqc_multiqc_generic.pbs)  
- Read trimming:  
  [`2_SHA_RNA-seq_trimmomatic_generic.sh`](./SHA_RNAseq/2_SHA_RNA-seq_trimmomatic_generic.sh)  
- Read alignment to reference genome:  
  [`3_SHA_RNA-seq_parallel-subread-PE_generic.sh`](./SHA_RNAseq/3_SHA_RNA-seq_parallel-subread-PE_generic.sh)  
- Feature counting:  
  [`4_SHA_RNA-seq_featureCounts_generic.pbs`](./SHA_RNAseq/4_SHA_RNA-seq_featureCounts_generic.pbs)  
- TPM calculation and quantile assignment:  
  [`5_SHA_RNA-seq_TPM_quantiles_generic.R`](./SHA_RNAseq/5_SHA_RNA-seq_TPM_quantiles_generic.R)  
