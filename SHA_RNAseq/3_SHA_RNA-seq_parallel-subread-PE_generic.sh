#!/bin/bash
#PBS -N Subread_mapping_paired_SHA_RNAseq
#PBS -l select=1:ncpus=8:mem=120gb
#PBS -l walltime=6:00:00
#PBS -j oe
#PBS -M email
#PBS -m ae

###Make sure to first build an index for your genome, cmd below
#subread-buildindex -o /path/to/index/basename /path/to/genome.fasta

###set the following variables
#Directory to folder containing your sequences
#Directory to output folder
#Path to your index
#File extension identifier for your R1 and R2 file, this will be 1P.fq.gz and 2P.fq.gz if they were trimmed with the parallel script in utils folder
#Your input sequence type, 1 for genomic DNA-seq and 0 for RNA-seq

indir="/path/to/trimmed_reads"
outdir="/path/to/subread_output"
index="/path/to/subread/index/basename"
R1identifier="1P.fq.gz"
R2identifier="2P.fq.gz"
seqtype="0"

###NOTE, make sure you see all of your raw R1 input files with the following ls command in your ${indir}
#ls *${R1identifier}

############################################################################################################################
module load subread/2.0.2 parallel

cd ${indir}

for f in $(ls *${R1identifier}); do echo ${f/${R1identifier}/}; done \
| parallel --jobs ${NCPUS} \
subread-align -t ${seqtype} -T ${NCPUS} --sortReadsByCoordinates -i ${index} \
-r {}${R1identifier} \
-R {}${R2identifier} \
-o ${outdir}/{}_subread.bam