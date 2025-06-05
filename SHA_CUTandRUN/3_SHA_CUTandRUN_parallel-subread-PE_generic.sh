#!/bin/bash
#PBS -N Subread_mapping_paired_SHA_CUTandRUN
#PBS -l select=1:ncpus=8:mem=120gb
#PBS -l walltime=6:00:00
#PBS -j oe
#PBS -M email
#PBS -m ae

#Set the following variables
indir="/path/to/trimmomatic_output"
outdir="/path/to/subread_align_output"
index="/path/to/subread/index/basename"
R1identifier="_1P.fq.gz"
R2identifier="_2P.fq.gz"
seqtype="1"

#NOTE, make sure you see all of your raw R1 input files with the following ls command in your ${indir}
#ls *${R1identifier}

############################################################################################################################
module load subread/2.0.2 parallel

#Make sure to first build an index for your genome, cmd below
#subread-buildindex -o /path/to/subread/index/basename /path/to/genome.fasta

cd ${indir}

for f in $(ls *${R1identifier}); do echo ${f/${R1identifier}/}; done \
| parallel --jobs ${NCPUS} \
subread-align -t ${seqtype} -T ${NCPUS} --sortReadsByCoordinates -i ${index} \
-r {}${R1identifier} \
-R {}${R2identifier} \
-o ${outdir}/{}_subread.bam
