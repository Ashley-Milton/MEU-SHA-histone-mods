#!/bin/bash
#PBS -N MACS3_intersect
#PBS -l select=1:ncpus=1:mem=120gb
#PBS -l walltime=6:00:00
#PBS -j oe
#PBS -M email
#PBS -m ae

module load bedtools2/2.30.0

export existing_gene="${topdir}/gene"

cd ${workingdir}
mkdir combine
export BDG_BASE=$(basename ${BDG/.bdg/})

#Replace content of bracket if you want to use all scaffolds, refer to 5_MEU_CUTandRUN_setup_intersect_directories_generic.sh for details
for i in $(awk '$1 !~ /^NW/' ${GTF} | awk '$3 == "gene"' | cut -f1 | sort | uniq); do
    bedtools intersect -a ${existing_gene}/gene_$i.bed -b ${BDG} -wa -wb | \
    awk -F '\t' '{if ($2 > $9 && $3 > $10) print $1"\t"$2"\t"$10"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11; else if ($2 < $9 && $3 < $10) print $1"\t"$9"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11; else if ($2 > $9 && $3 < $10) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11; else print $1"\t"$9"\t"$10"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11}' | \
    awk -F '\t' 'BEGIN {OFS = "\t"} {for (i=$2; i<$3; i++) print $1, i, i+1, $4, $5, $6, $7, $8}' | \
    awk -F '\t' '{if ($5 == "+") print $1","$2","$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$2-$6"\t"$8; else if ($5 == "-") print $1","$2","$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$7-$3"\t"$8}' > ${workingdir}/combine/combined_$i.tsv
done

#This will tar and compress the output files
tar cf - combine | pigz -p 4 > combine.tar.gz

#This will remove the combine directory, you can comment this out if you want to keep the uncompressed folder for immediate downstream processing
#rm -rf combine

