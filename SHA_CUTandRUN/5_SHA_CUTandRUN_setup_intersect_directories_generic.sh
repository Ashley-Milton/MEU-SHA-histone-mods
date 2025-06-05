#Set up path to files
export topdir="/path/to/TSS_MACS3_intersect"
export GTF="/path/to/genome_annotation.gtf"

#cd into your working directory
cd ${topdir}

#make sub-directories
mkdir gene
mkdir pileup
##put the MACS3 .bdg files you want to analyse in /pileup
cp /path/to/macs3_output/*pileup* ${topdir}/pileup
mkdir log
mkdir workingdir
##make a sub-directory for each condition
cd ${topdir}/pileup
for i in $(ls *.bdg); do
    mkdir ${topdir}/workingdir/${i/.bdg}
done

#Split your annotation gtf by chromosome and adding 10kb to each side of the gene's TSS
##the first awk command filters out the unplaced scaffolds (gets rid of NW and keeps the NC). if you want all scaffolds instead, then replace content inside the bracket with next line.
##(awk '$3 == "gene"' ${GTF} | cut -f1 | sort | uniq)
##NOTE: need to also replace content of bracket in 6_SHA_CUTandRUN_bedtools_intersect_pileup_gtf_generic.sh if you want to use all scaffolds
for i in $(awk '$1 !~ /^NW/' ${GTF} | awk '$3 == "gene"' | cut -f1 | sort | uniq); do
    awk -v i="$i" '$1 == i' ${GTF} | \
    awk '$3 == "gene"' | \
    awk -F '\t'  '{print $1"\t"$4-1"\t"$5"\t"$9"\t"$7}' | \
    awk -F '\t' '{if ($5 == "+") print $1"\t"$2-10000"\t"$2+10001"\t"$4"\t"$5"\t"$2"\t"$3; else if ($5 == "-") print $1"\t"$3-10001"\t"$3+10000"\t"$4"\t"$5"\t"$2"\t"$3 }' | \
    awk -F '\t' 'BEGIN {OFS = "\t"} {if ($2 < 0) $2 = 0; print $0}' > ${topdir}/gene/gene_$i.bed
done

#Submit your intersect jobs, point the path to your copy of the script
cd ${topdir}/pileup
for i in $(ls *.bdg); do
    qsub -o ${topdir}/log -v topdir=${topdir},workingdir=${topdir}/workingdir/${i/.bdg},BDG=${topdir}/pileup/${i},GTF=${GTF} /path/to/6_SHA_CUTandRUN_bedtools_intersect_pileup_gtf_generic.sh
done
