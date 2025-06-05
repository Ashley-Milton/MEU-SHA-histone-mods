#!/bin/bash
#PBS -N combine_tsv_SHA
#PBS -l ncpus=1,walltime=4:00:00,mem=100GB
#PBS -j oe
#PBS -M email
#PBS -m ae

cd /path/to/TSS_MACS3_intersect/workingdir

if [ "$Sex" = "Female" ]; then
    export Sex_Histone="F_TD_${Histone}"
    export Sex_HistoneDir=$(ls -d F_TD_${Histone}_*_treat_pileup | grep -E "F_TD_${Histone}_.*_treat_pileup")
elif [ "$Sex" = "Male" ]; then
    export Sex_Histone="M_TD_${Histone}"
    export Sex_HistoneDir=$(ls -d M_TD_${Histone}_*_treat_pileup | grep -E "M_TD_${Histone}_.*_treat_pileup")
fi

# Debugging output
echo "Sex: ${Sex}"
echo "Histone: ${Histone}"
echo "Sex_Histone: ${Sex_Histone}"
echo "Sex_HistoneDir: ${Sex_HistoneDir}"

# Check if Sex_HistoneDir is set and is a valid directory
if [ -z "$Sex_HistoneDir" ] || [ ! -d "$Sex_HistoneDir" ]; then
    echo "Error: Sex_HistoneDir is not set correctly or does not exist."
    exit 1
fi

echo ${Sex} > ${Sex_Histone}.running
echo ${Histone} >> ${Sex_Histone}.running
echo ${Sex_Histone} >> ${Sex_Histone}.running
echo ${Sex_HistoneDir} >> ${Sex_Histone}.running

cd ${Sex_HistoneDir}

cd combine

mkdir workingdir

awk -v OFS='\t' '{
    split($1, a, ",");
     if (a[1] == "chrX") 
        print $0, "X";
    else if (a[1] == "chrY") 
        print $0, "Y";
    else 
        print $0, "A";
}' *.tsv | awk -F'\t' -v OFS='\t' -v histone="$Histone" -v sex="$Sex" '{print $0"\t"histone"\t"sex}' | \
sed $'1iposition\tgene\tstrand\tgene_start\tgene_end\trelative_pos\theight\tc\tmod\tsex' > workingdir/${Sex_Histone}_info_headers.tsv

cut -f6,7,8,9,10 -d$'\t' workingdir/${Sex_Histone}_info_headers.tsv > workingdir/${Sex_Histone}_info_simplified.tsv

mv workingdir/*.tsv .

rm -rf workingdir

cd /path/to/TSS_MACS3_intersect/workingdir
mv ${Sex_Histone}.running ${Sex_Histone}.done
