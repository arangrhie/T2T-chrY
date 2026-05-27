#!/bin/sh

# Split ERRBASE, ERRSTRUCT and NGAP out

outdir=issues
mkdir -p $outdir chrY

beds=`ls /data/T2T-Y/globus/verkko-v2.2.1/annotations/seq_classes/2026-02_final-rev/patched/t2tv2/*.bed`

for bed in $beds
do
  out=`basename $bed`
  sample=`echo $out | awk -F '.' '{print $1}'`
  awk -F'\t' -v OFS='\t' '
    NR==FNR { color[$1]=$2; next }
    $(NF-2)=="ISSUE" {
      issue=$4
      print $1,$2,$3,issue,$5,".",$2,$3,color[issue]
    }' color_codes.map $bed > $outdir/$sample.issues.bed
  grep -v "random" $outdir/$sample.issues.bed > $outdir/$sample.issues.main.bed
  grep "random" $outdir/$sample.issues.bed > $outdir/$sample.issues.random.bed
  if [[ ! -s $outdir/$sample.issues.random.bed ]]; then
    rm $outdir/$sample.issues.random.bed
  fi
  awk -F'\t' -v OFS='\t' '{$1="chrY"; print $0}' $outdir/$sample.issues.main.bed > chrY/$sample.issues.main.bed
done

module load bedtools/2.31.1

for bed in $outdir/*.issues.main.bed
do
  # Merge adjacent issues of the same type (same color code) into a single issue
  # bedtools merge -i $bed > ${bed%.bed}.mrg.bed
  sample=`basename ${bed} | awk -F '.' '{print $1}'`
  
  # Make a non-overlapping version, by prioritizing NGAP > ERRSTRUCT > ERRBASE.
  awk '$4=="NGAP"' $outdir/$sample.issues.main.bed > $outdir/$sample.issues.NGAP.bed
  awk '$4=="ERRSTRUCT"' $outdir/$sample.issues.main.bed | \
    bedtools subtract -a - -b $outdir/$sample.issues.NGAP.bed \
    > $outdir/$sample.issues.ERRSTRUCT.bed
  
  cat $outdir/$sample.issues.NGAP.bed $outdir/$sample.issues.ERRSTRUCT.bed | \
    sort -k1,1V -k2,2n - > $outdir/$sample.issues.NGAP_ERRSTRUCT.bed
  awk '$4=="ERRBASE"' $outdir/$sample.issues.main.bed | \
    bedtools subtract -a - -b $outdir/$sample.issues.NGAP_ERRSTRUCT.bed \
    > $outdir/$sample.issues.ERRBASE.bed
  
  cat $outdir/$sample.issues.NGAP_ERRSTRUCT.bed $outdir/$sample.issues.ERRBASE.bed | \
    sort -k1,1V -k2,2n - > $outdir/$sample.issues.main.nonoverlap.bed

  # remove intermediate files
  rm $outdir/$sample.issues.NGAP.bed $outdir/$sample.issues.ERRSTRUCT.bed $outdir/$sample.issues.NGAP_ERRSTRUCT.bed $outdir/$sample.issues.ERRBASE.bed
done

