#!/bin/sh

# Collect regions from *-pathced.bed file, give the sequence class and its coloring based on color_codes.map

outdir=originY
mkdir -p $outdir chrY

beds=`ls /data/T2T-Y/globus/verkko-v2.2.1/annotations/seq_classes/2026-02_final-rev/patched/t2tv2/*.bed`

for bed in $beds
do
  out=`basename $bed`
  sample=`echo $out | awk -F '.' '{print $1}'`
  awk -F'\t' -v OFS='\t' '
    NR==FNR { color[$1]=$2; next }
    $(NF-2)=="LABEL" {
      group=$NF
      label=$(NF-1)

      # Apply sequence class remapping rules: PAR -> PAR1 or PAR2, otherwise use the group as the class
      seqclass=group
      if (group ~ /^PAR/)       seqclass=label
      if (group == "OTH")       seqclass="OTHER"
      col=color[seqclass]
      print $1,$2,$3,seqclass,$5,".",$2,$3,col
    }' color_codes.map $bed | sort -k1,1V -k2,2n > $outdir/$sample.seqclass.bed
  grep -v "random" $outdir/$sample.seqclass.bed > $outdir/$sample.seqclass.main.bed
  grep "random" $outdir/$sample.seqclass.bed > $outdir/$sample.seqclass.random.bed
  if [[ ! -s $outdir/$sample.seqclass.random.bed ]]; then
    rm $outdir/$sample.seqclass.random.bed
  fi
  awk -F'\t' -v OFS='\t' '{$1="chrY"; print $0}' $outdir/$sample.seqclass.main.bed > chrY/$sample.seqclass.main.bed
done
