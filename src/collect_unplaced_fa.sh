#!/bin/sh

module load samtools

mkdir -p unplaced # unplaced

for sample in $(cat samples.list)
do
  fa=`awk -v sample=$sample '$1==sample {print $2}' samples.wi_manual_update.map`
  fa=`echo $fa | sed 's/analysis\///g' | sed 's/refOriented.//g'`

  ls $fa
  if [[ ! -s $fa.fai ]]; then
    samtools faidx -@12 $fa
  fi

  tigs=`ls tigs/$sample.tigs.no_noise.txt`
  unplaced=`awk '$3!="SCAFFOLD" {print $2}' $tigs | tr '\n' ' '`
  if [[ $unplaced == "" ]]; then
    echo "$sample : no unplaced"
    continue
  fi

  tig_fa=unplaced/${sample}_chrY.unplaced.fa
  if [[ -s $tig_fa ]]; then
    rm $tig_fa
  fi

  echo "$sample : $unplaced"
  for tig in $unplaced
  do
    hasSeq=`grep -c $tig $fa.fai`
    # echo $hasSeq $sample has $tig
    if [[ $hasSeq == "1" ]]; then
      unptig=`echo $tig | sed 's/chrY_//g' | sed 's/haplotype1-//g' | sed 's/haplotype2-//g' | sed 's/mat-//g' | sed 's/pat-//g' | sed 's/unassigned-//g'`
      unptig="random$unptig"
      echo -e "${sample}\tchrY_${unptig}\t$tig" >> unplaced.map
      echo ">${sample}_chrY_${unptig}" >> $tig_fa
      echo "samtools faidx -@12 -n 60 $fa $tig"
      samtools faidx -@12 -n 60 $fa $tig | awk 'NR>1' >> $tig_fa
    fi
  done
  samtools faidx -@12 $tig_fa
  wc -l $tig_fa.fai
done
