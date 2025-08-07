#!/bin/sh

module load samtools/1.21

out=hgsvc_chrYv1
mkdir -p $out

nonY=hgsvc_nonYv1
mkdir -p $nonY

# Waiting for 27 - 29 in samples.wi_manual.map
for i in 28 29 # $(seq 1 27)
do
  ln=`sed -n ${i}p samples.wi_manual.map`
  sample=`echo $ln | awk '{print $1}'`
  asm=`echo $ln | awk '{print $2}'` # entire reforiented asm

  ## Scaffold(s)
# :<<'END' # Silence for debugging nonY part
  awk -v sample=$sample '$1==sample {print $2}' T2TScaffolds.csv > scaffolds.tmp
  numYs=`cat scaffolds.tmp | wc -l | awk '{print $1}'`
  echo "$sample numYs: $numYs"
  
  out_fa=$out/${sample}_chrY.fa
  if [[ -s $out_fa.gz ]]; then
    rm $out_fa.*
  fi

  # HG03456 is special, has 2 Y scaffolds (XYY)
  echo "Extracting scaffold(s)"
  if [[ $numYs -eq 2 ]]; then
    scaffolds=`cat scaffolds.tmp | tr '\n' ' '`
    j=1
    for scaff in $scaffolds
    do
      faHeader=">${sample}_chrY_$j"
      echo -e "${sample}\t${sample}_chrY_$j\t$scaff" >> scaff_unplaced.map
      echo "$faHeader" >> $out_fa
      samtools faidx -@12 -n 60 $asm $scaff | awk 'NR>1' >> $out_fa
      j=$((j + 1))
    done
  elif [[ $numYs -eq 1 ]]; then
    scaff=`cat scaffolds.tmp`
    faHeader=">${sample}_chrY"
    echo -e "${sample}\t${sample}_chrY\t$scaff" >> scaff_unplaced.map
    echo "$faHeader" > $out_fa
    samtools faidx -@12 -n 60 $asm $scaff | awk 'NR>1' >> $out_fa
  fi
  rm scaffolds.tmp
# END

  ## Unplaced
  echo
  echo "Extracting unplaced sequences"
  fa=`echo $asm | sed 's/analysis\///g' | sed 's/refOriented.//g'` # original assembly.fasta from verkko

  if [[ ! -s $fa.fai ]]; then
    echo "Indexing $fa"
    samtools faidx -@12 $fa
  fi

  tigs=`ls tigs/$sample.tigs.no_noise.txt`
  unplaced=`awk '$3!="SCAFFOLD" {print $2}' $tigs | tr '\n' ' '`
# :<<'END' # Silence for debugging nonY part
  if [[ $unplaced != "" ]]; then
    echo "$sample unplaced: $unplaced"
    for tig in $unplaced
    do
      hasSeq=`grep -c $tig $fa.fai`
      if [[ $hasSeq == "1" ]]; then
        unptig=`echo $tig | sed 's/chrY_//g' | sed 's/haplotype1-//g' | sed 's/haplotype2-//g' | sed 's/mat-//g' | sed 's/pat-//g' | sed 's/unassigned-//g'`
        unptig="random$unptig"
        echo -e "${sample}\t${sample}_chrY_${unptig}\t$tig" >> scaff_unplaced.map
        echo ">${sample}_chrY_${unptig}" >> $out_fa
        echo "samtools faidx -@12 -n 60 $fa $tig"
        samtools faidx -@12 -n 60 $fa $tig | awk 'NR>1' >> $out_fa
      fi
    done
  fi

  bgzip -@12 -i $out_fa
  samtools faidx -@12 $out_fa.gz
  wc -l $out_fa.gz.fai
  echo
# END

  ## nonY
  echo "$sample : Extracting nonY sequences"
  tigs=`ls tigs/$sample.tigs.csv`
  awk '{print $2}' $tigs > $sample.Y.tmp
  echo "chrY sequences:"
  cat $sample.Y.tmp
  grep -v -f $sample.Y.tmp $asm.fai | cut -f1 > $sample.nonY.list
  rm $sample.Y.tmp
  needRenaming=`grep chrY $sample.nonY.list`
  if [[ $needRenaming != "" ]]; then
    echo "The following needs renaming:"
    echo $needRenaming
    exit -1
  fi

  echo
  out_fa=$nonY/${sample}_nonY.fa
  echo "samtools faidx -@12 -r $sample.nonY.list -n 60 $asm > $out_fa"
  samtools faidx -@12 -r $sample.nonY.list -n 60 $asm > $out_fa
  rm $sample.nonY.list
  bgzip -@12 -i $out_fa
  samtools faidx -@12 $out_fa.gz
  wc -l $out_fa.gz.fai
  echo
done

