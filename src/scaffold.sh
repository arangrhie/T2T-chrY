#!/bin/sh

out=scaffolds
mkdir -p $out

module load samtools/1.21

prevSample=""
faHeader=""
scaffolds=chrY_scaffold.txt
samples_map=samples.wi_manual_update.map

LEN=`wc -l $scaffolds | awk '{print $1}'`

for i in $(seq 1 $LEN)
do
  ln=`sed -n ${i}p $scaffolds`
  sample=`echo $ln | awk '{print $1}'`
  seq=`echo $ln | awk '{print $2}'` # sequence
  dir=`echo $ln | awk '{print $3}'` # direction

  echo $sample $seq $dir

  if [[ $sample != $prevSample ]]; then
    if [[ $prevSample != "" ]]; then
      echo "Write $prevSample fa"
      echo "$faHeader" > $out/${prevSample}_chrY.fa
      cat tmp.seq | tr -d '\n' | sed -e "s/.\{60\}/&\n/g" >> $out/${prevSample}_chrY.fa # break lines every 60 bases
      echo "" >> $out/${prevSample}_chrY.fa # to add a newline at the end of the file
      echo "${prevSample}_chrY.fa"
    fi

    echo "Start a new fa sequence header and template"
    faHeader=">${sample}_chrY"
    fa=`awk -v sample=$sample '$1==sample {print $2}' $samples_map`
    rm tmp.seq
  else
    echo "Same sample: add gap"
    cat gap.200k.seq >> tmp.seq
  fi

  echo "Expand faSeq"
  if [[ $dir == "+" ]]; then
    set -x
    samtools faidx    -@12 $fa $seq | awk 'NR>1' >> tmp.seq
    set +x
  else
    set -x
    samtools faidx -i -@12 $fa $seq | awk 'NR>1' >> tmp.seq
    set +x
  fi

  prevSample=$sample;
done

echo "$faHeader" > $out/${prevSample}_chrY.fa
cat tmp.seq | tr -d '\n' | sed -e "s/.\{60\}/&\n/g" >> $out/${prevSample}_chrY.fa
echo "" >> $out/${prevSample}_chrY.fa # to add a newline at the end of the file
echo "${prevSample}_chrY.fa"
rm tmp.seq
