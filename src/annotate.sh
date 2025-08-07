#!/bin/sh

if [[ $# -lt 4 ]]; then
	echo "Usage: $0 <annotate_regions.txt> <sample> <asm_fa> <out>"
  echo "  annotate_regions.txt region to annotate:"
  echo "          <name of the target region> <tab> <chr:start-end> <tab> <color in final annotation>"
  echo "  sample  name of the sample to be annotated"
  echo "  asm_fa  fasta file to be annotated"
  echo "  out     output directory"
	exit 1
fi

target_regions=$1
sample=$2
asm_fa=$3
out=$4/$sample
# Fix chr to chrY so I can display them all on IGV under chrY
chr=chrY

len=7000 # set default to 7 kb?
if [[ -z $SLURM_CPUS_PER_TASK ]]; then
  SLURM_CPUS_PER_TASK=12
fi

module load samtools/1.21
module load minimap2/2.28
module load rustybam/0.1.33
module load bedtools/2.31.1
module load seqtk/1.4

mkdir -p $out

LEN=`wc -l $target_regions | awk '{print $1}'`

if [[ ! -s $out/${sample}_to_sat.paf ]]; then
  for i in $(seq 1 $LEN)
  do
    line=`sed -n ${i}p $target_regions`
    qry_name=`echo $line | awk '{print $1}'`
    region=`echo $line | awk '{print $2}'`
    if [[ $region == "." ]]; then
      continue
    fi

    echo "Mapping $asm_fa to $qry_name"
    set -x
    minimap2 \
      -x asm20 \
      --eqx \
      --MD \
      -t ${SLURM_CPUS_PER_TASK} \
      -c \
      $qry_name.fa \
      $asm_fa \
      >> $out/${sample}_to_sat.paf
    set +x
  done
fi

if [[ ! -s $out/sat_to_${sample}.paf ]]; then
  set -x
  rustybam invert $out/${sample}_to_sat.paf > $out/sat_to_${sample}.paf
  set +x
fi



#awk -v len=$len '{ if( $11 >= len ) print }' $out/sat_to_${sample}.paf > $out/sat_to_${sample}.filter.paf
cat $out/sat_to_${sample}.paf |\
  awk -v OFS='\t' '{if ($10<$11) {idy=100*$10/$11;} else {idy=100*$11/$10;} print $6,$8,$9,$1,idy,$5}' |\
  sort -k1,1V -k2,2n - > $out/sat_to_${sample}.bed
#echo "Final alignment file: $out/sat_to_${sample}.filter.bed, filtered at > $len bp"



if [[ ! -s $asm_fa.fai ]]; then
  echo "Indexing $asm_fa"
  samtools faidx -@$SLURM_CPUS_PER_TASK $asm_fa
fi
fai=`ls $asm_fa.fai`

seqtk telo -d 50000 $asm_fa > $out/${sample}.telo.bed
seqtk gap -l 1 $asm_fa > $out/${sample}.gap.bed

## Chromosomes of interest
chrs=`cat $out/sat_to_${sample}.bed | cut -f1 | sort -u | sort -k1,1V | tr '\n' ' '`
echo "chrs: $chrs"

for chr_hap in $chrs
do
  echo "chr_hap: $chr_hap"
  hap=`echo $chr_hap | sed "s/$sample\_//g"` # rename output to prevent having sample name twice for the final chrYs
  OUT="${out}/${sample}_${hap}.bed"
  TMP="${out}/${sample}_${hap}.tmp.bed"
  if [[ -s $OUT ]]; then
    rm $OUT
  fi
  
  # For each satellite
  for i in $(seq 1 $LEN)
  do
    line=`sed -n ${i}p $target_regions`
    sat=`echo $line | awk '{print $1}'`
    col=`echo $line | awk '{print $4}'`
    echo "Annotating $sat ..."

    if [[ $sat == "TEL" ]]; then
      awk -v chr_hap=$chr_hap '$1==chr_hap' ${out}/${sample}.telo.bed | \
        awk -v chr=$chr -v col=$col -v OFS="\t" '{print chr, $2, $3, "TEL", 100, ".", $2, $3, col}' >> $OUT
    
    elif [[ $sat == "GAP" ]]; then
      awk -v chr_hap=$chr_hap '$1==chr_hap' ${out}/${sample}.gap.bed | \
        awk -v chr=$chr -v col=$col  -v OFS="\t" '{print chr, $2, $3, "GAP", 100, ".", $2, $3, col}' >> $OUT
    elif [[ $sat == "SEQ" ]]; then
      # Fill in the rest as SEQ bg
      awk -v chr_hap=$chr_hap '$1==chr_hap' $fai | \
        awk -v chr=$chr -v col=$col -v OFS="\t" '{print chr, 0, $2, "SEQ", 100, ".", 0, $2, col}' > $TMP
    else
      awk -v chr_hap=$chr_hap -v sat=$sat '$1==chr_hap && $4==sat' ${out}/sat_to_${sample}.bed | \
        bedtools merge -s -d 500 -c 4,5,6 -o distinct,median,distinct -i - | \
        awk -v len=$len '$3-$2>len' |\
        awk -v chr=$chr -v col=$col -v OFS="\t" '{print chr, $2, $3, $4, $5, $6, $2, $3, col}' >> $OUT
    fi
  done

  bedtools subtract -a $TMP -b $OUT > ${OUT/.bed/_bg.bed}
  cat ${OUT/.bed/_bg.bed} $OUT | sort -k1,1V -k2,2n > $TMP
  mv $TMP $OUT
  rm ${OUT/.bed/_bg.bed}
done

# clean up
rm ${out}/${sample}.telo.bed ${out}/${sample}.gap.bed 
rm ${out}/sat_to_${sample}.*

