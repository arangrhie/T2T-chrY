#!/bin/sh

# This script reads a bed file and a fasta file, and annotates the fasta file according to the bed file by mapping the fasta file to chm13.

if [[ $# -lt 3 ]]; then
	echo "Usage: $0 <annotate_regions.txt> <asm_fa> <out>"
  echo "  annotate_regions.txt region to annotate:"
  echo "          <name of the target region> <tab> <chr:start-end> <tab> <color in final annotation>"
  echo "  asm_fa  fasta file to be annotated"
  echo "  out     output directory"
	exit 1
fi

target_regions=$1
asm_fa=$2
out=$3
asm_name=`basename $asm_fa`
asm_name=`echo $asm_name | sed 's/.gz$//g' | sed 's/\.fa$//g' | sed 's/\.fasta$//g'`
len=7000 # set default to 7 kb?
if [[ -z $SLURM_CPUS_PER_TASK ]]; then
  SLURM_CPUS_PER_TASK=12
fi
module load samtools/1.21
module load minimap2/2.28
module load rustybam/0.1.33
module load bedtools/2.31.1
module load seqtk

mkdir -p $out

LEN=`wc -l $target_regions | awk '{print $1}'`

if [[ ! -s $out/${asm_name}_to_sat.paf ]]; then
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
      >> $out/${asm_name}_to_sat.paf
    set +x
  done
fi

if [[ ! -s $out/sat_to_${asm_name}.paf ]]; then
  set -x
  rustybam invert $out/${asm_name}_to_sat.paf > $out/sat_to_${asm_name}.paf
  set +x
fi



#awk -v len=$len '{ if( $11 >= len ) print }' $out/sat_to_${asm_name}.paf > $out/sat_to_${asm_name}.filter.paf
cat $out/sat_to_${asm_name}.paf |\
  awk -v OFS='\t' '{if ($10<$11) {idy=100*$10/$11;} else {idy=100*$11/$10;} print $6,$8,$9,$1,idy,$5}' |\
  sort -k1,1V -k2,2n - > $out/sat_to_${asm_name}.bed
#echo "Final alignment file: $out/sat_to_${asm_name}.filter.bed, filtered at > $len bp"


# Fix chr to chrY so I can display them all on IGV under chrY
chr=chrY

if [[ ! -s $asm_fa.fai ]]; then
  echo "Indexing $asm_fa"
  samtools faidx -@$SLURM_CPUS_PER_TASK $asm_fa
fi
fai=`ls $asm_fa.fai`

seqtk telo -d 50000 $asm_fa > $out/${asm_name}.telo.bed
seqtk gap -l 1 $asm_fa > $out/${asm_name}.gap.bed

## Chromosomes of interest
chrs=`cat $out/sat_to_${asm_name}.bed | cut -f1 | sort -u | sort -k1,1V | tr '\n' ' '`
echo "chrs: $chrs"

for chr_hap in $chrs
do
  echo "chr_hap: $chr_hap"
  OUT="${out}/${chr_hap}.bed"
  TMP="${out}/${chr_hap}.tmp.bed"
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
      awk -v chr_hap=$chr_hap '$1==chr_hap' ${out}/${asm_name}.telo.bed | \
        awk -v chr=$chr -v col=$col -v OFS="\t" '{print chr, $2, $3, "TEL", 100, ".", $2, $3, col}' >> $OUT
    
    elif [[ $sat == "GAP" ]]; then
      awk -v chr_hap=$chr_hap '$1==chr_hap' ${out}/${asm_name}.gap.bed | \
        awk -v chr=$chr -v col=$col  -v OFS="\t" '{print chr, $2, $3, "GAP", 100, ".", $2, $3, col}' >> $OUT
    elif [[ $sat == "SEQ" ]]; then
      # Fill in the rest as SEQ bg
      awk -v chr_hap=$chr_hap '$1==chr_hap' $fai | \
        awk -v chr=$chr -v col=$col -v OFS="\t" '{print chr, 0, $2, "SEQ", 100, ".", 0, $2, col}' > $TMP
    else
      awk -v chr_hap=$chr_hap -v sat=$sat '$1==chr_hap && $4==sat' ${out}/sat_to_${asm_name}.bed | \
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
rm ${out}/${asm_name}.telo.bed ${out}/${asm_name}.gap.bed 
rm ${out}/sat_to_${asm_name}.*

