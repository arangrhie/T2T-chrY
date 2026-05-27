#!/bin/sh

if [[ $# -ne 1 ]]; then
	echo "Usage: $0 sample"
  echo "  sample: sample name, e.g. HG02486"
  echo "Pulls out input files from the following paths:"
	echo "  refcovered.sort.bed  bed file of T2T-Y covered regions by prior-Y sequences"
	echo "  seqclass.main.bed    bed file of sequence classes in the T2T-Y sequences"
	exit 1
fi

sample=$1
refcovered=`realpath /data/T2T-Y/y_comparison/GQC_AC/*/$sample/${sample}_oldY_vs_newY/${sample}_oldY_vs_${sample}_newY.refcovered.sort.bed`
seqclass=`realpath ../seqclass_ideogram/originY/$sample.seqclass.main.bed`

if ! [[ -s $refcovered && -s $seqclass ]]; then
  # echo "Error: input files not found or empty"
  # echo "  refcovered: $refcovered"
  # echo "  seqclass: $seqclass"
  exit 1
fi

out_dir=gaps_in_prior_wi_0
mkdir -p $out_dir

# Input bed file looks as following:
# HG02486_chrY    80      174060  HG02486#1#JAGYVM020000041.1.87.174063
# HG02486_chrY    174283  1088221 HG02486#1#JAGYVM020000041.1.174288.1088249
# HG02486_chrY    1089546 1192032 HG02486#1#JAGYVM020000041.1.1089590.1192079
# HG02486_chrY    1192097 1331830 HG02486#1#JAGYVM020000041.1.1192154.1331899
# HG02486_chrY    1332094 10467495        HG02486#1#JAGYVM020000041.1.1332138.10467588
# HG02486_chrY    10515861        22960221        HG02486#1#JAGYVM020000029.1.1.12444492
# HG02486_chrY    22960221        26447330        HG02486#1#JAGYVM020000029.1.12468031.15955135
# HG02486_chrY    26478407        32326428        HG02486#1#JAGYVM020000029.1.15955136.21803100
# HG02486_chrY    32329972        41653463        HG02486#1#JAGYVM020000029.1.21803101.31126591
# HG02486_chrY    41655883        47094410        HG02486#1#JAGYVM020000029.1.31126592.36565109

# Output file looks as following:
# CEN-DYZ3        1
# HET     2

# 1. Merge covered regions based on the 4th column (contig name); until first awk
# 2. Find where each contig start and end position belong in the newY sequence class; from bedtools intersect
# module load bedtools

cat $refcovered | grep -v random | \
  awk -v OFS='\t' '{
    qend=gensub(/^.*\./, "", 1, $NF);
    qstart=gensub(/^.*\./, "", 1, gensub(/(\.[^.]*){1}$/, "", 1, $NF));
    $4=gensub(/(\.[^.]*){2}$/, "", 1, $NF);
    print $0, qstart, qend}' | \
  awk -v OFS='\t' '{ if (NR == 1) { prev=$4; start=$2; end=$3; qstart=$5; qend=$6 }
    else { if ($4 == prev) { end=$3; qend=$6 }
	       else { print $1, start, end, prev, qstart, qend;
		          prev=$4; start=$2; end=$3; qstart=$5; qend=$6 }
	} } END { print $1, start, end, prev, qstart, qend }' | \
  bedtools intersect -loj -a $seqclass -b - -sortout | sort -k2,2n -k13,13V | \
    awk -v OFS='\t' -v sample=$sample '{
      if (FNR==NR) { gaps[$1]=0; next }
      if (FNR==1) {qry=$(NF-2); prevclass=$4;}
      if ($10 == ".") { gaps[$4]+=1; }
      else if (NR>1 && qry!="." && $(NF-2) != qry) { gaps[prevclass]+=1;
                                 if ($(NF-2)!="." && prevclass != $4) {gaps[$4]+=1};
                                 prevclass=$4; qry=$(NF-2);
                                }
      else { prevclass=$4; qry=$(NF-2); } }
    END { for (class in gaps) { print sample, class, gaps[class] } }' seqclass.list - \
> $out_dir/${sample}.gaps_in_prior.txt
