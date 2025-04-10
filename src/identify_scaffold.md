# Identify Y Scaffold

* Y Scaffold: Backbone Y chromosome assembly
* Y Unplaced: Sequences unable to place within the scaffold, due to ambiguity in the orientation
* refOriented: Initially identified in the verkko v2.2.1 freeze using mashmap

## Dependency
* aws cli/2.15.26
* samtools/1.21
* bedtools/2.31.1
* seqtk/1.4
* minimap2/2.28
* rustybam/0.1.33

## Idea
The Y chromosome is highly conserved for its structure except for the PAR, non ampliconic sequences, non heterochromatic sequences.

1. Mask the following regions in CHM13v2Y:
* PAR and HET from `[chm13v2.0_chrXY_sequence_class_v1.bed](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_chrXY_sequence_class_v1.bed)`
* IR, P1-5; entire `[chm13v2.0Y_inverted_repeats_v2.bed](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0Y_inverted_repeats_v2.bed)`
* Amplicons: `[chm13v2.0Y_amplicons_v2.bed](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0Y_amplicons_v2.bed)` - later used in the analysis. Could've included earlier.

2. Map the refOriented versions to the `masked_Y`

3. Find the proper orientation and order of the sequences belonging to Y

4. Check for presence of telomeres and gaps

5. Map the  refOriented versions missing telomeres to the non-masked Y to further identify missing PARs and telomeres



## 1. `chrY_masked`
Download files we need and generate a masked Y

```sh
aws s3 cp --no-sign-request s3://human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_chrXY_sequence_class_v1.bed .
aws s3 cp --no-sign-request s3://human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0Y_inverted_repeats_v2.bed
aws s3 cp --no-sign-request s3://human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0Y_amplicons_v2.bed

cat chm13v2.0Y_inverted_repeats_v2.bed | awk '{print $1"\t"$2"\t"$3}' > mask.bed
awk '$1="chrY" && ( $4~/^PAR/ || $4=="HET" ) {print $1"\t"$2"\t"$3}'  chm13v2.0_chrXY_sequence_class_v1.bed >> mask.bed
cat mask.bed | sort -k2,2n | bedtools merge -i - > mask.mrg.bed

seqtk seq -M mask.mrg.bed -n N chrY.fa > chrY_masked.fa
samtools faidx -@12 chrY_masked.fa
```

## 2. Alignments
Below is a snippet of the alignment. Because the qry sequence is much smaller, the assembly is aligned against the qry sequence then inverted for its coordinate. This is more sensitive for finding mutiple alignments with minimap2, which is designed to find one best hit of the querey sequence to the target.

```sh
ref=$1 # assembly fasta
ref_name=`basename $ref`
ref_name=`echo $ref_name | sed 's/.gz$//g' | sed 's/\.fa$//g' | sed 's/\.fasta$//g'`

qry=$2 # region of interest fasta
qry_name=`basename $qry`
qry_name=`echo $qry_name | sed 's/.gz$//g' | sed 's/\.fa$//g' | sed 's/\.fasta$//g'`

len=$3

minimap2 \
    -x asm10 \
    --eqx \
    --MD \
    -t ${SLURM_CPUS_PER_TASK} \
    -c \
    $qry \
    $ref \
    > ${ref_name}_to_${qry_name}.paf

rustybam invert ${ref_name}_to_${qry_name}.paf > ${qry_name}_to_${ref_name}.paf

awk -v len=$len '{ if( $11 >= len ) print }' ${qry_name}_to_${ref_name}.paf > ${qry_name}_to_${ref_name}.filter.paf

cat ${qry_name}_to_${ref_name}.filter.paf |\
  awk -v OFS='\t' -v qry_name=${qry_name} '{if ($10<$11) {idy=100*$10/$11;} else {idy=100*$11/$10;} print $6,$8,$9,qry_name,idy,$5}' |\
  sort -k1,1V -k2,2n - > ${qry_name}_to_${ref_name}.filter.bed
```

This script is available as `align.sh`, and was run through `align_job.sh` for a list of assemblies in slurm job arrays.

```sh
mkdir -p verkko_wi_manual_update # link the assembly.fasta files as $sample.fa
LEN=`wc -l samples.wi_manual_update.map | awk '{print $1}'`
for i in $(seq 1 $LEN)
do
  ln=`sed -n ${i}p samples.wi_manual_update.map`
  sample=`echo $ln | awk '{print $1}'`
  fa=`echo $ln | awk '{print $2}'`
  ln -sf $fa verkko_wi_manual_update/$sample.fa
done


# refOriented_fa_wi_Y.list : path to 105 assemblies with Ys
# split by 5 lines, starting from 10 to make the array jobs easier to track
split -d -a2 -l 5 --numeric-suffixes=10 refOriented_fa_wi_Y.list refOriented_fa_wi_Y.

qry=chrY_masked
list_prefix=refOriented_fa_wi_Y
len=10000 # ignore mappings < 10 kb, to exclude small hits to TEs

mkdir $qry && cd $qry

# Submit as a job array with 24 cores and 64g of memory
sh ~/codes/_submit_norm.sh 24 64g aln.$qry ../align_job.sh "../$list_prefix $qry.fa 10000" --array=10-30
```
Later, HG003 was added to the list of our samples and mapped through the same script. The initial HG02015 assembly had mis-scaffolded the Y (inverted the entire p to heterochromatic region). The `verkko-hi-c` mode was used instead after confirming the assembly was properly orienting the Y.

The same `align_job.sh` submission script was used for aligning to the unmasked Y.
```sh
qry=chrY
list_prefix=missing_PAR
len=10000
```

## 3. Filter alignments, identify Y scaffold
After manual inspection, the following regions in `mask.cov.bed` were excluded to collect the total uniquely mappable region.

```sh
cat mask.cov.bed
chrY    0       2458320
chrY    5877201 10389946
chrY    14891106        14926304
chrY    14929718        14964927
chrY    16781335        16811393
chrY    17066094        17332346
chrY    18362331        19776671
chrY    21881179        21945000
chrY    22270252        27124000
chrY    27243773        27254760
chrY    27449931        62072743
chrY    62122809        62460029
```
These include the following region not included in the original `mask.mrg.bed`:
* TSPY  AMPL, extend to IR3: ($3>5877201 && $4<10389946)
* P7    AMPL: ($3>16781335 && $4<16811393)
* P6    AMPL: ($3>17066093 && $4<17332347)
* P5    AMPL: ($3>18362331 && $4<19776672)
* X-DEG bSat: ($3>21881179 && $4<21945000)
* P1-3  AMPL: ($3>22270252 && $4<27124000)
* COMP MER5A1 ($3>27243773 && $4<27254760)

The final uniquely mappable regions add up to: `13819319` bp
```sh
# CHM13v2Y length: 62460029
mappable=`cat mask.cov.bed | awk 'BEGIN {sum=0;} { sum+=($3-$2) } END {print 62460029-sum}'`
echo $mappable
# 13819319
```

Alignments were chosen based on the following criteria:
* Exclude any mappings falling in the AMBL, bSat, COMP MER5A1 region
* Exclude any sequences larger than 90 Mbs in length (later inforced after breaking down the X-fused Y misassemblies)
* Alignment blocks are merged for its min and max coordinate along with its direction (`+` and `-`)
* Each sequence block was reported and evaluated if
  * it was originally identified as part of chrY
  * it was originally not identified as any chromosome (chr prefix) and the entire Y coverage was larger than 70%, or
  * the block mappable coverage is larger than 10% and the sequence is shorter than 30 Mbp, or
  * it maps on the PAR2 and the length is shorter than 30 Mbp and the entire Y coverage is larger than 0.8%

```sh
for sample in $(cat samples.list)
do
  cat chrY_masked/chrY_masked_to_$sample.*filter.paf | cut -f1-12 | sort -k3,3n | awk '$7<90000000' | \
    awk '! (($3>5877201 && $4<10389946) \
        || ($3>16781335 && $4<16811393) \
        || ($3>17066093 && $4<17332347) \
        || ($3>18362331 && $4<19776672) \
        || ($3>21881179 && $4<21945000) \
        || ($3>22270252 && $4<27124000) \
        || ($3>27243773 && $4<27254760) )' | \
    awk -v sample=$sample -v mappable=$mappable '\
     {  if ($5=="+") {
          seqBlock[$6"\t+"] = seqBlock[$6"\t+"] + $(NF-2); \
          seqLen[$6"\t+"] = $7
          if ( seqCoordMin[$6"\t+"] == 0 || $3 < seqCoordMin[$6"\t+"] ) seqCoordMin[$6"\t+"] = $3; \
          if ( seqCoordMax[$6"\t+"] < $4 ) seqCoordMax[$6"\t+"] = $4; \
        } \
        else { 
          seqBlock[$6"\t-"] = seqBlock[$6"\t-"] + $(NF-2); \
          seqLen[$6"\t-"] = $7
          if ( seqCoordMin[$6"\t-"] == 0 || $3 < seqCoordMin[$6"\t-"] ) seqCoordMin[$6"\t-"] = $3; \
          if ( seqCoordMax[$6"\t-"] < $4 ) seqCoordMax[$6"\t-"] = $4; \
        }  \
     } \
     END \
    { for (key in seqBlock) \
      { print sample "\t" key "\t" seqCoordMin[key] "\t" seqCoordMax[key] "\t" seqBlock[key] "\t" 100 * seqBlock[key] / mappable "\t" seqLen[key] "\t" 100 * seqLen[key] / $2 } }' | sort -k4,4n | \
      awk '{ if ($2 ~/^chrY/ ) { print $0 } \
        else { \
          if ( $2 !~/^chr/ && ( $(NF-2) > 70 || ( $NF>10 && $(NF-1)<30000000 ) || $4<62460029 && $5>62072743 && $(NF-1)<30000000 && $NF>0.8 )) { print $0 } \
        } \
      }' >> chrY_coverage.txt
done
```

At the end this generates `chrY_coverage.txt`, with  the following example:
```
#Sample	Seq	Direction	AlignBlockStart	AlignBlockEnd	BlockSize	MappableCov(%)	SeqSize	YSizeRatio(%)
HG003	chrY_haplotype1-0000002	+	2458320	62122809	14135503	102.288	59253811	94.8668
HG00321	chrY_haplotype2-0000065	+	2458320	27449931	13957698	101.001	37754623	60.4461
HG00321	chrY_haplotype2-0000041	+	62072743	62122808	50053	0.362196	7793185	12.4771
```
This way, we can identify that `chrY_haplotype2-0000065+` and `chrY_haplotype2-0000041+` can be placed in order as the larger Y scaffold.

Each sequence is then re-evaluated to filter out nested or duplicated sequence alignments.
```sh
cat chrY_coverage.txt | awk -v OFS="\t" 'BEGIN {genome=""; blockEnd=0;} \
  {
    if (genome==$1) {
      if ( blockEnd - $4 > 115000 || blockEnd > $5) {
        # print $0, "NESTED - ignore";
        next;
      } else if (seq==$2) {
        if ((blockEnd-blockStart) > $5-$4) {
          # print $0, "DupSeq - ignore";
          next;
        }
      } else {
        print ln, "SEQ"
      }
    } else {
      if (genome != "") print ln, "COMPLETE"
    }
    genome=$1
    seq=$2
    blockStart=$4
    blockEnd=$5
    ln=$0
  } END {print ln, "COMPLETE"}' >> chrY_scaffold.txt
```

The PARs were later identified and manually inserted from the non-masked chrY alignment (described later):

```sh
## HG01106	chrX_mat-0000007	+	PAR1
## HG01252	pat-0000707	+	PAR2
## HG02145	haplotype1-0000075	+	PAR1
```

At the end, we are adding a 200 kb gap and scaffold to `scaffolds/${sample}_chrY.fa`.

`scaffold.sh` takes `chrY_scaffold.txt` and `samples.wi_manual_update.map` as input and generates scaffolds under `scaffolds` folder. `gap.200k.seq` is solely made of 200k Ns, concatenated for inserting gaps.

```sh
./scaffold.sh
```

## 4. Check for presence of telomeres and gaps

```sh
for fa in $(ls scaffolds/*.fa)
do
  out=${fa/.fa/}
  seqtk telo -d 50000 $fa > $out.telo.bed
  seqtk gap -l1 $fa > $out.gap.bed
  cat $out.telo.bed | awk '{if ($2==0) t="p"; if ($3==$4) t=t"q";} END {print $1" "t}' >> scaffolds_telo.txt
done

```
Samples with no p-tels were collected in samples.PAR1missing.list, and likewise for q-tel missing in samples.PAR2missing.list.


## 5. Map the refOriented versions missing telomeres to the non-masked Y to further identify missing PARs and telomeres


### PAR1
```sh
PAR1=2458320
for sample in $(cat samples.PAR1missing.list)
do
  echo $sample
  cat chrY/chrY_to_${sample}.refOriented.filter.paf | cut -f1-12 | sort -k3,3n |\
    awk -v PAR1=$PAR1 '$3<PAR1 && $7 <20000000'
done

```

```
HG01106
chrY    62460029        5654    44233   +       chrX_mat-0000007        974310  12444   33208   19734   39609   60
chrY    62460029        63058   97156   +       chrX_mat-0000007        974310  35688   68330   29484   37256   60
chrY    62460029        106236  199623  +       chrX_mat-0000007        974310  76011   166629  88739   95266   60
chrY    62460029        200361  593915  +       chrX_mat-0000007        974310  165395  546167  373969  400357  60
chrY    62460029        593652  712945  +       chrX_mat-0000007        974310  546148  665061  117787  120419  60
chrY    62460029        713033  1016939 +       chrX_mat-0000007        974310  665154  974279  298812  314219  60
HG01358
HG01952
HG02145
chrY    62460029        331852  683802  +       haplotype1-0000075      724167  442     350838  341793  360553  60
chrY    62460029        674079  692013  -       haplotype1-0000074      18258   0       18145   17859   18220   60
chrY    62460029        683018  710919  +       unassigned-0000544      28726   0       28726   27834   28793   60
chrY    62460029        698057  712945  +       haplotype1-0000075      724167  396160  411037  14854   14911   60
chrY    62460029        713033  1021026 +       haplotype1-0000075      724167  411130  724167  301952  319078  60
chrY    62460029        1068185 1095490 +       chrY_haplotype1-0000002 8466197 0       27088   26948   27445   60
chrY    62460029        1095518 1443160 +       chrY_haplotype1-0000002 8466197 26916   378708  337363  362071  60
chrY    62460029        1931808 10111439        +       chrY_haplotype1-0000002 8466197 478736  8466165 7982615 8184445 60
HG03209
NA18608
NA20850
```

* For `HG01106`: `chrX_mat-0000007` is the only remaining contig. Another 150 Mbp sequence with chrX labeled exists, mapping to the same region.
* For `HG02145`: Choosing `haplotype1-0000075`. The other two are covering the same locus.


### PAR2

```sh
PAR2=62122810
for sample in $(cat samples.PAR2missing.list)
do
  echo $sample
  cat chrY/chrY_to_${sample}.refOriented.filter.paf | cut -f1-12 | sort -k3,3n |\
    awk -v PAR2=$PAR2 '$4>PAR2 && $7 <20000000'
done
```
```
HG01252
chrY    62460029        62302824        62453630        +       pat-0000707     160587  1       152293  150765  152333  60
chrY    62460029        62305622        62453630        +       pat-0000708     157790  1       149496  147968  149535  60
HG01433
HG02145
HG04157
chrY    62460029        61903453        62460019        +       chrY_haplotype2-0000084 15958177        13784906        14358011        556175  573496  34
NA21093
```
For `HG01252`: Choose `pat-0000707` over `pat-0000708` because it was slightly longer. Confirmed on the graph the nodes being connected to the Y.
