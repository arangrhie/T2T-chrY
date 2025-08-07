# HGSVC chrYs

* Y Scaffold: Backbone Y chromosome assembly
* Y Unplaced: Sequences unable to place within the scaffold, due to ambiguity in the orientation

Assemblies were generated using the same Verkko v2.2.1 versions, in hi-c mode or trio and thic (trio+hic) mode when parental reads were available.

The assembly graph was re-walked to generate a more contiguous path for the Ys that were manually resolvable, and are reflacted as 'manual' in the below table. Below shows which verkko mode was chosen for the final assembly scaffold and unplaced identification.

## Sample IDs and Verkko mode used

`samples.wi_manual.list`
| sample  | verkko mode |
| ------- | ----------- |
| HG00096 | manual      |
| HG00358 | manual      |
| HG00512 | manual      |
| HG00731 | manual      |
| HG01457 | manual      |
| HG01505 | thic        |
| HG01596 | manual      |
| HG01890 | hi-c        |
| HG02011 | manual      |
| HG02492 | thic      |
| HG02554 | hi-c        |
| HG02666 | manual      |
| HG02953 | hi-c        |
| HG03009 | manual      |
| HG03065 | manual      |
| HG03248 | manual      |
| HG03371 | thic        |
| HG03456 | manual      |
| HG03732 | thic        |
| NA18534 | hi-c        |
| NA18989 | manual      |
| NA19239 | manual      |
| NA19317 | manual      |
| NA19331 | hi-c        |
| NA19347 | manual      |
| NA19384 | manual      |
| NA19650 | thic        |
| NA19705 | manual      |
| NA20509 | manual      |

## Identify scaffolds

All Y assemblies have been checked for T2T status and orientation. The refOriented are better tuned to reflect this.
This time, we only need to greb the scaffold names with the more complete Y except for HG03371 and NA19384.

Collect path to the refOriented assembly:
```sh
rm samples.wi_manual.map

echo -e "Genome\tSeq\tLength\tnTels\tnGaps" > T2TScaffolds.csv

# Waiting for 3 23 - HG00512 NA19317 / 10 HG02492 - should be thic, not manual
for i in 3 23 # 10 # $(seq 1 29)
do
  ln=`sed -n ${i}p samples.wi_manual.list`
  sample=`echo $ln | awk '{print $1}'`
  mode=`echo $ln | awk '{print $2}'`
  if [[ $mode == "manual" ]]; then
    mode="manual-update"
  fi
  echo $i $sample

  asmdir="/data/Phillippy2/projects/hgsvc/assemblies-hgsvc/$sample/verkko-$mode"
  fasta=`ls $asmdir/analysis/assembly.refOriented.fasta`
  gfa=`ls $asmdir/analysis/assembly.homopolymer-compressed.add_telo_remove_rdna.noseq.gfa`
  paths=`ls $asmdir/assembly.paths.tsv`
  scfmap=`ls $asmdir/assembly.scfmap`
  mashmap=`ls $asmdir/analysis/assembly-ref.comp.mashmap`

  # Some are under re constructing the consensus, let's skip those
  if [[ -s $fasta ]]; then
    echo -e "$sample\t$fasta\t$gfa\t$paths\t$scfmap\t$mashmap" >> samples.wi_manual.map
    
    # Pick the larger chrY refOriented sequence
    if ! [[ -s $fasta.fai ]]; then
      module load samtools
      echo "index $fasta"
      samtools faidx -@12 $fasta
    fi
    scaff=`grep chrY $fasta.fai | awk '{print $1"\t"$2}' | sort -k2,2nr | head -n1`

    # Collect TEL and GAP status from the initial annotation files
    scaff_n=`echo $scaff | awk '{print $1}'`

    if [[ $mode == "manual-update" ]]; then mode=manu; fi
    if [[ $mode == "hi-c" ]]; then mode=hic; fi

    scaff_bed=`ls annotate_$mode/$sample/${sample}_$scaff_n.bed`
    numTel=`grep -c TEL $scaff_bed`
    numGap=`grep -c GAP $scaff_bed`
    echo -e "$sample\t$scaff\t$numTel\t$numGap\t$mode" >> T2TScaffolds.csv
  fi
done
```
* HG03456 has 2 Y chromosomes (XYY karyotype). The second Y (chrY_haplotype1-0000015) was manually added to the T2TScaffolds.csv.
* HG03371 and NA19384 are missing the qTel and the pTel, also disconnected on the graph.

## Identify unplaced

```sh
module load python/3.10
scaffold=T2TScaffolds.csv

samples=samples.wi_manual.map
num_samples=`wc -l $samples | awk '{print $1}'`
echo "Number of samples: $num_samples"
# 29

out=tigs
mkdir -p $out

for i in 28 29 # 27 # $(seq 1 $num_samples)
do
  ln=`sed -n ${i}p $samples`
  sample=`echo $ln  | awk '{print $1}'`
  gfa=`echo $ln     | awk '{print $3}'`
  paths=`echo $ln   | awk '{print $4}'`
  scfmap=`echo $ln  | awk '{print $5}'`
  mashmap=`echo $ln | awk '{print $6}'`
  asm_path=`dirname $mashmap`
  rdna=$asm_path/rdna.nodes
  echo $sample
  python $tools/T2T-chrY/src/CollectYPathAndNodes.py -i $gfa -o $out/$sample.tigs.csv \
    -s $sample -m $mashmap -r $rdna -p $paths -f $scfmap -c $scaffold
  python $tools/T2T-chrY/src/RemoveNoisyNodes.py $gfa $out/$sample.tigs.csv
  echo
done
```

* Note: Dropping 1 unplaced in HG00512, haplotype2-0000185 which is a dummy node in the PAR1 connecting Y and X PAR1.


## Extract final Y sequences

```sh
./extract_fa.sh
```
* Final sequences are under `hgsvc_chrYv1`, chrY scaffold and unplaced named ${sample}_chrY.fa.gz and indexed.
* nonY sequences are under `hgsvc_nonYv1` and named ${sample}_nonY.fa.gz and indexed.

### scaffolds
Unlike the HPRC versions, scaffolds are already properly oriented for the Y scaffold. Thus, we only need to extract the sequence and rename. The new name will be `${sample}_chrY`.

### unplaced
Using those under `tigs/$sample.tigs.no_noise.txt`, renaming to match `${sample}_chrY_randomXXXXXXX`.

### nonY
Make sure we don't include the noisy nodes excluded from the Y.



## Annotate

```sh
list=hgsvc_chrYv1.list
split -l5 -d -a2 --numeric-suffixes=10 $list $list.
sh ~/codes/_submit_quick.sh 24 32g annotate_v1 annotate_falist.sh "$list annotate_v1" "--array=10-15"

# The 2 that came in later
ls $PWD/hgsvc_chrYv1/HG00512_chrY.fa.gz >  $list.16
ls $PWD/hgsvc_chrYv1/NA19317_chrY.fa.gz >> $list.16
sh ~/codes/_submit_quick.sh 24 32g annotate_v1 annotate_falist.sh "$list annotate_v1" "--array=16"
61827404
```


## Upload
```sh
module load aws

cd hgsvc_chrYv1/
aws s3 sync . s3://human-pangenomics/T2T/scratch/chrY/hgsvc/v1/chrY_assemblies/

cd ../hgsvc_nonYv1/
aws s3 sync . s3://human-pangenomics/T2T/scratch/chrY/hgsvc/v1/nonY_assemblies/

cd ../hgsvc_chrYv1_annotation/
# for bed in $(ls ../annotate_v1/*/*.bed); do ln -s $bed; done
ln -s ../annotate_v1/HG00512/HG00512_chrY.bed
ln -s ../annotate_v1/NA19317/NA19317_chrY.bed
aws s3 sync . s3://human-pangenomics/T2T/scratch/chrY/hgsvc/v1/chrY_annotation/
```