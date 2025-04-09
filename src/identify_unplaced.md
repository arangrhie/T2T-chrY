# Identify Y unplaced sequences

## Dependency
* samtools/1.21
* bedtools/2.31.1
* seqtk/1.4
* python/3.10

## Idea
Y usually comes out as its own component, or connected to X, with a few smaller components disconnected from the Y.

1. Identify Y chromosome components using the nodes used in the scaffold

2. Filter out noisy nodes

3. Collect and rename to 'random' unplaced sequences

## 1. Identify Y chromosome tigs

Several input files are used in this step, most are output from Verkko v2.2.1 and the analysis script, which runs mashmap on each nodes to the compressed, PAR / ACRO masked T2T-CHM13v2.


### scaffolds
Collect only the Genome and Seq parts to find the path in the graph for the unplaced part.
```sh
echo -e "Genome\tSeq" > T2TScaffolds.csv
cut -f1-2 chrY_scaffold.txt >> T2TScaffolds.csv
```

### paths.tsv
For the manual_update versions, analysis step wasn't run and the paths.tsv was not generated.
They were re-generated for convenience using `convert_path.py`, which takes `6-layoutContigs/consensus_paths.txt` as input.

### Traverse the graph, find all 'unplaced' patahs and remove 'noisy' nodes
Using the unitigs and paths identified in T2TScaffolds.csv, we can use them as 'reliable Y' and pull out all nodes connected to this component (sub-graph).
In addition, paths containing nodes with best hit to Y but not connected to non-Y chromosomes or rDNA are also added as a potential `UNPLACED` or `UNPLACED_Cmpnt` sequence.
This step is done using `CollectYPathAndNodes.py`, taking the following files as input:
* gfa: `analysis/assembly.homopolymer-compressed.add_telo_remove_rdna.noseq.gfa` or the default `assembly.homopolymer-compressed.noseq.gfa` - they are the same regardless of verkko versions (e.g. hi-c, thic, or trio)
* paths: `assembly.paths.tsv`
* scfmap: `assembly.scfmap`
* mashmap: `analysis/assembly-ref.comp.mashmap` - they are the same regardless of verkko versions (e.g. hi-c, thic, or trio), using the same node tigs as input
* rdna: `analysis/rdna.nodes` - they are the same regardless of verkko versions (e.g. hi-c, thic, or trio), using the same node tigs as input
* scaffold: `T2TScaffolds.csv`

The output is written to `tigs/${sample}.tigs.csv` (tab telimited), with the following in each column:
* SampleID
* Seq: Sequence name as in the assembly.fasta
* Resolution: `SCAFFOLD` / `UNPLACED` (connected to the component with `SCAFFOLD`) / `UNPLACED_Cmpnt` (Not connected to the component with `SCAFFOLD`)
* Scfmap: path name as in scfmap and paths.tsv
* Path: nodes and gaps

The graph is traversed one more time to check the connectivity of the singleton nodes identified as `UNPLACED` or `UNPLACED_Cmpnt` and its' unique sequence length and hifi coverage. This is implemented in `RemoveNoisyNodes.py`.
* Avg. node coverage: The node coverage is collected for the nodes in the SCAFFOLD Path, if the node coverage is between 10x and 100x
* All incoming and outgoing edge combinations are checked for its overlapping bps and extracted from the node length to get the unique bps in the unplaced node (length - overlaps)
* A node is considered 'noise' if
  * the node has only 1 incoming and 1 outgoing edge or is a singleton with less than 2 edges and
  * node coverage is over 0 (0 coverage indicates the same node is used in multiple places in the entire graph) and
  * node coverage is less than half the avg. node coverage
  * unique sequence length is less than 13 kb (`padding`, 20 kb hifi read length / est. homopolymer compression factor 1.5) or half the 13 kb if the node had only 1 edge

The output is written as `tigs/${sample}.tigs.no_noise.txt` with the same columns as in tigs.csv, excluding all unplaced sequences.

Now, the second column contains all the sequences that are used in the final Y chromosome assembly.

```sh
scaffold=T2TScaffolds.csv

samples=samples.wi_manual_update.map
num_samples=`wc -l $samples | awk '{print $1}'`
echo "Number of samples: $num_samples"
# 106

out=tigs
mkdir -p $out

for i in $(seq 1 $num_samples)
do
  ln=`sed -n ${i}p $samples`
  sample=`echo $ln | awk -F' ' '{print $1}'`
  gfa=`echo $ln | awk -F' ' '{print $3}'`
  paths=`echo $ln | awk -F' ' '{print $4}'`
  scfmap=`echo $ln | awk -F' ' '{print $5}'`
  mashmap=`echo $ln | awk -F' ' '{print $6}'`
  asm_path=`dirname $mashmap`
  rdna=$asm_path/rdna.nodes
  echo $sample
  python CollectYPathAndNodes.py -i $gfa -o $out/$sample.tigs.csv \
    -s $sample -m $mashmap -r $rdna -p $paths -f $scfmap -c $scaffold
  python RemoveNoisyNodes.py $gfa $out/$sample.tigs.csv
  echo
done
```

Any UNPLACED or UNPLACED_Cmpnt are pulled out to make the unplaced fa sequences.

## 2. Unplaced fa sequences

Each Seq name contains a unique identifier, regardless of its haplotype designation in Verkko.
Since we know the Y is all paternal, we ignore the automatic designation and simply are renaming the unplaced sequences to `${sample}_chrY_randomXXXXXXX`, with the `XXXXXXX` coming  from the unique identifier.

This is done with `collect_unplaced_fa.sh`.
Verkko has its own criteria and removes isolated singleton contigs if it meets certain criteria, so the 'unplaced' sequence may or may not be present in `assembly.fasta`. This is checked one more time and then placed into `unplaced/${sample}_chrY.unplaced.fa`.

## 3. Concatenate Scaffold and Unplaced
Lastly, the sequences can be simply concatenated to make the final `${sample}_chrY.fa`.

```sh
for sample in samples.list
do
  cat scaffolds/${sample}_chrY.fa unplaced/${sample}_chrY.unplaced.fa > T2T-chrYv2/${sample}_chrY.fa
  bgzip -i -@12 T2T-chrYv2/${sample}_chrY.fa
  samtools faidx -@12 T2T-chrYv2/${sample}_chrY.fa.gz
done
```

## 4. Annotate
As a sanity check, the assemblies were annotated with preliminary chosen sequences in `annotatea_regions.txt` using the colorcodes in the file.
The reference sequence is fixed to chrY for visualization convenience, with each bed file created for each sequence within the fa file (including all unplaced).

```sh
for fa in $(ls T2T-chrYv2/*.fa.gz)
do
  ./annotate.sh annotate_regions.txt $fa annotate_v2
done
```

## 5. Note on samples with `verkko_manual_update`
During the process of identifying unplaced contigs, we discovered T2T contigs with unplaced sequences - meaning there were collapsed sequences not used in the SCAFFOLD path. Some were identified as hi-c scaffolding errors, some were rukki errors. These mistakes are under investigation by the Verkko developers and may be included as fixes in a future version of the Verkko release.

Also preliminary annotation of the ampliconic genes identified some samples deemed to be interesting for futher investigation - carrying more copies of certain genes than others.

The sample's assembly graph (sometimes going up to the hifi graph) were manually investigated and the paths were re-built. If a hairpin-like connection was present with both walks with equal support, a conservative approach was chosen to prevent mis-orienting and the nodes in each path were put in the unplaced sequences. Some of the initial T2T contigs were no longer T2T contigs and had gaps for the unplaced sequences.

