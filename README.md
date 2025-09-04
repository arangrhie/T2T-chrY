# T2T-chrY
HPRC / HGSVC / CEPH Verkko r2 Y chromosome assemblies

## What's available
* chrY_assemblies: Each sequence is named `${SAMPLE}_chrY` for the Y scaffold (or contig) and `${SAMPLE}_chrY_randomXXXXXXX` for the unplaced sequences.
* nonY_assemblies: Non Y sequence set, after removing all Y sequences (in `chrY_assemblies`) as well as small contigs that were in the assembly graph connected to the Y scaffold or had best hits to chrY with `mashmap` but were filtered out due to lack of coverage support or uniqueness of the sequences. Verkko `exemplar` sequences identified as `MT`, `rDNA`, and `EBV` (if any) are also included for the HPRC / HGSVC samples. For validation purposes, chrY_assemblies + nonY_assemblies should be concatenated and used for read mapping.
* chrY_annotation: Rough annotation of the ampliconic palindrome sequences and satellite repeats based on `minimap2`.
Note that all 'sequence' are set to chrY in the bed files to load them on the same chrY view for visualization purposes.

## Data access
Below can be downloaded with aws cli:
```sh
mkdir -p hprc hgsvc ceph
aws s3 cp --recursive s3://human-pangenomics/T2T/scratch/chrY/v2/  hprc/
aws s3 cp --recursive s3://human-pangenomics/T2T/scratch/hgsvc/v1/ hgsvc/
aws s3 cp --recursive s3://human-pangenomics/T2T/scratch/ceph/v1/  ceph/
```

* HPRC: 106 genomes
  * [chrY_assemblies](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/scratch/chrY/v2/chrY_assemblies)
  * [nonY_assemblies](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/scratch/chrY/v2/nonY_assemblies)
  * [chrY_annotation](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/scratch/chrY/v2/chrY_annotation)

* HGSVC: 29 genomes
  * [chrY_assemblies](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/scratch/chrY/hgsvc/v1/chrY_assemblies/)
  * [nonY_assemblies](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/scratch/chrY/hgsvc/v1/nonY_assemblies)
  * [chrY_annotation](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/scratch/chrY/hgsvc/v1/chrY_annotation)

* CEPH: 8 genomes; 2 will be shared later via dbGAP
  * [chrY_assemblies](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/scratch/chrY/ceph/v1/chrY_assemblies)
  * [nonY_assemblies](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/scratch/chrY/ceph/v1/nonY_assemblies)
  * chrY_annotation; TBA

Below are details from each cohort.

## HPRC

### Samples manually updated
In total, 11 samples are not using the same assembly as in the initial [Verkko v2.2.1 release 2](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=submissions/6807247E-4F71-45D8-AECE-9E5813BA1D9F--verkko-v2.2.1-release2_asms/).

### Version used and its resolution

| sample | verkko mode | resolution |
| ------ | ----------- | ---------- |
| HG01167 | manual | Initial Y walked through the shortest path, ignoring one unplaced sequence. Resolved as one T2T contig. |
| HG01192 | manual | Initial Y had 2 unplaced sequences with no gaps. Resolved with 6 gaps in a T2T scaffold and 6 unplaced contigs which couldn't unambiguously walked through. |
| HG02015 | hi-c | Initial Y was mis-orienting the entire p arm to hetero chromatin. Using Hi-C version instead. Likely excluding from analysis as the assembly is too fragmented. |
| HG02486 | manual | Initial Y walked through the shortest path, ignoring one unplaced sequence. Resolved as one T2T scaffold with a gap for the missing edge. |
| HG02514 | manual | Initial Y walked through the shortest path, ignoring three unplaced sequence. Resolved as one T2T contig. |
| HG03139 | manual | Initial Y walked through the shortest path, ignoring one unplaced sequence. Resolved as one T2T scaffold with a gap for the missing edge. |
| HG03688 | manual | Initial Y walked through the shortest path, ignoring two unplaced sequence. Resolved as one T2T scaffold with a gap for the missing edge. |
| HG04187 | manual | Initial Y too fragmented, manually checked and walked through. |
| NA18608 | manual | Initial Y mis-scaffolded with X, manually fixed and merged the contig manually for the missing edge. Resolved as one T2T contig. |
| NA18971 | manual | Initial Y too fragmented, manually checked and walked through. Resolved as one T2T scaffold with two gaps and one unplaced sequence. |
| NA18974 | manual | Initial Y too fragmented, manually checked and walked through. Resolved as one T2T scaffold with two gaps and ten unplaced sequences, all in the centromere. |

### Y scaffolds renamed and its original sequences used to build the scaffold

| NewName | scaffolded sequences in order, wi orientation | Assembly ver. |
| ------- | --------------------------------------------- | ------------- |
| HG00126_chrY | chrY_haplotype1-0000017+ | refOriented |
| HG00140_chrY | chrY_haplotype1-0000004+ | refOriented |
| HG00280_chrY | chrY_haplotype1-0000024+ | refOriented |
| HG00290_chrY | chrY_haplotype1-0000013+ | refOriented |
| HG002_chrY | chrY_pat-0000702+ | refOriented |
| HG00321_chrY | chrY_haplotype2-0000065+_chrY_haplotype2-0000041+ | refOriented |
| HG00329_chrY | chrY_haplotype1-0000012+ | refOriented |
| HG005_chrY | chrY_pat-0001218+_chrY_pat-0001233+ | refOriented |
| HG00609_chrY | chrY_haplotype1-0000011+ | refOriented |
| HG00621_chrY | chrY_pat-0000611- | refOriented |
| HG00642_chrY | chrY_haplotype1-0000004+ | refOriented |
| HG00658_chrY | chrY_pat-0000544+_chrY_pat-0000537+ | refOriented |
| HG00673_chrY | chrY_haplotype1-0000012+_chrY_haplotype2-0000096+ | refOriented |
| HG00706_chrY | chrY_haplotype1-0000020+ | refOriented |
| HG00738_chrY | chrY_pat-0000564+ | refOriented |
| HG01074_chrY | chrY_haplotype1-0000019+ | refOriented |
| HG01099_chrY | chrY_haplotype1-0000025+ | refOriented |
| HG01106_chrY | chrX_mat-0000007+_chrY_pat-0000759+ | refOriented |
| HG01109_chrY | chrY_pat-0000969+_chrY_pat-0000951+ | refOriented |
| HG01243_chrY | chrY_haplotype1-0000012+_chrY_haplotype1-0000021+ | refOriented |
| HG01252_chrY | chrY_pat-0000532+_pat-0000707+ | refOriented |
| HG01255_chrY | chrY_haplotype1-0000016+ | refOriented |
| HG01258_chrY | chrY_pat-0000774+ | refOriented |
| HG01261_chrY | chrY_haplotype1-0000012+ | refOriented |
| HG01358_chrY | chrY_pat-0000434+ | refOriented |
| HG01433_chrY | chrY_pat-0000521+_chrY_pat-0000511+_chrY_pat-0000498+ | refOriented |
| HG01530_chrY | chrY_haplotype1-0000004+ | refOriented |
| HG01928_chrY | chrY_pat-0000542+ | refOriented |
| HG01934_chrY | chrY_haplotype1-0000010+ | refOriented |
| HG01943_chrY | chrY_haplotype1-0000026- | refOriented |
| HG01952_chrY | chrY_haplotype1-0000014+ | refOriented |
| HG02015_chrY | chrY_haplotype1-0000001+_chrY_haplotype2-0000121+ | refOriented, hi-c |
| HG02027_chrY | chrY_pat-0000390+ | refOriented |
| HG02040_chrY | chrY_haplotype1-0000017- | refOriented |
| HG02055_chrY | chrY_pat-0000629+_chrY_pat-0000631+ | refOriented |
| HG02056_chrY | chrY_pat-0000451+ | refOriented |
| HG02071_chrY | chrY_pat-0000661- | refOriented |
| HG02074_chrY | chrY_pat-0000739- | refOriented |
| HG02083_chrY | chrY_pat-0000362+ | refOriented |
| HG02129_chrY | chrY_haplotype1-0000005+_chrY_haplotype2-0000101+ | refOriented |
| HG02132_chrY | chrY_haplotype2-0000165+_chrY_haplotype2-0000173+ | refOriented |
| HG02135_chrY | chrY_pat-0000670+ | refOriented |
| HG02145_chrY | haplotype1-0000075+_chrY_haplotype1-0000002+_chrY_haplotype1-0000157+_chrY_haplotype1-0000004+_chrY_haplotype1-0000024+_chrY_haplotype2-0000213- | refOriented |
| HG02258_chrY | chrY_haplotype1-0000016+ | refOriented |
| HG02391_chrY | chrY_haplotype1-0000014+ | refOriented |
| HG02392_chrY | chrY_haplotype1-0000002+ | refOriented |
| HG02572_chrY | chrY_haplotype1-0000021+ | refOriented |
| HG02602_chrY | chrY_pat-0000505+_chrY_pat-0000506- | refOriented |
| HG02647_chrY | chrY_pat-0000655- | refOriented |
| HG02668_chrY | chrY_pat-0000504-_chrY_pat-0000487+ | refOriented |
| HG02698_chrY | chrY_haplotype1-0000017+ | refOriented |
| HG02717_chrY | chrY_mat-0000024+ | refOriented |
| HG02735_chrY | chrY_haplotype1-0000003+ | refOriented |
| HG02738_chrY | chrY_pat-0000697+ | refOriented |
| HG02965_chrY | chrY_haplotype1-0000013+_chrY_haplotype2-0000167+ | refOriented |
| HG02984_chrY | chrY_haplotype1-0000013- | refOriented |
| HG03017_chrY | chrY_haplotype1-0000011- | refOriented |
| HG03050_chrY | chrY_pat-0000469+ | refOriented |
| HG03098_chrY | chrY_haplotype1-0000004+_chrY_haplotype1-0000006+ | refOriented |
| HG03130_chrY | chrY_haplotype1-0000006+_chrY_haplotype2-0000110+ | refOriented |
| HG03209_chrY | chrY_haplotype1-0000002+ | refOriented |
| HG03225_chrY | chrY_haplotype1-0000015- | refOriented |
| HG03471_chrY | chrY_haplotype1-0000006+ | refOriented |
| HG03492_chrY | chrY_pat-0000829-_chrY_pat-0000828+ | refOriented |
| HG03521_chrY | chrY_haplotype1-0000022+ | refOriented |
| HG03579_chrY | chrY_haplotype1-0000023+ | refOriented |
| HG03654_chrY | chrY_haplotype1-0000012+ | refOriented |
| HG03710_chrY | chrY_haplotype1-0000009+ | refOriented |
| HG03742_chrY | chrY_haplotype1-0000022+ | refOriented |
| HG03942_chrY | chrY_pat-0000811-_chrY_pat-0000813+ | refOriented |
| HG04115_chrY | chrY_pat-0001036+ | refOriented |
| HG04157_chrY | chrY_haplotype1-0000029+_chrY_haplotype1-0000023+_chrY_haplotype2-0000084+ | refOriented |
| HG04160_chrY | chrY_haplotype1-0000022+ | refOriented |
| HG04199_chrY | chrY_haplotype1-0000001+ | refOriented |
| HG04204_chrY | chrY_haplotype1-0000022+ | refOriented |
| HG04228_chrY | chrY_haplotype1-0000011+ | refOriented |
| NA18522_chrY | chrY_haplotype1-0000020+ | refOriented |
| NA18612_chrY | chrY_haplotype2-0000088-_chrY_haplotype1-0000005+ | refOriented |
| NA18620_chrY | chrY_haplotype1-0000001+_chrY_unassigned-0000173+ | refOriented |
| NA18747_chrY | chrY_haplotype1-0000010+ | refOriented |
| NA18879_chrY | chrY_haplotype1-0000021+ | refOriented |
| NA18952_chrY | chrY_haplotype1-0000008- | refOriented |
| NA18983_chrY | chrY_haplotype1-0000006+ | refOriented |
| NA19043_chrY | chrY_haplotype1-0000017+ | refOriented |
| NA19443_chrY | chrY_haplotype1-0000025+_chrY_haplotype2-0000191+ | refOriented |
| NA19682_chrY | chrY_haplotype1-0000008- | refOriented |
| NA19700_chrY | chrY_haplotype1-0000009+ | refOriented |
| NA20346_chrY | chrY_haplotype1-0000017+ | refOriented |
| NA20752_chrY | chrY_haplotype1-0000011- | refOriented |
| NA20805_chrY | chrY_haplotype1-0000003+ | refOriented |
| NA20809_chrY | chrY_haplotype1-0000030+_chrY_haplotype2-0000094+ | refOriented |
| NA20850_chrY | chrY_haplotype1-0000018+ | refOriented |
| NA20870_chrY | chrY_haplotype1-0000014- | refOriented |
| NA20905_chrY | chrY_haplotype1-0000006+ | refOriented |
| NA21093_chrY | chrY_haplotype1-0000021+ | refOriented |
| HG01167_chrY | haplotype1-0000001+ | manual |
| HG01192_chrY | haplotype1-0000011+ | manual |
| HG02486_chrY | pat-0000295+ | manual |
| HG02514_chrY | haplotype1-0000011- | manual |
| HG03139_chrY | haplotype1-0000013- | manual |
| HG03688_chrY | haplotype1-0000013+ | manual |
| HG04187_chrY | pat-0000905+ | manual |
| NA18608_chrY | haplotype1-0000026+ | manual+patching |
| NA18971_chrY | haplotype1-0000013+ | manual |
| NA18974_chrY | haplotype1-0000003- | manual |
| HG003_chrY | chrY_haplotype1-0000002+ | refOriented, not part of HPRC release |

### Samples likely to exclude from analysis

| Set  | Sample  | Reason  |
| ---- | ------- | ------- |
| HPRC | HG02015 | Assembly too fragmented. 18 gaps in the scaffold with 84 unplaced sequences. |
| HPRC | HG02145 | Y is 8x vs. 35x in the rest of the haplotype. Cellline is likely loosing Y in some cells. Assembly too fragmented, missing PAR/telomeres. 9 gaps in the scaffold with 30 unplaced sequences. |

## HGSVC
See [hgsvc](hgsvc/hgsvc.md).

## CEPH
TBA