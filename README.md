# T2T-chrY
HPRC / HGSVC Verkko r2 Y chromosome assemblies

## HPRC r2 Verkko Assemblies with the curated Y
* [chrY_assemblies](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/scratch/chrY/v2/chrY_assemblies): Each sequence is named `${SAMPLE}_chrY` for the Y scaffold (or contig) and `${SAMPLE}_chrY_randomXXXXXXX` for the unplaced sequences.
* [nonY_assemblies](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/scratch/chrY/v2/nonY_assemblies): Non Y sequence set, after removing all Y sequences (in `chrY_assemblies`) as well as small contigs that were in the assembly graph connected to the Y scaffold or had best hits to chrY with `mashmap` but were filtered out due to lack of coverage support or uniqueness of the sequences. Verkko `exemplar` sequences identified as `MT`, `rDNA`, and `EBV` (if any) are also included. For validation purposes, chrY_assemblies + nonY_assemblies should be concatenated and used for read mapping.
* [chrY_annotation](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/scratch/chrY/v2/chrY_annotation): Rough annotation of the ampliconic palindrome sequences and satellite repeats based on `minimap2`.
Note that all 'sequence' are set to chrY in the bed files to load them on the same chrY view for visualization purposes.
