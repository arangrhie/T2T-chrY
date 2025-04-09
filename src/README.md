# Curating Y assemblies

Initial Y assemblies were identified based on `mashmap`, by mapping the `assembly.fasta` against a masked version of the T2T-CHM13v2. The reference was masked for PARs and the acrocentric short arms. Chromosomal sequences were identified based on the longest matched sequences and oriented accordingly ([freeze2 Verkko assemblies](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=submissions/6807247E-4F71-45D8-AECE-9E5813BA1D9F--verkko-v2.2.1-release2_asms/), under `${sample}/verkko-${mode}/analysis/${sample}.assembly.refOriented.fasta.gz`).

The ampliconic palindromic sequences and the heterochromatic sequences are challenging to map and compare against each other, and were sometimes erroneously chosen for the inverted orientation compared to T2T-CHM13 (HG002Y). 
In addition, pieces of Y sequences (often in PAR1 or heterochromatin) were not connected due to the ambiguity within the XPAR1 and YPAR1 path or not enough Hi-C reads to be confidently map in the heterochromatin sequences.

Thus, refOriented sequences were mapped to a further masked version of the Y to extend chrY scaffolds and recover missing telomeres and PAR sequences. See [identify_scaffold.md](TODO).

Once the Y scaffolds were identified, each node in the graph connected to the Y scaffold were visited and pulled out as a putative unplaced Y sequence. This was done to rescue centromeric and ampliconic repeats that Verkko was unable to resolve automatically. <br>
Components not connected to the Y scaffolds but with strong matches to chrY were also visited to validate and put in the unplaced bin. See [identify_unplaced.md](identify_unplaced.md).

During the process, the assembly path of [10 genomes](verkko_manual_update.samples.list) were manually re-evaluated and corrected if there was a mis-assembly.

Finally, the assemblies were renamed to
* `${sample}_chrY`: scaffold, with 200kb gap if two sequences were brought together and oriented accordingly from p to q arm
* `${sample}_chrY_randomXXXXXXX`: unplaced, orientation follows as in the original Verkko assembly.fasta file (at random)

Preliminary annotation of the sequences were used to validate the final assembly orientation and composition.


