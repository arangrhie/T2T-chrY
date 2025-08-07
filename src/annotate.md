# Preliminary region annotation of the Y region



## Extracting target regions to annotate
For each region defined in `annotate_region.txt`, a target region sequence is extracted in its desired orientation and renamed using `generate_targets.sh`.
* `annotate_region.txt`
```
HSat1A	chr13:434802-563385	-	0,222,96
HSat1B	chr13:10165987-10219675	+	27,153,139
HSat3	chr13:5000117-5152650	-	51,81,137
HSat3Y	chrY:58654600-58660873	-	51,0,153
aSat	chr13:15547594-17498291	-	153,0,0
bSat	chr13:4973578-4985238	+	250,153,255
ACRO	chr13:5338347-5358151	+	0,204,204
CER	chr13:5369401-5424523	-	0,204,204
CEN	chrY:10699757-10745416	-	153,0,0
SST1	chr13:12301367-12440010	-	0,204,204
DJarm	chr13:5555305-5643740	+	204,0,204
DJflank	chr13:5650824-5770548	+	153,0,153
rDNA	chr13:5770549-5821272	+	153,102,255
P5-AZFb	chrY:18878839-18989510	+	255,153,0
teal2	chrY:23213314-23328994	+	0,128,128
P3-AZFc	chrY:23327300-23566501	+	60,84,164
blue2	chrY:23328993-23582346	+	60,84,164
red2	chrY:24037533-24167130	+	231,35,35
green2	chrY:24167430-24482113	+	11,130,67
yellow1	chrY:24482112-25159839	+	255,153,0
gray1	chrY:25389078-25503792	+	125,112,96
SEQ		.	.	204,204,204
TEL		.	.	255,102,51
GAP		.	.	0,0,0
```

## Annotate
Annotation was performed by `annotate.sh` on each chrY assembly, submited via `annotate_falist.sh`.

An input fasta file is then mapped to each target using `minimap2`, and its alignments is reversed for reference and target using `rustybam`. The resulting output is then merged back into one bed file using the color codes in `annotate_region.txt`.

### HPRC assembly annotation
Under `/data/Phillippy/projects/hprc-assemblies/assign_Y/`
```sh
ls T2T-chrYv2/*.fa > T2T-chrYv2_fa.list

split -l 5 -d -a 2 --numeric-suffixes=10 T2T-chrYv2_fa.list T2T-chrYv2_fa.list.
sh ~/codes/_submit_norm.sh 12 12g annotate_chrYv2 annotate_falist.sh "T2T-chrYv2_fa.list annotate_v2" "--array=10-30"
```

### HGSVC assembly annotation
Under `/data/Phillippy2/projects/hgsvc/annotation`.
* b1 - batch 1, 21 hi-c scaffolded refOriented fasta files.

```sh
cd ../
ls $PWD/assemblies-hgsvc/*/verkko-hi-c/analysis/*.refOriented.fasta > annotation/refOriented_b1.list
cd annotation

list=refOriented_b1.list

split -l5 -d -a2 --numeric-suffixes=10 $list $list.
sh ~/codes/_submit_norm.sh 24 32g annotate_v1 annotate_falist.sh "$list annotate_v1" "--array=10-14"
```
*later renamed `annotate_v1` to `annotate_hic`

* b2 - batch 2, 8 hi-c scaffolded refOriented fasta files.
```sh
for sample in $(cat samples_b2.list)
do
  ls /data/Phillippy2/projects/hgsvc/assemblies-hgsvc/$sample/verkko-hi-c/analysis/*.refOriented.fasta >> refOriented_b2.list
done

list=refOriented_b2.list

split -l1 -d -a2 --numeric-suffixes=10 $list $list.
sh ~/codes/_submit_norm.sh 24 32g annotate_hic annotate_falist.sh "$list annotate_hic" "--array=10-16" # 59269253
sh ~/codes/_submit_norm.sh 24 32g annotate_hic annotate_falist.sh "$list annotate_hic" "--array=17"
```


