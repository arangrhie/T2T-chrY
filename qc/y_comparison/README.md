# GQC Comparison using the new Ys as the reference

## GQC

## Output
* `*/*/*_oldY_vs_newY/*_oldY_vs_*_newY.refcovered.sort.bed`: Used for assembly issue comparison under `../asm_issues/`

* Copying over the bed files
```sh
for bed in $(ls /data/T2T-Y/y_comparison/GQC_AC/*/*/*_oldY_vs_newY/*_oldY_vs_*_newY.refcovered.sort.bed)
do
  dir_path=`dirname $bed`
  new_dir=`echo $dir_path | sed 's/\/data\/T2T-Y\/y_comparison\///g'`
  mkdir -p $new_dir
  set -x
  cp $bed $new_dir/
  set +x
done
```
