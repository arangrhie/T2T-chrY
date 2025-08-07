#!/bin/sh

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <fa.list> <output_dir>"
  echo "Example: $0 T2T-chrYv1_fa.list annotate_v2"
  exit 1
fi

i=$SLURM_ARRAY_TASK_ID
if [[ -z $i ]]; then
  fa_list=$1
else
  fa_list=$1.$i
fi

out=$2

# for fa in $(cat T2T-chrYv1_fa.list.$i)
# for fa in $(cat T2T-chrYv1_fa.scaffold.list unplaced_fa.list)
# for fa in $(cat manual_update_fa.list)
for fa in $(cat $fa_list)
do
  # sample=`echo $fa | awk -F '/' '{print $(NF-3)}'` # for naive refOriented verkko output fa paths
  # sample=`echo $fa | awk -F '/' '{print $NF}' | sed 's/_chrY.fa.gz//g'` # for renamed chrY fa paths
  sample=`echo $fa | awk -F '/' '{print $(NF-1)}'`
	$tools/T2T-chrY/src/annotate.sh $tools/T2T-chrY/src/annotate_regions.txt $sample $fa $out
done

