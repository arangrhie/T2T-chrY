#!/bin/sh

echo "Usage: $0"
echo "Count the number of issues, sequence class length and issues length"
echo "Output: in issue_counts_lengths_by_seqclass.nonoverlap.txt"
echo -e "  Sample\tClass\tIssueCount\tClassLength\tIssueLength"

module load bedtools/2.31.1
tmp=issue_counts_by_seqclass.txt
# out=issue_counts_lengths_by_seqclass.txt
out=issue_counts_lengths_by_seqclass.nonoverlap.txt

echo -e "Sample\tSeqClass\tNumIssues\tClassLength\tIssueLength\tIssueType" > $out

for sample in $(cut -f1 sample_fa.map)
do
  echo $sample
  seqclass=../seqclass_ideogram/originY/$sample.seqclass.main.bed
  issues=../seqclass_ideogram/issues/$sample.issues.main.nonoverlap.bed

  for issuetype in ERRBASE ERRSTRUCT NGAP
  do
    echo "  issue: $issuetype"
    # sample <tab> sequence_class <tab> number_of_issues <tab> total_length_of_class
    awk -F '\t' -v issuetype=$issuetype '$4==issuetype' $issues | \
      bedtools intersect -c -a $seqclass -b - | \
        awk -F'\t' -v sample=$sample '{class_cnt[$4]+=$NF; class_len[$4]+=($3-$2); } END {
          for (class in class_cnt) {
            print sample"\t"class"\t"class_cnt[class]"\t"class_len[class];
          } }' | sort -k2,2V > $tmp
    
    # sample <tab> sequence_class <tab> number_of_issues <tab> total_length_of_class <tab> total_length_of_issues
    awk -F '\t' -v issuetype=$issuetype '$4==issuetype' $issues | \
      bedtools intersect -wao -a $seqclass -b - | \
        awk -F'\t' -v sample=$sample -v issuetype=$issuetype '{issues_len[$4]+=$NF} END {
          for (class in issues_len) {
            print class"\t"issues_len[class]"\t"issuetype;
          } }' | sort -k1,1V | join -1 2 -2 1 $tmp - | awk '{print $2"\t"$1"\t"$3"\t"$4"\t"$5"\t"$6}' \
        >> $out
    rm $tmp
  done
done

