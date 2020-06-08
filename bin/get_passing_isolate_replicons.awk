#!/usr/bin/awk -f
# Record samples and associated replicons that:
#   - pass mapping filter
#   - have more than one SNP
NR != 1 && $10 == "p" && $7 > 0 {
  sub(/_mapping_stats.tsv/, "", FILENAME); print $1"\t"FILENAME
}
