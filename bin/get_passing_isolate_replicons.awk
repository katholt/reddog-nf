#!/usr/bin/awk -f
# Output isolates and associated replicons that pass the mapping filter
NR != 1 && $10 == "p" {
  sub(/_mapping_stats.tsv/, "", FILENAME); print $1"\t"FILENAME
}
