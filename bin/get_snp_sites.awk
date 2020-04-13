#!/bin/awk -f
BEGIN {
  OFS="\t"
  print "#CHROM", "POS"
}

# Print variant position - regexes do the following:
#   1. Skip commented lines
#   2. Skip INDELs (info field starting with 'I')
#   3. Skip non-SNPs (alt field containing ',')
#   4. Skip when ref and alt are the same
/^[^#]/ && $8 !~ /^I/ && $4 !~ /,/ && $4 != $5 {
  print $1, $2
}
