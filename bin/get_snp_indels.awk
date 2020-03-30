#!/bin/awk
BEGIN {
    OFS="\t"
    snps = indels = 0
  } $1 !~ /^#/ {
    replicons[$1]++
    if ($8 ~ /^I/) {
      indels++
    } else {
      snps++
    }
  } END {
    if (length(replicons) > 1) exit 1
    if (length(replicons) == 0) exit 0
    for (replicon in replicons) printf "%s", replicon
    print "", snps, indels
}
