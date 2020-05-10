#!/bin/awk -f
BEGIN {
  OFS="\t"
}

$1 !~ /^#/ {
  replicons[$1]
  if ($8 ~ /^I/) {
    indels[$1]++
  } else {
    snps[$1]++
  }
}

END {
  if (length(replicons) == 0) exit 0
  for (replicon in replicons) {
    # Get INDELs/SNPs counts, if none set to zero
    replicon_snps = replicon_indels = 0
    if (replicon in indels) {
      replicon_indels = indels[replicon]
    }
    if (replicon in snps) {
      replicon_snps = snps[replicon]
    }
    print replicon, replicon_snps, replicon_indels
  }
}
