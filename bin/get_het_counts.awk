#!/bin/awk -f
BEGIN {
  OFS="\t"
}

$1 !~ /^#/ {
  replicon_hets[$1]++
}

END {
  if (length(replicon_hets) == 0) exit 0
  for (replicon in replicon_hets)  {
    print replicon, replicon_hets[replicon]
  }
}
