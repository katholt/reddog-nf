#!/bin/awk -f
BEGIN {
  OFS="\t"
}

$1 !~ /^#/ {
  replicons[$1]++
  hets++
}

END {
  if (length(replicons) > 1) exit 1
  if (length(replicons) == 0) exit 0
  for (replicon in replicons) printf "%s", replicon
  print "", hets
}
