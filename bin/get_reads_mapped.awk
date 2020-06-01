#!/bin/awk -f
BEGIN {
  OFS="\t"
}

{
  mapped_replicon[$3]++
  mapped_total++
}

END {
  print "total_mapped", mapped_total
  for (replicon in mapped_replicon) {
    print replicon, mapped_replicon[replicon]
  }
}
