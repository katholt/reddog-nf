#!/bin/awk -f
BEGIN {
  OFS="\t"
}

{
  # Count mapped reads for each replicon - ummapped flag not set
  if (and($2, 0x4) != 0x4) {
    mapped_replicon[$3]++
    mapped_total++
  }
    reads_total++
}

END {
  print "total_reads", reads_total
  print "total_mapped", mapped_total
  for (replicon in mapped_replicon) {
    print replicon, mapped_replicon[replicon]
  }
}
