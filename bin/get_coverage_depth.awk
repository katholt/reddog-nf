#!/bin/awk -f
BEGIN {
  OFS="\t"
}

{
  depths[$1] += $4
  size[$1] += 1
  if ($4 > 0) {
    bases_mapped[$1] += 1
  }
}

END {
  for (replicon in depths) {
    if (bases_mapped[replicon] == 0) {
      average_depth = 0.0
    } else {
      average_depth = depths[replicon] / bases_mapped[replicon]
    }
    printf "%s\t%.2f\t%.2f\t%d\n", replicon, average_depth, bases_mapped[replicon] / size[replicon] * 100, size[replicon]
  }
}
