#!/bin/awk -f
BEGIN {
  OFS="\t"
}

{
  depths[$1] += $2
  size[$1] += 1
  if ($2 > 0) {
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
    print replicon, average_depth, bases_mapped[replicon] / size[replicon] * 100, size[replicon]
  }
}
