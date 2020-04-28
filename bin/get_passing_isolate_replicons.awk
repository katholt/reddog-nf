#!/usr/bin/awk -f
# Record samples and associated replicons that:
#   - pass mapping filter
#   - have more than one SNP
! /isolate/ && $11 == "p" && $8 > 0 {
  replicons_pass[$2][$1]
}

END {
  for (sample in replicons_pass) {
    if (length(replicons_pass[sample]) > 0) {
      printf "%s", sample
      for (replicon in replicons_pass[sample]) {
        printf "\t%s", replicon
      }
      printf "\n"
    }
  }
}
