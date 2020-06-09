#!/bin/awk -f
BEGIN {
  # Check arguments
  if (mapped_fp == "") {
    print "error: -v mapped_fp=<filepath> is required" > "/dev/stderr"
    exit 1
  }
  if (unmapped_fp == "") {
    print "error: -v unmapped_fp=<filepath> is required" > "/dev/stderr"
    exit 1
  }
  # Force creation of outputs
  system("touch " mapped_fp)
  system("touch " unmapped_fp)
  # Create pipe processes
  mapped_ps = "samtools sort > " mapped_fp
  unmapped_ps = "samtools view -OBAM > " unmapped_fp
}

# Header
/^@(HD|SQ|RG|PG)(\t[A-Za-z][A-Za-z0-9]:[ -~]+)+$/ {
  print $0 | mapped_ps
  print $0 | unmapped_ps
}

# Alignment entries
{
  # Check FLAG for unmapped bit
  if (and($2, 0x4) == 0x4) {
    print | unmapped_ps
  } else {
    print | mapped_ps
  }
}
