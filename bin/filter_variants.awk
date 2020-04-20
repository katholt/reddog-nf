#!/bin/awk -f
BEGIN {
  # Check arguments
  if (snp_fp == "") {
		print "error: -v snp_fp=<filepath> is required" > "/dev/stderr"
    exit 1
  }
  if (het_fp == "") {
		print "error: -v het_fp=<filepath> is required" > "/dev/stderr"
    exit 1
  }
}

# Print header
/^#/ {
  print $0 > snp_fp
  print $0 > het_fp
}

# Print passing SNPs and hets to appropriate output
# Variants with quality >= 30
! /^#/ && $6 >= 30 {
  # Check allele frequency and genotype allele (0 refers to the reference allele)
  # Check that we have only one alternative allele
  # Write to appropriate output
  if (($8 ~ /AF1=1/ || $10 !~ /^0/) && $5 !~ /,/) {
    print $0 > snp_fp
  } else if ($8 !~ /^IND/) {
    print $0 > het_fp
  }
}
