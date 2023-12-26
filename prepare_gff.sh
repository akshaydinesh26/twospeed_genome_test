#!/usr/bin/bash

# obtain gene line from gff file.
gff=$1;
out=$2;
#create bed with chr start end strand name
grep -v "^#" $gff | awk -F"[\t;]" 'BEGIN{OFS=",";} $3=="gene" {print $1,$4,$5,$7,$9}' | sed "s/ID=//g" > $out;