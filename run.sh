#!/bin/bash

### user data
if [ "$#" -ne 4 ]; then
  echo -e "Illegal number of parameters, usage:"
  echo -e "\trun.sh raw_read1 raw_read2 mapping_table output_prefix"
  exit 0
else
  raw_read1=$1
  raw_read2=$2
  mapping_table=$3
  grpbase=$4
fi

### run
exec ./data2usearch.sh ${raw_read1} ${raw_read2} ${mapping_table} ${grpbase}
Rscript ./usearch2phyloseq.R ${mapping_table} ${grpbase}
