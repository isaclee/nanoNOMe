#!/bin/bash
root=/mnt/d/Data/analysis
mdir=$root/mcall
mbeddir=$root/mbed
mods="gpc"
scr=/home/ubuntu/Code/nanopore-methylation-utilities/mtsv2bedGraph.py

for mod in $mods; do
  dir=$mdir
  outdir=$mbeddir
  [ -e $outdir ]||mkdir $outdir
  mtsvs=$(find $dir -name "*$mod*tsv.gz")
  for mtsv in $mtsvs; do
    base=$(basename "${mtsv%%.*}")
    echo $base
    mbed=$outdir/$base.$mod.meth.bed.gz
    com="gunzip -c $mtsv | head -n1000 | $scr -m $mod --nome |\
      sort -k1,1 -k2,2n | bgzip > $mbed"
#    eval $com
#    exit #testing
  done
done
