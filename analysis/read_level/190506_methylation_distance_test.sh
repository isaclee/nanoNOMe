#!/bin/bash
root="/kyber/Data/Nanopore/projects/nanonome/analysis"
scr="/home/isac/Code/nanoNOMe/scripts/methylation_distance.py"
outdir=$root/plots/mdist

mbed="$root/data/nanonome/pooled/mbed/GM12878_nanoNOMe.pooled.gpc.meth.bed.gz"

# let's first start with ctcf region
reg="$root/data/gm12878/GM12878_CTCF.noTSS.center.2000bp.bed"
shuf -n100 $reg | $scr -v -i $mbed
