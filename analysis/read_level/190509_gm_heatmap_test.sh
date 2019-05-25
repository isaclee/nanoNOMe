#!/bin/bash
srcdir=$(dirname $(readlink -f $0))
root=/kyber/Data/Nanopore/projects/nanonome/analysis
plotdir=$root/plots/heatmap
[ -e $plotdir ]||mkdir $plotdir
methdir=$root/data/nanonome/pooled/mbed
reg=ctcf_test
bed=$root/data/gm12878/GM12878_CTCF.noTSS.center.2000bp.select.bed

for cell in GM12878;do
  echo $cell
  for mod in cpg gpc;do
    prefix=${cell}_nanoNOMe.pooled.$mod
    pres="$pres $prefix"
  done
done
echo $pres

meth=$methdir/{}.meth.bed.gz
plotpath=$plotdir/{}.$reg.heatmap_oddsratio.pdf
log=$plotdir/{}.$reg.heatmap_oddsratio.log
# for testing, use only the first -n entry
com="python -u $srcdir/../../scripts/readlevelHeatmap_oddsRatio.py -v \
  -r $bed -i $meth -o $plotpath 2> $log"
#com="python -u $srcdir/../../script/readlevelHeatmap.py -v -w 20 -f 0.25 \
#  -i $meth -r $bed -o $plotpath 2> $log"
echo $com
parallel "$com" ::: $pres
