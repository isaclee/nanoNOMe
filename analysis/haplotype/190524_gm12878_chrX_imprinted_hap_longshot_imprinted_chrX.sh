#!/bin/bash
root=/kyber/Data/Nanopore/projects/nanonome/analysis/data
name=190516_GM12878_nanoNOMe
bam=$root/nanonome/pooled/bam/$name.pooled.bam
subdir=$root/nanonome/subset
regname=chrXandimprint
bamsub=$subdir/$name.$regname.bam
ref=/mithril/Data/NGS/Reference/hg38_noalt/hg38_noalt.fa
outdir=$root/hap/$regname
[ -e $outdir ]||mkdir $outdir

# merge regs
allreg=$root/hg38/hg38_chrX_and_imprintedgenes_slop10kb.bed
if [ ! -e $allreg ];then
  regs=$root/hg38/hg38_imprinted_genes.slop10kb.bed
  chrxend=$(grep chrX $ref.fai | awk '{ print $2 }')
  # do bedtools merge to remove overlapping regions
  echo "chrX,0,$chrxend,chrX,.,.,." |\
    tr "," "\t" | cat - $regs |\
    sort -k1,1 -k2,2n |\
    bedtools merge -i stdin > $allreg
fi

# first subset bam
if [ ! -e $bamsub ];then
  echo subsetting bam
  [ -e $subdir ]||mkdir $subdir
  samtools view -@ 8 -hb $bam -L $allreg > $bamsub
  samtools index $bamsub
fi

# longshot
# max cov : knowing that I have ~ 40x cov, let's do mean + 5 * sqrt(40) = ~70
# using Q = 5, E = 0.25, y = 5
if [ "$1" == "phase" ];then
  # since they recommend separating by chrom, let's just parallelize by entry in bed file (distinct non-overlapping regions)
  regs=$(awk '{ print $1":"$2+1"-"$3 }' $allreg | tr "\n" " ")
  reg="{}"
  outpre=$outdir/$name.$regname.longshot.$reg
  log=$outpre.log
  vcf=$outpre.vcf
  arg="-r $reg -Q 5 -E 0.25 -y 5"
  com="longshot -C 70 -F -p $outpre $arg \
  --bam $bamsub --ref $ref --out $vcf &> $log"
  parallel $com ::: $regs
fi

