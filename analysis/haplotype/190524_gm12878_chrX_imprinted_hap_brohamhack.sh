#!/bin/bash
root=/kyber/Data/Nanopore/projects/nanonome/analysis/data
name=190516_GM12878_nanoNOMe
bam=$root/nanonome/pooled/bam/$name.pooled.bam
subdir=$root/nanonome/subset
regname=chrXandimprint
regs=$root/hg38/hg38_imprinted_genes.slop10kb.bed
bamsub=$subdir/$name.$regname.bam
ref=/mithril/Data/NGS/Reference/hg38_noalt/hg38_noalt.fa
outdir=$root/hap

# merge regs
allreg=$root/hg38/hg38_chrX_and_imprintedgenes_slop10kb.bed
if [ ! -e $allreg ];then
  chrxend=$(grep chrX $ref.fai | awk '{ print $2 }')
  cp $regs $allreg
  echo "chrX,0,$chrxend" | tr "," "\t" >> $allreg
fi

# first subset bam
if [ ! -e $bamsub ];then
  echo subsetting bam
  [ -e $subdir ]||mkdir $subdir
  samtools view -@ 8 -hb $bam -L $allreg > $bamsub
  samtools index $bamsub
fi

# takes too long - let's try just peg10 region
regname=peg10
reg=chr7:94651325-94674695
outpre=$outdir/$name.$regname.longshot
vcf=$outpre.vcf
# longshot
if [ ! -e $vcf ];then
  # max cov : knowing that I have ~ 40x cov, let's do mean + 5 * sqrt(40) = ~70
  echo longshot
  #outpre=$outdir/$name.$regname.longshot
  log=$outpre.log
  if [ "$regname" == "peg10" ];then
    arg="-r $reg"
  else 
    arg=""
  fi
  longshot -C 70 -F -p $outpre $arg \
    --bam $bamsub --ref $ref --out $vcf &> $log
fi

# igv - mcall still going


