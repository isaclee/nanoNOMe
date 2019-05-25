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


# compare to hapcut2
old_dir=/dilithium/Data/Nanopore/projects/nomeseq/analysis/pooled/hap/hapcut_ase
old_mat=$old_dir/peg10.maternal.readID.txt
old_pat=$old_dir/peg10.paternal.readID.txt

mat_num=$(wc -l $old_mat | awk '{ print $1 }')
pat_num=$(wc -l $old_pat | awk '{ print $1 }')
# compare to annotation 
variants=/dilithium/Data/Nanopore/projects/nomeseq/analysis/database/gm12878/variants/NA12878.vcf.gz

# test ranges of combinations of -Q, -E, and -y
regnames="PEG10 TP53 H19 C9orf85 IGF2R GAB1"
q_range="1 5 10 20"
e_range="0.1 0.2 0.25 0.3"
y_range="1 5 10"
vcfreport=$outdir/$name.longshot_test.txt
phasereport=$outdir/$name.phasereport.txt
# longshot
if [ "$1" == "phase" ];then
  # max cov : knowing that I have ~ 40x cov, let's do mean + 5 * sqrt(40) = ~70
  echo "region,Q,E,y,annotation,longshot,overlap,sensitivity(%),precision(%)" > $vcfreport
  echo "region,Q,E,y,totreads,hap1,hap2,hap1mat(peg10),hap2mat(peg10),hap1pat(peg10),hap2pat(peg10)" > $phasereport
  for regname in $regnames; do
    echo "reg : $regname"
    # get reg
    reg=$(grep $regname $regs | awk '{ print $1":"$2+1"-"$3 }')
    # test ranges of combinations of -Q, -E, and -y
    for q in $q_range; do
      echo "Q : $q"
      for e in $e_range; do
        echo "E : $e"
        for y in $y_range; do
          echo "y : $y"
          outpre=$outdir/$name.$regname.Q$q.E$e.y$y.longshot
          log=$outpre.log
          vcf=$outpre.vcf
          arg="-r $reg -Q $q -E $e -y $y"
          com="longshot -C 70 -F -p $outpre $arg \
            --bam $bamsub --ref $ref --out $vcf &> $log"
          eval $com
          regvar=$outdir/NA12878.variants.$regname.vcf
          tabix $variants $reg > $regvar
          anno=$(cat $regvar | wc -l)
          detect=$(grep -v "##" $vcf | wc -l)
          overlap=$(awk 'NR==FNR{c[$2]++;next};c[$2]' $regvar $vcf | wc -l)
          echo "$regname,$q,$e,$y,$anno,$detect,$overlap,$(($overlap*100/$anno)),$(($overlap*100/$detect))" >> $vcfreport
          hap1fp=$outpre.hap1.readID.txt
          hap2fp=$outpre.hap2.readID.txt
          totnum=$(samtools view $bamsub $reg | cut -f1 | uniq | wc -l)
          samtools view $outpre.hap1.bam | cut -f1 | uniq > $hap1fp
          samtools view $outpre.hap2.bam | cut -f1 | uniq > $hap2fp
          hap1num=$(wc -l $hap1fp | awk '{ print $1 }')
          hap2num=$(wc -l $hap2fp | awk '{ print $1 }')
          if [ "$regname" == "PEG10" ];then
            hap1mat=$(grep -f $hap1fp $old_mat | wc -l)
            hap2mat=$(grep -f $hap2fp $old_mat | wc -l)
            hap1pat=$(grep -f $hap1fp $old_pat | wc -l)
            hap2pat=$(grep -f $hap2fp $old_pat | wc -l)
            echo "$regname,$q,$e,$y,$totnum,$hap1num,$hap2num,$hap1mat,$hap2mat,$hap1pat,$hap2pat" >> $phasereport
          else
            echo "$regname,$q,$e,$y,$totnum,$hap2num,$hap2num" >> $phasereport
          fi
        done
      done
    done
  done
fi

## conclusion : using Q = 5, E = 0.25, y = 5
