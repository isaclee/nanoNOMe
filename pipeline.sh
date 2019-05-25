#!/bin/bash 
set -eo pipefail

ROOT=/kyber/Data/Nanopore/projects/nanonome/analysis
STEP="all"
while :
do
  case "$1" in
    -p | --path ) # root path
      ROOT="$2"
      shift 2
      ;;
    -s | --step ) # step to perform
      STEP="$2"
      shift 2
      ;;
    * ) break
      ;;
  esac
done

echo "analysis root directory : $ROOT"
[ -e $ROOT/tmp ]||mkdir $ROOT/tmp
echo "performing $STEP step(s)"
CODEROOT=$(readlink -f $PWD/analysis)

############################################################
#
# Step 00. Download data and install software
# 
############################################################

if [[ $STEP =~ all|step0|hg38Annotation ]];then
  DBDIR=$ROOT/annotations/hg38
  [ -e $DBDIR ]||mkdir -p $DBDIR
  # genes
  genedb=$DBDIR/hg38_genes.gtf.gz
	[ -e $genedb ]||\
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz \
      -O $genedb
  # cgi
  cgidb=$DBDIR/hg38_cgi.txt.gz
	[ -e $cgidb ]||\
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExtUnmasked.txt.gz \
      -O $cgidb
  # repetitive elements
	repdb=$DBDIR/hg38_repeats.txt.gz
	[ -e $repdb ]||\
		wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/nestedRepeats.txt.gz \
			-O $repdb
  # liftover
  chain=$DBDIR/hg19ToHg38.over.chain
	if [ ! -e $chain ];then
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz \
      -O $chain.gz
    gunzip $chain.gz
  fi
  # parsing -
  # genes
  $CODEROOT/annotations/parse_genes.sh $genedb
  # repeats
  $CODEROOT/annotations/parse_repeats.sh $repdb
  # cgi
  gunzip -c $cgidb | awk 'OFS="\t"{ print $2,$3,$4 }' > $DBDIR/hg38_cgi.bed
  # shuffled bed for random sampling
  cp annotations/hg38_shuffle.bed $DBDIR
  
fi

if [[ $STEP =~ all|step0|gmAnnotation ]];then
  DBDIR=$ROOT/annotations/gm12878
  [ -e $DBDIR ]||mkdir -p $DBDIR
  # ctcf
  # download chip data
  CHIP=$DBDIR/GM12878_CTCF_ChIP.bed
  if [ ! -e $CHIP ];then
    wget https://www.encodeproject.org/files/ENCFF356LIU/@@download/ENCFF356LIU.bed.gz \
      -O $CHIP.gz 
    gunzip $CHIP.gz
  fi
  awk 'OFS="\t"{ print $1,$4-1,$5,".",".",$7 }' annotations/CTCF_bindingsites.txt \
    > annotations/CTCF_bindingsites_hg19.bed
  CTCF=annotations/CTCF_bindingsites_hg38.bed
  # first liftover
  liftOver annotations/CTCF_bindingsites_hg19.bed \
    $ROOT/annotations/hg38/hg19ToHg38.over.chain \
    $CTCF $CTCF.unmapped
  # parse to get GM12878
  $CODEROOT/annotations/parse_GM12878_CTCF.sh $DBDIR $CTCF
  #DNAse
  DNASE=$DBDIR/GM12878_DNAse.bed
  if [ ! -e $DNASE ];then
    wget https://www.encodeproject.org/files/ENCFF590NWM/@@download/ENCFF590NWM.bed.gz \
      -O $DNASE.gz 
    gunzip $DNASE.gz
  fi

  $CODEROOT/annotations/parse_GM12878_DNAse.sh $DBDIR $DNASE
  #RNA-seq
  RNA=$DBDIR/GM12878_RNAseq
  [ -e ${RNA}_1.tsv ]||\
    wget https://www.encodeproject.org/files/ENCFF212CQQ/@@download/ENCFF212CQQ.tsv \
    -O ${RNA}_1.tsv
  [ -e ${RNA}_2.tsv ]||\
    wget https://www.encodeproject.org/files/ENCFF350QZU/@@download/ENCFF350QZU.tsv \
    -O ${RNA}_2.tsv
  Rscript $CODEROOT/annotations/parse_GM12878_RNAseq.R $ROOT 2> /dev/null
fi

if [[ $STEP =~ all|step0|gmData ]];then
  DBDIR=$ROOT/data/gm12878
  [ -e $DBDIR ]||mkdir -p $DBDIR
  # mnase
  PRE=$DBDIR/GM12878_MNase
  [ -e $PRE.bigWig ]||\
    wget https://www.encodeproject.org/files/ENCFF000VME/@@download/ENCFF000VME.bigWig \
    -O $PRE.bigWig
  if [ ! -e $PRE.bedGraph.gz ];then
    bigWigToBedGraph $PRE.bigWig $PRE.bed 
    bgzip $PRE.bedGraph 
    tabix -p bed $PRE.bedGraph.gz
  fi
  MNASE=${PRE}_hg38.bedGraph
  if [ ! -e $MNASE.gz ];then
    liftOver $PRE.bed.gz $ROOT/annotations/hg38/hg19ToHg38.over.chain \
      $MNASE $MNASE.unmapped
    sort -k1,1 -k2,2n -T $ROOT/tmp $MNASE.sorted && bgzip $MNASE.sorted
    mv $MNASE.sorted $MNASE.gz
    tabix -p bed $MNASE.gz
  fi
fi

if [[ $STEP =~ all|step0|bcanAnnotation ]];then
  DBDIR=$ROOT/annotations/breastcancer
  [ -e $DBDIR ]||mkdir -p $DBDIR
  echo "bcan genes"
  $CODEROOT/annotations/parse_bcangenes.sh $ROOT
  cp annotations/bcan_10a_vs_231_repeats.bed $ROOT/annotations/breastcancer/
  # enhancer
  echo "enhancers"
  PRE=$DBDIR/MCF10A_enhancer
  [ -e ${PRE}_hg19.txt ]||\
    wget http://www.enhanceratlas.org/data/AllEPs/MCF10A_EP.txt \
    -O ${PRE}_hg19.txt
  $CODEROOT/annotations/parse_bcan_enhancers.sh $DBDIR
  # expression data
  echo "rnaseq data"
  EXP=$DBDIR/bcan_rnaseq.txt
  wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75168/suppl/GSE75168_MCF10A_MCF7_MDA-MB-231_HTSeq_Counts.txt.gz \
    -O $EXP.gz
  gunzip $EXP.gz
  Rscript $CODEROOT/annotations/181031_bcan_rnaseq.R
fi

############################################################
#
# Step 0b. Basecalling, alignment, methylation calling,
#         and converting to methylation bed format
#
############################################################


############################################################
#
# Step 1. methylation model data analysis
#
############################################################

if [[ $STEP =~ all|step1|events ]];then
  echo "event comparison for select kmers"
  Rscript $CODEROOT/model_training/ecoli_event_kmer_comparison.R 2> /dev/null
fi
if [[ $STEP =~ all|step1|roc ]];then
  echo "roc curve"
  # cpg
  python scripts/test_methylmodel.py roc -v \
    --methylated $ROOT/data/na12879_methylation/NA12878_CpG.cpg.meth.subset.tsv \
    --unmethylated $ROOT/data/na12878_methylation/NA12878_unmethylated.cpg.meth.subset.tsv \
    -o $ROOT/plots/model/NA12878_CpG_ROC.pdf -n 100000
  # gpc
  python scripts/test_methylmodel.py roc -v --motif GC \
    --methylated $ROOT/data/na12878_methylation/NA12878_GpC.gpc.meth.subset.tsv \
    --unmethylated $ROOT/data/na12878_methylation/NA12878_unmethylated.gpc.meth.subset.tsv \
    -o $ROOT/plots/model/NA12878_GpC_ROC.pdf -n 100000
fi

############################################################
#
# Step 2. Figure 2 metaplots, expression plots, repeats
#
############################################################

# metaplots
if [[ $STEP =~ all|step2|metaplots ]];then
  echo "metaplots (ctcf,dnase-seq)"
  DBDIR=$ROOT/annotations/gm12878
  OUTDIR=$ROOT/plots/metaplots
  [ -e $OUTDIR ]||mkdir -p $OUTDIR
  # methylation
  CELLS="GM12878 GM12878_BSseq_1"
  DIR=$ROOT/pooled/methylation/mfreq_all
  for CELL in $CELLS ; do
    cpg=$DIR/$CELL.cpg.methfreq.txt.gz
    gpc=$DIR/$CELL.gpc.methfreq.txt.gz
    # CTCF
    Rscript script/nanonome_plots.R metaplotByDistance -c $cpg -g $gpc \
      -r $DBDIR/GM12878_CTCF_hg38.center.2000bp.bed \
      -o $OUTDIR/$CELL.CTCF.metaplot.pdf 2> /dev/null
    # DNAse
    Rscript script/nanonome_plots.R metaplotByDistance -c $cpg -g $gpc \
      -r $DBDIR/GM12878_DNAse_top.center.2000bp.bed \
      -o $OUTDIR/$CELL.DNAse.metaplot.pdf 2> /dev/null
  done
  # mnase
  analysis/frequency/180921_mnase_metaplot.sh $ROOT
fi

# expression plots
if [[ $STEP =~ all|step2|gm12878exp ]];then
  echo "GM12878 tss chromatin by expression plots"
  Rscript $CODEROOT/frequency/181015_gm12878_tss_chromatin_by_exp.R $ROOT 2> /dev/null
fi
# coverage
if [[ $STEP =~ all|step2|gm12878coverage ]];then
  echo "GM12878 BS-seq vs nanoNOMe coverage comparison"
  # first global plots
  Rscript $CODEROOT/frequency/180707_bsseqVSnano_coveragecomp.R $ROOT 2> /dev/null
fi
# specific regions, read-level
if [[ $STEP =~ all|step2|gm12878readlevel ]];then
  echo "GM12878 BS-seq vs nanoNOMe read-level comparison"
  $CODEROOT/igv/gm12878_igv_plotting.sh $ROOT
fi

############################################################
#
# Step 3. Figure 3 Breast cancer comparison - 
#         promoter region comparison,
#         enhancer promter chromatin state correlation (heatmap?)
#
############################################################

# bcan promoters, read-level
if [[ $STEP =~ all|step3|bcangenes ]];then
  echo "bresat cancer sample gene promoters comparisons"
#  $CODEROOT/frequency/180910_bcan_promoter.R $ROOT 2> /dev/null
#  $CODEROOT/igv/bcan_igv_plotting.sh $ROOT bcanpromoter
  echo "co-methylation heatmap"
  $CODEROOT/read_level/181024_bcan_heatmap.sh $ROOT rna
fi
# bcan enhancer-promoter - do all transcripts instead of just bcan genes


############################################################
#
# Step 4. Figure 4 SV analysis - 
#         exclusive SV falling in bcan gene with 
#         interesting chromatin pattern
#
############################################################
if [[ $STEP =~ all|step4|breastcancer_expression ]];then
  echo "bresat cancer sample expression comparisons"
#  Rscript $CODEROOT/bcan/181031_bcan_rnaseq.R $ROOT 3> /dev/null
  echo "igv"
#  $CODEROOT/bcan/bcan_igv_plotting.sh $ROOT MCF7
#  $CODEROOT/bcan/bcan_igv_plotting.sh $ROOT MDAMB231
  echo "co-occurrence"
  $CODEROOT/bcan/bcan_heatmap.sh $ROOT MCF7
  $CODEROOT/bcan/bcan_heatmap.sh $ROOT MDAMB231

fi

if [[ $STEP =~ all|step4|sv ]];then
  echo "bresat cancer sample sv comparisons"
  $CODEROOT/sv/181020_bcanSV.sh
  echo "igv"
  $CODEROOT/igv/bcan_igv_plotting.sh $ROOT bcansv
fi
