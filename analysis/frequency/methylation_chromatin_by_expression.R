#!/usr/bin/Rscript
library(getopt)
library(optparse)
library(tidyverse)
library(GenomicRanges)
srcdir=dirname(get_Rscript_filename())
if (is.na(srcdir)){ srcdir="." }
if (rstudioapi::isAvailable()){ 
    srcdir = dirname(rstudioapi::getSourceEditorContext()$path)}
source(file.path(srcdir,"../../scripts/methylation_plot_utils.R"))

root = commandArgs(trailingOnly=TRUE)[1]
# default path if not provided
if (is.na(root)){ root="/kyber/Data/Nanopore/projects/nanonome/analysis" }
setwd(root)
# data
cpgpath="data/nanonome/pooled/mfreq/GM12878_nanoNOMe.pooled.cpg.mfreq.txt.gz"
gpcpath="data/nanonome/pooled/mfreq/GM12878_nanoNOMe.pooled.gpc.mfreq.txt.gz"
pd=tibble(type=c("cpg","gpc"),filepath=c(cpgpath,gpcpath))
# annotations
tss.fp="data/hg38/hg38_genes.TSS.200bp.bed"
cgi.fp="data/hg38/hg38_cgi.bed"
# rnaseq database
exp.fp=c("data/gm12878/GM12878_rnaseq_1.tsv",
         "data/gm12878/GM12878_rnaseq_2.tsv")
# output pdf name
plotpath = "plots/GM12878_methylation_chromatin_by_expression.pdf"

#read in tss db
tss.gr=load_db(tss.fp)
tss.gr$id = sapply(strsplit(tss.gr$id,"[.]"),"[[",1)

# expression data parsing
qtiles=seq(0,1,.25)
exp.list=lapply(exp.fp,read_tsv)
exp.list=lapply(seq_along(exp.list),function(i){
    x=exp.list[[i]]
    id=sapply(strsplit(x$gene_id,"[.]"),"[[",1)
    y=tibble(id=id,transcripts=x$"transcript_id(s)",fpkm=x$FPKM_ci_upper_bound)
    y=y[which(y$id %in% tss.gr$id),]
    qs=quantile(y$fpkm,probs=c(0.25,0.5,0.75))
    y%>%mutate(qtile=ifelse(fpkm>=qs[1],2,1),
               qtile=ifelse(fpkm>=qs[2],3,qtile),
               qtile=ifelse(fpkm>=qs[3],4,qtile),
               rep=i)
})
exp=do.call(rbind,exp.list)
# using only genes that are consistent b/w the two replicates, taking average
exp.qtile=exp%>%select(-fpkm)%>%
    spread(rep,qtile)%>%filter(`1`==`2`)%>%mutate(qtile=`1`)
exp.fpkm=exp[which(exp$id %in% exp.qtile$id),] %>%
    group_by(id,transcripts,qtile)%>%summarize(fpkm=mean(fpkm))

# merge with tss.gr
tss.gr=tss.gr[which(tss.gr$id %in% exp.fpkm$id)]
tss.gr$qtile=exp.fpkm$qtile[match(tss.gr$id,exp.fpkm$id)]
tss.gr$fpkm=exp.fpkm$fpkm[match(tss.gr$id,exp.fpkm$id)]

## export high and low expression genes
if (FALSE){
    tss.high = tss.gr[which(tss.gr$qtile==4)]
    tss.low = tss.gr[which(tss.gr$qtile==1)]
    bed.high = GRangesTobed(tss.high)
    bed.low = GRangesTobed(tss.low)
    names.high = bed.high[,"id"]
    names.low = bed.low[,"id"]
    outdir=file.path(root,"data/gm12878")
    highpath = file.path(outdir,"GM12878_rnaseq_highexpression.geneid.txt")
    lowpath = file.path(outdir,"GM12878_rnaseq_lowexpression.geneid.txt")
    write_tsv(names.high,highpath,col_names=F)
    write_tsv(names.low,lowpath,col_names=F)
}

# read in data around TSS
dat.tss=lapply(seq(dim(pd)[1]),function(i){
    tabix_mfreq(pd$filepath[i],tss.fp)})

regmeth.list=lapply(seq_along(dat.tss),function(i){
    getRegionMeth(dat.tss[[i]],tss.gr)%>%
        mutate(calltype=pd$type[i])
})
regmeth=do.call(rbind,regmeth.list)

dat.spread = regmeth%>%mutate(avgcov = totcov/numsites) %>%
    filter( avgcov > 5 & totcov > 10 & numsites > 5 & avgcov < 100)%>%
    select(-totcov,-numsites,-avgcov)%>%
    spread(calltype,freq)%>%na.omit() 
dat.spread = bind_cols(dat.spread,
                       as.data.frame(tss.gr[dat.spread$feature.index]))

plotter <- function(dat.plt,plotpath) {
    # plot
    g = ggplot(dat.plt,aes(y=cpg,x=gpc))+theme_bw()+
        facet_wrap(~qtile)+lims(x=c(0,1),y=c(0,1))+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(color="black"),
              strip.background = element_blank(),
              plot.title = element_text(size=10,hjust=0),
              panel.border = element_rect(colour = "black"))+
        labs(y="CpG methylation",x="GpC accessibility",
             title="Chromatin state by expression")+
        coord_fixed()
    # scatter plot
    g.scatter = g + geom_point(size=0.2,alpha=0.5) +
        geom_rug(size=0.05,alpha=0.5)
    g.2d = g + stat_density_2d(aes(fill=..level..),geom="polygon")
    # heatmap
    hist = bin_pairwise_methylation(dat.plt)
    g.heat = g +
        geom_tile(inherit.aes=F,data=hist,
                  mapping=aes(y=y,x=x,fill=count))+
        scale_fill_distiller(palette="Spectral")+
        theme(legend.position="right")
    # scatter plot in one
    g.one = ggplot(dat.plt,aes(y=cpg,x=gpc))+theme_bw()+
        geom_point(aes(color=factor(qtile)),size=0.5,alpha=0.5)+
        lims(y=c(0,1))+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(color="black"),
              legend.position="right",
              plot.title = element_text(size=10,hjust=0))+
        labs(y="CpG methylation",x="GpC accessibility",
             title="Chromatin state by expression",
             color="quartile")

    # density
    dat.gather = gather(dat.plt,calltype,freq,-qtile)
    g.density = g +
        geom_density(inherit.aes=F,data=dat.gather,
                     mapping=aes(x=freq,y=..scaled..))+
        facet_grid(calltype~qtile) +
        labs(y="Density",x="Methylation frequency",
             title="Methylation density by expression")
    # print
    pdf(plotpath,useDingbats=F,height=5,width=8)
    print(g.one)
    print(g.2d)
    print(g.density)
    print(g.scatter)
    print(g.heat)
#    print(g.density)
    dev.off()
}
dat.plt=dat.spread %>%
    select(c("gpc","cpg","qtile"))

# subsample by min qtile
minnum = min(table(dat.plt$qtile))
dat.plt=dat.plt%>%group_by(qtile)%>%
    do(sample_n(.,minnum))

plotter(dat.plt,plotpath)

# cgi 
# cgi db
cgi.cnames=c("idx","chrom","start","end")
cgi = GRanges(read_tsv(cgi.fp,cgi.cnames))
tss.cgi = tss.gr[overlapsAny(tss.gr,cgi)]
cgi.id = tss.cgi$id

dat.plt=dat.spread[match(tss.cgi$id,tss.gr$id),] %>%
    na.omit()%>%
    select(c("gpc","cpg","qtile"))

# subsample by min qtile
minnum = min(table(dat.plt$qtile))
dat.plt=dat.plt%>%group_by(qtile)%>%
    do(sample_n(.,minnum))

plotter(dat.plt,plotpath)
# non-cgi
dat.plt = dat.spread[-match(tss.cgi$id,tss.gr$id),] %>%
    na.omit() %>%
    select(c("gpc","cpg","qtile"))
# subsample by min qtile
minnum = min(table(dat.plt$qtile))
dat.plt=dat.plt%>%group_by(qtile)%>%
    do(sample_n(.,minnum))

plotter(dat.plt,plotpath)
