#!/usr/bin/env Rscript
library(getopt)
srcdir=dirname(get_Rscript_filename())
if (is.na(srcdir)){ srcdir="." } # if executing within emacs
if (isRStudio <- Sys.getenv("RSTUDIO") == "1"){ 
srcdir=dirname(rstudioapi::getSourceEditorContext()$path)} # if in rstudio
source(file.path(srcdir,"../../scripts/ggplot_theme.R"))
source(file.path(srcdir,"../../scripts/methylation_plot_utils.R"))
library(BSgenome.Hsapiens.UCSC.hg38)
library(parallel)
cores = detectCores()-4

## data
root="/kyber/Data/Nanopore/projects/nanonome/analysis"
subdir=file.path(root,"data/subset/repeats/mfreq")
outdir=file.path(root,"plots")

fps = system(paste0("find ",subdir," -name \"*mfreq.txt.gz\""),intern=T)
bases = basename(fps)
cells = sapply(strsplit(bases,"_"),"[[",1)
assays = sapply(strsplit(sapply(strsplit(bases,"_"),"[[",2),"[.]"),"[[",1)
dbtypes = sapply(strsplit(bases,"[.]"),tail,4)[1,]
calltypes = sapply(strsplit(bases,"[.]"),tail,5)[1,]
pd = tibble(cell = cells, assay = assays, dbtype=dbtypes,
          calltype = calltypes, fp = fps)

### looking at cpg only for now
#pd = pd[which(pd$calltype=="cpg"),]

## region
regfp = file.path(root,"data/hg38/hg38_repeats.bed")
regs = load_db(regfp,extracols="type")
## remove alternative chroms
chroms = levels(seqnames(regs))
chroms.noalt = chroms[-grep("_",chroms)]
regs = regs[which(seqnames(regs) %in% chroms.noalt)]

rtypes = unique(dbtypes)
rtype = "DNA"
sum.list = list()
## doing one region at a time to save ram
for (rtype in rtypes){
cat(paste0(rtype,"\n"))
regs.rtype = regs[which(regs$type==rtype)]
pd.rtype = pd[pd$dbtype==rtype,]
## reading in data in regions
dat.list = mclapply(mc.cores=cores,seq_along(pd.rtype$fp),function(i){
tabix_mfreq(pd.rtype$fp[i]) %>%
  mutate(cell=pd.rtype$cell[i],
         dbtype=pd.rtype$dbtype[i],
         calltype=pd.rtype$calltype[i],
         assay=pd.rtype$assay[i])
})
dat.all = do.call(rbind,dat.list) %>%
  filter(cov > 2) # at least 3 reads
  
rm(dat.list);gc()

## true number of cpgs per region
regs.fix = regs.rtype
strand(regs.fix)[which(strand(regs.fix)=="*")] = "+" # fix unstranded entries
start(regs.fix) = start(regs.fix)+1 # from bed to coord
regs.cpg = regs.fix
regs.gpc = regs.fix
end(regs.cpg) = end(regs.cpg)+1
start(regs.gpc) = start(regs.gpc)-1
genome = BSgenome.Hsapiens.UCSC.hg38
#  start(regs.cpg) = start(regs.cpg)-1
seqs.cpg = getSeq(genome,regs.cpg)
seqs.gpc = getSeq(genome,regs.gpc)
cpgcounts = vcountPattern("CG",seqs.cpg)
gpccounts = vcountPattern("GC",seqs.gpc)

## granges overlap
dat.gr = GRanges(dat.all)
ovl = findOverlaps(dat.gr,regs.fix,ignore.strand=T,select="all")
dat.ovl = dat.all[queryHits(ovl),]
dat.ovl$regidx = subjectHits(ovl)
rm(dat.all);gc()


## summarize per region
dat.sum = dat.ovl %>%
  group_by(cell,calltype,dbtype,assay,regidx) %>%
  summarize(n = n(),
            cov = sum(cov),
            avgcov = cov/n,
            meth = mean(meth)/cov) %>%
  filter(n>2) # at least 3 sites
## get the number of cpg loci that were measured out of all available per region
cpgidx = which(dat.sum$calltype == "cpg")
gpcidx = which(dat.sum$calltype == "gpc")
dat.sum$nn = dat.sum$n#/cpgcounts[dat.sum$regidx]
dat.sum$nn[cpgidx] = dat.sum$n[cpgidx]/cpgcounts[dat.sum[cpgidx,]$regidx]
dat.sum$nn[gpcidx] = dat.sum$n[gpcidx]/gpccounts[dat.sum[gpcidx,]$regidx]

summary(dat.sum$nn[cpgidx])
summary(dat.sum$nn[gpcidx])
summary(dat.sum$nn)
dat.sum[which(dat.sum$nn==1.5),]
dat.ovl[dat.ovl$regidx==36463,]
regs.fix[36463]



sum.list[[rtype]] = dat.sum
rm(dat.ovl);gc()
}
lapply(sum.list,function(x){
summary(x$nn)
})
dat.sum = do.call(rbind,sum.list)
## get some numbers
bins = seq(0,1,0.1)
dat.spread = dat.sum %>% ungroup() %>%
  select(cell,calltype,dbtype,assay,regidx,nn) %>%
  spread(assay,nn) %>%
  replace_na(list(BSseq=0,nanoNOMe=0)) %>%
  mutate(del = nanoNOMe-BSseq,
         category = cut(BSseq,bins,include.lowest=T))

dat.correct = dat.spread %>%
  gather(assay,n,-cell,-calltype,-dbtype,-regidx,-del,-category) %>%
  ungroup()

dat.nome = dat.correct[which(dat.correct$assay=="nanoNOMe"),]

dat.spread[which(dat.spread$calltype=="gpc"),]
dat.sum[which(dat.sum$regidx==8),]

dat.correct %>% 
  group_by(cell,calltype,dbtype,assay) %>%
  filter(n>0) %>%
  summarize(n())

dat.nome %>%
  group_by(cell,calltype,dbtype,assay,category) %>%
  summarize(mean(n))
  dat.spread[which(dat.spread$BSseq<0.1),] %>%
  group_by(calltype,dbtype)%>%
  summarize(n = mean(nanoNOMe))

## plotting
plotpath = file.path(outdir,"gm12878_repeat_covcomparison.pdf")
pdf(plotpath,width=4,height=3)
## box plot
ggplot(dat.correct,aes(x=dbtype,y=n,color=assay)) +
facet_wrap(~calltype)+
geom_boxplot() 
## violin
ggplot(dat.correct,aes(x=dbtype,y=n,color=assay)) +
  facet_wrap(.~calltype)+
  geom_violin() 
## scatter
ggplot(dat.spread,aes(x=BSseq,y=nanoNOMe))+
  facet_wrap(.~calltype) +
  geom_point(alpha=0.3) + coord_fixed() 
## density
ggplot(dat.correct,aes(x=n))+
  facet_wrap(assay~calltype,scale="free",nrow=2) +
  geom_density()
## 2d density
ggplot(dat.spread,aes(x=BSseq,y=nanoNOMe)) +
  facet_wrap(dbtype~calltype) +
  stat_density_2d(geom = "raster", aes(fill = stat(density)), contour = FALSE) + coord_fixed() 
  
dev.off()
