configfile: 
  "snakemake_config.yml"
workdir:
  config['workdir']
smkdir = config['codedir']+"/snakemake/"
include: 
  smkdir+"downloaded_data_parse.smk"
include:
  smkdir+"model_validation.smk"
include: 
  smkdir+"frequency_analysis.smk"
include:
  smkdir+"nanonome_data_parse.smk"
include:
  smkdir+"sv_analysis.smk"

samples_tb = pd.read_csv(config['codedir']+"/sample_info.csv")

rule parse_downloaded_data:
  input:
    "data/gm12878/GM12878_CTCF.noTSS.center.2000bp.bed",
    "data/hg38/hg38_genes.TSS.slop2000bp.bed",
    "data/hg38/hg38_cgi.bed"

rule model_validation:
  input:
    "plots/NA12878_ROC_CpG_methylation_cpg_model.pdf",
    "plots/NA12878_ROC_GpC_methylation_gpc_model.pdf"

rule frequency_analysis:
  input:
    "plots/GM12878_nanoNOMe.pooled.metaplot.CTCF.noTSS.center.2000bp.pdf",
    "plots/GM12878_BSseq.metaplot.CTCF.noTSS.center.2000bp.pdf",
    "plots/GM12878_MNase.metaplot.CTCF.noTSS.center.2000bp.pdf"

rule nanonome_data:
  input:
    expand("analysis/mbed/{sample}.{mod}.meth.bed.gz",
      sample=samples_tb['sample'],
        mod=["cpg","gpc"])

rule sv_analysis:
  input:
    expand("data/bcan/MCF10A_vs_{sample}_SVcomparison_{svtype}_survivor.flank200bp.bed",
      sample=["MCF7","MDAMB231"],
      svtype=["del","ins"]),
    expand("data/bcan/MCF10A_vs_{sample}_SVcomparison_{svtype}_survivor.shuffle.flank200bp.bed",
      sample=["MCF7","MDAMB231"],
      svtype=["del","ins"]),
    expand("data/bcan/MCF10A_vs_{sample}_SVcomparison_tra_survivor.TRAregion200bp.bed",
      sample=["MCF7","MDAMB231"]),
    expand("data/bcan/MCF10A_vs_{sample}_SVcomparison_tra_survivor.shuffle.TRAregion200bp.bed",
      sample=["MCF7","MDAMB231"])


    
