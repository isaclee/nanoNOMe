#!/usr/bin/snakemake --snakefile
import pandas as pd
from snakemake.utils import validate

rule mtsv_to_mbed:
	input:
		"{dir}/mcall/{sample}.{mod}.meth.tsv.gz"
	params:
		config['utildir']
	output:
		"{dir}/mbed/{sample}.{mod}.meth.bed.gz"
	shell:
		"gunzip -c {input} | "
		"{params}/mtsv2bedGraph.py -m {wildcards.mod} --nome | "
		"sort -T {wildcards.dir}/mbed -k1,1 -k2,2n | bgzip "
		"> {output} && "
		"tabix -p bed {output}"

rule mbed_to_mfreq:
	input:
		"{dir}/mbed/{sample}.{mod}.meth.bed.gz"
	params:
		config['codedir']
	output:
		mfreq="{dir}/mfreq/{sample}.{mod}.mfreq.txt.gz",
		log="{dir}/mfreq/{sample}.{mod}.mfreq.log"
	shell:
		"python -u {params}/scripts/parseMethylbed.py frequency -v "
		"-i {input} -m {wildcards.mod} 2> {output.log} | "
		"bgzip > {output.mfreq} && "
		"tabix -b 2 -e 2 {output.mfreq}"

rule mfreq_to_wig:
	input:
		"{dir}/mfreq/{sample}.{mod}.mfreq.txt.gz"
	params:
		config['codedir']
	output:
		methwig=temp("{dir}/bigwig/{sample}.{mod}.methylation.wig"),
		covwig=temp("{dir}/bigwig/{sample}.{mod}.methcoverage.wig"),
		log="{dir}/bigwig/{sample}.{mod}.wig.log"
	shell:
		"python {params}/scripts/makeWig.py -v -i {input} "
		"-o {output.methwig} -c {output.covwig} &> {output.log}"

rule wig_to_bigwig:
	input:
		"{dir}/bigwig/{sample}.{mod}.{type}.wig",
		"data/hg38/hg38_genomesize.txt"
	output:
		"{dir}/bigwig/{sample}.{mod}.{type}.bw"
	shell:
		"wigToBigWig {input} {output}"
