#!/usr/bin/env python3
# readlevel analysis of NOMe-seq : plot heatmap of co-methylation
import sys
import os
import argparse
from collections import namedtuple
import pysam
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None
from methylbed_utils import bed_to_coord,coord_to_bed,MethRead,tabix
from itertools import combinations_with_replacement
from bisect import bisect
import re
import time
start_time = time.time()

comb = list(combinations_with_replacement([0,1,2,3],2))
states = ["CG:M","CG:U","GC:M","GC:U"]
comb_states = [ states[i]+"_"+states[j] for (i,j) in comb ]

def parseArgs():
    parser = argparse.ArgumentParser( description='plot heatmap of methylation co-occurrence')
    parser.add_argument('-v','--verbose',action='store_true',default=False,
            help="verbose output")
    parser.add_argument('-c', '--cpg', type=os.path.abspath, required=True,
            help="CpG read-level methylation bed file")
    parser.add_argument('-g', '--gpc', type=os.path.abspath, required=True,
            help="GpC read-level methylation (accessibility) bed file")
    parser.add_argument('-f','--fasta',type=os.path.abspath,required=False,
            help="genome fasta file - if supplied, will output distance in genomic context")
    parser.add_argument('-r','--regions',type=argparse.FileType('r'),
            required=False,default=sys.stdin, help="regions in bed format (default stdin)")
    parser.add_argument('-o', '--out', type=argparse.FileType('w'), required=False,
            default=sys.stdout,help="output file path (default stdout)")
    args = parser.parse_args()
    return args

def make_tabix_dict(tabixlist) :
    outdict = dict()
    for line in tabixlist :
        read = MethRead(line)
        outdict[read.qname] = read
    return outdict

def get_methcontext_positions(read,methstate,start,end) :
    context_ind = np.where(read.callarray[:,1]==methstate)[0]
    context_pos = read.callarray[context_ind,0]
    pos_range = context_pos[np.where((context_pos>=start) & (context_pos<=end))[0]]
    return pos_range

def calculate_distances(pos_one,pos_two) :
    if len(pos_one) == 0 or len(pos_two) == 0 : return list()
    dist_list = list()
    pos_list = [pos_one,pos_two]
    list_idx = np.argsort([len(x) for x in pos_list])
    reference = pos_list[list_idx[1]]
    for x in pos_list[list_idx[0]] :
        reference_exclude = reference[reference != x ] # not counting distance to itself
        idx = bisect(reference_exclude,x)
        if len(reference_exclude) == 0 : # if itself is the only entry, move on
            continue
        elif idx == 0 : # if x is less than the all of reference, just take diff from lowest ref number
            min_dist = reference_exclude[0] - x
        elif idx == len(reference_exclude) : # if x is larger
            min_dist = x - reference_exclude[-1]
        else : 
            dists = [ x - reference_exclude[idx-1], reference_exclude[idx] - x ]
            min_dist = np.max(dists) # maximum of two distances
        dist_list.append(min_dist)
    return dist_list

def get_mdistance(cpgpath,gpcpath,coords,verbose=False) :
    # load data
    cpg_data = tabix(cpgpath,coords)
    gpc_data = tabix(gpcpath,coords)
    chrom,start,end = coord_to_bed(coords) 
    # make dict from data
    cpg_dict = make_tabix_dict(cpg_data)
    gpc_dict = make_tabix_dict(gpc_data)
    # filter out reads that don't span the entire region within 10bp
    keys_include = list()
    for qname in gpc_dict:
        read = gpc_dict[qname]
        if ( read.start-10 <= start and
                read.end+10 >= end ) :
            keys_include.append(qname)
    # filter out reads that don't have both cpg and gpc record
    for key in set(keys_include) : 
        if key not in cpg_dict.keys() :
            keys_include.pop(key)
    if len(keys_include) == 0 :  # no reads
        if verbose : 
            print("no reads in {}".format(coords),file=sys.stderr)
        return
    if verbose : 
        print("calculating distances in {} reads for {}".format(
            len(keys_include),coords),file=sys.stderr)
    dist_all = np.empty((len(comb),0)).tolist()
    distlist = list()
    context_list = list()
    for key in keys_include:
        cpg_methpos = get_methcontext_positions(cpg_dict[key],1,start,end)
        cpg_unmethpos = get_methcontext_positions(cpg_dict[key],0,start,end)
        gpc_methpos = get_methcontext_positions(gpc_dict[key],1,start,end)
        gpc_unmethpos = get_methcontext_positions(gpc_dict[key],0,start,end)
        methpos_list = [cpg_methpos,cpg_unmethpos,gpc_methpos,gpc_unmethpos]
        distance_list = [ calculate_distances(methpos_list[i],methpos_list[j]) 
                for (i,j) in comb ]
        dist_all = [ dist_all[i]+distance_list[i] for i in range(len(comb)) ]
    return dist_all

if __name__=="__main__":
    args=parseArgs()
    if args.fasta : 
        comb_genome = list(combinations_with_replacement([0,1],2))
        states_genome = ["CG","GC"]
        genome_states = [ states_genome[i]+"_"+states_genome[j] for (i,j) in comb_genome ]
        fasta = pysam.FastaFile(args.fasta)
        genome_all = np.empty((len(comb_genome),0)).tolist()
    glist = list()
    print("distance\tcount\trelationship\tregion",file=args.out)
    dists_all = np.empty((len(comb),0)).tolist()
    num_region = 0
    for reg in args.regions :
        coords = bed_to_coord(reg)
        distances = get_mdistance(args.cpg,args.gpc,coords,args.verbose)
        if args.verbose : 
            print("outputting counts of distances for {}".format(coords),file=sys.stderr)
        counts_list = [ np.unique(x,return_counts=True) for x in distances ]
        for i in range(len(comb)) :
            counts = counts_list[i]
            out = [ "{}\t{}\t{}\t{}".format(
                counts[0][j],counts[1][j],comb_states[i],coords)
                for j in range(len(counts[0]))]
            [ print(x,file=args.out) for x in out ]
            dists_all[i] = dists_all[i] + distances[i]
        num_region += 1
        # distance in the genome - CG vs CG, CG vs GC, GC vs GC
        if args.fasta : 
            chrom,start,end = coord_to_bed(coords) 
            seq = fasta.fetch(reference=chrom,start=start,end=end).upper()
            cpg_pos =  np.array([m.start() for m in re.finditer("CG",seq) ])
            gpc_pos =  np.array([m.start() for m in re.finditer("GC",seq) ])
            genome_pos = [cpg_pos,gpc_pos]
            genome_distlist = [calculate_distances(genome_pos[i],genome_pos[j])
                    for (i,j) in comb_genome]
            for i in range(len(comb_genome)) :
                counts = np.unique(genome_distlist[i],return_counts=True)
                out = [ "{}\t{}\t{}\t{}".format(
                    counts[0][j],counts[1][j],"Genome - "+genome_states[i],coords)
                    for j in range(len(counts[0]))]
                [ print(x,file=args.out) for x in out ]
                genome_all[i] = genome_all[i] + genome_distlist[i]
    if args.verbose : 
        print("outputting counts for all regions",file=sys.stderr)
    counts_list = [ np.unique(x,return_counts=True) for x in dists_all ]
    for i in range(len(comb)) :
        counts = counts_list[i]
        out = [ "{}\t{}\t{}\t{}".format(
            counts[0][j],counts[1][j],comb_states[i],"all")
            for j in range(len(counts[0]))]
        [ print(x,file=args.out) for x in out ]
    if args.fasta :
        for i in range(len(comb_genome)) :
            counts = np.unique(genome_all[i],return_counts=True)
            out = [ "{}\t{}\t{}\t{}".format(
                counts[0][j],counts[1][j],"Genome - "+genome_states[i],"all")
                for j in range(len(counts[0]))]
            [ print(x,file=args.out) for x in out ]

        fasta.close()
    if args.verbose : print("time elapsed : {} seconds for {} regions".format(time.time()-start_time,num_region),file=sys.stderr)
