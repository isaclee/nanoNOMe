#!/usr/bin/env python3
# readlevel analysis of NOMe-seq : plot heatmap of co-methylation
import sys
import os
import argparse
from collections import namedtuple
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None
from methylbed_utils import bed_to_coord,coord_to_bed,MethRead,tabix
import time
start_time = time.time()

def parseArgs():
    parser = argparse.ArgumentParser( description='plot heatmap of methylation co-occurrence')
    parser.add_argument('-v','--verbose',action='store_true',default=False,
            help="verbose output")
    parser.add_argument('-i', '--input', type=os.path.abspath, required=True,
            help="read-level methylation bed file")
    parser.add_argument('-r','--regions',type=argparse.FileType('r'),
            required=False,default=sys.stdin, help="regions in bed format (default stdin)")
    parser.add_argument('-o', '--out', type=argparse.FileType('w'), required=False,
            default=sys.stdout,help="output file path (default stdout)")
    parser.add_argument('-c','--cov',type=float,required=False,default=10,
            help="coverage threshold")
    args = parser.parse_args()
    try : 
        args.out
    except NameError :
        args.out = "heatmap.pdf"
    return args

def calculate_distances(nlist,start,end) :
    range_idx = np.where((nlist>=start) & (nlist<=end))[0]
    if len(range_idx) == 0 or len(nlist) <= 1 : return 0,list()
    dlist = np.zeros(len(range_idx))
    i = 0
    if range_idx[0] == 0 :
        dlist[0] = nlist[1]-nlist[0]
        range_idx = range_idx[1:]
        i = 1
    if len(range_idx) == 0 : return 1,dlist
    if range_idx[-1] == len(nlist)-1:
        dlist[-1] = nlist[-1]-nlist[-2]
        range_idx = range_idx[:-1]
    if len(range_idx) == 0 : return 1,dlist
    for j in range_idx :
        diff1 = nlist[j]-nlist[j-1]
        diff2 = nlist[j+1]-nlist[j]
        dlist[i] = max(diff1,diff2)
        i += 1
    return 1,dlist

def get_mdistance(datapath,coords,covthr=10,verbose=False) :
    data = tabix(datapath,coords)
    chrom,start,end = coord_to_bed(coords)
    rlist = list()
    for line in data :
        read = MethRead(line)
        if ( read.start <= start and
                read.end >= end ) :
            rlist.append(read)
    if len(rlist) < covthr : 
        if verbose : 
            print("not enough coverage for {}".format(coords),file=sys.stderr)
        return
    if verbose : 
        print("calculating distances in {} reads for {}".format(len(rlist),coords),file=sys.stderr)
    distlist = list()
    numreads = 0
    for read in rlist :
        methind = np.where(read.callarray[:,1]==1)[0]
        methpos = read.callarray[methind,0]
        i,dists = calculate_distances(methpos,start,end)
        numreads += i
        [ distlist.append(x) for x in dists ]
    if numreads < covthr :
        if verbose : 
            print("not enough reads with methylation data for {}".format(coords),file=sys.stderr)
        return
    if verbose : 
        print("calculated distances in {} reads for {}".format(len(rlist),coords),file=sys.stderr)
    return distlist

if __name__=="__main__":
    args=parseArgs()
    glist = list()
    print("distance\tcount\tregion",file=args.out)
    dists_all = list()
    num_region = 0
    for reg in args.regions :
        coords = bed_to_coord(reg)
        distances = get_mdistance(args.input,coords,args.cov,args.verbose)
        if distances is None : continue
        if args.verbose : 
            print("outputting counts of distances for {}".format(coords),file=sys.stderr)
        counts = np.unique(distances,return_counts=True)
        outlist = ["{}\t{}\t{}".format(counts[0][i],counts[1][i],coords)
                for i in range(len(counts[0]))]
        [ print(x,file=args.out) for x in outlist ]
        [dists_all.append(x) for x in distances]
        num_region += 1
    if args.verbose : 
        print("outputting counts for all regions",file=sys.stderr)
    counts = np.unique(dists_all,return_counts=True)
    outlist = ["{}\t{}\t{}".format(counts[0][i],counts[1][i],"all")
            for i in range(len(counts[0]))]
    [ print(x,file=args.out) for x in outlist ]
    if args.verbose : print("time elapsed : {} seconds for {} regions".format(time.time()-start_time,num_region),file=sys.stderr)
