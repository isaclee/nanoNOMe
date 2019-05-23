#!/usr/bin/env python3
# readlevel analysis of NOMe-seq : histogram of distances between accessible marks
import sys
import os
import argparse
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
#    parser.add_argument('-o', '--out', type=str, required=False,
#            default="heatmap.pdf",help="output file path (default heatmap.pdf)")
#    parser.add_argument('-w','--window',type=int,required=False,default=10,
#            help="binning window for calculating methylation")
    parser.add_argument('-c','--cov',type=float,required=False,default=10,
            help="coverage threshold")
    args = parser.parse_args()
    try : 
        args.out
    except NameError :
        args.out = "heatmap.pdf"
    return args

class HeatmapRegion :
    def __init__(self,regline):
        self.regline = regline
        self.regfields = regline.strip().split("\t")
        self.coord = bed_to_coord(regline)
        self.chrom,self.start,self.end = coord_to_bed(self.coord)
        self.id = self.regfields[3]
        if len(self.regfields) > 6 :
            self.name = self.regfields[6]
        else :
            self.name = self.id
        self.center = self.start + np.floor((self.end-self.start)/2)
        self.metharrays=[]
        self.totreads=0
    def binning(self,num,win) :
        return round(num/win)*win
    def initDict(self,win) :
        start = self.binning(self.start,win)
        end = self.binning(self.end,win)
        self.keys = list(range(start,end,win))
        self.dist = [ x - self.binning(self.center,win) for x in self.keys ]
        self.centerkey = self.binning(self.center,win)
    def addread(self,read):
        self.metharrays.append(read.callarray)
        self.totreads+=1
    def makeMatrix(self,win) :
        self.initDict(win)
        self.matrix = np.zeros((len(self.keys),len(self.keys)))
        self.coverage = np.zeros(len(self.keys))
        self.meth = np.zeros(len(self.keys))
        for metharray in self.metharrays:
            methlist = list()
            unmethlist = list()
            for pos,call in metharray :
                key=self.binning(pos,win)
                if call == 1 :
                    methlist.append(key)
                elif call == 0 :
                    unmethlist.append(key)
            methind = list(np.intersect1d(self.keys,methlist,return_indices=True)[1])
            self.coverage[methind] += 1
            self.meth[methind] += 1
            unmethind = list(np.intersect1d(self.keys,unmethlist,return_indices=True)[1])
            self.coverage[unmethind] += 1
            for y in methind :
                yind = methind.index(y)
                for x in methind[0:yind] :
                    self.matrix[y,x]+=1
            for y in unmethind :
                yind = unmethind.index(y)
                for x in unmethind[yind+1:] :
                    self.matrix[y,x]+=1
        # resolve coverage
        cov_idx = np.nonzero(self.coverage)[0]
        for y in cov_idx :
            for x in cov_idx :
                if self.matrix[y,x] == 0 :
                    self.matrix[y,x] = 0.1
        return 
    def plot(self,thr,window) :
        # title
        title="{} {}:{}-{} ({})".format(self.name,
                self.chrom,
                self.start,
                self.end,
                self.regfields[5])
        # plot range
        dist_min = min(self.dist)
        dist_max = max(self.dist)
        ymax = (dist_max-dist_min)/2
        # plot heatmap
        mat_flat = np.array(np.ndarray.flatten(self.matrix))
        xlist = np.array(self.dist*len(self.dist))
        ylist = np.array(np.repeat(self.dist,len(self.dist)))
        mat = np.transpose([xlist,ylist,mat_flat])
        df_all = pd.DataFrame(mat,columns=["x","y","z"])
        # rotate graph
        df_all['newx'] = (df_all['x']+df_all['y'])/2
        df_all['newy'] = (df_all['y']-df_all['x'])/2
        # split meth vs unmeth
        df_meth = df_all.loc[df_all['y']>df_all['x']]
        df_unmeth = df_all.loc[df_all['y']<df_all['x']]
        # filter out non data points
        df_meth = df_meth.loc[df_meth['z']>0]
        df_unmeth = df_unmeth.loc[df_unmeth['z']>0]
        # normalize by max count
        df_meth['z'] = df_meth['z']/np.max(df_meth['z'])
        df_unmeth['z'] = df_unmeth['z']/np.max(df_unmeth['z'])
        # thresholding
        df_meth['z'][df_meth['z']<thr] = 0
        df_unmeth['z'][df_unmeth['z']<thr] = 0
        # flip y axis for unmeth matrix
        df_unmeth['newy'] = -df_unmeth['newy']
        # label and merge
        df_meth['lab'] = 'Methylation'
        df_unmeth['lab'] = 'Unmethylation'
        df = df_meth.append(df_unmeth)
#        df['z'] = 1-df['z'] # flip scale for fill color
        if(len(df) == 0) : return
        # ggplot
        g = (ggplot(df) +
                facet_grid(['lab','.'])+
                geom_tile(aes(x='newx',y='newy',fill='z'),size=0.1) +
                lims(x=(dist_min,dist_max),
                    y=(0,ymax)) +
                scale_fill_distiller(type='seq',palette="YlOrRd",
                    limits=(0,1),
                    breaks=(0,1),
                    labels=('Low','High'),
                    name="Co-occurrence") +
                labs(x="Distance to center",y=None,title=title) +
                theme_bw() +
                theme(panel_grid=element_blank(),
                    axis_text=element_text(color="black"),
                    axis_text_y=element_blank(),
                    axis_ticks_major_y=element_blank(),
                    figure_size=(3,3))
                )
        # plot average
        cov_idx = np.nonzero(self.coverage)[0]
        df = pd.DataFrame(np.transpose([np.array(self.dist)[cov_idx],
            self.meth[cov_idx]/self.coverage[cov_idx]]), 
            columns=["x","y"])
        g_avg = (ggplot(df) +
                geom_line(aes(x='x',y='y'))+
                lims(x=(dist_min,dist_max),y=(0,1))+
                labs(x="Distance to center",y="Meth Freq",
                    title=title)+
                theme_bw()+
                theme(panel_grid=element_blank(),
                    axis_text=element_text(color="black"),
                    figure_size=(3,1))
                )
        return (g_avg,g)


def makekey(region):
    return ".".join([str(x) for x in region.items])

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
    print("distance\tfrequency\tregion",file=args.out)
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
