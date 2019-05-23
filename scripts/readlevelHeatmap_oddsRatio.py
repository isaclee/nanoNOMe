#!/usr/bin/python
# readlevel analysis of NOMe-seq : plot heatmap of co-methylation
import scipy.stats as ss
import sys
import os
import argparse
import gzip
from collections import namedtuple
from plotnine import *
import bisect
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from methylbed_utils import bed_to_coord,coord_to_bed,MethRead,tabix
blues=["#6BAED6","#4292C6","#2171B5","#08519C","#08306B","#041938"]
import time
start_time = time.time()

def parseArgs():
    parser = argparse.ArgumentParser( description='plot heatmap of methylation co-occurrence')
    parser.add_argument('-v','--verbose',action='store_true',default=False,
            help="verbose output")
    parser.add_argument('-i', '--input', type=os.path.abspath, required=True,
            help="read-level methylation bed file")
    parser.add_argument('-r','--regions',type=argparse.FileType('r'),
            required=False,default=sys.stdin, help="windows in bed format (default stdin)")
    parser.add_argument('-o', '--out', type=str, required=False,
            default="heatmap.pdf",help="output file path (default heatmap.pdf)")
    parser.add_argument('-w','--window',type=int,required=False,default=20,
            help="binning window for calculating methylation")
    parser.add_argument('-f','--frequency',type=float,required=False,default=0.1,
            help="frequency threshold")
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
        if len(self.regfields) > 6 :
            self.title= "{} {} ({})".format(self.regfields[6],
                    self.coord,self.regfields[5])
        elif len(self.regfields) > 4 :
            self.title= "{} {} ({})".format(self.regfields[3],
                    self.coord,self.regfields[5])
        elif len(self.regfields) > 3 : 
            self.title= "{} ({})".format(self.coord,self.regfields[5])
        else :
            self.title= self.coord

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
        dim = len(self.keys)
        c=1e6
        self.ratiomatrix = np.ones((dim,dim))/c # initialize with 1 (add constant to all numbers for log calculation)
        self.methmatrix = np.ones((dim,dim))/c # initialize with 1 (add constant to all numbers for log calculation)
        self.unmethmatrix = np.ones((dim,dim))/c  # initialize with 1 (add constant to all numbers for log calculation)
        self.onemethmatrix = np.ones((dim,dim))/c  # initialize with 1 (add constant to all numbers for log calculation)
        self.twomethmatrix = np.ones((dim,dim))/c  # initialize with 1 (add constant to all numbers for log calculation)
        self.coverage = np.zeros(dim)
        self.covmat = np.zeros((dim,dim))
        self.meth = np.zeros(dim)
        # for each read
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
            # for each pair of co-occurred marks, populate the matrices
            for y in methind :
                yind = methind.index(y)
                for x in methind[0:yind] :
                    self.methmatrix[y,x]+=1
                    self.covmat[y,x]+=1
                yind = bisect.bisect_right(unmethind,y)
                for x in unmethind[0:yind] :
                    self.onemethmatrix[y,x]+=1
                    self.covmat[y,x]+=1
            for y in unmethind :
                yind = unmethind.index(y)
                for x in unmethind[0:yind] :
                    self.unmethmatrix[y,x]+=1
                    self.covmat[y,x]+=1
                yind = bisect.bisect_right(methind,y)
                for x in methind[0:yind] :
                    self.twomethmatrix[y,x]+=1
                    self.covmat[y,x]+=1
        # odds ratio
        for y in range(dim) :
            for x in range(y) :
                count_meth = self.methmatrix[y,x]
                count_unmeth = self.unmethmatrix[y,x]
                count_onemeth = self.onemethmatrix[y,x]
                count_twometh = self.twomethmatrix[y,x]
                ratio = np.log(count_meth*count_unmeth/(count_onemeth*count_twometh))
                self.ratiomatrix[y,x] = ratio
#        print(self.matrix,file=sys.stderr)
        return 
    def plot(self,thr,window) :
        # plot range
        dist_min = min(self.dist)
        dist_max = max(self.dist)
        ymax = (dist_max-dist_min)/2
        # plot heatmap
        mat_flat = np.array(np.ndarray.flatten(self.ratiomatrix))
        xlist = np.array(self.dist*len(self.dist))
        ylist = np.array(np.repeat(self.dist,len(self.dist)))
        mat = np.transpose([xlist,ylist,mat_flat])
        df_all = pd.DataFrame(mat,columns=["x","y","z"])
        # coverage df
        covmat_flat = np.array(np.ndarray.flatten(self.covmat))
        covmat = np.transpose([xlist,ylist,covmat_flat])
        df_cov = pd.DataFrame(covmat,columns=["x","y","z"])
        # filter out non data points
        df_all = df_all.loc[df_cov['z']>5]
        # thresholding
        df_all['newz'] = 0
        df_all['newz'][df_all['z']>10] = 1
        df_all['newz'][df_all['z']<-10] = -1
        # quit if no data
        if df_all.shape[0] == 0 : return
        # flip y axis
        df_all['newy'] = -df_all['y']
        df_plot = df_all.loc[df_all['y']>df_all['x']]
#        # flip scale for fill color
#        df_plot['z'] = -df_plot['z']
        # ggplot
        g = (ggplot(df_plot) +
                geom_tile(aes(x='x',y='y',fill='newz'),size=0.1) +
                lims(x=(dist_min,dist_max),
                    y=(dist_min,dist_max)) +
                scale_fill_distiller(type='div',palette="RdBu", 
                    limits=(-1,1),
                    breaks=(-1,0,1),
                    labels=('(,-10)','[-10,10]','(10,)'), 
#                    name="Co-occurrence") +
                    name="Odds Ratio") +
                labs(x="Distance to center",y=None,title=self.title) +
                theme_bw() +
                theme(panel_grid=element_blank(),
                    axis_text=element_text(color="black"),
                    axis_text_y=element_blank(),
                    axis_ticks_major_y=element_blank(),
                    figure_size=(3,6),
                    aspect_ratio=1)
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
                    title=self.title)+
                theme_bw()+
                theme(panel_grid=element_blank(),
                    axis_text=element_text(color="black"),
                    figure_size=(3,1))
                )
        return (g_avg,g)


def makekey(region):
    return ".".join([str(x) for x in region.items])

def readlevelHeatmap(datapath,reg,thr=0.1,window=10,verbose=False) :
    heat = HeatmapRegion(reg)
    if verbose : print(heat.coord,file=sys.stderr)
    data = tabix(datapath,heat.coord)
    for line in data :
        read = MethRead(line)
        if ( read.start <= heat.start and
                read.end >= heat.end ) :
            heat.addread(read)
    if verbose : 
        print("{} reads covering the entire region out of total {} in the region".format(
            heat.totreads,len(data)),file=sys.stderr)
    heat.makeMatrix(window)
    g = heat.plot(thr,window)
    if verbose : print(heat.title,file=sys.stderr)
    return g

if __name__=="__main__":
    args=parseArgs()
    glist = list()
    for reg in args.regions :
        g = readlevelHeatmap(args.input,reg,args.frequency,args.window,args.verbose)
        if g is not None : 
            glist.append(g)
    if args.verbose : print("generating plots",file=sys.stderr)
    save_as_pdf_pages([g for x in glist for g in x ],filename=args.out)
    if args.verbose : print("time elapsed : {} seconds".format(time.time()-start_time),file=sys.stderr)
