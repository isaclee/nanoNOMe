#! /usr/bin/env python3
# much of this is taken from the following script with adjustments for current version of nanopolish
# https://github.com/jts/methylation-analysis/blob/master/calculate_call_accuracy.py
import os
import sys
import argparse
import gzip
import numpy as np
import pandas as pd
import random
from plotnine import *
from collections import namedtuple

def parseArgs() :
    # dir of source code
    srcpath=sys.argv[0]
    srcdir=os.path.dirname(os.path.abspath(srcpath))
    # parser
    parser = argparse.ArgumentParser(description='test models using testing set')
    subparsers = parser.add_subparsers()
    # parent parser
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('-o','--out',type=os.path.abspath,required=False,
            help="output file path for plots")
    parent_parser.add_argument('-v','--verbose',action='store_true',default=False,
            help="verbose output")
    parent_parser.add_argument('--methylated',type=os.path.abspath,required=True,
            help="methylated data (gzipped)")
    parent_parser.add_argument('--unmethylated',type=os.path.abspath,required=True,
            help="unmethylated data (gzipped)")
    parent_parser.add_argument('--model',type=str,required=False,
            default="cpg",help="methylation model (cpg,gpc,etc, default cpg)")
    # parser for frequency
    parser_roc= subparsers.add_parser('roc',parents=[parent_parser],
            help = 'plot roc curve')
    parser_roc.add_argument('-n','--num',type=int,required=False,default=100000,
            help="number of data points")
    parser_roc.add_argument('--table',action='store_true',default=False,
            help="output a table of log likelihoods instead of plots")
    parser_roc.add_argument('--roctable',action='store_true',default=False,
            help="output a table of ROC instead of plots")
    parser_roc.add_argument('--thrtable',action='store_true',default=False,
            help="output a table of thresholds instead of plots")
    parser_roc.set_defaults(func=getROC)
    # parse args
    args = parser.parse_args()
    args.srcdir=srcdir
    return args

class MethylQuery :
    def __init__(self,fields) :
        self.fields = fields
        self.chrom = fields[0]
        self.ratio = float(fields[5])
        self.numsites = int(fields[9])
        self.sequence = fields[10]

def filter_query(fp,motif,num) :
    sites = list()
    with gzip.open(fp,'rt') as fh :
        for line in fh :
            if "chrom" in line :
                continue
            query = MethylQuery(line.strip().split("\t"))
            if query.numsites != 1 :
                continue
            assert(query.sequence.count(motif) == 1)
            sites.append(query.ratio)
        random.shuffle(sites)
    return sites[0:num]

def getROC(args) :
    if args.verbose : print("generating ROC curves to {}".format(args.out),file=sys.stderr)
    # reading in the data
    if args.model == "cpg" :
        motif = "CG"
    elif args.model == "gpc" :
        motif = "GC"
    if args.verbose : print("reading in unmethylated data",file=sys.stderr)
    unmeth_sites = filter_query(args.unmethylated,motif,args.num)
    if args.verbose : print("reading in methylated data",file=sys.stderr)
    meth_sites = filter_query(args.methylated,motif,args.num)
    if args.table :
        vals = unmeth_sites + meth_sites
        labs = len(unmeth_sites)*["unmeth"] + len(meth_sites)*["meth"]
        df = pd.DataFrame({'values':vals,'truth':labs})
        df.to_csv(sys.stdout,index=False)
        return
    # thresholds
    meth_array = np.array(meth_sites)
    unmeth_array = np.array(unmeth_sites)
    # ROC
    increments = 0.5
    thresholds = np.arange(-20,20,increments)
    tp = list()
    fp = list()
    for t in thresholds :
        if args.verbose : print("t :"+str(t),file=sys.stderr)
        tp.append(sum(meth_array>t))
        fp.append(sum(unmeth_array>t))
    roc_array = np.column_stack((thresholds,tp,fp,
        [len(unmeth_array)-x for x in fp],
        [len(meth_array)-x for x in tp],
        [float(x)/len(meth_array) for x in tp],
        [float(x)/len(unmeth_array) for x in fp]))
    roc_df = pd.DataFrame(roc_array,
            columns=['thr','tp','fp','tn','fn','recall','fallout']) 
    if args.roctable :
        roc_df.to_csv(sys.stdout,index=False)
        return
    # get auc
    sortidx = np.argsort(roc_df['fallout'])
    fallout_sorted = np.append(np.array(roc_df['fallout'])[sortidx],1)
    recall_sorted = np.array(roc_df['recall'])[sortidx]
    recall_sorted = np.append(recall_sorted,recall_sorted[-1]) # x=1 value for extending roc to x=1
    fallout_del = np.diff(fallout_sorted)
    auc_list = recall_sorted[1:] * fallout_del
    auc = np.sum(auc_list)
    # call thr
    meth_array_abs = abs(meth_array)
    unmeth_array_abs = abs(unmeth_array)
    all_array_abs = np.concatenate((meth_array_abs,unmeth_array_abs))
    call_thresholds = np.arange(0,5,0.1)
    called = list()
    correct = list()
    for t in call_thresholds :
        if args.verbose : print("t :"+str(t),file=sys.stderr)
        called.append(sum(all_array_abs>t))
        correct.append(sum(meth_array>t)+sum(unmeth_array<-t))
    called = np.array(called)
    correct = np.array(correct)
    callthr_df = pd.DataFrame(
            {'thr' : np.append(call_thresholds,call_thresholds), 
                'data' : np.append(correct/called,[x/len(all_array_abs) for x in called]),
                'type' : ["correct"]*len(call_thresholds)+["called"]*len(call_thresholds)})
    if args.thrtable :
        callthr_df.to_csv(sys.stdout,index=False)
        return
    if not args.out :
        print("output file path required",file=sys.stderr)
        return
    # lines at thr = 2.5
    xint = 2.5
    yidx = np.argwhere(callthr_df['thr'] == 2.5)
    yint1 = callthr_df['data'][yidx[0]]
    yint2 = callthr_df['data'][yidx[1]]
#    thr = 0.95
#    intidx = callthr_df.index[callthr_df['data']>thr][0]
#    xint = callthr_df['thr'][intidx]
#    yint = callthr_df['data'][intidx]
    # plot
    g_roc = (ggplot(roc_df,aes(x='fallout',y='recall')) +
            geom_line() + lims(x=(0,1),y=(0,1)) +
            labs(title = "AUC = {}".format(round(auc,3)),
                x="False positive rate",
                y="True postivie rate") + 
            theme_bw() + 
            theme(panel_grid=element_blank(),
                axis_text=element_text(color="black"),
                figure_size=(3,3))
            ) 
    g_callthr = (ggplot(callthr_df,aes(x='thr',y='data'))+
            geom_line(aes(color='type',group='type'))+
            geom_vline(xintercept=xint,linetype="dotted")+
            scale_x_continuous(breaks=(0,xint,2.5,5))+
            geom_hline(yintercept=yint1,linetype="dotted")+
            geom_hline(yintercept=yint2,linetype="dotted")+
            scale_y_continuous(breaks=(0,0.25,0.5,yint1,yint2),limits=(0,1))+
            labs(x='Call threshold',y='Fraction')+
            theme_bw()+
            theme(panel_grid=element_blank(),
                axis_text=element_text(color="black"),
                figure_size=(3,3))
            )
    save_as_pdf_pages([g_roc,g_callthr],filename=args.out)

if __name__=="__main__":
    args=parseArgs()
    args.func(args)

