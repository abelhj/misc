#!/usr/bin/env python

import networkx as nx
import numpy as np
import re, math, argparse, gzip
import cProfile, pstats, resource
import sys, code, pickle
from collections import Counter, defaultdict
from utils2 import *

                                        
def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', '--infile', type=str, default=None)
  parser.add_argument('-s', '--singletons', type=str, default=None)
  parser.add_argument('-t', '--temp_prefix', type=str, default=None)
  parser.add_argument('-f', '--infile2', type=str, default=None)

  args = parser.parse_args()
  min_counts_strict=5
  usage_denom=1024*1024

  homvar={}
  with gzip.open(args.infile2, 'rt') as fp:
    line=fp.readline().strip()
    while line:
      homvar[line]=1
      line=fp.readline().strip()



  with gzip.open(args.infile, 'rt') as fp, gzip.open(args.temp_prefix+'.het.txt.gz', 'wt') as fhet, gzip.open(args.temp_prefix+'.mixed.txt.gz', 'wt') as fmix, gzip.open(args.temp_prefix+'.hom.txt.gz', 'wt') as fhom:
    line=fp.readline().strip()
    ct=0
    while line:
      if ct%1000==0:
        sys.stderr.write(str(ct)+'\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom)+'\n')
      ll=re.split('[\t]', line)
      [loc1, loc2]=ll[0:2]      
      if not (loc1 in homvar or loc2 in homvar):
        tp='het-het'
        print(line, file=fhet)
      elif not (loc1 in homvar and loc2 in homvar):
        tp='het-hom'
        print(line, file=fmix)
      else:
        tp='hom-hom'
        cts=list(map(int, ll[2:6]))
        mns=list(map(float, ll[6:10]))
        if loc2<loc1:
          [loc1, loc2]=[loc2, loc1]
          cts=[cts[0], cts[2], cts[1], cts[3]]
          mns=[mns[0], mns[2], mns[1], mns[3]]
        [passf, orient, nn, dist]=pre_filter_strict_pass(cts, mns, 0.95, tp)
        if passf:
          print(loc1+'\t'+loc2+'\t'+str(int(dist))+'\t'+orient+'\t'+str(nn), file=fhom)
      line=fp.readline().strip()
      ct+=1
  
  Ghom=nx.Graph()
  ct=0
  with gzip.open(args.temp_prefix+'.hom.txt.gz', 'rt') as fp:
    line=fp.readline().strip()
    while line:
      if ct%1000==0:
        sys.stderr.write(str(ct)+'\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom)+'\n')
      [loc1, loc2, dist, orient, nn]=re.split('[\t]', line)
      Ghom.add_edge(loc1, loc2, dist=int(dist))
      line=fp.readline().strip()
      ct+=1
  
  sys.stderr.write('graph in memory\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom)+'\n')  
  gg=list(Ghom.subgraph(cc) for cc in sorted(nx.connected_components(Ghom), key=len, reverse=True))
  sys.stderr.write('comp decomp\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom)+'\n')
  
  with gzip.open(args.temp_prefix+'.comps.hom.txt.gz', 'wt') as fp:
    for ii in range(len(gg)):
      for edge in gg[ii].edges(data=True):
        print(str(ii)+'\t'+edge[0]+'\t'+edge[1]+'\t'+str(edge[2]['dist']), file=fp)

  Ghom=None; homvar=None; gg=None
  oldcomp=-1; treeid=0; ct=0; gg=nx.Graph()

  #code.interact(local=locals())

  with gzip.open(args.temp_prefix+'.comps.hom.txt.gz', 'rt') as fp, gzip.open(args.temp_prefix+'.hom.l2c.txt.gz', 'wt') as l2c, gzip.open(args.temp_prefix+'.hom.c2t.txt.gz', 'wt') as c2t:
    line=fp.readline().strip()
    while line:
      if ct%1000==0:
        sys.stderr.write(str(ct)+'\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom)+'\n')
      [comp, loc1, loc2, dist]=re.split('[\t]', line)
      #print(str(treeid)+'\t'+str(comp))
      comp=int(comp); dist=int(dist)
      #if comp>2:
      #  code.interact(local=locals())  
      if comp==oldcomp or oldcomp==-1:
        gg.add_edge(loc1, loc2, dist=dist)
        if oldcomp==-1:
          oldcomp=comp
      else:
        hom_bridges=remove_bridges(gg, min_counts_strict, 'hom-hom')
        gg1=list(gg.subgraph(cc) for cc in sorted(nx.connected_components(gg), key=len, reverse=True))
        for jj in range(len(gg1)):
          if gg1[jj].number_of_nodes()>2:
            tr=nx.minimum_spanning_tree(gg1[jj], weight='dist')
            for edge in tr.edges(data=True):
              print(str(treeid)+'\t'+edge[0]+'\t'+edge[1]+'\t'+str(edge[2]['dist']),  file=c2t)
            for node in tr.nodes():
              print(node+'\t'+str(treeid), file=l2c)
            treeid+=1
        gg=nx.Graph()
        gg.add_edge(loc1, loc2, dist=dist)
        oldcomp=comp
      line=fp.readline().strip()
      ct+=1

                                             

if __name__ == "__main__":
    main()
