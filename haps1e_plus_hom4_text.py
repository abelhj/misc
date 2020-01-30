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
  Ghom=nx.Graph()

  homvar={}
  with gzip.open(args.infile2, 'rt') as fp:
    line=fp.readline().strip()
    while line:
      homvar[line]=1
      line=fp.readline().strip()

  with gzip.open(args.infile, 'rt') as fp, gzip.open(args.temp_prefix+'.het.txt.gz', 'wt') as fhet, gzip.open(args.temp_prefix+'.mixed.txt.gz', 'wt') as fmix:
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
          Ghom.add_edge(loc1, loc2, orient=orient, dist=int(dist), wt=nn)
      line=fp.readline().strip()
      ct+=1

  hom_bridges=remove_bridges(Ghom, min_counts_strict, 'hom-hom')
  loc2comphom={}; comp2treehom={};

  gg=list(Ghom.subgraph(cc) for cc in sorted(nx.connected_components(Ghom), key=len, reverse=True))
  for ii in range(len(gg)):
    print(str(ii))
    tr=nx.minimum_spanning_tree(gg[ii], weight='dist')
    comp2treehom[ii]=tr
    for node in tr.nodes():
      loc2comphom[node]=ii
  
  Ghom=None
  with gzip.open(args.temp_prefix+'.hom.l2c.txt.gz', 'wt') as fp:
    for loc in loc2comphom.keys():
      print(loc+'\t'+str(loc2comphom[loc]), file=fp)

  with gzip.open(args.temp_prefix+'.hom.c2t.txt.gz', 'wt') as fp:
    for comp in comp2treehom.keys():
      tr=comp2treehom[comp]
      for edge in tr.edges(data=True):
        print(str(comp)+'\t'+edge[0]+'\t'+edge[1]+'\t'+edge[2]['orient']+'\t'+str(edge[2]['dist'])+'\t'+str(edge[2]['wt']),  file=fp)
                                           
                                             

if __name__ == "__main__":
    main()
