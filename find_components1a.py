#!/usr/bin/env python

import gzip
import networkx as nx
import re
from collections import defaultdict
import math
import argparse
import cProfile , pstats , resource
import sys
import code
#code.interact(local=locals())

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', '--infile', type=str, default=None)
  parser.add_argument('-o', '--outfile', type=str, default=None)
  args = parser.parse_args()

  
  usage_denom=1024

  Gloc=nx.Graph()
  infile=args.infile

  with gzip.open(infile, 'rb') as fp:
    line=fp.readline().strip().decode()
    ct=0
    while line :
      if ct%1000==0:
        sys.stderr.write(str(ct)+'\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom)+'\n')
      ll=re.split('[\t]', line)
      cts=list(map(int, ll[2:6]))
      mns=list(map(float, ll[6:10]))
      sds=list(map(float, ll[10:14])) 
      [loc1, loc2]=ll[0:2]
      if pre_filter_strict_pass(cts):
        Gloc.add_edge(loc1, loc2);
        Gloc.edges[loc1, loc2].update({'conf': True, 'cts' : cts, 'mns': mns, 'sds' : sds, 'line': line})
      elif pre_filter_loose_pass(cts):
        Gloc.add_edge(loc1, loc2);
        Gloc.edges[loc1, loc2].update({'conf': False, 'cts' : cts, 'mns': mns, 'sds' : sds, 'line': line})
      line=fp.readline().strip().decode()
      ct+=1

  sys.stderr.write('finished loading graph\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom)+'\n')
  comp_loose=list(Gloc.subgraph(c) for c in sorted(nx.connected_components(Gloc), key=len, reverse=True))
  code.interact(local=locals())
  for ii in range(len(comp_loose)):
    gg_loose=comp_loose[ii]
    selected_edges = [(u,v) for u,v,e in gg_loose.edges(data=True) if e['conf'] == True]
    gg_conf=gg_loose.edge_subgraph(selected_edges)
    gg=list(gg_conf.subgraph(cc) for cc in sorted(nx.connected_components(gg_conf), key=len, reverse=True))
    for jj in range(len(gg)):
      edgelist=list(gg[jj].edges())
      for edge in edgelist:
        curedge=gg[jj].edges[edge]
        print(str(ii)+'\t'+str(jj)+'\t'+str(curedge['line'])+'\t')
                                                         
#code.interact(local=locals())

def pre_filter_strict_pass(cts, minfrac=0.95):
  [aa, bb, cc, dd]=cts
  nn=aa+bb+cc+dd
  passf=False
  bal=-1
  if nn<2:
    passf=False
  elif (aa+dd)*1.0/nn>minfrac:
    bal=aa*1.0/(aa+dd)
  elif (bb+cc)*1.0/nn>minfrac:
    bal=bb*1.0/(bb+cc)
  if bal>0.1 and bal<0.9:
    passf=True
  return passf

def pre_filter_loose_pass(cts, minfrac=0.90):
  [aa, bb, cc, dd]=cts
  nn=aa+bb+cc+dd
  passf=False
  if nn < 2:
    passf=False
  if  (aa+dd)*1.0/nn>minfrac:
    passf=True
  elif (bb+cc)*1.0/nn>minfrac:
    passf=True
  return passf
                                                                            


if __name__ == "__main__":
    main()
