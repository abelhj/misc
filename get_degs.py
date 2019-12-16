#!/usr/bin/env python

import gzip
import networkx as nx
import re
from collections import defaultdict
import math
import argparse
import cProfile , pstats , resource
import sys

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
      dists=list(map(float, ll[6:10]))
      [loc1, loc2]=ll[0:2]
      if pre_filter_strict_pass(cts):
        Gloc.add_edge(loc1, loc2);
        Gloc.edges[loc1, loc2].update({'conf': True, 'line': line, 'dists': dists})
      elif pre_filter_loose_pass(cts):
        Gloc.add_edge(loc1, loc2);
        Gloc.edges[loc1, loc2].update({'conf': False, 'line': line, 'dists': dists})
      line=fp.readline().strip().decode()
      ct+=1

  sys.stderr.write('finished loading graph\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom)+'\n')
  comp_loose=list(Gloc.subgraph(c) for c in sorted(nx.connected_components(Gloc), key=len, reverse=True))
  for ii in range(len(comp_loose)):
    gg_loose=comp_loose[ii]
    selected_edges = [(u,v) for u,v,e in gg_loose.edges(data=True) if e['conf'] == True]
    gg_conf=gg_loose.edge_subgraph(selected_edges)
    gg=list(gg_conf.subgraph(cc) for cc in sorted(nx.connected_components(gg_conf), key=len, reverse=True))
    for jj in range(len(gg)):
      nodelist=list(gg[jj].nodes())
      for node in nodelist:
        print(node+'\t'+str(gg_loose.degree(node))+'\t'+str(gg[jj].degree(node)))
      #edgelist=list(gg[jj].edges())
      #for edge in edgelist:
      #curedge=gg[jj].edges[edge]
      #print(str(ii)+'\t'+str(jj)+'\t'+curedge['line']+'\t'+str(curedge['conf']))
                                                         


                       


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

def pre_filter_loose_pass(cts, minfrac=0.95):
  [aa, bb, cc, dd]=cts
  nn=aa+bb+cc+dd
  passf=False
  if aa==nn or bb==nn or cc==nn or dd==nn:
    return passf
  if pre_filter_strict_pass(aa, bb, cc, dd, minfrac-0.05):
    passf=True
  else:
    bal=-1
    if (aa+dd)*1.0/nn>minfrac:
      bal=aa*1.0/(aa+dd)
    elif (bb+cc)*1.0/nn>minfrac:
      bal=bb*1.0/(bb+cc)
    if bal>0.01 and bal<0.99:
      passf=True
  return passf
                                                                            

      
  


if __name__ == "__main__":
    main()
