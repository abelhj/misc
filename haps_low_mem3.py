#!/usr/bin/env python

#from __future__ import print_function
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
  parser.add_argument('-d', '--debugfile', type=str, default=None)
  args = parser.parse_args()

  #cp = cProfile.Profile()
  #cp.enable()
  usage_denom=1024
  print('start\tMemory usage info (Mb):\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom))


  Gloc=nx.Graph()
  infile=args.infile
  linect=0
  with gzip.open(infile, 'rb') as fp:
    line=fp.readline().strip()
    while line:
      if(linect % 1000 == 0):
        print("line "+str(linect)+'\tMemory (Mb):\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom))
      [id, loci, alleles, poses]=re.split('[\t]', line.decode())
      locus=re.split(',', loci)
      allele=re.split(',', alleles)
      pos=list(map(int, re.split(',', poses)))
      if len(pos)>2:
        for ii in range(len(pos)-1):
          for jj in range(1, len(pos)):
            dist=abs(pos[ii]-pos[jj])
            add_edge(locus[ii], locus[jj], allele[ii], allele[jj], dist , Gloc)
      line=fp.readline().strip()
      linect+=1


  edgelist=list(Gloc.edges) 
  with open(args.debugfile, 'w') as outd:
    for edge in edgelist:
      curedge=Gloc.edges[edge]
      mn=curedge['dists']
      cts=curedge['counts']
      for ii in range(4):
        if cts[ii]>0:
          mn[ii]=round(mn[ii]*1/cts[ii],2)
        else:
          mn[ii]=-1
      mnstr='\t'.join(map(str, mn))
      ctstr='\t'.join(map(str, cts))
      print(edge[0]+'\t'+edge[1]+'\t'+ctstr+'\t'+mnstr, file=outd)

          

  #cp.disable()
  #cp.print_stats()


def add_edge(loc1, loc2, all1, all2, dist, gg):
  if not gg.has_edge(loc1, loc2):
    gg.add_edge(loc1, loc2)
    curedge=gg.edges[loc1, loc2]
    curedge.update({'l1': loc1, 'l2': loc2, 'counts':[0,0,0,0], 'dists':[0, 0, 0, 0]})
    if curedge['l1']==loc1:
      if all1=='a' and all2=='a':
        curedge['counts'][0]+=1
        curedge['dists'][0]+=dist
      elif all1=='a' and all2=='r':
        curedge['counts'][1]+=1
        curedge['dists'][1]+=dist
      elif all1=='r' and all2=='a':
        curedge['counts'][2]+=1
        curedge['dists'][2]+=dist
      elif all1=='r' and all2=='r':
        curedge['counts'][3]+=1
        curedge['dists'][3]+=dist
    elif curedge['l2']==loc1:
      if all1=='a' and all2=='a':
        curedge['counts'][0]+=1
        curedge['dists'][0]+=dist
      elif all1=='r' and all2=='a':
        curedge['counts'][1]+=1
        curedge['dists'][1]+=dist
      elif all1=='a' and all2=='r':
        curedge['counts'][2]+=1
        curedge['dists'][2]+=dist
      elif all1=='r' and all2=='r':
        curedge['counts'][3]+=1
        curedge['dists'][3]+=dist
  





if __name__ == "__main__":
    main()
