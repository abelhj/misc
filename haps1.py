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
  elif  (aa+dd)*1.0/nn>minfrac or (bb+cc)*1.0/nn>minfrac:
    passf=True
  return passf


def phase_conf_component(ggsub):
    maxnode=list(ggsub.nodes)[0]
    maxval=-1
    for node in ggsub.nodes:
      deg=ggsub.degree(node)
      if deg>maxval:
        maxval=deg
        maxnode=node
    ggsub.nodes[maxnode]['phased_all']=[0,1]
    starts=[]
    starts.append(maxnode)
    #phased_nodes=set()
    #phased_nodes.add(maxnode)
    [h0, h1]=phase_from_node(ggsub, maxnode)
    return([h0, h1])

def get_phased_allele(gg, node, hap):
  refalt=gg.nodes[node]['phased_all'][hap]
  return(node+'_'+str(refalt))

def phase_from_node(gg1, startnode):
  hap0=[]
  hap1=[]
  if len(gg1)==1:
    hap0.append(get_phased_allele(gg1, startnode, 0))
    hap1.append(get_phased_allele(gg1, startnode, 1))
  else:
    edgelist=list(nx.bfs_edges(gg1, startnode))
    edge0=edgelist[0]
    curnode=edge0[0]
    hap0.append(get_phased_allele(gg1, curnode, 0))
    hap1.append(get_phased_allele(gg1, curnode, 1))
    for ii in range(len(edgelist)):
      edge0=edgelist[ii]
      curnode=edge0[0]
      edata=gg1.edges[edge0]
      if 'phased_all' in gg1.nodes[curnode]:
        alleles=gg1.nodes[curnode]['phased_all']
        if edata['order'][0]==curnode:
          nextnode=edata['order'][1]
        else:
          nextnode=edata['order'][0]
        oddsratio=(edata['cts'][0]+1.0)*(edata['cts'][3]+1.0)/((edata['cts'][1]+1.0)*(edata['cts'][2]+1.0))
        if oddsratio>1:
          gg1.nodes[nextnode]['phased_all']=alleles
        else:
          gg1.nodes[nextnode]['phased_all']=[alleles[1], alleles[0]]
        #phased_nodes.add(nextnode)
        hap0.append(get_phased_allele(gg1, nextnode, 0))
        hap1.append(get_phased_allele(gg1, nextnode, 0))
  return [hap0, hap1]
                                                                                                                                                                                                        

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
        Gloc.edges[loc1, loc2].update({'order': [loc1, loc2], 'conf': True, 'cts' : cts, 'mns': mns, 'sds' : sds})
      elif pre_filter_loose_pass(cts):
        Gloc.add_edge(loc1, loc2);
        Gloc.edges[loc1, loc2].update({'order': [loc1, loc2], 'conf': False, 'cts' : cts, 'mns': mns, 'sds' : sds})
      line=fp.readline().strip().decode()
      ct+=1

  sys.stderr.write('finished loading graph\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom)+'\n')
  comp_loose=list(Gloc.subgraph(c) for c in sorted(nx.connected_components(Gloc), key=len, reverse=True))
  for ii in range(len(comp_loose)):
    gg_loose=comp_loose[ii]
    selected_edges = [(u,v) for u,v,e in gg_loose.edges(data=True) if  e['conf'] == True]
    gg_conf=gg_loose.edge_subgraph(selected_edges)
    brlist=list(nx.bridges(gg_conf))
    for edge in brlist:
      if sum(gg_loose.edges[edge]['cts'])<3 or min(gg_loose.edges[edge]['cts'])<1:
        gg_loose.edges[edge].update({'conf': False})
    selected_edges = [(u,v) for u,v,e in gg_loose.edges(data=True) if  e['conf'] == True]
    gg_conf=gg_loose.edge_subgraph(selected_edges)
    gg=list(gg_conf.subgraph(cc) for cc in sorted(nx.connected_components(gg_conf), key=len, reverse=True))
    subgraphs1=[]
    subgraphs2=[]
    for jj in range(len(gg)):
      [h1, h2]=phase_conf_component(gg[jj])
      subgraphs1.append(h1)
      subgraphs2.append(h2)
    allhaps=[subgraphs1, subgraphs2]
                                                         

                                             

if __name__ == "__main__":
    main()