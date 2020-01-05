#!/usr/bin/env python

from __future__ import print_function
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
  parser.add_argument('-e', '--edgefile', type=str, default=None)
  parser.add_argument('-p', '--hapfile', type=str, default=None)
  parser.add_argument('--strict', default=False, action='store_true')
  args = parser.parse_args()

  cp = cProfile.Profile()
  cp.enable()
  usage_denom=1024*1000
  print('before')
  print('Memory usage info (Mb):\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom))

  Gall=nx.Graph()
  Gloc=nx.Graph()
  infile=args.infile

  #load allele graph
  reads=defaultdict(dict)
  readid={}
  locid={}
  id2loc={}
  [readct, locct, allct]=[0,0,0]
  with gzip.open(infile, 'rb') as fp:
    line=fp.readline().strip()
    while line :
      [allele, readstr, readpos, strand]=re.split('[\t]', line)
      [loc, ref, alt, refalt]=re.split('_', allele)
      locstr=loc+'_'+ref+'_'+alt
      readpos=int(readpos)
      if not readstr in readid:
        readid[readstr]=readct
        readct+=1
      read=readid[readstr]
      print(str(read))
      if not locstr in locid:
        lid='loc'+str(locct)
        locid[locstr]=lid
        id2loc[lid]=locstr
        locct+=1
      locus=locid[locstr]
      allele=locus+'_'+refalt
      reads[read][allele]=readpos
      
      if not Gloc.has_node(locus):
        Gloc.add_node(locus)
      if not Gall.has_node(allele):
        Gall.add_node(allele, refalt=refalt)
      for r1 in reads[read]:
        if r1!=allele:
          if not Gall.has_edge(allele, r1):
            Gall.add_edge(allele, r1, count=0, dist=0, dist_sq=0, reads=[])
          curedge=Gall.edges[allele, r1]
          rr=curedge['reads']
          rr.append(read)
          dist=abs(reads[read][r1]-reads[read][allele])
          Gall.edges[allele, r1].update({'count': curedge['count']+1, 'dist': curedge['dist']+dist, 'dist_sq': curedge['dist_sq']+dist*dist, 'reads' : rr})
      line=fp.readline().strip()

  del reads
  print('allele graph')
  print('Memory usage info (Mb):\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom))

  #load locus graph
  for edge in Gall.edges:
    curedge=Gall.edges[edge]
    mndist=curedge['dist']*1.0/curedge['count']
    sddist=math.sqrt(curedge['dist_sq']*1.0/curedge['count']-mndist*mndist)
    Gall.edges[edge].update({'mean_dist' : mndist, 'sd_dist': sddist})
    add_counts_edge(Gloc, Gall, edge[0], edge[1], curedge['count'])
  print('locus graph')
  print('Memory usage info (Mb):\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom))

  #create strict locus graph and prune edges
  Gconf=Gloc.copy()
  edgelist=list(Gloc.edges)
  for edge in edgelist:
    curedge=Gloc.edges[edge]
    [aa, bb, cc, dd]=[curedge['r', 'r'], curedge['r', 'a'], curedge['a', 'r'], curedge['a', 'a']]
    if (aa+bb==0 and cc*dd>0) or (cc+dd==0 and aa*bb>0) or (aa+cc==0 and bb*dd>0) or (bb+dd==0 and aa*cc>0): # could add addl paralog filter here
      remove_allele_edges(edge[0], edge[1], Gall)
      Gloc.remove_edge(edge[0], edge[1])
      Gconf.remove_edge(edge[0], edge[1])
    else:
      oddsratio=(aa+1.0)*(dd+1.0)/(bb+1.0)/(cc+1.0)
      if (oddsratio>2 or oddsratio<0.5) and (aa*dd>0 or bb*cc>0):
        pass
      else:
        Gconf.remove_edge(edge[0], edge[1])

  print('strict locus  graph')
  print('Memory usage info (Mb):\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom))


  with open(args.edgefile, 'w') as outf, open(args.hapfile, 'w') as outh:
    comp_loose=list(Gloc.subgraph(c) for c in sorted(nx.connected_components(Gloc), key=len, reverse=True))
    #iterate over connected components of loose locus graph, phase and merge if possible
    for ii in range(len(comp_loose)):
      gg_loose=comp_loose[ii]
      gg_conf=Gconf.subgraph(list(gg_loose.nodes()))
      gg=list(gg_conf.subgraph(cc).copy() for cc in sorted(nx.connected_components(gg_conf), key=len, reverse=True))
      subgraphs1=[]
      subgraphs2=[]
      for jj in range(len(gg)):
        ggsub1=gg[jj]
        [h1, h2]=phase_conf_component(ggsub1)
        subgraphs1.append(h1)
        subgraphs2.append(h2)
      if(len(gg)>1):
        allhaps=merge_subgraphs(subgraphs1, subgraphs2, Gall)
      else:
        allhaps=[list([h1]), list([h2])]
      print(str(ii)+'\t'+str(jj)+'\t'+str(len(allhaps[0]))+'*')
      for jj in range(len(allhaps[0])):
        for hapid in range(2):
          minforest=nx.minimum_spanning_tree(Gall.subgraph(allhaps[hapid][jj]), weight='mean_dist')
          mintree=list(minforest.subgraph(cc) for cc in nx.connected_components(minforest))
          for treeid in range(len(mintree)):
            treelist=list(mintree[treeid].edges)
            for tredge in treelist:
              curedge=mintree[treeid].edges[tredge]
              [id1, refalt1]=re.split('_', tredge[0])
              [id2, refalt2]=re.split('_', tredge[1])
              edgestr=str(round(curedge['mean_dist'], 3))+';'+str(round(curedge['sd_dist'], 3))+';'+str(curedge['count'])
              print(str(ii)+'_'+str(jj)+'_'+str(hapid+1)+'_'+str(treeid)+'\t'+edgestr+'\t'+id2loc[id1]+'_'+refalt1+'\t'+id2loc[id2]+'_'+refalt2, file=outf)
            terminal_nodes=[]
            for node1 in list(mintree[treeid].nodes):
              if mintree[treeid].degree(node1)==1:
                terminal_nodes.append(node1)
            for aa in range(len(terminal_nodes)-1):
              for bb in range(aa+1, len(terminal_nodes)):
                sp=nx.shortest_path(mintree[treeid], terminal_nodes[aa], terminal_nodes[bb])
                for nodeii in range(len(sp)):
                  node1=sp[nodeii]
                  [id1, refalt1]=re.split('_', node1)
                  if nodeii==0:
                    outstr=str(mintree[treeid].degree(node1))+'\t.'
                  else:
                    prevedge=mintree[treeid].edges[sp[nodeii], sp[nodeii-1]]
                    outstr=str(mintree[treeid].degree(node1))+'\t'+str(prevedge['count'])+'_'+str(prevedge['dist'])+'_'+str(prevedge['dist_sq'])
                  print(str(ii)+'_'+str(jj)+'_'+str(hapid+1)+'_'+str(treeid)+'\t'+str(aa)+'_'+str(bb)+'\t'+outstr+'\t'+id2loc[id1]+'_'+refalt, file=outh)        
      print('-----------------------------------------------------------------')
      print('Memory usage info (Mb):\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom))

  cp.disable()
  cp.print_stats()


def add_counts_edge(gg, ggall, node1, node2, count):
  locus1=re.split('_', node1)[0]
  locus2=re.split('_', node2)[0]
  if not gg.has_edge(locus1, locus2):
    gg.add_edge(locus1, locus2)
    gg.edges[locus1, locus2].update({'first': locus1, 'second': locus2, ('a', 'a') : 0, ('a', 'r') : 0, ('r', 'a') : 0, ('r', 'r') : 0})
  curedge=gg.edges[locus1, locus2]
  if curedge['first']==locus1:
    curedge[(ggall.nodes[node1]['refalt'], ggall.nodes[node2]['refalt'])]=count
  else:
    curedge[(ggall.nodes[node2]['refalt'], ggall.nodes[node1]['refalt'])]=count

def remove_allele_edges(node1, node2, ggall):
  all=['r', 'a']
  for aa in all:
    for bb in all:
      if ggall.has_edge(node1+'_'+aa, node2+'_'+bb):
        ggall.remove_edge(node1+'_'+aa, node2+'_'+bb)


def phase_conf_component(ggsub,  flip=False, strict=True):
  maxnode=list(ggsub.nodes)[0]
  maxval=-1
  for node in ggsub.nodes:
    deg=ggsub.degree(node)
    if deg>maxval:
      maxval=deg
      maxnode=node
  ggsub.nodes[maxnode]['phased_all']={1: 'r', 2: 'a'}
  if flip:
    ggsub.nodes[maxnode]['phased_all']={1: 'a', 2: 'r'}
  starts=[]
  starts.append(maxnode)
  phased_nodes=set()
  phased_nodes.add(maxnode)
  [h1, h2]=phase_from_node(ggsub, maxnode, phased_nodes,  strict)
  return([h1, h2])

def get_phased_allele(gg, node, hap):
  refalt=gg.nodes[node]['phased_all'][hap]
  return(node+'_'+refalt)


def phase_from_node(gg1, startnode, phased_nodes,  strict=True):
  hap1=[]
  hap2=[]
  if len(gg1)==1:
    hap1.append(get_phased_allele(gg1, startnode, 1))
    hap2.append(get_phased_allele(gg1, startnode, 2))
  else:
    edgelist=list(nx.bfs_edges(gg1, startnode))
    edge0=edgelist[0]
    curnode=edge0[0]
    hap1.append(get_phased_allele(gg1, curnode, 1))
    hap2.append(get_phased_allele(gg1, curnode, 2))
    for ii in range(len(edgelist)):
      edge0=edgelist[ii]
      curnode=edge0[0]
      edata=gg1.edges[edge0]
      if 'phased_all' in gg1.nodes[curnode]:
        cur=[gg1.nodes[curnode]['phased_all'][1], gg1.nodes[curnode]['phased_all'][2]]
        if edata['first']==curnode:
          nextnode=edata['second']
          [aa, bb, cc, dd]=[edata[(cur[0], 'r')], edata[(cur[0], 'a')], edata[(cur[1], 'r')], edata[(cur[1], 'a')]]
        else:
          nextnode=edata['first']
          [aa, bb, cc, dd]=[edata[('r', cur[0])], edata[('a', cur[0])], edata[('r', cur[1])], edata[('a', cur[1])]]
        oddsratio=(aa+1.0)*(dd+1.0)/(bb+1.0)/(cc+1.0)
        if not strict or (strict and  (aa*dd>0 or bb*cc>0)): #see both haps
          if oddsratio>1:
            gg1.nodes[nextnode]['phased_all']={1: 'r', 2: 'a'}
          else:
            gg1.nodes[nextnode]['phased_all']={1: 'a', 2: 'r'}
          phased_nodes.add(nextnode)
          hap1.append(get_phased_allele(gg1, nextnode, 1))
          hap2.append(get_phased_allele(gg1, nextnode, 2))
  return [hap1, hap2]
  


def check_join(G, n1i, n2i, n1j, n2j):
  [aa, bb, cc, dd]=[0,0,0,0]
  reads1=set()
  reads2=set()
  for node in n1i:
    if G.has_node(node):
      nb=set(list(G.neighbors(node)))
      nb1=nb.intersection(n1j)
      nb2=nb.intersection(n2j)
      for node2 in nb1:
        curedge=G.edges[node, node2]
        for read in curedge['reads']:
          reads1.add(read)
      for node2 in nb2:
        curedge=G.edges[node, node2]
        for read in curedge['reads']:
          reads2.add(read)
  aa=len(reads1)
  bb=len(reads2)
  reads1=set()
  reads2=set()
  for node in n2i:
    if G.has_node(node):
      nb=set(list(G.neighbors(node)))
      nb1=nb.intersection(n1j)
      nb2=nb.intersection(n2j)
      for node2 in nb1:
        curedge=G.edges[node, node2]
        for read in curedge['reads']:
          reads1.add(read)
      for node2 in nb2:
        curedge=G.edges[node, node2]
        for read in curedge['reads']:
          reads2.add(read)
  cc=len(reads1)
  dd=len(reads2)
  return([aa, bb, cc, dd])

def merge_subgraphs(sg1, sg2, G):
  done=False
  haps1=[]
  haps2=[]
  sg1_temp=list(sg1)
  sg2_temp=list(sg2)
  max1=sg1_temp[0]
  max2=sg2_temp[0]
  while not done:
    if len(sg1_temp)==1:
      done=True
      haps1.append(max1)
      haps2.append(max2)
      break
    adds=[]
    for jj in range(1, len(sg1_temp)):
      [aa, bb, cc, dd]=check_join(G, max1, max2, sg1_temp[jj], sg2_temp[jj])
      atot=aa+bb+cc+dd
      print('attempt to merge subgraph '+str(jj)+' '+str([aa, bb, cc, dd]))
      if atot>1 and (aa+dd==atot or bb+cc==atot):
        print('merge')
        adds.append(jj)
        if  aa+dd==atot:
          max1.extend(sg1_temp[jj])
          max2.extend(sg2_temp[jj])
        elif bb+cc==atot:
          max1.extend(sg2_temp[jj])
          max2.extend(sg1_temp[jj])
    if len(adds)>0:
      adds.sort(reverse=True)
      for jj in adds:
        sg1_temp.pop(jj)
        sg2_temp.pop(jj)
    else:
      haps1.append(sg1_temp.pop(0))
      haps2.append(sg2_temp.pop(0))
      max1=sg1_temp[0]
      max2=sg2_temp[0]
  return [haps1, haps2]



if __name__ == "__main__":
    main()
