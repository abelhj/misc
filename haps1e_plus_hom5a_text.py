#!/usr/bin/env python

import networkx as nx
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
  usage_denom=1024
  Gloc=nx.Graph()
  Ghom=nx.Graph()

  homvar={}
  with gzip.open(args.infile2, 'rt') as fp:
    line=fp.readline().strip()
    while line:
      homvar[line]=1
      line=fp.readline().strip()

  sz=get_obj_size(homvar)/usage_denom
  sys.stderr.write('homvar\t'+str(sz))

  Gloc=nx.Graph()
  Gloconf=nx.Graph()

  with gzip.open(args.temp_prefix+'.het.txt.gz', 'rt') as fp:
    line=fp.readline().strip()
    ct=0
    while line:
      if ct%1000==0:
        sys.stderr.write(str(ct)+'\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom)+'\n')
      ll=re.split('[\t]', line)
      [loc1, loc2]=ll[0:2]
      tp='het-het'
      cts=list(map(int, ll[2:6]))
      mns=list(map(float, ll[6:10]))
      if loc2<loc1:
        [loc1, loc2]=[loc2, loc1]
        cts=[cts[0], cts[2], cts[1], cts[3]]
        mns=[mns[0], mns[2], mns[1], mns[3]]
      [passf, orient, nn, dist]=pre_filter_strict_pass(cts, mns, 0.95, tp)
      if passf:
        Gloc.add_edge(loc1, loc2,  conf=True, cts=cts, mns=mns, wt=nn, orient=orient, dist=dist)
      else:
        Gloc.add_node(loc1)
        Gloc.add_node(loc2)
      line=fp.readline().strip()
      ct+=1

  loc2comp={}; comp2tree={};



  with gzip.open(args.singletons, 'rt') as fp:
    line=fp.readline().strip()
    while line:
      if not Gloc.has_node(line):
        Gloc.add_node(line)
      line=fp.readline().strip()

  treeid=0
  gg=list(Gloc.subgraph(cc) for cc in sorted(nx.connected_components(Gloc), key=len, reverse=True))
  for ii in range(len(gg)):
    gg_temp=gg[ii].copy()
    print(str(ii)+'\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom/usage_denom))
    het_bridges=remove_bridges(gg_temp, min_counts_strict, 'het-hom')
    gg1=list(gg_temp.subgraph(cc).copy() for cc in sorted(nx.connected_components(gg_temp), key=len, reverse=True))
    for jj in range(len(gg1)):
      if gg1[jj].number_of_nodes()>2:
        tr=nx.minimum_spanning_tree(gg1[jj], weight='dist')
      else:
        tr=gg1[jj]
      comp2tree[treeid]=tr
      for node in tr.nodes():
        loc2comp[node]=treeid
      treeid+=1


  with gzip.open(args.temp_prefix+'.het.l2c.txt.gz', 'wt') as fp:
    for loc in loc2comp.keys():
      print(loc+'\t'+str(loc2comp[loc]), file=fp)

  with gzip.open(args.temp_prefix+'.het.c2t.txt.gz', 'wt') as fp:
    for comp in comp2tree.keys():
      tr=comp2tree[comp]
      for edge in tr.edges(data=True):
        print(str(comp)+'\t'+edge[0]+'\t'+edge[1]+'\t'+edge[2]['orient']+'\t'+str(edge[2]['dist'])+'\t'+str(edge[2]['wt']),  file=fp)

  loc2comphom={}
  comp2treehom={}


  with gzip.open(args.temp_prefix+'.hom.l2c.txt.gz', 'rt') as fp:
    line=fp.readline().strip()
    while line:
      [loc, comp]=re.split('[\t]', line)
      loc2comphom[loc]=int(comp)
      line=fp.readline().strip()

  oldcomp=-1
  with gzip.open(args.temp_prefix+'.hom.c2t.txt.gz', 'rt') as fp:
    line=fp.readline().strip()
    tr=nx.Graph()
    while line:
      [comp, loc1, loc2, dist]=re.split('[\t]', line)
      comp=int(comp); dist=int(dist)
      if comp==oldcomp or oldcomp<0:
        tr.add_edge(loc1, loc2, dist=dist)
        if oldcomp<0:
          oldcomp=comp
      else:
        comp2treehom[oldcomp]=tr.copy()
        tr=nx.Graph()
        oldcomp=comp
        tr.add_edge(loc1, loc2, dist=dist)
      line=fp.readline().strip()



  oldcomp=-1
  with gzip.open(args.temp_prefix+'.het.c2t.txt.gz', 'rt') as fp:
    line=fp.readline().strip()
    tr=nx.Graph()
    while line:
      [comp, loc1, loc2, orient, dist, wt]=re.split('[\t]', line)
      comp=int(comp); dist=int(float(dist)); wt=int(wt)
      if comp==oldcomp or oldcomp<0:
        tr.add_edge(loc1, loc2, dist=dist, wt=wt, orient=orient)
        if oldcomp<0:
          oldcomp=comp
      else:
        comp2tree[oldcomp]=tr.copy()
        tr=nx.Graph()
        oldcomp=comp
        tr.add_edge(loc1, loc2, dist=dist, wt=wt, orient=orient)
      line=fp.readline().strip()

  with gzip.open(args.temp_prefix+'.het.l2c.txt.gz', 'rt') as fp:
    line=fp.readline().strip()
    while line:
      [loc, comp]=re.split('[\t]', line)
      comp=int(comp)
      loc2comp[loc]=comp
      if not comp in comp2tree:
        tr=nx.Graph()
        tr.add_node(loc)
        comp2tree[comp]=tr
      line=fp.readline().strip()

  
  with gzip.open(args.temp_prefix+'hom.bed.gz', 'wt') as outb:
    for comp in comp2treehom.keys():
      id='hom_'+str(comp)
      mintree=comp2treehom[comp]
      str1=comp2bed(mintree.nodes, id, '0,150,0')
      if str1 is not None:
        print(str1, file=outb)

  with gzip.open(args.temp_prefix+'het.bed.gz', 'wt') as outb:
    for comp in comp2tree.keys(): 
      id='het_'+str(comp)
      print(id)
      mintree=comp2tree[comp]
      str1=comp2bed(mintree.nodes, id, '150,0,0')
      if str1 is not None:
        print(str1, file=outb)
  
#  code.interact(local=locals()) 

  Gmix=nx.Graph(); superg=nx.Graph()

  with gzip.open(args.temp_prefix+'.mixed.txt.gz', 'rt') as fp:
    line=fp.readline().strip()
    ct=0
    while line:
      if ct%1000==0:
        sys.stderr.write(str(ct)+'\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom)+'\n')
      ll=re.split('[\t]', line); cts=list(map(int, ll[2:6])); mns=list(map(float, ll[6:10]))
      [loc1, loc2]=ll[0:2] 
      tp='het-hom'
      [passf, orient, nn, dist]=pre_filter_strict_pass(cts, mns, 0.95, tp)
      if passf:
        if not loc2 in homvar :
          [loc1, loc2]=[loc2, loc1]; cts=[cts[0], cts[2], cts[1], cts[3]]; mns=[mns[0], mns[2], mns[1], mns[3]]
        if loc1 in loc2comp and loc2 in loc2comphom:
          Gmix.add_edge(loc1, loc2,  conf=True, cts=cts, mns=mns, wt=nn, orient=orient, dist=dist)
          hetcomp=loc2comp[loc1]; homcomp=loc2comphom[loc2];
          node1='het_'+str(hetcomp); node2='hom_'+str(homcomp)
          if not superg.has_edge(node1, node2):
            superg.add_edge(node1, node2, dist=[], wt=0, ct=0, sum_dist=0)
          superg.edges[node1, node2]['wt']+=nn
          superg.edges[node1, node2]['dist'].append(dist)
          superg.edges[node1, node2]['ct']+=1
          superg.edges[node1, node2]['sum_dist']+=dist
      line=fp.readline().strip()
      ct+=1

  gg=list(superg.subgraph(cc) for cc in sorted(nx.connected_components(superg), key=len, reverse=True))
  comp2treemixed={};  loc2compmixed={}

  code.interact(local=locals())

  for ii in range(len(gg)):
    id='mixed_'+str(ii)
    tocomp=[]
    for node in gg[ii].nodes():
      [tp, id]=re.split('_', node)
      print(node)
      if tp=='het':
        tr=comp2tree[int(id)]
      else:
        tr=comp2treehom[int(id)]
      for node in tr.nodes():
        tr.nodes[node]['tp']=tp
      tocomp.append(tr)
    Gcomp=nx.compose_all(tocomp)
    Gsub=Gmix.subgraph(Gcomp.nodes())
    Gcomp1=nx.compose(Gsub, Gcomp)
    mintree=nx.minimum_spanning_tree(Gcomp1, weight='dist')
    comp2treemixed[ii]=mintree
    for node in mintree.nodes():
      loc2compmixed[node]=ii

#  with open(args.temp_prefix+'.mixed.p', 'wb') as f:
#    ll=[loc2compmixed, comp2treemixed]
#    pickle.dump(ll, f)


  with gzip.open(args.temp_prefix+'mixed.bed.gz', 'wt') as outb:
    for comp in comp2treemixed.keys(): 
      id='mixed_'+str(comp)
      print(id)
      mintree=comp2treemixed[comp]
      str1=comp2bed(mintree.nodes, id, '0,0,150')
      if str1 is not None:
        print(str1, file=outb)

  code.interact(local=locals())
  for node in superg.nodes():
    deg=superg.degree(node)
    [tp, id]=re.split('_', node)
    if tp=='hom':
      print(node+'\t'+str(deg))



                                           
                                             

if __name__ == "__main__":
    main()
