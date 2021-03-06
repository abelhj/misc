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
  parser.add_argument('-c', '--compfile', type=str, default=None)
  parser.add_argument('-b', '--bedfile', type=str, default=None)
  parser.add_argument('-e', '--edgefile', type=str, default=None)
  parser.add_argument('-p', '--hapfile', type=str, default=None)

  args = parser.parse_args()
  min_counts_strict=5
  cp = cProfile.Profile()
  cp.enable()
  usage_denom=1024
  Gloc=nx.Graph()
  Ghom=nx.Graph()

  homvar={}
  with gzip.open(args.infile2, 'rt') as fp:
    line=fp.readline().strip()
    while line:
      homvar[line]=1
      line=fp.readline().strip()

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
        [passf, orient, nn, dist]=pre_filter_loose_pass(cts, mns, 2, 0.95, tp)
        if passf:
          Gloconf.add_edge(loc1, loc2,  conf=False, cts=cts, mns=mns, wt=nn, orient=orient, dist=dist)
      line=fp.readline().strip()
      ct+=1
     
  het_bridges=remove_bridges(Gloc, min_counts_strict, 'het-hom')
  loc2comp={}; comp2tree={};

  with gzip.open(args.singletons, 'rt') as fp:
    line=fp.readline().strip()
    while line:
      if not Gloc.has_node(line):
        Gloc.add_node(line)
      line=fp.readline().strip()

  gg=list(Gloc.subgraph(cc) for cc in sorted(nx.connected_components(Gloc), key=len, reverse=True))
  for ii in range(len(gg)):
    print(str(ii)+'\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom/usage_denom))
    if gg[ii].number_of_nodes()>2:
      tr=nx.minimum_spanning_tree(gg[ii], weight='dist')
    else:
      tr=gg[ii]
    comp2tree[ii]=tr
    for node in tr.nodes():
      loc2comp[node]=ii

  ll=[loc2comp, comp2tree]
  with open(args.temp_prefix+'.het.p', 'wb') as f:
    pickle.dump(ll, f)
  
  with open(args.temp_prefix+'.hom.p', 'rb') as f:
    [loc2comphom, comp2treehom]=pickle.load(f)

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
  comp2treemixed={}
  loc2compmixed={}

  for ii in range(len(gg)):
    id='mixed_'+str(ii)
    tocomp=[]
    for node in gg[ii].nodes():
      [tp, id]=re.split('_', node)
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

  with open(args.temp_prefix+'.mixed.p', 'wb') as f:
    ll=[loc2compmixed, comp2treemixed]
    pickle.dump(ll, f)


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

  for ii in range(18):
    id='mixed_'+str(ii)
    tocomp=[]
    for node in ['hom_228', 'het_374']:
      [tp, id]=re.split('_', node)
      if tp=='het':
        tr=comp2tree[int(id)]
      else:
        tr=comp2treehom[int(id)]
      for node in tr.nodes():
        tr.nodes[node]['tp']=tp
      tocomp.append(tr)
    Gcomp=nx.compose_all(tocomp); Gsub=Gmix.subgraph(Gcomp.nodes()); Gcomp1=nx.compose(Gsub, Gcomp)
    mintree=nx.minimum_spanning_tree(Gcomp1, weight='dist')
    comp2treemixed[ii]=mintree
    for node in mintree.nodes():
      loc2compmixed[node]=ii
    


hets=list(superg.neighbors('hom_228'))
ghom=comp2treehom[228]

g1=comp2treehom[228]
ghet1=comp2tree[282]
ghet2=comp2tree[374]
tocomp=[g1, ghet1, ghet2]
Gcomp=nx.compose_all(tocomp); Gsub=Gmix.subgraph(Gcomp.nodes()); Gcomp1=nx.compose(Gsub, Gcomp)
[start, stop]=nx.periphery(Gcomp1)
p1=nx.periphery(g1)
phet1=nx.periphery(ghet1)
phet2=nx.periphery(ghet2)                                      
pp=[p1, phet1, phet2]
dists={}
for px in pp:
  for ii in [0, 1]:
    dist=nx.shortest_path_length(Gcomp1, start, px[ii], weight='dist')
    dists[px[ii]]=dist

endpts=[k  for k, v in sorted(dists.items(), key=lambda item: item[1])]
for px in pp:
  d0=nx.shortest_path_length(Gcomp1, start, px[0], weight='dist')
  d1=nx.shortest_path_length(Gcomp1, start, px[1], weight='dist')
  print(str(d0)+'\t'+str(d1))
  if d1<d0:
    px.reverse()

homtree=comp2treehom[10]
list(superg.neighbors('hom_10'))
['het_272', 'het_335', 'het_88305', 'het_88343', 'het_88344', 'het_390']
hettree={}
pp={}
trees=[homtree]
pp['hom_10']=nx.periphery(homtree)
for nb in list(superg.neighbors('hom_10')):
  [pre, num]=re.split( '_', nb)
  print(num)
  hettree[nb]=comp2tree[int(num)]
  pp[nb]=nx.periphery(hettree[nb])
  trees.append(hettree[nb])
Gcomp=nx.compose_all(trees); Gsub=Gmix.subgraph(Gcomp.nodes()); Gcomp1=nx.compose(Gsub, Gcomp)
tr=nx.minimum_spanning_tree(Gcomp1, weight='dist')
[start, end]=nx.periphery(tr)
dists={}
for px in pp.keys():
  if len(pp[px])>1:
    d0=nx.shortest_path_length(Gcomp1, start, pp[px][0], weight='dist')
    d1=nx.shortest_path_length(Gcomp1, start, pp[px][1], weight='dist')
    print(str(d0)+'\t'+str(d1))
    dists[pp[px][0]]=d0
    dists[pp[px][1]]=d1
    if d1<d0:
      pp[px].reverse()



endpts=[k  for k, v in sorted(dists.items(), key=lambda item: item[1])]
hethom=[]
for pt in endpts:
  if pt in loc2comphom:
    hethom.append('hom')
  elif pt in loc2comp:
    hethom.append('het')
  else:
    print('error')
ll=[]
heton=True
for tt in hethom:
  if tt=='het':
    ll.append(heton)
    heton=not heton
  else:
    ll.append('hom')

ll1=[]
has_started=False
inhet=False
chunks=[]
for tt in ll:
  if tt=='hom':
    has_started=True
  elif not tt:
    
done=False
ff=0
while !done:


chunk=['chr10:50275417_C_T', 'chr10:50328897_C_T']
 

for px in pp:
  d0=nx.shortest_path_length(Gcomp1, start, px[0], weight='dist')
  d1=nx.shortest_path_length(Gcomp1, start, px[1], weight='dist')
  print(str(d0)+'\t'+str(d1))
  if d1<d0:
    px.reverse()


if __name__ == "__main__":
    main()
