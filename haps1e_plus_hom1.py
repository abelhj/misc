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
from collections import Counter
from utils import *
import pickle


                                        
def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', '--infile', type=str, default=None)
  parser.add_argument('-s', '--singletons', type=str, default=None)
  parser.add_argument('-f', '--infile2', type=str, default=None)
  parser.add_argument('-c', '--compfile', type=str, default=None)
  parser.add_argument('-b', '--bedfile', type=str, default=None)
  parser.add_argument('-e', '--edgefile', type=str, default=None)
  parser.add_argument('-p', '--hapfile', type=str, default=None)
  args = parser.parse_args()
  min_counts_strict=5

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

  with gzip.open(args.infile, 'rt') as fp:
    line=fp.readline().strip()
    ct=0
    while line:
      if ct%1000==0:
        sys.stderr.write(str(ct)+'\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom)+'\n')
      ll=re.split('[\t]', line)
      cts=list(map(int, ll[2:6]))
      mns=list(map(float, ll[6:10]))
      tot=sum(cts)
      [loc1, loc2]=ll[0:2]
      tp='het-hom'
      if loc1 in homvar and loc2 in homvar:
        tp='hom-hom'
      elif not (loc1 in homvar or loc2 in homvar):
        tp='het-het'
      if loc2<loc1:
        [loc1, loc2]=[loc2, loc1]
        cts=[cts[0], cts[2], cts[1], cts[3]]
        mns=[mns[0], mns[2], mns[1], mns[3]]
      conf=False
      if pre_filter_strict_pass(cts, 0.95, tp):
        conf=True
        if tp=='het-het':
          Gloc.add_edge(loc1, loc2,  conf=conf, cts=cts, mns=mns, wt=tot)
        elif tp=='hom-hom' and conf:
          Ghom.add_edge(loc1, loc2,  conf=conf, cts=cts, mns=mns, wt=tot)
      line=fp.readline().strip()
      ct+=1
   

  remove_bridges(Gloc, min_counts_strict)
  remove_bridges(Ghom, min_counts_strict)

  for edge in Ghom.edges(data=True):
    nn=sum(edge[2]['cts'])
    for ii in range(4):
      if edge[2]['cts'][ii]>0.9*nn:
        edge[2]['dist']=edge[2]['mns'][ii]

  for edge in Gloc.edges(data=True):
    nn=sum(edge[2]['cts'])
    if edge[2]['cts'][0]+edge[2]['cts'][3]>0.8*nn:
      edge[2]['orient']='outer'
      edge[2]['dist']=0.5*(edge[2]['mns'][0]+edge[2]['mns'][3])
    if edge[2]['cts'][1]+edge[2]['cts'][2]>0.8*nn:
      edge[2]['orient']='inner'
      edge[2]['dist']=0.5*(edge[2]['mns'][1]+edge[2]['mns'][2])

  loc2comphom={}; comp2treehom={}; loc2comp={}; comp2tree={}

  gg=list(Ghom.subgraph(cc) for cc in sorted(nx.connected_components(Ghom), key=len, reverse=True))
  for ii in range(len(gg)):
    print(str(ii))
    tr=nx.minimum_spanning_tree(gg[ii], weight='dist')
    comp2treehom[ii]=tr
    for node in tr.nodes():
      loc2comphom[node]=ii

  Ghom=None
  gg=None

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


  with gzip.open(args.bedfile, 'wt') as outb:
    for comp in comp2treehom.keys():
      id='hom_'+str(comp)
      mintree=comp2treehom[comp]
      #nodes=[int(re.split('[:_]', node)[1]) for node in mintree.nodes]
      #chrs=[re.split('[:_]', node)[0] for node in mintree.nodes]
      nodes=[int(re.split('[:_]', node)[1]) for node in mintree.nodes if not 'ctg' in node]
      chrs=[re.split('[:_]', node)[0] for node in mintree.nodes if not 'ctg' in node]
      chr=Counter(chrs).most_common(1)[0][0]
      nodes.sort()
      dists=[]
      ones=[]
      for kk in range(len(nodes)):
        dists.append(nodes[kk]-nodes[0]+1)
        ones.append(1)
        dstr=','.join(map(str, dists))
        onestr=','.join(map(str,ones))
      print(chr+'\t'+str(nodes[0])+'\t'+str(max(nodes))+'\t'+id+'\t100\t.\t'+str(nodes[0])+'\t'+str(max(nodes))+'\t150,150,0\t'+str(len(dists))+'\t'+onestr+'\t'+dstr, file=outb)

  with gzip.open('../het.chr10.bed.gz', 'wt') as outb:
    for comp in comp2tree.keys(): 
      id='het_'+str(comp)
      mintree=comp2tree[comp]
      nodes=[int(re.split('[:_]', node)[1]) for node in mintree.nodes if not 'ctg' in node]
      chrs=[re.split('[:_]', node)[0] for node in mintree.nodes if not 'ctg' in node]
      chr=Counter(chrs).most_common(1)[0][0]
      nodes.sort()
      dists=[]
      ones=[]
      for kk in range(len(nodes)):
        dists.append(nodes[kk]-nodes[0]+1)
        ones.append(1)
        dstr=','.join(map(str, dists))
        onestr=','.join(map(str,ones))
      print(chr+'\t'+str(nodes[0])+'\t'+str(max(nodes))+'\t'+id+'\t100\t.\t'+str(nodes[0])+'\t'+str(max(nodes))+'\t0,150,0\t'+str(len(dists))+'\t'+onestr+'\t'+dstr, file=outb)

  ll=[loc2comp, comp2tree, loc2comphom, comp2treehom]
  with open("../hom.p", 'wb') as f:
    pickle.dump(ll, f)

  code.interact(local=locals())


    
  Gmix=nx.Graph(); superg=nx.Graph()

  with gzip.open(args.infile, 'rt') as fp:
    line=fp.readline().strip()
    ct=0
    while line:
      if ct%1000==0:
        sys.stderr.write(str(ct)+'\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom)+'\n')
      ll=re.split('[\t]', line)
      cts=list(map(int, ll[2:6]))
      mns=list(map(float, ll[6:10]))
      tot=sum(cts)
      [loc1, loc2]=ll[0:2] 
      if (loc1 in homvar and loc2 not in homvar) or (loc2 in homvar and loc1 not in homvar):
        tp='het-hom'
        if pre_filter_strict_pass(cts, 0.95, tp):
          conf=True
          if not loc2 in homvar :
            [loc1, loc2]=[loc2, loc1]
            cts=[cts[0], cts[2], cts[1], cts[3]]
            mns=[mns[0], mns[2], mns[1], mns[3]]
          if loc1 in loc2comp and loc2 in loc2comphom:
            hetcomp=loc2comp[loc1]; homcomp=loc2comphom[loc2]
            Gmix.add_edge(loc1, loc2,  conf=conf, cts=cts, mns=mns, wt=tot)
      line=fp.readline().strip()
      ct+=1


  superg=nx.Graph()
  for edge in Gmix.edges(data=True):
    loc1=edge[0]; loc2=edge[1];
    if loc1 in loc2comp and loc2 in loc2comphom:
      hetcomp=loc2comp[loc1]; homcomp=loc2comphom[loc2];
      node1='het_'+str(hetcomp); node2='hom_'+str(homcomp)
      num=0; denom=0
      for ii in range(4):
        if edge[2]['cts'][ii]>0:
          num+=edge[2]['mns'][ii]
          denom+=1
      mn_dist=int(1.0*num/denom)
      if not superg.has_edge(node1, node2):
        superg.add_edge(node1, node2, dist=[], wt=0, ct=0, sum_dist=0)
      superg.edges[node1, node2]['wt']+=edge[2]['wt']
      superg.edges[node1, node2]['dist'].append(mn_dist)
      superg.edges[node1, node2]['wt']+=1
      superg.edges[node1, node2]['sum_dist']+=mn_dist

  for edge in Gmix.edges(data=True):
    num=0; denom=0
    for ii in range(4):
      if edge[2]['cts'][ii]>0:
        num+=edge[2]['mns'][ii]
        denom+=1
    edge[2]['dist']=1.0*num/denom

  gg=list(superg.subgraph(cc) for cc in sorted(nx.connected_components(superg), key=len, reverse=True))


with gzip.open(args.bedfile, 'wt') as outb:
    for comp in comp2treehom.keys():
      id='hom_'+str(comp)
      mintree=comp2treehom[comp]
      #nodes=[int(re.split('[:_]', node)[1]) for node in mintree.nodes]                                                                                                              
      #chrs=[re.split('[:_]', node)[0] for node in mintree.nodes]                                                                                                                    
      nodes=[int(re.split('[:_]', node)[1]) for node in mintree.nodes if not 'ctg' in node]
      chrs=[re.split('[:_]', node)[0] for node in mintree.nodes if not 'ctg' in node]
      chr=Counter(chrs).most_common(1)[0][0]
      nodes.sort()
      dists=[]
      ones=[]
      for kk in range(len(nodes)):
        dists.append(nodes[kk]-nodes[0]+1)
        ones.append(1)
        dstr=','.join(map(str, dists))
        onestr=','.join(map(str,ones))
      print(chr+'\t'+str(nodes[0])+'\t'+str(max(nodes))+'\t'+id+'\t100\t.\t'+str(nodes[0])+'\t'+str(max(nodes))+'\t150,150,0\t'+str(len(dists))+'\t'+onestr+'\t'+dstr, file=outb)

with gzip.open(args.bedfile, 'wt') as outb:
    for comp in comp2treehom.keys():
      id='hom_'+str(comp)
      mintree=comp2treehom[comp]
      #nodes=[int(re.split('[:_]', node)[1]) for node in mintree.nodes]                                                                                                              
      #chrs=[re.split('[:_]', node)[0] for node in mintree.nodes]                                                                                                                    

Gloose=nx.Graph()


with gzip.open(args.infile, 'rt') as fp:
    line=fp.readline().strip()
    ct=0
    while line:
      if ct%1000==0:
        sys.stderr.write(str(ct)+'\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom)+'\n')
      ll=re.split('[\t]', line)
      cts=list(map(int, ll[2:6])); mns=list(map(float, ll[6:10]))
      tot=sum(cts)
      [loc1, loc2]=ll[0:2]
      tp='het-hom'
      if not (loc1 in homvar or loc2 in homvar):
        tp='het-het'
        if loc2<loc1:
          [loc1, loc2]=[loc2, loc1]
          cts=[cts[0], cts[2], cts[1], cts[3]]; mns=[mns[0], mns[2], mns[1], mns[3]]
        conf=False
        if pre_filter_loose_pass(cts, 0.90, 1, tp) and not Gloc.has_edge(loc1, loc2):
          Gloose.add_edge(loc1, loc2,  conf=conf, cts=cts, mns=mns, wt=tot)
      line=fp.readline().strip()
      ct+=1

comp2supercomp={}
for ii in range(len(gg)):
  for node in gg[ii].nodes():
    comp2supercomp[node]=ii

edges=[]
for edge in Gloose.edges():
  if edge[0] in loc2comp and edge[1] in loc2comp and not loc2comp[edge[0]]==loc2comp[edge[1]]:
    comp0='het_'+str(loc2comp[edge[0]])
    comp1='het_'+str(loc2comp[edge[1]])
    if comp0 in comp2supercomp and comp1 in comp2supercomp:
      s0=comp2supercomp[comp0]
      s1=comp2supercomp[comp1]
      if s0==s1:
        edges.append(edge)
        #print(str(edge))
        spl=nx.shortest_path_length(supercomp2tree[s0], source=edge[0], target=edge[1], weight='dist')
        print(comp0+'\t'+comp1+'\t'+str(spl)+'\t'+str( Gloose.edges[edge[0], edge[1]]))


supercomp2tree={}
with gzip.open('../mixed.chr10.bed.gz', 'wt') as outb:
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
    supercomp2tree[ii]=mintree
    nodes=[int(re.split('[:_]', node)[1]) for node in mintree.nodes if not 'ctg' in node]
    chrs=[re.split('[:_]', node)[0] for node in mintree.nodes if not 'ctg' in node]
    chr=Counter(chrs).most_common(1)[0][0]
    nodes.sort()
    dists=[]; ones=[]
    for kk in range(len(nodes)):
      dists.append(nodes[kk]-nodes[0]+1)
      ones.append(1)
    dstr=','.join(map(str, dists))
    onestr=','.join(map(str,ones))
    print(chr+'\t'+str(nodes[0])+'\t'+str(max(nodes))+'\t'+id+'\t100\t.\t'+str(nodes[0])+'\t'+str(max(nodes))+'\t150,150,0\t'+str(len(dists))+'\t'+onestr+'\t'+dstr, file=outb)




      with open('../graph.temp.txt', 'w') as fp:
        for edge in tr.edges():
          print(edge[0]+'\t'+edge[1], file=fp)

  with open('../nodes.temp.txt', 'w') as fp:
    for node in tr.nodes:
      tp=tr.nodes[node]['tp']
      if tp=='het':
        comp=loc2comp[node]
      else:
        comp=loc2comphom[node]
      print(node+'\t'+tp+'\t'+tp+'_'+str(comp), file=fp)
        
      bfs=list(nx.bfs_edges(tr, list(tr.nodes())[0]))
      with open('../graph.temp.txt', 'w') as fp:
        for ii in range(len(bfs)):
          edata=tr.edges[bfs[ii]]
          [prevnode, curnode]=bfs[ii]
          print(prevnode+'\t'+curnode, file=fp)
  


#for edge in superg.edges(data=True):
#  edge[2]['min_dist']=min(superg.edges[node1, node2]['dist'])
#  
#  for ii in comp2tree.keys():
#    phase_conf_component_tree(comp2tree[ii])        



  
  
#
#
#    haps=phase_conf_component_tree(gg[jj])
#    bad_edges=[(u,v,e) for u,v,e in Gloc.subgraph(gg[jj].nodes).edges(data=True) if  e['conf'] == False]
#    for kk in range(len(bad_edges)):
#      edge=bad_edges[kk]
#      if gg[jj].has_node(edge[0]) and gg[jj].has_node(edge[1]):
#        if check_edge(edge, gg[jj]):
#          Gloc.edges[edge[0], edge[1]].update({'conf': True})
#    ec=check_phase(gg[jj])
#    if ec>0:
#      sys.stderr.write('phase error\n')
#      sys.exit(1)
#    for loc in list(gg[jj].nodes):
#      Gloc.nodes[loc]['phased_all']=gg_conf.nodes[loc]['phased_all']
#      loc2comp[loc]=compid
#    tr=nx.minimum_spanning_tree(gg[jj], weight='wt')
#    comp2haps[compid]=tr
#    compid+=1
#
#
#
#  superg=nx.Graph()
#
##hic
##  with gzip.open(args.infile2, 'rt') as fp:
##    line=fp.readline().strip()
##    while line:
##      ll=re.split('[\t]', line)
##      cts=tuple(map(int, ll[2:6]))
##      locs=list(map(str, ll[0:2]))
##      if (cts[0]==0 and cts[3]==0) or (cts[1]==0 and cts[2]==0):
##        if (locs[0] in loc2comp) and (locs[1] in loc2comp):
##          add_edge_to_hic(superg, locs, cts, loc2comp, Gloc)
##          #code.interact(local=locals())
##      line=fp.readline().strip()
#
#  with gzip.open(args.infile, 'rt') as fp:
#    line=fp.readline().strip()
#    while line:
#      ll=re.split('[\t]', line)
#      cts=tuple(map(int, ll[2:6]))
#      locs=list(map(str, ll[0:2]))
#      mns=list(map(float, ll[6:10]))
#      if (cts[0]==0 and cts[3]==0) or (cts[1]==0 and cts[2]==0):
#        if (locs[0] in loc2comp) and (locs[1] in loc2comp):
#          add_edge_to_supergraph_hifi_ont(superg, locs, cts, mns, loc2comp, Gloc, 'hifi')
#      line=fp.readline().strip()
#
#
#  with gzip.open(args.infile2, 'rt') as fp:
#    line=fp.readline().strip()
#    while line:
#      ll=re.split('[\t]', line)
#      cts=tuple(map(int, ll[2:6]))
#      locs=list(map(str, ll[0:2]))
#      mns=list(map(float, ll[6:10]))
#      if (cts[0]==0 and cts[3]==0) or (cts[1]==0 and cts[2]==0):
#        if (locs[0] in loc2comp) and (locs[1] in loc2comp):
#          add_edge_to_supergraph_hifi_ont(superg, locs, cts, mns, loc2comp, Gloc, 'ont')
#      line=fp.readline().strip()
#
#  code.interact(local=locals())
#
#  for edge in superg.edges(data=True):
#    for dt in ['hifi', 'ont']:
#      curedge=edge[2][dt]
#      for ii in range(4):
#        if curedge['cts'][ii]>0:
#          curedge['dists'][ii]=curedge['dists'][ii]/curedge['cts'][ii]
#    mm1=min(edge[2]['hifi']['dists'])
#    mm2=min(edge[2]['ont']['dists'])
#    if mm1>0:
#      edge[2]['min_dist']=mm1
#    else:
#      edge[2]['min_dist']=mm2
#        
#
#  #selected_edges3 =  [(u,v) for u,v,e in superg.edges(data=True) if  ( e['pairct'] > 1 and ((e['cts'][0]==0 and e['cts'][3]==0) or (e['cts'][1]==0 and e['cts'][2]==0)))]
#  #selected_edges4 =  [(u,v) for u,v,e in superg.edges(data=True) if  ( e['pairct'] > 1 and ((e['cts'][0]==0 and e['cts'][3]==0 and e['cts'][1]*e['cts'][2]>0) or (e['cts'][1]==0 and e['cts'][2]==0 and e['cts'][0]*e['cts'][3]>0)))]
#  
#  seb = [(u,v) for u,v,e in superg.edges(data=True) if  e['hifi']['pairct']>0 and e['ont']['pairct']>0 and ((e['hifi']['outer']==e['hifi']['pairct'] and e['ont']['outer']>0.95*e['ont']['pairct']) or (e['hifi']['inner']==e['hifi']['pairct'] and e['ont']['inner']>0.95*e['ont']['pairct'])) and len(e['ont']['loci'][1])>1 and len(e['ont']['loci'][0])>1]
#  seb = [(u,v) for u,v,e in superg.edges(data=True) if  e['hifi']['pairct']>=0 and e['ont']['pairct']>0 and ((e['hifi']['outer']==e['hifi']['pairct'] and e['ont']['outer']>0.95*e['ont']['pairct']) or (e['hifi']['inner']==e['hifi']['pairct'] and e['ont']['inner']>0.95*e['ont']['pairct'])) and len(e['ont']['loci'][1])>3 and len(e['ont']['loci'][0])>3]
#
#
#  seb_loose = [(u,v) for u,v,e in superg.edges(data=True) if  e['hifi']['pairct']>=0 and e['ont']['pairct']>0 and ((e['hifi']['outer']>=0.95*e['hifi']['pairct'] and e['ont']['outer']>0.9*e['ont']['pairct']) or (e['hifi']['inner']>=0.95*e['hifi']['pairct'] and e['ont']['inner']>0.9*e['ont']['pairct'])) and len(e['ont']['loci'][1])>3 and len(e['ont']['loci'][0])>3]
#
#  superg2=superg.edge_subgraph(seb_loose).copy()
#
#  #superg2=superg.edge_subgraph(seb)
#
#  for edge in superg2.edges(data=True):
#    scale=10
#    edge[2]['cts']=[0,0,0,0]
#    edge[2]['conf']=True
#    for ii in range(4):
#      edge[2]['cts'][ii]=edge[2]['ont']['cts'][ii]+scale*edge[2]['hifi']['cts'][ii]
#  
#    
#  
#  super_comps=list(superg2.subgraph(c).copy() for c in sorted(nx.connected_components(superg2), key=len, reverse=True))
#  for ii in range(len(super_comps)):
#    hh=phase_conf_component_tree(super_comps[ii])
#    ct=check_phase(super_comps[ii])
#    print(str(ii)+'\t'+str(ct))
#
#
#  tr=nx.minimum_spanning_tree(superg2, weight='min_dist')
#  #seb = [(u,v,e) for u,v,e in superg.edges(data=True) if  e['hifi']['pairct']>0 and len(e['hifi']['loci'][0])>1 and len(e['hifi']['loci'][1])>1 and e['ont']['pairct']==0)]
#  
##  selected_edges_ont = [(u,v) for u,v,e in superg.edges(data=True) if  ( e['pairct'] > 2 and (e['inner']*1.0/e['pairct']>0.85 or e['outer']*1.0/e['pairct']>0.85))]
##  selected_edges_hifi = [(u,v) for u,v,e in superg.edges(data=True) if  ( e['pairct'] > 2 and ((e['inner']*1.0/e['pairct']>0.85 and e['cts'][1]>0.1*e['inner'] and e['cts'][1]<0.9*e['inner']) or ( e['outer']*1.0/e['pairct']>0.85 and e['cts'][0]>0.1*e['outer']and e['cts'][0]<0.9*e['outer'])))] 
##
##  for edge in selected_edges_ont:
##      node1=edge[0]
##      node2=edge[1]
##      wt=(superg.edges[edge]['pairct']+0.0)/(comp2haps[node1].number_of_nodes()*comp2haps[node2].number_of_nodes())
##      superg.edges[edge]['conf']=True
##      superg.edges[edge]['wt']=wt
##      #nn1=str(node1)+'_'+str(comp2haps[node1].number_of_nodes())
##      #nn2=str(node2)+'_'+str(comp2haps[node2].number_of_nodes())
##
##  superg2=superg.edge_subgraph(selected_edges_ont).copy()
##  super_comps=super_comps=list(superg2.subgraph(c) for c in sorted(nx.connected_components(superg2), key=len, reverse=True))
##  #for edge in selected_edges3:
##  #  print(str(superg.edges[edge]['cts']))
##  
##  with gzip.open(args.edgefile, 'wt') as outf:
##    for node in superg.nodes():
##      nn=str(node)+'_'+str(comp2haps[node].number_of_nodes())
##      print(nn, file=outf)
##    for edge in selected_edges5:
##      node1=edge[0]
##      node2=edge[1]
##      superg.edges[edge]['conf']=True
##      nn1=str(node1)+'_'+str(comp2haps[node1].number_of_nodes())
##      nn2=str(node2)+'_'+str(comp2haps[node2].number_of_nodes())
##      print(nn1+'\t'+nn2,file=outf)
##
##
##  superg2=superg.edge_subgraph(selected_edges4).copy()
##  #for  edge in selected_edges4:
##  #  superg2.add_edge(edge[0], edge[1])
##  for node in superg.nodes():
##    superg2.add_node(node)
##
##  tr=nx.maximum_spanning_tree(superg2, weight='pairct')
##  super_comps=list(tr.subgraph(c) for c in sorted(nx.connected_components(tr), key=len, reverse=True))
##  for ii in range(len(super_comps)):
##    gg=super_comps[ii]
##    haps=phase_conf_component_tree(gg)
##    for node in gg.nodes():
##      superg2.nodes[node]['phased_all']=gg.nodes[node]['phased_all']
##    check_phase(gg)
##
##    
##  sys.exit(0)
##
##  with gzip.open(args.edgefile, 'wt') as outf, gzip.open(args.hapfile, 'wt') as outh, gzip.open(args.bedfile, 'wt') as outb, gzip.open(args.compfile, 'wt') as outc:
##
##
###        
##
##    #add_to_supergraph1(superg, loc2comp, Gloc, args.infile, id='hifi')
##    add_to_supergraph1(superg, loc2comp, Gloc, args.infile2, id='hic')
##    sys.stderr.write('post hash\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom)+'\n')
##    code.interact(local=locals()) 
##    edgelist=list(superg.edges(data=True))
##    for ii in range(len(edgelist)):
##      score_edge(edgelist[ii])
##      
##    selected_edges2 =  [(u,v) for u,v,e in superg.edges(data=True) if  e['conf1'] == 2]
##    sconf=superg.edge_subgraph(selected_edges2).copy()
##    orphans = [ n for n in superg.nodes() if not sconf.has_node(n)]
##    for node in orphans:
##      sconf.add_node(node)
##    super_comps=list(sconf.subgraph(c) for c in sorted(nx.connected_components(sconf), key=len, reverse=True))
##    
##    superg2=nx.Graph()
##    comp2super={}
##    for ii in range(len(super_comps)):
##      gg_loose=super_comps[ii]
##      haps=phase_conf_component_tree(gg_loose)
##      for node in super_comps[ii].nodes():
##        comp2super[node]=ii
##
##    selected_edges =  [(u,v) for u,v,e in superg.edges(data=True) if  e['conf1'] == 1]
##    sconf2=superg.edge_subgraph(selected_edges).copy()
##    orphans = [ n for n in superg.nodes() if not sconf2.has_node(n)]
##    for node in orphans:
##      sconf2.add_node(node)
##    edgelist=list(sconf2.edges(data=True))
##    for edge in edgelist:
##      add_edge_to_supergraph1(edge, superg2, comp2super, sconf)
##
##
##    edgelist=list(superg2.edges(data=True))
##    for ii in range(len(edgelist)):
##      score_edge(edgelist[ii])
##
##    selected_edges3 =  [(u,v) for u,v,e in superg2.edges(data=True) if  e['conf1'] > 0 ]
##
##    for edge in selected_edges2:
##      [loc1, loc2]=superg.edges[edge]['loc']
##      if loc2<loc1:
##        [loc1, loc2]=[loc2, loc1]
##      [ph1, ph2]=[Gloc.nodes[loc1]['phased_all'], Gloc.nodes[loc2]['phased_all']]
##      cts=superg.edges[edge]['cts']
##      if ph1==[1,0] and ph2==[1,0]:
##        pass
##      elif ph1==[0,1] and ph2==[0,1]:
##        cts.reverse()
##      elif ph1==[1,0] and ph2==[0,1]:
##        cts=[cts[1], cts[0], cts[3], cts[2]]
##      elif ph1==[0,1] and ph2==[1,0]:
##        cts=[cts[2], cts[3], cts[0], cts[1]]
##      Gloc.add_edge(loc1, loc2, cts=cts, wt=1/sum(cts), conf=True)
##
##    for edge in selected_edges3:
##      [loc1, loc2]=superg2.edges[edge]['loc']
##      if loc2<loc1:
##        [loc1, loc2]=[loc2, loc1]
##      [ph1, ph2]=[Gloc.nodes[loc1]['phased_all'], Gloc.nodes[loc2]['phased_all']]
##      cts=superg2.edges[edge]['cts']
##      if ph1==[1,0] and ph2==[1,0]:
##        pass
##      elif ph1==[0,1] and ph2==[0,1]:
##        cts.reverse()
##      elif ph1==[1,0] and ph2==[0,1]:
##        cts=[cts[1], cts[0], cts[3], cts[2]]
##      elif ph1==[0,1] and ph2==[1,0]:
##        cts=[cts[2], cts[3], cts[0], cts[1]]
##      Gloc.add_edge(loc1, loc2, cts=cts, wt=1/sum(cts), conf=True)
##      
##
##
##    selected_edges = [(u,v) for u,v,e in Gloc.edges(data=True) if  e['conf'] == True]
##    gg_conf=Gloc.edge_subgraph(selected_edges).copy()
##    orphans=[ n for n in Gloc.nodes() if not gg_conf.has_node(n)]
##    for node in orphans:
##      gg_conf.add_node(node)
##    gg=list(gg_conf.subgraph(cc) for cc in sorted(nx.connected_components(gg_conf), key=len, reverse=True))
##    for jj in range(len(gg)):
##        haps=phase_conf_component_tree(gg[jj])
##        Gall=nx.Graph()
##        for edge in list(gg[jj].edges()):
##          [loc0, loc1]=edge
##          if loc1<loc2:
##            [loc0, loc1]=[loc1, loc0]
##          curedge=gg[jj].edges[edge]
##          add_allele_edges(Gall, loc0, loc1, curedge['cts'], curedge['mns'])
##          print(str(ii)+'\t'+str(jj)+'\t'+edge[0]+'\t'+edge[1]+'\t'+str(curedge['cts']), file=outc)
##        if gg[jj].number_of_edges()==0:
##          singleton=list(gg[jj].nodes())[0]
##          Gall.add_node(singleton+'_0')
##          Gall.add_node(singleton+'_1')
##        for hapid in range(2):
##          minforest=nx.minimum_spanning_tree(Gall.subgraph(haps[hapid]), weight='mean_dist')
##          mintree=list(minforest.subgraph(cc) for cc in nx.connected_components(minforest))
##          for treeid in range(len(mintree)):
##            id=(str(ii)+'_'+str(jj)+'_'+str(hapid)+'_'+str(treeid))
##            treelist=list(mintree[treeid].edges)
##            nodes=[int(re.split('[:_]', node)[1]) for node in mintree[treeid].nodes]
##            chrs=[re.split('[:_]', node)[0] for node in mintree[treeid].nodes]
##            chr=Counter(chrs).most_common(1)[0][0]
##            nodes.sort()
##            dists=[]
##            ones=[]
##            for kk in range(len(nodes)):
##              dists.append(nodes[kk]-nodes[0]+1)
##              ones.append(1)
##              dstr=','.join(map(str, dists))
##              onestr=','.join(map(str,ones))
##            print(chr+'\t'+str(nodes[0])+'\t'+str(max(nodes))+'\t'+id+'\t100\t.\t'+str(nodes[0])+'\t'+str(max(nodes))+'\t150,150,0\t'+str(len(dists))+'\t'+onestr+'\t'+dstr, file=outb)
##            for tredge in treelist:
##              curedge=mintree[treeid].edges[tredge]
##              edgestr=str(round(curedge['mean_dist'], 3))+';'+str(round(curedge['sd_dist'], 3))+';'+str(curedge['ct'])
##              print(id+'\t'+edgestr+'\t'+tredge[0]+'\t'+tredge[1], file=outf)
##            terminal_nodes=[]
##            if len(list(mintree[treeid].nodes))>1:
##              for node1 in list(mintree[treeid].nodes):
##                if mintree[treeid].degree(node1)==1:
##                  terminal_nodes.append(node1)
##              for aa in range(len(terminal_nodes)-1):
##                for bb in range(aa+1, len(terminal_nodes)):
##                  sp=nx.shortest_path(mintree[treeid], terminal_nodes[aa], terminal_nodes[bb])
##                  for nodeii in range(len(sp)):
##                    node1=sp[nodeii]
##                    if nodeii==0:
##                      outstr=str(mintree[treeid].degree(node1))+'\t.'
##                    else:
##                      prevedge=mintree[treeid].edges[sp[nodeii], sp[nodeii-1]]
##                      outstr=str(mintree[treeid].degree(node1))+'\t'+str(prevedge['ct'])+'_'+str(prevedge['mean_dist'])+'_'+str(prevedge['sd_dist'])
##                    print(id+'\t'+str(aa)+'_'+str(bb)+'\t'+outstr+'\t'+node1, file=outh)
##            else:
##              node1=list(mintree[treeid].nodes)[0]
##              print(id+'\t0_0\t0\t0_0_0\t'+node1, file=outh)
#                  
#    
#      
##  cp.disable()
##y  cp.print_stats()
                                            

                                           
                                             

if __name__ == "__main__":
    main()
