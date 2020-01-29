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


def pre_filter_strict_pass(cts, mns, minfrac=0.95, tp='het-het'):
  nn=sum(cts)
  passf=False; bal=-1; orient=''; dist=-1
  if nn>=2:
    thresh=minfrac*nn
    if tp=='het-het':
      ors=[[0,3], [1,2]]
      for ort in ors:
        [aa, bb]=ort
        if cts[aa]+cts[bb]>thresh:
          orient=str(aa)+'_'+str(bb)
          bal=cts[aa]*1.0/(cts[aa]+cts[bb])
          dist=0.5*(mns[aa]+mns[bb])
          break
      if bal>0.1 and bal<0.9:
        passf=True
    elif tp=='het-hom':
      ors=[[0,1], [0,2], [2,3], [1,3]]
      for ort in ors:
        [aa, bb]=ort
        if cts[aa]+cts[bb]>thresh:
          orient=str(aa)+'_'+str(bb)
          bal=cts[aa]*1.0/(cts[aa]+cts[bb])
          dist=0.5*(mns[aa]+mns[bb])
          break
      if bal>0.1 and bal<0.9:
        passf=True
    elif tp=='hom-hom':
      mm=max(cts)
      for ii in range(4):
        if cts[ii]==mm:
          orient=str(ii)
          dist=mns[ii]
      if mm>thresh:
        passf=True
  return [passf, orient, nn, dist]


def pre_filter_loose_pass(cts, mns, minfrac=0.90, minct=2, tp='het-het'):
  nn=sum(cts)
  passf=False; orient=''; dist=-1
  if nn>=minct:
    thresh=minfrac*nn
    if tp=='het-het':
      ors=[[0,3], [1,2]]
      for ort in ors:
        [aa, bb]=ort
        if cts[aa]+cts[bb]>thresh:
          orient=str(aa)+'_'+str(bb)
          dist=0.5*(mns[aa]+mns[bb])
          passf=True
          break
    elif tp=='het-hom':
      ors=[[0,1], [0,2], [2,3], [1,3]]
      for ort in ors:
        [aa, bb]=ort
        if cts[aa]+cts[bb]>thresh:
          orient=str(aa)+'_'+str(bb)
          dist=0.5*(mns[aa]+mns[bb])
          passf=True
          break
    elif tp=='hom-hom':
      mm=max(cts)
      for ii in range(4):
        if cts[ii]==mm:
          orient=str(ii)
          dist=mns[ii]
      if mm>thresh:
        passf=True
  return [passf, orient, nn, dist]


def comp2bed(trnodes, id, col):  
  nodes=[int(re.split('[:_]', node)[1]) for node in trnodes if not 'ctg' in node]
  if len(nodes)>0:
    chrs=[re.split('[:_]', node)[0] for node in trnodes if not 'ctg' in node]
    chr=Counter(chrs).most_common(1)[0][0]
    nodes.sort()
    dists=[]
    ones=[]
    for kk in range(len(nodes)):
      dists.append(nodes[kk]-nodes[0]+1)
      ones.append(1)
    dstr=','.join(map(str, dists))
    onestr=','.join(map(str,ones))
    outstr=chr+'\t'+str(nodes[0])+'\t'+str(max(nodes))+'\t'+id+'\t100\t.\t'+str(nodes[0])+'\t'+str(max(nodes))+'\t'+col+'\t'+str(len(dists))+'\t'+onestr+'\t'+dstr
    return outstr


def odds_ratio(ll):
  [aa, bb, cc, dd]=ll
  return 1.0*(aa+1)*(dd+1)/((bb+1)*(cc+1))

def check_edge(edge, gg):
  oddsratio=odds_ratio(edge[2]['cts'])
  ph0=gg.nodes[edge[0]]['phased_all']
  ph1=gg.nodes[edge[1]]['phased_all']
  if (oddsratio>1 and ph0==ph1) or (oddsratio<1 and ph1!=ph0):
    return True
  else:
    return False

def check_phase(ggsub):
  edgelist=list(ggsub.edges(data=True))
  error_count=0
  for ii in range(len(edgelist)):
    edata=edgelist[ii]
    if edata[2]['conf']:
      if not check_edge(edata, ggsub):
        sys.stderr.write('Error')
        error_count+=1
  return error_count


def check_phase1(ggsub):
  edgelist=list(ggsub.edges(data=True))
  error_count=0
  bad_edges=[]
  for ii in range(len(edgelist)):
    edata=edgelist[ii]
    if edata[2]['conf']:
      if not check_edge(edata, ggsub):
        sys.stderr.write('Error')
        error_count+=1
        bad_edges.append(edata)
  return bad_edges

def phase_conf_component_tree(ggsub):
  tr=nx.minimum_spanning_tree(ggsub, weight='wt')
  startnode=''
  for node in tr.nodes:
    if tr.degree(node)<=1:
      startnode=node
      break
  tr.nodes[startnode]['phased_all']=[0,1]
  ggsub.nodes[startnode]['phased_all']=[0,1]
  [h0, h1] = phase_from_node_tree(ggsub, tr, startnode)
  return([h0, h1])


def phase_conf_component_tree2(ggsub):
  tr=nx.minimum_spanning_tree(ggsub, weight='wt')
  startnode=''
  for node in tr.nodes:
    if tr.degree(node)<=1:
      startnode=node
      break
  tr.nodes[startnode]['phased_all']=[0,1]
  ggsub.nodes[startnode]['phased_all']=[0,1]
  [h0, h1] = phase_from_node_tree(ggsub, tr, startnode)
  return([h0, h1], tr, startnode)

def phase_from_node_tree(gg1, tree, startnode):
  hap0=[]
  hap1=[]
  hap0.append(get_phased_allele(gg1, startnode, 0))
  hap1.append(get_phased_allele(gg1, startnode, 1))
  if len(gg1)>1:
    bfs=list(nx.bfs_edges(tree, startnode))
    for ii in range(len(bfs)):
      edata=gg1.edges[bfs[ii]]
      [prevnode, curnode]=bfs[ii]
      if 'phased_all' in gg1.nodes[prevnode]:
        alleles=gg1.nodes[prevnode]['phased_all']
        oddsratio=odds_ratio(edata['cts'])
        if oddsratio>1:
          gg1.nodes[curnode]['phased_all']=alleles
        else:
          gg1.nodes[curnode]['phased_all']=[alleles[1], alleles[0]] 
        hap0.append(get_phased_allele(gg1, curnode, 0))
        hap1.append(get_phased_allele(gg1, curnode, 1))
  return([hap0, hap1])      

def get_phased_allele(gg, node, hap):
  refalt=gg.nodes[node]['phased_all'][hap]
  return(str(node)+'_'+str(refalt))

def add_allele_edges(gg, loc1, loc2, cts, mns):
    [aa, bb, cc, dd]=cts
    oddsratio=(aa+1.0)*(dd+1.0)/((bb+1.0)*(cc+1.0))
    nn=sum(cts)
    if oddsratio>1:
      if aa>0:
        gg.add_edge(loc1+'_0', loc2+'_0', ct=aa, mean_dist=mns[0])
      if dd>0:
        gg.add_edge(loc1+'_1', loc2+'_1', ct=dd, mean_dist=mns[3])
    elif oddsratio<1:
      if bb>0:
        gg.add_edge(loc1+'_0', loc2+'_1', ct=bb, mean_dist=mns[1])
      if cc>0:
        gg.add_edge(loc1+'_1', loc2+'_0', ct=cc, mean_dist=mns[2])

def add_to_supergraph1(superg, comphash, Gconf, infile, id):

  with gzip.open(infile, 'rt') as fp:
    line=fp.readline().strip()
    while line:
      ll=re.split('[\t]', line)
      cts=tuple(map(int, ll[2:6]))
      mns=tuple(map(float, ll[6:10]))
      locs=list(map(str, ll[0:2]))
      if pre_filter_loose_pass(cts):
        add_edge_to_supergraph(superg, locs, cts, mns, comphash, Gconf)
      line=fp.readline().strip()


def add_edge_to_hic(superg, locs, cts, comphash, Gconf):
  [loc1, loc2]=locs
  [c1, c2]=[comphash[loc1], comphash[loc2]]
  if c1!=c2:
    if c2<c1:
      [loc1, loc2]=[loc2, loc1]
      [c1, c2]=[c2, c1]
      cts=(cts[0], cts[2], cts[1], cts[3])
    [ph1, ph2]=Gconf.nodes[loc1]['phased_all'], Gconf.nodes[loc2]['phased_all']
    if ph1==[1,0] and ph2==[1,0]:
      trans=1
    elif ph1==[0,1] and ph2==[0,1]:
      trans=2
    elif ph1==[1,0] and ph2==[0,1]:
      trans=3
    elif ph1==[0,1] and ph2==[1,0]:
      trans=3
    cts=transform(cts, trans)

    if not superg.has_edge(c1, c2):
      superg.add_edge(c1, c2, cts=[0,0,0,0], cts_all=[], loc_pairs=[], evtype='hi', pairct=0)
    curedge=superg.edges[c1, c2]
    curedge['pairct']+=1
    curedge['cts_all'].append(cts)
    curedge['loc_pairs'].append([loc1, loc2])
    for ii in range(4):
      curedge['cts'][ii]+=cts[ii]

def add_edge_to_supergraph_hifi_ont(superg, locs, cts, mns, comphash, Gconf, dt):
  [loc1, loc2]=locs
  [c1, c2]=[comphash[loc1], comphash[loc2]]
  if c1!=c2:
    if c2<c1:
      [loc1, loc2]=[loc2, loc1]
      [c1, c2]=[c2, c1]
      cts=(cts[0], cts[2], cts[1], cts[3])
      mns=(mns[0], mns[2], mns[1], mns[3])
    [ph1, ph2]=Gconf.nodes[loc1]['phased_all'], Gconf.nodes[loc2]['phased_all']
    if ph1==[1,0] and ph2==[1,0]:
      trans=1
    elif ph1==[0,1] and ph2==[0,1]:
      trans=2
    elif ph1==[1,0] and ph2==[0,1]:
      trans=3
    elif ph1==[0,1] and ph2==[1,0]:
      trans=3
    cts=transform(cts, trans)
    mns=transform(mns, trans)
    if not superg.has_edge(c1, c2):
      superg.add_edge(c1, c2, 
       hifi={'cts' : [0,0,0,0], 'dists' : [0, 0, 0, 0], 'outer': 0, 'inner' : 0, 'max_outer' : 0, 'max_inner':0, 'cts_all': [], 'loci' : [set(), set()], 'pairct' : 0}, 
       ont={'cts' : [0,0,0,0], 'dists' : [0, 0, 0, 0], 'outer': 0, 'inner' : 0, 'max_outer' : 0, 'max_inner':0, 'cts_all': [], 'loci' : [set(), set()], 'pairct' : 0})
    curedge=superg.edges[c1, c2][dt]
    curedge['pairct']+=1
    curedge['cts_all'].append(cts)
    curedge['loci'][0].add(loc1)
    curedge['loci'][1].add(loc2)
    if cts[1]==0 and cts[2]==0:
      curedge['outer']+=1
      if cts[0]>curedge['max_outer']:
        curedge['max_outer']=cts[0]
      if cts[3]>curedge['max_outer']:
        curedge['max_outer']=cts[3]
    elif cts[0]==0 and cts[3]==0:
      curedge['inner']+=1
      if cts[1]>curedge['max_inner']:
        curedge['max_inner']=cts[1]
      if cts[2]>curedge['max_inner']:
        curedge['max_inner']=cts[2]
    for ii in range(4):
      if(cts[ii]>0):
        curedge['cts'][ii]+=1
        curedge['dists'][ii]+=abs(mns[ii])


def add_edge_to_supergraph_hifi_ont_hic(superg, locs, cts, comphash, Gconf, dt):
  [loc1, loc2]=locs
  [c1, c2]=[comphash[loc1], comphash[loc2]]
  if c1!=c2:
    if c2<c1:
      [loc1, loc2]=[loc2, loc1]
      [c1, c2]=[c2, c1]
      cts=(cts[0], cts[2], cts[1], cts[3])
    [ph1, ph2]=Gconf.nodes[loc1]['phased_all'], Gconf.nodes[loc2]['phased_all']
    if ph1==[1,0] and ph2==[1,0]:
      trans=1
    elif ph1==[0,1] and ph2==[0,1]:
      trans=2
    elif ph1==[1,0] and ph2==[0,1]:
      trans=3
    elif ph1==[0,1] and ph2==[1,0]:
      trans=3
    cts=transform(cts, trans)
    if not superg.has_edge(c1, c2):
      superg.add_edge(c1, c2, hifi={'cts' : [0,0,0,0], 'outer': 0, 'inner' : 0, 'max_outer' : 0, 'max_inner':0, 'cts_all': [], 'loci' : [set(), set()], 'pairct' : 0}, 
        ont={'cts' : [0,0,0,0], 'outer': 0, 'inner' : 0, 'max_outer' : 0, 'max_inner':0, 'cts_all': [], 'loci' : [set(), set()], 'pairct' : 0},
        hic={'cts' : [0,0,0,0], 'outer': 0, 'inner' : 0, 'max_outer' : 0, 'max_inner':0, 'cts_all': [], 'loci' : [set(), set()], 'pairct' : 0})
    curedge=superg.edges[c1, c2][dt]
    curedge['pairct']+=1
    curedge['cts_all'].append(cts)
    curedge['loci'][0].add(loc1)
    curedge['loci'][1].add(loc2)
    if cts[1]==0 and cts[2]==0:
      curedge['outer']+=1
      if cts[0]>curedge['max_outer']:
        curedge['max_outer']=cts[0]
      if cts[3]>curedge['max_outer']:
        curedge['max_outer']=cts[3]
    elif cts[0]==0 and cts[3]==0:
      curedge['inner']+=1
      if cts[1]>curedge['max_inner']:
        curedge['max_inner']=cts[1]
      if cts[2]>curedge['max_inner']:
        curedge['max_inner']=cts[2]
    for ii in range(4):
      if(cts[ii]>0):
        curedge['cts'][ii]+=1


def add_edge_to_supergraph(superg, locs, cts, mns, comphash, Gconf):
  [loc1, loc2]=locs
  [c1, c2]=[comphash[loc1], comphash[loc2]]
  if c1!=c2:
    if c2<c1:
      [loc1, loc2]=[loc2, loc1]
      [c1, c2]=[c2, c1]
      cts=(cts[0], cts[2], cts[1], cts[3])
      mns=(mns[0], mns[2], mns[1], mns[3])
    [ph1, ph2]=Gconf.nodes[loc1]['phased_all'], Gconf.nodes[loc2]['phased_all']
    if ph1==[1,0] and ph2==[1,0]:
      trans=1
    elif ph1==[0,1] and ph2==[0,1]:
      trans=2
    elif ph1==[1,0] and ph2==[0,1]:
      trans=3
    elif ph1==[0,1] and ph2==[1,0]:
      trans=3
    dist=max(mns)
    cts=transform(cts, trans)
    mns=transform(mns, trans)
    if not superg.has_edge(c1, c2):
      superg.add_edge(c1, c2, cts_all={}, min_dist={}, min_loc_pairs={}, evtype={}, pairct=0)
    curedge=superg.edges[c1, c2]
    curedge['pairct']+=1
    if not cts in curedge['cts_all']:
      curedge['cts_all'][cts]=1
      curedge['min_dist'][cts]=dist
      curedge['min_loc_pairs'][cts]=[loc1, loc2]
      curedge['evtype'][cts]=id
      curedge['mns_dict'][cts]=mns
    else:
      curedge['cts_all'][cts]+=1
      if dist<curedge['min_dist'][cts]:
        curedge['min_dist'][cts]=dist
        curedge['min_loc_pairs'][cts]=[loc1, loc2]
        curedge['evtype'][cts]=id
        curedge['mns_dict'][cts]=mns

def add_edge_to_supergraph1(edge, superg, comphash, sconf):

  [loc1, loc2]=[edge[0], edge[1]]
  [c1, c2]=[comphash[loc1], comphash[loc2]]
  if c1!=c2:
    if c2<c1:
      [loc1, loc2]=[loc2, loc1]
      [c1, c2]=[c2, c1]
    [ph1, ph2]=[sconf.nodes[loc1]['phased_all'], sconf.nodes[loc2]['phased_all']]
    if ph1==[1,0] and ph2==[1,0]:
      trans=1
    elif ph1==[0,1] and ph2==[0,1]:
      trans=2
    elif ph1==[1,0] and ph2==[0,1]:
      trans=3
    else:
      trans=4
    if not superg.has_edge(c1, c2):
      superg.add_edge(c1, c2, pairct=0, cts_all={}, min_dist={}, min_loc_pairs={}, evtype={})
    ee=superg.edges[c1, c2]
    ee['pairct']+=edge[2]['pairct']
    for cts in edge[2]['cts_all'].keys():
      ctstrans=transform(cts, trans)
      if not ctstrans in ee['cts_all']:
        ee['cts_all'][ctstrans]=edge[2]['cts_all'][cts]
        ee['min_dist'][ctstrans]=edge[2]['min_dist'][cts]
        ee['min_loc_pairs'][ctstrans]=edge[2]['min_loc_pairs'][cts]
        ee['evtype'][ctstrans]=edge[2]['evtype'][cts]
        ee['mns_dict'][ctstrans]=edge[2]['mns_dict'][cts]
      else:
        ee['cts_all'][ctstrans]=ee['cts_all'][ctstrans]+edge[2]['cts_all'][cts]
        if edge[2]['min_dist'][cts]<ee['min_dist'][ctstrans]:
          ee['min_dist'][ctstrans]=edge[2]['min_dist'][cts]
          ee['min_loc_pairs'][ctstrans]=edge[2]['min_loc_pairs'][cts]
          ee['mns_dict'][ctstrans]=edge[2]['mns_dict'][cts]

def remove_bridges(GG, min_counts_strict, tp):
  bridges=[]
  brlist=list(nx.bridges(GG))
  for edge in brlist:
    curedge=GG.edges[edge]
    if curedge['wt']<min_counts_strict:
      if tp=='hom-hom' or abs(curedge['mns'][0]-curedge['mns'][3])>350 or abs(curedge['mns'][1]-curedge['mns'][2])>350:
        bridges.append(curedge)
        GG.remove_edge(edge[0], edge[1])
  return bridges

def transform(cts, trans=1):
  
  if trans==2:
    cts=(cts[3], cts[2], cts[1], cts[0])
  elif trans==3:
    cts=(cts[1], cts[0], cts[3], cts[2])
  elif trans==4:
    cts=(cts[2], cts[3], cts[0], cts[1])
  return cts

def score_edge(edge, min_frac=0.95):
  [inner, outer, tot, maxct_inner, maxct_outer, maxa, maxb, maxc, maxd, aa, bb, cc, dd]=[0,0,0,0,0,0,0,0,0,0,0,0,0]
  [mindist_outer, mindist_inner]=[1e9, 1e9]
  [edge_inner, edge_outer]=[None, None]
  for cts in edge[2]['cts_all'].keys():
    mult=edge[2]['cts_all'][cts]
    if cts[1]+cts[2]==0:
      outer+=mult
      [aa, dd]=[aa+cts[0]*mult, dd+cts[3]*mult]
      if cts[0]+cts[3]>=maxct_outer:
        maxct_outer=cts[0]+cts[3]
        if edge[2]['min_dist'][cts]<mindist_outer:
          [mindist_outer, edge_outer]=[edge[2]['min_dist'][cts], edge[2]['min_loc_pairs'][cts]]
      if cts[0]>maxa:
        maxa=cts[0]
      if cts[3]>maxd:
        maxd=cts[3]
    elif cts[0]+cts[3]==0:
      inner+=mult
      [bb, cc]=[bb+cts[1]*mult, cc+cts[2]*mult]
      if cts[1]+cts[2]>=maxct_inner:
        maxct_inner=cts[1]+cts[2]
        if edge[2]['min_dist'][cts]<mindist_inner:
          [mindist_inner, edge_inner]=[edge[2]['min_dist'][cts], edge[2]['min_loc_pairs'][cts]]
      if cts[1]>maxb:
        maxb=cts[1]
      if cts[2]>maxc:
        maxc=cts[2]
    tot+=mult
  [frac_inner, frac_outer]=[inner/(0.0+tot), outer/(0.0+tot)]
  edge[2]['cts_io']=[inner, outer, tot]
  edge[2]['cts_max']=[maxa, maxb, maxc, maxd]
  edge[2]['cts_tot']=[aa, bb, cc, dd]
  conf=0
  tot=aa+bb+cc+dd
  if frac_inner>min_frac and bb+cc>min_frac*tot:
    edge[2]['dir']='inner'
    edge[2]['cts']=[0, maxb, maxc, 0]
    edge[2]['dist']=mindist_inner
    edge[2]['loc']=edge_inner
    conf=1
    if bb>(1-min_frac)*(bb+cc) and bb<min_frac*(bb+cc) and maxb+maxc>2:
      conf=2
  elif frac_outer>min_frac and aa+dd>min_frac*tot:
    edge[2]['dir']='outer'
    edge[2]['cts']=[maxa, 0, 0, maxd]
    edge[2]['dist']=mindist_outer
    edge[2]['loc']=edge_outer
    conf=1
    if aa>(1-min_frac)*(aa+dd) and aa<min_frac*(aa+dd) and maxa+maxd>2:
      conf=2
  edge[2]['conf1']=conf


