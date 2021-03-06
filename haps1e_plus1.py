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


def pre_filter_strict_pass(cts, minfrac=0.95):
  [aa, bb, cc, dd]=cts
  nn=aa+bb+cc+dd
  passf=False
  bal=-1
  if nn<2:
    return False
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
    return False
  elif  (aa+dd)*1.0/nn>minfrac or (bb+cc)*1.0/nn>minfrac:
    passf=True
  return passf

def phase_conf_component_shortest(ggsub):
    maxnode=list(ggsub.nodes)[0]
    maxval=-1
    for node in ggsub.nodes:
      deg=ggsub.degree(node)
      if deg>maxval:
        maxval=deg
        maxnode=node
    ggsub.nodes[maxnode]['phased_all']=[0,1]
    [h0, h1]=phase_from_node_shortest(ggsub, maxnode)
    return([h0, h1])

def phase_from_node_shortest(gg1, startnode):
  hap0=[]
  hap1=[]
  if len(gg1)==1:
    hap0.append(get_phased_allele(gg1, startnode, 0))
    hap1.append(get_phased_allele(gg1, startnode, 1))
  else:
    dp=nx.dijkstra_predecessor_and_distance(gg1, startnode, weight='tot')
    kk=list(dp[0].keys())
    for ii in range(len(kk)):
      if ii==0:
        curnode=startnode
        hap0.append(get_phased_allele(gg1, curnode, 0))
        hap1.append(get_phased_allele(gg1, curnode, 1))
      else:
        curnode=kk[ii]
        prevnode=dp[0][kk[ii]][0]
        edata=gg1.edges[curnode, prevnode]
        if 'phased_all' in gg1.nodes[prevnode]:
          alleles=gg1.nodes[prevnode]['phased_all']
          oddsratio=(edata['cts'][0]+1.0)*(edata['cts'][3]+1.0)/((edata['cts'][1]+1.0)*(edata['cts'][2]+1.0))
          if oddsratio>1:
            gg1.nodes[curnode]['phased_all']=alleles
          else:
            gg1.nodes[curnode]['phased_all']=[alleles[1], alleles[0]] 
          hap0.append(get_phased_allele(gg1, curnode, 0))
          hap1.append(get_phased_allele(gg1, curnode, 1))
  return [hap0, hap1]


def get_phased_allele(gg, node, hap):
  refalt=gg.nodes[node]['phased_all'][hap]
  return(node+'_'+str(refalt))


def add_allele_edges(gg, loc1, loc2, cts, mns, sds):
    [aa, bb, cc, dd]=cts
    oddsratio=(aa+1.0)*(dd+1.0)/((bb+1.0)*(cc+1.0))
    nn=sum(cts)
    if oddsratio>1:
      if aa>0:
        gg.add_edge(loc1+'_0', loc2+'_0', ct=aa, mean_dist=mns[0], sd_dist=sds[0])
      if dd>0:
        gg.add_edge(loc1+'_1', loc2+'_1', ct=dd, mean_dist=mns[3], sd_dist=sds[3])
    elif oddsratio<1:
      if bb>0:
        gg.add_edge(loc1+'_0', loc2+'_1', ct=bb, mean_dist=mns[1], sd_dist=sds[1])
      if cc>0:
        gg.add_edge(loc1+'_1', loc2+'_0', ct=cc, mean_dist=mns[2], sd_dist=sds[2])
                                        
def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', '--infile', type=str, default=None)
  parser.add_argument('-s', '--singletons', type=str, default=None)
  parser.add_argument('-f', '--infile2', type=str, default=None)
  parser.add_argument('-c', '--compfile', type=str, default=None)
  args = parser.parse_args()
  min_counts_strict=5

  cp = cProfile.Profile()
  cp.enable()
  usage_denom=1024
  Gloc=nx.Graph()
  infile=args.infile

  with gzip.open(infile, 'rt') as fp:
    line=fp.readline().strip()
    ct=0
    while line :
      if ct%1000==0:
        sys.stderr.write(str(ct)+'\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom)+'\n')
      ll=re.split('[\t]', line)
      cts=list(map(int, ll[2:6]))
      mns=list(map(float, ll[6:10]))
      tot=sum(cts)
      [loc1, loc2]=ll[0:2]
      if pre_filter_strict_pass(cts):
        Gloc.add_edge(loc1, loc2, order=[loc1, loc2], conf=True, cts=cts, mns=mns, wt=1.0/tot)
      elif pre_filter_loose_pass(cts):
        Gloc.add_edge(loc1, loc2,  order=[loc1, loc2], conf=False, cts=cts, mns=mns, wt=1.0/tot)
      line=fp.readline().strip()
      ct+=1


  with gzip.open(args.singletons, 'rt') as fp:
    snp=fp.readline().strip()
    if not Gloc.has_node(snp):
      Gloc.add_node(snp)


  supercomphash={}
  comphash={}
  compid=0
  sys.stderr.write('finished loading graph\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom)+'\n')
  with gzip.open(args.compfile, 'wt') as outc:
    comp_loose=list(Gloc.subgraph(c) for c in sorted(nx.connected_components(Gloc), key=len, reverse=True))
    for ii in range(len(comp_loose)):
      gg_loose=comp_loose[ii]
      selected_edges = [(u,v) for u,v,e in gg_loose.edges(data=True) if  e['conf'] == True]
      gg_conf=gg_loose.edge_subgraph(selected_edges)
      brlist=list(nx.bridges(gg_conf))
      for edge in brlist:
        curedge=gg_loose.edges[edge]
        if sum(curedge['cts'])<min_counts_strict:
          if abs(curedge['mns'][0]-curedge['mns'][3])>350 or abs(curedge['mns'][1]-curedge['mns'][2])>350:
            gg_loose.edges[edge].update({'conf': False})
      selected_edges = [(u,v) for u,v,e in gg_loose.edges(data=True) if  e['conf'] == True]
      bad_edges = [(u,v,e) for u,v,e in gg_loose.edges(data=True) if  e['conf'] == False]
      gg_conf=gg_loose.edge_subgraph(selected_edges)
      orphans=[ n for n in gg_loose.nodes() if not gg_conf.has_node(n)]
      gg=list(gg_conf.subgraph(cc) for cc in sorted(nx.connected_components(gg_conf), key=len, reverse=True))
      code.interact(local=locals())
      for jj in range(len(gg)):
        haps=phase_conf_component_shortest(gg[jj])
        for kk in range(len(bad_edges)):
          edge=bad_edges[kk]
          if gg[jj].has_node(edge[0]) and gg[jj].has_node(edge[1]):
            oddsratio=(edge[2]['cts'][0]+1.0)*(edge[2]['cts'][3]+1.0)/((edge[2]['cts'][1]+1.0)*(edge[2]['cts'][2]+1.0))
            if (oddsratio>1 and gg[jj].nodes[edge[0]]['phased_all']==gg[jj].nodes[edge[1]]['phased_all']) or (oddsratio<1 and not gg[jj].nodes[edge[0]]['phased_all']==gg[jj].nodes[edge[1]]['phased_all']):
              Gloc.edges[edge[0], edge[1]]['conf']=True
        for edge in list(gg[jj].edges()):
          oddsratio=(edge['cts'][0]+1.0)*(edge['cts'][3]+1.0)/((edge['cts'][1]+1.0)*(edge['cts'][2]+1.0))
          if oddsratio>1 and not  edge[0]['phased_all']==edge[1]['phased_all']:
            sys.stderr.write('error\n')
          elif oddsratio<1 and edge[0]['phased_all']==edge[1]['phased_all']:
            sys.stderr.write('error\n')
        for nn in list(gg[jj].nodes()):
          supercomphash[nn]=ii
          comphash[nn]=compid
        compid+=1
      for jj in range(len(orphans)):
        nn=orphans[jj]
        supercomphash[nn]=ii
        comphash[nn]=compid
        compid+=1


  code.interact(local=locals())
  superg=nx.Graph()
  sys.stderr.write('post hash\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom)+'\n')

  with gzip.open(args.infile2, 'rt') as fp:
    line=fp.readline().strip()
    ct=0
    while line:
      if ct%1000==0:
        sys.stderr.write(str(ct)+'\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom)+'\n')
        ll=re.split('[\t]', line)
        cts=list(map(int, ll[2:6]))
        tot=sum(cts)
        [loc1, loc2]=ll[0:2]
        if pre_filter_loose_pass(cts):
          conf=0
          if pre_filter_strict_pass(cts):
            conf=1
          if loc1 in comphash and loc2 in comphash:
            c1=comphash[loc1]
            c2=comphash[loc2]
            if c1!=c2:
              if c2<c1:
                [loc1, loc2]=[loc2, loc1]
                [c1, c2]=[c2, c1]
                cts.reverse()
              ph1=Gloc.nodes[loc1]['phased_all']
              ph2=Gloc.nodes[loc2]['phased_all']
              if ph1==[1,0] and ph2==[1,0]:
                pass
              elif ph1==[0,1] and ph2==[0,1]:
                cts.reverse()
              elif ph1==[1,0] and ph2==[0,1]:
                cts=[cts[1], cts[0], cts[3], cts[2]]
              elif ph1==[0,1] and ph2==[1,0]:
                cts=[cts[2], cts[3], cts[0], cts[1]]
              if not superg.has_edge(c1, c2):
                cc=[0,0]
                cc[conf]+=1
                superg.add_edge(c1, c2, cts_all=[cts], loc_pairs=[[loc1, loc2]], pairct=1, conf=cc)
              else:
                curedge=superg.edges[c1, c2]
                curedge['loc_pairs'].append([loc1, loc2])
                curedge['pairct']+=1
                curedge['cts_all'].append(cts)
                curedge['conf'][conf]+=1
      line=fp.readline().strip()
      ct+=1

  code.interact(local=locals())                                                                                                                             
                            
  cp.disable()
  cp.print_stats()
                                           
                                             

if __name__ == "__main__":
    main()
