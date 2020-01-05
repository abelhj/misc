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



                                        
def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', '--infile', type=str, default=None)
  parser.add_argument('-b', '--bedfile', type=str, default=None)
  parser.add_argument('-e', '--edgefile', type=str, default=None)
  parser.add_argument('-p', '--hapfile', type=str, default=None)
  parser.add_argument('-c', '--compfile', type=str, default=None)
  parser.add_argument('-s', '--singletons', type=str, default=None)
  args = parser.parse_args()
  min_counts_strict=5

  cp = cProfile.Profile()
  cp.enable()
  usage_denom=1024
  Gloc=nx.Graph()

  with gzip.open(args.infile, 'rt') as fp:
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
      if loc2<loc1:
        [loc1, loc2]=[loc2, loc1]
        cts=[cts[0], cts[2], cts[1], cts[3]]
        mns=[mns[0], mns[2], mns[1], mns[3]]
      if pre_filter_strict_pass(cts):
        Gloc.add_edge(loc1, loc2,  conf=True, cts=cts, mns=mns,  wt=1.0/tot)
      elif pre_filter_loose_pass(cts):
        Gloc.add_edge(loc1, loc2,  conf=False, cts=cts, mns=mns, wt=1.0/tot)
      else:
        Gloc.add_node(loc1)
        Gloc.add_node(loc2)
      line=fp.readline().strip()
      ct+=1

  with gzip.open(args.singletons, 'rt') as fp:
    line=fp.readline().strip()
    while line:
      if not Gloc.has_node(line):
        Gloc.add_node(line)
      line=fp.readline().strip()

  remove_bridges(Gloc, min_counts_strict)
  sys.stderr.write('finished loading graph\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom)+'\n')
  with gzip.open(args.edgefile, 'wt') as outf, gzip.open(args.hapfile, 'wt') as outh, gzip.open(args.bedfile, 'wt') as outb, gzip.open(args.compfile, 'wt') as outc:
    comp_loose=list(Gloc.subgraph(c) for c in sorted(nx.connected_components(Gloc), key=len, reverse=True))
    for ii in range(len(comp_loose)):
      gg_loose=comp_loose[ii]
      selected_edges = [(u,v) for u,v,e in gg_loose.edges(data=True) if  e['conf'] == True]
      bad_edges = [(u,v,e) for u,v,e in gg_loose.edges(data=True) if  e['conf'] == False]
      gg_conf=gg_loose.edge_subgraph(selected_edges).copy()
      orphans=[ n for n in gg_loose.nodes() if not gg_conf.has_node(n)]
      for node in orphans:
        gg_conf.add_node(node)
      gg=list(gg_conf.subgraph(cc).copy() for cc in sorted(nx.connected_components(gg_conf), key=len, reverse=True))
      for jj in range(len(gg)):
        haps=phase_conf_component_tree(gg[jj])
        for kk in range(len(bad_edges)):
          edge=bad_edges[kk]
          if gg[jj].has_node(edge[0]) and gg[jj].has_node(edge[1]):
            oddsratio=(edge[2]['cts'][0]+1.0)*(edge[2]['cts'][3]+1.0)/((edge[2]['cts'][1]+1.0)*(edge[2]['cts'][2]+1.0))
            if (oddsratio>1 and gg[jj].nodes[edge[0]]['phased_all']==gg[jj].nodes[edge[1]]['phased_all']) or (oddsratio<1 and not gg[jj].nodes[edge[0]]['phased_all']==gg[jj].nodes[edge[1]]['phased_all']):
              gg[jj].add_edge(edge[0], edge[1],  conf=True, mns=edge[2]['mns'], cts=edge[2]['cts'], wt=edge[2]['wt'])
        ec=check_phase(gg[jj])
        if ec>0:
          sys.stderr.write('phase error '+str(jj)+'\n')
          code.interact(local=locals())
         
        Gall=nx.Graph()
        for edge in list(gg[jj].edges()):
          curedge=gg[jj].edges[edge]
          [loc0, loc1]=edge
          if loc1<loc0:
            [loc0, loc1]=[loc1, loc0]
          add_allele_edges(Gall, loc0, loc1, curedge['cts'], curedge['mns'])
          print(str(ii)+'\t'+str(jj)+'\t'+edge[0]+'\t'+edge[1]+'\t'+str(curedge['cts'])+'\t'+str(curedge['mns']), file=outc)
        if gg[jj].number_of_edges()==0:
          singleton=list(gg[jj].nodes())[0]
          Gall.add_node(singleton+'_0')
          Gall.add_node(singleton+'_1')
        for hapid in range(2):
          minforest=nx.minimum_spanning_tree(Gall.subgraph(haps[hapid]), weight='mean_dist')
          mintree=list(minforest.subgraph(cc) for cc in nx.connected_components(minforest))
          for treeid in range(len(mintree)):
            id=(str(ii)+'_'+str(jj)+'_'+str(hapid)+'_'+str(treeid))
            treelist=list(mintree[treeid].edges)
            nodes=[int(re.split('[:_]', node)[1]) for node in mintree[treeid].nodes]
            chrs=[re.split('[:_]', node)[0] for node in mintree[treeid].nodes]
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
            for tredge in treelist:
              curedge=mintree[treeid].edges[tredge]
              edgestr=str(round(curedge['mean_dist'], 3))+';'+str(curedge['ct'])
              print(id+'\t'+edgestr+'\t'+tredge[0]+'\t'+tredge[1], file=outf)
            terminal_nodes=[]
            if len(list(mintree[treeid].nodes))>1:
              for node1 in list(mintree[treeid].nodes):
                if mintree[treeid].degree(node1)==1:
                  terminal_nodes.append(node1)
              for aa in range(len(terminal_nodes)-1):
                for bb in range(aa+1, len(terminal_nodes)):
                  sp=nx.shortest_path(mintree[treeid], terminal_nodes[aa], terminal_nodes[bb])
                  for nodeii in range(len(sp)):
                    node1=sp[nodeii]
                    if nodeii==0:
                      outstr=str(mintree[treeid].degree(node1))+'\t.'
                    else:
                      prevedge=mintree[treeid].edges[sp[nodeii], sp[nodeii-1]]
                      outstr=str(mintree[treeid].degree(node1))+'\t'+str(prevedge['ct'])+'_'+str(prevedge['mean_dist'])
                    print(id+'\t'+str(aa)+'_'+str(bb)+'\t'+outstr+'\t'+node1, file=outh)
            else:
              node1=list(mintree[treeid].nodes)[0]
              print(id+'\t0_0\t0\t0_0_0\t'+node1, file=outh)

                  
      sys.stderr.write(str(ct)+'\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom)+'\n')
      
  cp.disable()
  cp.print_stats()
                                           
                                             

if __name__ == "__main__":
    main()
