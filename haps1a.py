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
from statistics import mode
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
        hap0.append(get_phased_allele(gg1, nextnode, 0))
        hap1.append(get_phased_allele(gg1, nextnode, 0))
  return [hap0, hap1]

def add_allele_edges(gg, loc1, loc2, cts, mns, sds):
  [aa, bb, cc, dd]=cts
  #code.interact(local=locals())  
  nn=sum(cts)
  if (aa+dd)*1.0/nn>0.5:
    gg.add_edge(loc1+'_0', loc2+'_0', ct=aa, mean_dist=mns[0], sd_dist=sds[0])
    gg.add_edge(loc1+'_1', loc2+'_1', ct=dd, mean_dist=mns[3], sd_dist=sds[3])
  else:
    gg.add_edge(loc1+'_0', loc2+'_1', ct=bb, mean_dist=mns[2], sd_dist=sds[2])
    gg.add_edge(loc1+'_1', loc2+'_0', ct=cc, mean_dist=mns[1], sd_dist=sds[1])
    

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', '--infile', type=str, default=None)
  parser.add_argument('-b', '--bedfile', type=str, default=None)
  parser.add_argument('-e', '--edgefile', type=str, default=None)
  parser.add_argument('-p', '--hapfile', type=str, default=None)
  args = parser.parse_args()

  
  usage_denom=1024
  Gloc=nx.Graph()
  #Gall=nx.Graph()
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
        Gloc.add_edge(loc1, loc2, order=[loc1, loc2], conf=True, cts=cts, mns=mns, sds=sds)
      elif pre_filter_loose_pass(cts):
        Gloc.add_edge(loc1, loc2,  order=[loc1, loc2], conf=False, cts=cts, mns=mns, sds=sds)
      line=fp.readline().strip().decode()
      ct+=1

  sys.stderr.write('finished loading graph\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom)+'\n')
  with gzip.open(args.edgefile, 'wt') as outf, gzip.open(args.hapfile, 'wt') as outh, gzip.open(args.bedfile, 'wt') as outb:
    comp_loose=list(Gloc.subgraph(c) for c in sorted(nx.connected_components(Gloc), key=len, reverse=True))
    for ii in range(len(comp_loose)):
      gg_loose=comp_loose[ii]
      selected_edges = [(u,v) for u,v,e in gg_loose.edges(data=True) if  e['conf'] == True]
      bad_edges = [(u,v) for u,v,e in gg_loose.edges(data=True) if  e['conf'] == False]
      gg_conf = gg_loose.edge_subgraph(selected_edges)
      brlist = list(nx.bridges(gg_conf))
      for edge in brlist:
        curedge=gg_loose.edges[edge]
        if sum(curedge['cts'])<3 or min(curedge['cts'])<1:
          if abs(curedge['mns'][0]-curedge['mns'][0])>350 or abs(curedge['mns'][1]-curedge['mns'][2])>350:
            gg_loose.edges[edge].update({'conf': False})
      selected_edges = [(u,v) for u,v,e in gg_loose.edges(data=True) if  e['conf'] == True]
      bad_edges = [(u,v) for u,v,e in gg_loose.edges(data=True) if  e['conf'] == False]
      gg_conf=gg_loose.edge_subgraph(selected_edges)
      changed=True
      min_counts_strict=5
      while changed:
        changed = False
        for edge in bad_edges:
          if gg_conf.has_node(edge[0]) and gg_conf.has_node(edge[1]):
            if edge[1] in nx.node_connected_component(gg_conf, edge[0]):
              gg_loose.edges[edge].update({'conf': True})
              selected_edges.append(edge)
              changed=True
            else:
              cts=gg_loose.edges[edge]['cts']
              nn=sum(cts)
              if nn>min_counts_strict and (cts[0]+cts[3]==nn or  cts[1]+cts[2]==nn):
                gg_loose.edges[edge].update({'conf': True})
                selected_edges.append(edge)
                changed=True
          elif (gg_conf.has_node(edge[0]) and not gg_conf.has_node(edge[1])) or (gg_conf.has_node(edge[1]) and not gg_conf.has_node(edge[0])):
            cts=gg_loose.edges[edge]['cts']
            nn=sum(cts)
            if nn>min_counts_strict and (cts[0]+cts[3]==nn or cts[1]+cts[2]==nn):
              gg_loose.edges[edge].update({'conf': True})
              selected_edges.append(edge)
              changed=True
        bad_edges = [(u,v) for u,v,e in gg_loose.edges(data=True) if  e['conf'] == False]
        gg_conf=gg_loose.edge_subgraph(selected_edges)
      gg=list(gg_conf.subgraph(cc) for cc in sorted(nx.connected_components(gg_conf), key=len, reverse=True))
      Gall=nx.Graph()
      for jj in range(len(gg)):
        for edge in list(gg[jj].edges()):
          curedge=gg[jj].edges[edge]
          add_allele_edges(Gall, edge[0], edge[1], curedge['cts'], curedge['mns'], curedge['sds']) 
          print(str(ii)+'\t'+str(jj)+'\t'+edge[0]+'\t'+edge[1]+'\t'+str(curedge['cts'])+'\t'+str(curedge['mns'])+'\t'+str(curedge['sds']))
      for jj in range(len(gg)):
        haps=phase_conf_component(gg[jj]))
        for hapid in range(2):
          minforest=nx.minimum_spanning_tree(Gall.subgraph(haps[hapid]), weight='mean_dist')
          mintree=list(minforest.subgraph(cc) for cc in nx.connected_components(minforest))
          for treeid in range(len(mintree)):
            id=(str(ii)+'_'+str(jj)+'_'+str(hapid)+'_'+str(treeid))
            treelist=list(mintree[treeid].edges)
            nodes=[int(re.split('[:_]', node)[1]) for node in mintree[treeid].nodes]
            chr=mode([re.split('[:_]', node)[0] for node in mintree[treeid].nodes])
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
              edgestr=str(round(curedge['mean_dist'], 3))+';'+str(round(curedge['sd_dist'], 3))+';'+str(curedge['ct'])
              print(id+'\t'+edgestr+'\t'+tredge[0]+'\t'+tredge[1], file=outf)
            terminal_nodes=[]
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
                    outstr=str(mintree[treeid].degree(node1))+'\t'+str(prevedge['ct'])+'_'+str(prevedge['mean_dist'])+'_'+str(prevedge['sd_dist'])
                  print(id+'\t'+str(aa)+'_'+str(bb)+'\t'+outstr+'\t'+node1, file=outh)
      sys.stderr.write(str(ct)+'\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom)+'\n')
      #print('Memory usage info (Mb):\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom))

#  cp.disable()
#  cp.print_stats()
                                           
                                             

if __name__ == "__main__":
    main()
