from __future__ import print_function
import gzip
import networkx as nx
import re
from collections import defaultdict
import math
import argparse

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', '--infile', type=str, default=None)
  parser.add_argument('-e', '--edgefile', type=str, default=None)
  parser.add_argument('--strict', default=False, action='store_true')
  args = parser.parse_args()

  Gall=nx.Graph()
  Gloc=nx.Graph()
  infile=args.infile


  #load allele graph
  reads=defaultdict(dict)
  with gzip.open(infile, 'rb') as fp:
    line=fp.readline().strip()
    while line:
      [chr, posall, read, ll, readpos, strand]=re.split('[;:\t]', line)
      [pos, ref, alt, refalt]=re.split('_', posall)
      readpos=int(readpos)
      allele=chr+':'+posall
      locus=chr+':'+pos
      reads[read][allele]=readpos
      if not locus in Gloc.nodes:
        add_locus_node(Gloc, locus, posall)
      if not allele in Gall.nodes:
        add_allele_node(Gall, allele, locus)
      for r1 in reads[read].keys():
        if r1!=allele:
          edge0=[allele, r1]
          if not edge0 in Gall.edges:
            Gall.add_edge(allele, r1)
            Gall.edges[edge0]['count']=0
            Gall.edges[edge0]['dist']=0
            Gall.edges[edge0]['dist_sq']=0
          Gall.edges[edge0]['count']=Gall.edges[edge0]['count']+1
          dist=abs(reads[read][r1]-reads[read][allele])
          Gall.edges[edge0]['dist']=Gall.edges[edge0]['dist']+dist
          Gall.edges[edge0]['dist_sq']=Gall.edges[edge0]['dist_sq']+dist*dist
      line=fp.readline().strip()

  #load locus graph
  for edge in Gall.edges:
    mndist=Gall.edges[edge]['dist']/Gall.edges[edge]['count']
    Gall.edges[edge]['mean_dist']=mndist
    Gall.edges[edge]['sd_dist']=math.sqrt(Gall.edges[edge]['dist_sq']/Gall.edges[edge]['count']-mndist*mndist)
    node1=Gall.nodes[edge[0]]
    node2=Gall.nodes[edge[1]]
    add_counts_edge(Gloc, node1, node2, Gall.edges[edge]['count'])

  #create strict locus graph and prune edges
  Gconf=Gloc.copy()
  bad_edges=[]
  edgelist=list(Gloc.edges)
  for edge in edgelist:
    curedge=Gloc.edges[edge]
    [aa, bb, cc, dd]=[curedge['r', 'r'], curedge['r', 'a'], curedge['a', 'r'], curedge['a', 'a']]
    if (aa+bb==0 and cc*dd>0) or (cc+dd==0 and aa*bb>0) or (aa+cc==0 and bb*dd>0) or (bb+dd==0 and aa*cc>0):
      bad_edges.append(edge)
      remove_allele_edges(Gloc, Gall, edge)
      Gloc.remove_edge(edge[0], edge[1])
      Gconf.remove_edge(edge[0], edge[1])
    else:
      oddsratio=(aa+1.0)*(dd+1.0)/(bb+1.0)/(cc+1.0)
      if (oddsratio>2 or oddsratio<0.5) and (aa*dd>0 or bb*cc>0):
        pass
      else:
        Gconf.remove_edge(edge[0], edge[1])


  with open(args.edgefile, 'w') as outf:
    comp_loose=list(Gloc.subgraph(c).copy() for c in nx.connected_components(Gloc))
    #iterate over connected components of loose locus graph, phase and merge if possible
    for ii in range(len(comp_loose)):
      gg_loose=comp_loose[ii].copy()
      gg_conf=Gconf.subgraph(list(gg_loose.nodes())).copy()
      gg=list(gg_conf.subgraph(cc).copy() for cc in nx.connected_components(gg_conf))
      for jj in range(len(gg)):
        ggsub1=gg[jj]
        [h1, h2]=phase_conf_component(ggsub1)
        print(str(ii)+'\t'+str(jj)+'\t'+str(gg_conf.number_of_nodes())+'\t'+str(ggsub1.number_of_nodes())+'\t'+str(len(h1)))
      if len(gg)>1:
        merge_all(gg_conf, ii, Gloc)
      gg=list(gg_conf.subgraph(cc).copy() for cc in nx.connected_components(gg_conf))
      if len(gg)>1 and args.strict==False:
        merge_all(gg_conf, ii, Gloc, strict=False)
      #list of connected graphs, post-merging
      gg=list(gg_conf.subgraph(cc).copy() for cc in nx.connected_components(gg_conf))
      for jj in range(len(gg)):
        ggsub1=gg[jj]
        haps=phase_conf_component(ggsub1, strict=False)
        print(str(ii)+'\t'+str(jj)+'\t'+str(gg_conf.number_of_nodes())+'\t'+str(ggsub1.number_of_nodes())+'\t'+str(len(haps[0])))
        #subgraphs of allele graph corresponding to allele on each haplotype
        for hapid in range(2):
          Gallsub=Gall.subgraph(haps[hapid]).copy()
          minforest=nx.minimum_spanning_tree(Gallsub, weight='mean_dist')
          mintree=list(minforest.subgraph(cc).copy() for cc in nx.connected_components(minforest))
          for treeid in range(len(mintree)):
            treelist=list(mintree[treeid].edges)
            for tredge in treelist :
              edgestr=str(mintree[treeid].edges[tredge]['mean_dist'])+';'+str(round(mintree[treeid].edges[tredge]['sd_dist'], 3))+';'+str(mintree[treeid].edges[tredge]['count'])
              print(str(ii)+'_'+str(jj)+'_'+str(hapid+1)+'_'+str(treeid)+'\t'+edgestr+'\t'+tredge[0]+'\t'+tredge[1], file=outf)
      print('-----------------------------------------------------------------')



def add_locus_node(gg, locus, pos):
  gg.add_node(locus)
  splpos=re.split('_', pos)
  gg.nodes[locus]['alleles']={'r': splpos[1], 'a': splpos[2]}

def add_allele_node(gg, allele, locus):
  gg.add_node(allele)
  spl=re.split('_', allele)
  gg.nodes[allele]['refalt']=spl[3]
  gg.nodes[allele]['locus']=locus

def add_counts_edge(gg, node1, node2, count):
  if not [node1['locus'], node2['locus']] in gg.edges:
    gg.add_edge(node1['locus'], node2['locus'])
    curedge=gg.edges[node1['locus'], node2['locus']]
    curedge['first']=node1['locus']
    curedge['second']=node2['locus']
    curedge[('a', 'a')]=0
    curedge[('a', 'r')]=0
    curedge[('r', 'a')]=0
    curedge[('r', 'r')]=0
  curedge=gg.edges[node1['locus'], node2['locus']]
  if curedge['first']==node1['locus']:
    curedge[(node1['refalt'], node2['refalt'])]=count
  else:
    curedge[(node2['refalt'], node1['refalt'])]=count


def  remove_allele_edges(ggloc, ggall, edge):
  node0=[]
  node1=[]
  node0.append(edge[0]+'_'+ggloc.nodes[edge[0]]['alleles']['r']+'_'+ggloc.nodes[edge[0]]['alleles']['a']+'_'+'a')
  node0.append(edge[0]+'_'+ggloc.nodes[edge[0]]['alleles']['r']+'_'+ggloc.nodes[edge[0]]['alleles']['a']+'_'+'r')
  node1.append(edge[1]+'_'+ggloc.nodes[edge[1]]['alleles']['r']+'_'+ggloc.nodes[edge[1]]['alleles']['a']+'_'+'a')
  node1.append(edge[1]+'_'+ggloc.nodes[edge[1]]['alleles']['r']+'_'+ggloc.nodes[edge[1]]['alleles']['a']+'_'+'r')
  for ii in range(2):
    for jj in range(2):
      if [node0[ii], node1[jj]] in ggall.edges:
        ggall.remove_edge(node0[ii], node1[jj])


def phase_conf_component(ggsub, flip=False, strict=True):
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
  bad_nodes=set()
  phased_nodes=set()
  phased_nodes.add(maxnode)
  [h1, h2]=phase_from_node(ggsub, maxnode, phased_nodes, bad_nodes, strict)
  return([h1, h2])

def get_phased_allele(gg, node, hap):
  return(node+'_'+gg.nodes[node]['alleles']['r']+'_'+gg.nodes[node]['alleles']['a']+'_'+gg.nodes[node]['phased_all'][hap])

def phase_from_node(gg1, startnode, phased_nodes, bad_nodes, strict=True):
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
      if 'phased_all' in gg1.nodes[curnode].keys():
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
      else:
        bad_nodes.add(nextnode)
  return [hap1, hap2]
  


def merge_all(gg_conf, ii,  G, strict=True):
  done=False
  time0=True
  comps=[]
  gg_cur=gg_conf.copy()
  while not done:
    gg=list(gg_cur.subgraph(cc).copy() for cc in sorted(nx.connected_components(gg_cur), key=len, reverse=True))
    for jj in range(len(gg)):
      ggsub1=gg[jj]
      [h1, h2]=phase_conf_component(ggsub1, strict=False)
    if len(gg)==1:
      done=True
      break  
    ggmax=gg[0]
    adds=0
    for jj in range(1, len(gg)):
      [aa, bb, cc, dd]=check_join(G, ggmax, gg[jj])
      atot=aa+bb+cc+dd
      oddsratio=(aa+1.0)*(dd+1.0)/(cc+1.0)/(bb+1.0)
      if atot>0 and  (aa+dd==atot or (not strict and oddsratio>5)):
        print('OR='+str(round(oddsratio, 3))+'\tcounterev='+str(atot-aa-dd)+'\tjoin\t'+str(strict)+'\t'+str(adds)+'\t'+str(jj)+' of '+str(len(gg)))
        be=get_bridge_edge(G, ggmax, gg[jj], 'join')
        add_bridge(G, gg_conf, be)
        add_bridge(G, gg_cur, be)
        gg_temp=list(gg_cur.subgraph(cc).copy() for cc in sorted(nx.connected_components(gg_cur), key=len, reverse=True))
        ggmax=gg_temp[0]
        [h1, h2]=phase_conf_component(ggmax, strict=False)
        adds=adds+1
      elif atot>0 and (bb+cc==atot or (not strict and oddsratio<0.2)):
        print('OR='+str(round(oddsratio, 3))+'\tcounterev='+str(atot-bb-cc)+'\tflip\t'+str(strict)+'\t'+str(adds)+'\t'+str(jj)+' of '+str(len(gg)))
        be=get_bridge_edge(G, ggmax, gg[jj], 'flip')
        add_bridge(G, gg_conf, be)
        add_bridge(G, gg_cur, be)
        gg_temp=list(gg_cur.subgraph(cc).copy() for cc in sorted(nx.connected_components(gg_cur), key=len, reverse=True))
        ggmax=gg_temp[0]
        [h1, h2]=phase_conf_component(ggmax, strict=False)
        adds=adds+1
    if adds==0:
      comps.append(ggmax)
      for maxnode in list(ggmax.nodes):
        gg_cur.remove_node(maxnode)



def check_join(G, g1, g2):
  n1=g1.nodes()
  n2=g2.nodes()
  aa=[0,0,0,0]
  for node in n1:
    nb=set(list(G.neighbors(node))).intersection(n2)
    if len(nb)>0:
      for node2 in nb:
        curedge=G.edges[node, node2]
        atemp=[0,0,0,0]
        if curedge['first']==node:
          first=g1.nodes[node]['phased_all']
          second=g2.nodes[node2]['phased_all']
        else:
          first=g2.nodes[node2]['phased_all']
          second=g1.nodes[node]['phased_all']
        atemp=[curedge[(first[1], second[1])], curedge[(first[1], second[2])], curedge[(first[2], second[1])], curedge[(first[2], second[2])]]
        print(str(atemp))
        for ii in range(4):
          aa[ii]=aa[ii]+atemp[ii]    
  return aa

def get_bridge_edge(G, g1, g2, move):
  n1=g1.nodes()
  n2=g2.nodes()
  aa=[0,0,0,0]
  for node in n1:
    nb=set(list(G.neighbors(node))).intersection(n2)
    if len(nb)>0:
      for node2 in nb:
        curedge=G.edges[node, node2]
        atemp=[0,0,0,0]
        if curedge['first']==node:
          first=g1.nodes[node]['phased_all']
          second=g2.nodes[node2]['phased_all']
        else:
          first=g2.nodes[node2]['phased_all']
          second=g1.nodes[node]['phased_all']
        atemp=[curedge[(first[1], second[1])], curedge[(first[1], second[2])], curedge[(first[2], second[1])], curedge[(first[2], second[2])]]
        tot=atemp[0]+atemp[1]+atemp[2]+atemp[3]
        if move=='flip' and atemp[1]+atemp[2]==tot:
          return [node, node2]
        elif move=='join' and atemp[0]+atemp[3]==tot:
          return [node, node2]

def add_bridge(G, gg0, be):
  gg0.add_edge(be[0], be[1])
  for key in G.edges[be[0], be[1]].keys():
    gg0.edges[be[0], be[1]][key]= G.edges[be[0], be[1]][key]

if __name__ == "__main__":
    main()
