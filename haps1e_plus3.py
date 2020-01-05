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

def check_phase(ggsub):
  edgelist=list(ggsub.edges(data=True))
  error_count=0
  for ii in range(len(edgelist)):
    edata=edgelist[ii]
    oddsratio=(edata[2]['cts'][0]+1.0)*(edata[2]['cts'][3]+1.0)/((edata[2]['cts'][1]+1.0)*(edata[2]['cts'][2]+1.0))
    ph0=ggsub.nodes[edata[0]]['phased_all']
    ph1=ggsub.nodes[edata[1]]['phased_all']
    if (oddsratio>1 and ph1!=ph0) or (oddsratio<1 and ph1==ph0):
      sys.stderr.write('Error')
      error_count+=1
  return error_count
    

def phase_conf_component_tree(ggsub):
  tr=nx.minimum_spanning_tree(ggsub, weight='wt')
  startnode=''
  for node in tr.nodes:
    if tr.degree(node)<=1:
      startnode=node
      break
  #sys.stderr.write(startnode)
  tr.nodes[startnode]['phased_all']=[0,1]
  ggsub.nodes[startnode]['phased_all']=[0,1]
  [h0, h1] = phase_from_node_tree(ggsub, tr, startnode)
  return([h0, h1])
  

def phase_from_node_tree(gg1, tree, startnode):
  hap0=[]
  hap1=[]
  hap0.append(get_phased_allele(gg1, startnode, 0))
  hap1.append(get_phased_allele(gg1, startnode, 1))
  if len(gg1)>1:
    bfs=list(nx.bfs_edges(tree, startnode))
    for ii in range(len(bfs)):
      edata=gg1.edges[bfs[ii]]
      prevnode=bfs[ii][0]
      curnode=bfs[ii][1]
      if 'phased_all' in gg1.nodes[prevnode]:
        alleles=gg1.nodes[prevnode]['phased_all']
        oddsratio=(edata['cts'][0]+1.0)*(edata['cts'][3]+1.0)/((edata['cts'][1]+1.0)*(edata['cts'][2]+1.0))
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

def add_to_supergraph(superg, comphash, Gconf, infile, id):

  with gzip.open(infile, 'rt') as fp:
    line=fp.readline().strip()
    ct=0
    while line:
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
            ph1=Gconf.nodes[loc1]['phased_all']
            ph2=Gconf.nodes[loc2]['phased_all']
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
              superg.add_edge(c1, c2, cts_all=[cts], loc_pairs=[[loc1, loc2]], evtype=[id], pairct=1, conf=cc)
            else:
              curedge=superg.edges[c1, c2]
              curedge['loc_pairs'].append([loc1, loc2])
              curedge['pairct']+=1
              curedge['cts_all'].append(cts)
              curedge['conf'][conf]+=1
              curedge['evtype'].append(id)
      line=fp.readline().strip()
      ct+=1


                                        
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
  #code.interact(local=locals())

  with gzip.open(infile, 'rt') as fp:
    line=fp.readline().strip()
    ct=0
    while line :
      if ct%1000==0:
        sys.stderr.write(str(ct)+'\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom)+'\n')
      ll=re.split('[\t]', line)
      cts=list(map(int, ll[2:6]))
      mns=list(map(float, ll[6:10]))
      sds=list(map(float, ll[10:14]))
      tot=sum(cts)
      [loc1, loc2]=ll[0:2]
      if pre_filter_strict_pass(cts):
        Gloc.add_edge(loc1, loc2, order=[loc1, loc2], conf=True, cts=cts, mns=mns, wt=1.0/tot)
      elif pre_filter_loose_pass(cts):
        Gloc.add_edge(loc1, loc2,  order=[loc1, loc2], conf=False, cts=cts, mns=mns, wt=1.0/tot)
      line=fp.readline().strip()
      ct+=1

  loc2comp={}
  compid=0

  sys.stderr.write('finished loading graph\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom)+'\n')
  with gzip.open(args.compfile, 'wt') as outc:
    comp_loose=list(Gloc.subgraph(c) for c in sorted(nx.connected_components(Gloc), key=len, reverse=True))
    for ii in range(10):
    #for ii in range(len(comp_loose)):
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
      gg_conf=gg_loose.edge_subgraph(selected_edges).copy()
      orphans=[ n for n in gg_loose.nodes() if not gg_conf.has_node(n)]
      for node in orphans:
        gg_conf.add_node(node)
      gg=list(gg_conf.subgraph(cc) for cc in sorted(nx.connected_components(gg_conf), key=len, reverse=True))
      for jj in range(len(gg)):
        haps=phase_conf_component_tree(gg[jj])
        for kk in range(len(bad_edges)):
          edge=bad_edges[kk]
          if gg[jj].has_node(edge[0]) and gg[jj].has_node(edge[1]):
            oddsratio=(edge[2]['cts'][0]+1.0)*(edge[2]['cts'][3]+1.0)/((edge[2]['cts'][1]+1.0)*(edge[2]['cts'][2]+1.0))
            if (oddsratio>1 and gg[jj].nodes[edge[0]]['phased_all']==gg[jj].nodes[edge[1]]['phased_all']) or (oddsratio<1 and not gg[jj].nodes[edge[0]]['phased_all']==gg[jj].nodes[edge[1]]['phased_all']):
              gg_loose.edges[edge[0], edge[1]].update({'conf': True})
      for node in gg_loose.nodes():
        gg_loose.nodes[node]['phased_all']=gg_conf.nodes[node]['phased_all']
      bad_edges = [(u,v) for u,v,e in gg_loose.edges(data=True) if  e['conf'] == False]
      g1=gg_loose.copy()
      g1.remove_edges_from(bad_edges)
      ec=check_phase(g1)
      #if ec>0:
        #code.interact(local=locals())
      sys.stderr.write(str(ii)+'\t'+str(ec)+'\n')
      gg=list(g1.subgraph(cc) for cc in sorted(nx.connected_components(g1), key=len, reverse=True))
      for jj in range(len(gg)):
        for loc in list(gg[jj].nodes()):
          Gloc.nodes[loc]['phased_all']=gg[jj].nodes[loc]['phased_all']
          loc2comp[loc]=compid
        compid+=1
        
  superg=nx.Graph()
  add_to_supergraph(superg, loc2comp, Gloc, args.infile2, id='ont')
  sys.stderr.write('post hash\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom)+'\n')
  add_to_supergraph(superg, loc2comp, Gloc, args.infile, id='hifi')
  sys.stderr.write('post hash\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom)+'\n')

  super_comps=list(superg.subgraph(c) for c in sorted(nx.connected_components(superg), key=len, reverse=True))
  code.interact(local=locals())
  edgelist=list(superg.edges(data=True))
  sys.exit(0)                  
#  for ii in range(len(edgelist)):
#    edge=edgelist[ii]
#    [inner, outer, tot, aa, bb, cc, dd]=[0,0,0,0,0,0,0]
#    for cts in edge[2]['cts_all']:
#      if cts[1]+cts[2]==0:
#        outer+=1
#        if cts[0]>aa:
#          aa=cts[0]
#        if cts[3]>dd
#          dd=cts[3]
#      elif cts[0]+cts[3]==0:
#        inner+=1
#        if cts[1]>bb:
#          bb=cts[1]
#        if cts[2]>cc
#          cts[2]=cc
#      tot+=1
#    frac_inner=inner/(0.0+tot)
#    frac_outer=outer/(0.0+tot)
#    conf=0
#    if edge[2]['conf'][1]>0:
#      if frac_inner>0.95:
#        edge[2]['cts']=[aa, bb, cc, dd]
#        edge[2]['conf1']=1
#      elif frac_outer>0.95:
#        edge[2]['cts'][aa, bb, cc, dd]
#        edge[2]['conf1']=1
#      else:
#        edge[2]['conf1']=0
#    else:
#       edge[2]['conf1']=0
#
#  selected_edges = [(u,v) for u,v,e in superg.edges(data=True) if  e['conf1'] == 1]
#  bad_edges = [(u,v) for u,v,e in superg.edges(data=True) if  e['conf1'] == 0]
#  superg_conf=superg.edge_subgraph(selected_edges).copy()
#
#  for edge in superg.edges(data=True):
#    if 'hifi' in edge[2]['evtype']:
#      sys.stderr.write(str(edge)+'\n')
#  
#  super_conf_comps=list(superg_conf.subgraph(c) for c in sorted(nx.connected_components(superg_conf), key=len, reverse=True))
#
#  ii=0
#  superg_conf=super_comps[ii].edge_subgraph(selected_edges).copy()
#  super_conf_comps=list(superg_conf.subgraph(c) for c in sorted(nx.connected_components(superg_conf), key=len, reverse=True))
#  code.interact(local=locals()) 
#

#
#
#              gg[jj].add_edge(edge[0], edge[1], order=edge[2]['order'], conf=True, mns=edge[2]['mns'], cts=edge[2]['cts'], wt=edge[2]['wt'])
#
#
#        Gall=nx.Graph()
#        for edge in list(gg[jj].edges()):
#          curedge=gg[jj].edges[edge]
#          if edge[0]==curedge['order'][0]:
#            [loc0, loc1]=edge
#          else:
#            [loc1, loc0]=edge
#          add_allele_edges(Gall, loc0, loc1, curedge['cts'], curedge['mns'])
#          print(str(ii)+'\t'+str(jj)+'\t'+edge[0]+'\t'+edge[1]+'\t'+str(curedge['cts'])+'\t'+str(curedge['mns'])), file=outc)
#        if gg[jj].number_of_edges()==0:
#          singleton=list(gg[jj].nodes())[0]
#          Gall.add_node(singleton+'_0')
#          Gall.add_node(singleton+'_1')
#        for hapid in range(2):
#          minforest=nx.minimum_spanning_tree(Gall.subgraph(haps[hapid]), weight='mean_dist')
#          mintree=list(minforest.subgraph(cc) for cc in nx.connected_components(minforest))
#          for treeid in range(len(mintree)):
#            id=(str(ii)+'_'+str(jj)+'_'+str(hapid)+'_'+str(treeid))
#            treelist=list(mintree[treeid].edges)
#            nodes=[int(re.split('[:_]', node)[1]) for node in mintree[treeid].nodes]
#            chrs=[re.split('[:_]', node)[0] for node in mintree[treeid].nodes]
#            chr=Counter(chrs).most_common(1)[0][0]
#            nodes.sort()
#            dists=[]
#            ones=[]
#            for kk in range(len(nodes)):
#              dists.append(nodes[kk]-nodes[0]+1)
#              ones.append(1)
#              dstr=','.join(map(str, dists))
#              onestr=','.join(map(str,ones))
#            print(chr+'\t'+str(nodes[0])+'\t'+str(max(nodes))+'\t'+id+'\t100\t.\t'+str(nodes[0])+'\t'+str(max(nodes))+'\t150,150,0\t'+str(len(dists))+'\t'+onestr+'\t'+dstr, file=outb)
#            for tredge in treelist:
#              curedge=mintree[treeid].edges[tredge]
#              edgestr=str(round(curedge['mean_dist'], 3))+';'+str(round(curedge['sd_dist'], 3))+';'+str(curedge['ct'])
#              print(id+'\t'+edgestr+'\t'+tredge[0]+'\t'+tredge[1], file=outf)
#            terminal_nodes=[]
#            if len(list(mintree[treeid].nodes))>1:
#              for node1 in list(mintree[treeid].nodes):
#                if mintree[treeid].degree(node1)==1:
#                  terminal_nodes.append(node1)
#              for aa in range(len(terminal_nodes)-1):
#                for bb in range(aa+1, len(terminal_nodes)):
#                  sp=nx.shortest_path(mintree[treeid], terminal_nodes[aa], terminal_nodes[bb])
#                  for nodeii in range(len(sp)):
#                    node1=sp[nodeii]
#                    if nodeii==0:
#                      outstr=str(mintree[treeid].degree(node1))+'\t.'
#                    else:
#                      prevedge=mintree[treeid].edges[sp[nodeii], sp[nodeii-1]]
#                      outstr=str(mintree[treeid].degree(node1))+'\t'+str(prevedge['ct'])+'_'+str(prevedge['mean_dist'])+'_'+str(prevedge['sd_dist'])
#                    print(id+'\t'+str(aa)+'_'+str(bb)+'\t'+outstr+'\t'+node1, file=outh)
#            else:
#              node1=list(mintree[treeid].nodes)[0]
#              print(id+'\t0_0\t0\t0_0_0\t'+node1, file=outh)
#                  
#      sys.stderr.write(str(ct)+'\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom)+'\n')
      
                                           
                                             

if __name__ == "__main__":
    main()
