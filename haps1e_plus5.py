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

def odds_ratio([aa, bb, cc, dd]):
  return 1.0*(aa+1)*(dd+1)/((bb+1)*(cc+1))

def check_edge(edge, gg):
  oddsratio=odds_ratio(edge[2]['cts'])
  ph0=gg.nodes[edge[0]['phased_all']
  ph1=gg.nodes[edge[1]['phased_all']
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
      tot=sum(cts)
      [loc1, loc2]=ll[0:2]
      if pre_filter_loose_pass(cts):
        if loc1 in comphash and loc2 in comphash:
          [c1, c2]=[comphash[loc1], comphash[loc2]]
          if c1!=c2:
            if c2<c1:
              [loc1, loc2]=[loc2, loc1]
              [c1, c2]=[c2, c1]
              cts=(cts[0], cts[2], cts[1], cts[3])
              mns=(mns[0], mns[2], mns[1], mns[3])
            [ph1, ph2]=Gconf.nodes[loc1]['phased_all'], Gconf.nodes[loc2]['phased_all']
            if ph1==[1,0] and ph2==[1,0]:
              pass
            elif ph1==[0,1] and ph2==[0,1]:
              cts=(cts[3], cts[2], cts[1], cts[0])
              mns=(mns[3], mns[2], mns[1], mns[0])
            elif ph1==[1,0] and ph2==[0,1]:
              cts=(cts[1], cts[0], cts[3], cts[2])
              mns=(mns[1], mns[0], mns[3], mns[2])
            elif ph1==[0,1] and ph2==[1,0]:
              cts=(cts[2], cts[3], cts[0], cts[1])
              mns=(mns[2], mns[3], mns[0], mns[1])
            dist=max(mns)
            if not superg.has_edge(c1, c2):
              superg.add_edge(c1, c2, cts_all={cts : 1}, min_dist={cts : dist}, mns_dict={cts: mns}, min_loc_pairs={cts : [loc1, loc2]}, evtype={cts : id}, pairct=1)
            else:
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
      line=fp.readline().strip()

def remove_bridges(GG):

  selected_edges = [(u,v) for u,v,e in GG.edges(data=True) if  e['conf'] == True]
  gg_conf=GG.edge_subgraph(selected_edges)
  brlist=list(nx.bridges(gg_conf))
  for edge in brlist:
    curedge=GG.edges[edge]
    if sum(curedge['cts'])<min_counts_strict:
      if abs(curedge['mns'][0]-curedge['mns'][3])>350 or abs(curedge['mns'][1]-curedge['mns'][2])>350:
        GG.edges[edge].update({'conf': False})


def transform(cts, trans=1):
  
  if trans==2:
    cts=(cts[3], cts[2], cts[1], cts[0])
  elif trans==3:
    cts=(cts[1], cts[0], cts[3], cts[2])
  elif trans==4:
    cts=(cts[2], cts[3], cts[0], cts[1])
  return cts

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
      superg.add_edge(c1, c2, pairct=edge[2]['pairct'], cts_all={}, min_dist={}, min_loc_pairs={}, evtype={})
    else:
      superg.edges[c1, c2]['pairct']+=edge[2]['pairct']
    ee=superg.edges[c1, c2]
    for cts in edge[2]['cts_all'].keys():
      ctstrans=transform(cts, trans)
      if not ctstrans in ee['cts_all']:
        ee['cts_all'][ctstrans]=edge[2]['cts_all'][cts]
        ee['min_dist'][ctstrans]=edge[2]['min_dist'][cts]
        ee['min_loc_pairs'][ctstrans]=edge[2]['min_loc_pairs'][cts]
        ee['evtype'][ctstrans]=edge[2]['evtype'][cts]
      else:
        ee['cts_all'][ctstrans]=ee['cts_all'][ctstrans]+edge[2]['cts_all'][cts]
        if edge[2]['min_dist'][cts]<ee['min_dist'][ctstrans]:
          ee['min_dist'][ctstrans]=edge[2]['min_dist'][cts]
          ee['min_loc_pairs'][ctstrans]=edge[2]['min_loc_pairs'][cts]


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
        Gloc.add_edge(loc1, loc2,  conf=True, cts=cts, mns=mns, wt=1.0/tot)
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

  remove_bridges(Gloc)


  loc2comp={}
  compid=0
  comp2haps={}

  
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
      gg=list(gg_conf.subgraph(cc) for cc in sorted(nx.connected_components(gg_conf), key=len, reverse=True))
      for jj in range(len(gg)):
        haps=phase_conf_component_tree(gg[jj])
        for kk in range(len(bad_edges)):
          edge=bad_edges[kk]
          if gg[jj].has_node(edge[0]) and gg[jj].has_node(edge[1]):
            if check_edge(edge, gg[jj]):
              gg_loose.edges[edge[0], edge[1]].update({'conf': True})
        ec=check_phase(gg[jj])
        if ec>0:
          sys.stderr.write('phase error\n')
          sys.exit(1)
        for loc in list(gg[jj].nodes):
          gg_loose.nodes[loc]['phased_all']=gg_conf.nodes[loc]['phased_all']
          loc2comp[loc]=compid
        compid+=1
#        
    superg=nx.Graph()
    add_to_supergraph1(superg, loc2comp, Gloc, args.infile, id='hifi')
    add_to_supergraph1(superg, loc2comp, Gloc, args.infile2, id='ont')
    sys.stderr.write('post hash\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom)+'\n')
    
    edgelist=list(superg.edges(data=True))
    for ii in range(len(edgelist)):
      score_edge(edgelist[ii])
      
    selected_edges2 =  [(u,v) for u,v,e in superg.edges(data=True) if  e['conf1'] == 2]
    sconf=superg.edge_subgraph(selected_edges2).copy()
    orphans = [ n for n in superg.nodes() if not sconf.has_node(n)]
    for node in orphans:
      sconf.add_node(node)
    super_comps=list(sconf.subgraph(c) for c in sorted(nx.connected_components(sconf), key=len, reverse=True))
    
    superg2=nx.Graph()
    comp2super={}
    for ii in range(len(super_comps)):
      gg_loose=super_comps[ii]
      haps=phase_conf_component_tree(gg_loose)
      for node in super_comps[ii].nodes():
        comp2super[node]=ii

    selected_edges =  [(u,v) for u,v,e in superg.edges(data=True) if  e['conf1'] == 1]
    sconf2=superg.edge_subgraph(selected_edges).copy()
    orphans = [ n for n in superg.nodes() if not sconf2.has_node(n)]
    for node in orphans:
      sconf2.add_node(node)
    edgelist=list(sconf2.edges(data=True))
    for edge in edgelist:
      add_edge_to_supergraph1(edge, superg2, comp2super, sconf)


    edgelist=list(superg2.edges(data=True))
    for ii in range(len(edgelist)):
      score_edge(edgelist[ii])

    selected_edges3 =  [(u,v) for u,v,e in superg2.edges(data=True) if  e['conf1'] > 0 ]

    for edge in selected_edges2:
      [loc1, loc2]=superg.edges[edge]['loc']
      if loc2<loc1:
        [loc1, loc2]=[loc2, loc1]
      [ph1, ph2]=[Gloc.nodes[loc1]['phased_all'], Gloc.nodes[loc2]['phased_all']]
      cts=superg.edges[edge]['cts']
      if ph1==[1,0] and ph2==[1,0]:
        pass
      elif ph1==[0,1] and ph2==[0,1]:
        cts.reverse()
      elif ph1==[1,0] and ph2==[0,1]:
        cts=[cts[1], cts[0], cts[3], cts[2]]
      elif ph1==[0,1] and ph2==[1,0]:
        cts=[cts[2], cts[3], cts[0], cts[1]]
      Gloc.add_edge(loc1, loc2, cts=cts, wt=1/sum(cts), conf=True)

    for edge in selected_edges3:
      [loc1, loc2]=superg2.edges[edge]['loc']
      if loc2<loc1:
        [loc1, loc2]=[loc2, loc1]
      [ph1, ph2]=[Gloc.nodes[loc1]['phased_all'], Gloc.nodes[loc2]['phased_all']]
      cts=superg2.edges[edge]['cts']
      if ph1==[1,0] and ph2==[1,0]:
        pass
      elif ph1==[0,1] and ph2==[0,1]:
        cts.reverse()
      elif ph1==[1,0] and ph2==[0,1]:
        cts=[cts[1], cts[0], cts[3], cts[2]]
      elif ph1==[0,1] and ph2==[1,0]:
        cts=[cts[2], cts[3], cts[0], cts[1]]
      Gloc.add_edge(loc1, loc2, cts=cts, wt=1/sum(cts), conf=True)
      


    selected_edges = [(u,v) for u,v,e in Gloc.edges(data=True) if  e['conf'] == True]
    gg_conf=Gloc.edge_subgraph(selected_edges).copy()
    orphans=[ n for n in Gloc.nodes() if not gg_conf.has_node(n)]
    for node in orphans:
      gg_conf.add_node(node)
    gg=list(gg_conf.subgraph(cc) for cc in sorted(nx.connected_components(gg_conf), key=len, reverse=True))
    for jj in range(len(gg)):
        haps=phase_conf_component_tree(gg[jj])
        Gall=nx.Graph()
        for edge in list(gg[jj].edges()):
          [loc0, loc1]=edge
          if loc1<loc2:
            [loc0, loc1]=[loc1, loc0]
          curedge=gg[jj].edges[edge]
          add_allele_edges(Gall, loc0, loc1, curedge['cts'], curedge['mns'])
          print(str(ii)+'\t'+str(jj)+'\t'+edge[0]+'\t'+edge[1]+'\t'+str(curedge['cts']), file=outc)
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
              edgestr=str(round(curedge['mean_dist'], 3))+';'+str(round(curedge['sd_dist'], 3))+';'+str(curedge['ct'])
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
                      outstr=str(mintree[treeid].degree(node1))+'\t'+str(prevedge['ct'])+'_'+str(prevedge['mean_dist'])+'_'+str(prevedge['sd_dist'])
                    print(id+'\t'+str(aa)+'_'+str(bb)+'\t'+outstr+'\t'+node1, file=outh)
            else:
              node1=list(mintree[treeid].nodes)[0]
              print(id+'\t0_0\t0\t0_0_0\t'+node1, file=outh)
                  
    
      
  cp.disable()
  cp.print_stats()
                                           

                                           
                                             

if __name__ == "__main__":
    main()
