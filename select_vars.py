#!/usr/bin/env python

import gzip
import re
from collections import defaultdict
import math
import argparse
import cProfile , pstats , resource
import sys
import code
from collections import Counter
import numpy as np

                                        
def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', '--infile', type=str, default=None)
  parser.add_argument('-o', '--snpfile', type=str, default=None)
  args = parser.parse_args()
  usage_denom=1024
  
  #tokeep_homref={}
  tokeep={}

  oldread=''
  keep_pos=[]
  maybe_pos={}

  with gzip.open(args.infile, 'rt') as fp:
    line=fp.readline().strip()
    ct=0
    while line :
      if ct%1000==0:
        sys.stderr.write(str(ct)+'\t'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/usage_denom)+'\n')
      [snp, read, pos, snp0, gt]=re.split('[\t]', line)
      pos=int(pos)
      if read!=oldread:
        if oldread!='':
          keep_pos.sort()
          for ii in range(1, len(keep_pos)):
            dist=keep_pos[ii]-keep_pos[ii-1]
            if dist>1000:
              pp=[]
              start=keep_pos[ii-1]; stop=keep_pos[ii]
              positions=[start, stop]
              for pos in maybe_pos.keys():
                if pos>start and pos<stop:
                  pp.append(pos)
              pp=np.array(pp)
              while len(pp)>0 and dist>1000:
                choices=pp[pp<start+1000]
                if len(choices)>0:
                  start=max(choices)
                else:
                  start=min(pp)
                positions.append(start); pp=pp[pp>start]
                tokeep[maybe_pos[start]]=1
                dist=stop-start
              positions.sort()
              #a = input('').split(" ")[0]
        keep_pos=[]
        maybe_pos={}
        oldread=read
      else:
        if gt in ['1', '2'] or snp in tokeep:
          tokeep[snp]=1
          keep_pos.append(pos)
        else:
          maybe_pos[pos]=snp
      line=fp.readline().strip()
      ct+=1

  with gzip.open(args.snpfile, 'wt') as fo:
    for snp in tokeep.keys():
      print(snp, file=fo)


if __name__ == "__main__":
    main()
