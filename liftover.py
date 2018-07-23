#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

from collections import defaultdict
from bisect import bisect_left, bisect_right

class GapLessBlock:
   '''
   Coordinates of alignment without gaps between genomes.
   Blocks have the 4 following attributes:
      s1: block start in first genome
      e1: block end in first genome
      s2: block start in second genome
      e2: block end in second genome

   The class implements comparisons with __lt__ and __gt__.
   The < operator compares numbers to e1, and the > operator
   compares numbers to e2. Since 'bisect_left' uses __lt__
   and 'bisect_right' uses __gt__, this gives a way to look
   for a block in genome 1 or genome 2 by changing the bisection
   method.
   '''
   def __init__(self, s1, e1, s2, e2):
      self.s1 = s1
      self.e1 = e1
      self.s2 = s2
      self.e2 = e2

   def __lt__(self, other):
      return self.e1 < other

   def __gt__(self, other):
      return self.e2 > other

   def boundaries(self):
      return (self.s1, self.e1, self.s2, self.e2)


class SeqPair:
   def __init__(self, seqname1, seqname2):
      self.seqname1 = seqname1
      self.seqname2 = seqname2
      self.blocks = list()

   def append(self, block):
      self.blocks.append(block)


class ChainIndex:

   def __init__(self):
      self.seqnames = defaultdict(set)

   @staticmethod
   def construct_from_chain_file(f):
      cidx = ChainIndex()
      s1 = s2 = float('nan')
      for line in f:
         if line[0] == '#' or line.rstrip() == '':
            continue
         if line.startswith('chain'): 
            #chain 110776 Y 3667352 + 1664227 1665401 A6s:Y 3667231 + 1747115 1748289 9
            items = line.split()
            s1 = int(items[5])
            s2 = int(items[10])
            # Create new 'SeqPair' for this chromosome.
            pair = SeqPair(items[2], items[7])
            # Create two seqname entries that point to it.
            cidx.seqnames[items[2]].add(pair)
            cidx.seqnames[items[7]] = cidx.seqnames[items[2]]
            continue
         try:
            sz,n1,n2 = line.split()
         except ValueError:
            n1 = n2 = 0
            sz = int(line.rstrip())
         # Update end pointers.
         e1 = s1+int(sz)
         e2 = s2+int(sz)
         # Store block.
         pair.append(GapLessBlock(s1,e1,s2,e2))
         # Update start pointers.
         s1 = e1+int(n1)+1
         s2 = e2+int(n2)+1
      return cidx


   def liftover(self, seqname, pos):
      for pair in self.seqnames[seqname]:
         if seqname == pair.seqname1:
            # We need to search genome 1.
            i = bisect_left(pair.blocks, pos)
            block = pair.blocks[i]
            s1,e1, s2,e2 = block.boundaries()
            if s1 <= pos <= e1:
               return (pair.seqname2, s2 + (pos-s1))
         elif seqname == pair.seqname2:
            # We need to search genome 2.
            i = bisect_right(pair.blocks, pos)
            block = pair.blocks[i]
            s1,e1, s2,e2 = block.boundaries()
            if s2 <= pos <= e2:
               return (pair.seqname1, s1 + (pos-s2))
      # Absent
      return None


if __name__ == '__main__':
   with open(sys.argv[1]) as f:
      idx = ChainIndex.construct_from_chain_file(f)
   print idx.liftover('3R', 1000000)
   print idx.liftover('A6s:3R', 1000098)
