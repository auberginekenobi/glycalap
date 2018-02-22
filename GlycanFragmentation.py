# jtsorren

import scipy 
from scipy.stats import poisson
import random
import GlycanParser

def GlycanFragmentation(n, glycan, lmbda=2.5):
  '''
  Glycan Fragmentation will take a glycan, acyclic undirected graph, and
  simulate fragments of the graph. The mean fragment length is defined as 2.5.
  The below method will return fragmented acyclic undirected graphs that follow
  the structure of object Glycan.
  '''
  fragments = []
  for i in range(n): # n defines the number of fragments to generate
    size = poisson.rvs(lmbda) 
    while size > len(glycan.names):
      size = poisson.rvs(lmbda)
    start = random.randint(0,len(glycan.names) - 1) # selection of starting node
    j = 0
    names, bonds = [], []
    prev = None
    while j < size: # loop until fragment is desired size
      d = glycan.bonds[start]
      try:
        d.pop(-1) # if node has root remove it
      except KeyError:
        pass
      nextNode = random.choice(list(d.keys())) # randomly select next node
      if nextNode == prev:
        if len(list(d.keys())) == 1:
          break # if the next node can only return to previous node then fragment is maximized 
        else:
          keys = list(d.keys())
          keys.remove(prev) # do not allow return to previous node
          nextNode = random.choice(keys)
      prev = start 
      start = nextNode 
      name = glycan.names[start] # add the next node
      bond = {start: d[start]} # bond leading to node above
      names.append(name), bonds.append(bond)
      j += 1
    fragment = GlycanParser.Glycan(names, bonds)
    if len(fragment.bonds) == 0:
      pass
    else:
      fragment.bonds.remove(fragment.bonds[0]) # remove bond leading to the first node
      newKeys = list(range(1, len(fragment.names)))
      fragment.bonds.append({})
      for key in newKeys: # remove positional information
        fragment.bonds[key - 1] = {key: list(fragment.bonds[key - 1].values())[0]}
      for key in newKeys: # add backward edges
        fragment.bonds[key] = {**fragment.bonds[key], **{key - 1: (fragment.bonds[key - 1][key][0], fragment.bonds[key - 1][key][2], fragment.bonds[key - 1][key][1])}}        
      fragments.append(fragment)
  return fragments
  
test = GlycanParser.Glycan('alpha-L-Fucp-(1->2)-beta-D-Galp-(1->4)[alpha-L-Fucp-(1->3)]-beta-D-GlcpNAc-(1->4)-beta-D-Glcp(1->')

for fragment in GlycanFragmentation(10, test):
  print(fragment.names, fragment.bonds)