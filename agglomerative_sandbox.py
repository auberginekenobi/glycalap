# -*- coding: utf-8 -*-
'''
Agglomerative sandbox
for playing around with the developing agglomerative_assemble
author Owen Chapman auberginekenobi
2/25/2018
'''

from GlycanParser import *
from GlycanFragmentation import *
from overlap_combine import *
from agglomerative_assemble import *


VIM="alpha-D-NeupAc-(2->3)-beta-D-Galp-(1->4)-beta-D-GlcpNAc(1->"
a=Glycan(VIM)
b=Glycan(VIM)
c=Glycan(VIM)
l = [a,b,c]
score, alns = overlap_score(a,b)
print(score)
print("%%%%%")
print(alns[0])
print(alns[1])
print('%%%%%')
print(a)
overlap_combine(a,b,alns[3])
print(a==b)
