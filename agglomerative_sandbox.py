# -*- coding: utf-8 -*-
'''
Agglomerative sandbox
for playing around with the developing agglomerative_assemble
author Owen Chapman auberginekenobi
2/25/2018
'''

from GlycanParser import *
from overlap_combine import *
from agglomerative_assemble import *


sialyl_lewis_a = "alpha-D-NeupAc-(2->3)-beta-D-Galp-(1->3)[alpha-L-Fucp-(1->4)]-beta-D-GlcpNAc(1->"
lewis_a = "beta-D-Galp-(1->3)[alpha-L-Fucp-(1->4)]-?-D-GlcpNAc(1->"
notVIM = "alpha-D-NeupAc-(2->3)-beta-D-Galp-(1->3)-beta-D-GlcpNAc(1->"
a=Glycan(sialyl_lewis_a)
b=Glycan(lewis_a)
c=Glycan(notVIM)
l = [a,b,c]
print(agglomerative_assemble(l))
