# -*- coding: utf-8 -*-
'''
function overlap_score
author Owen Chapman auberginekenobi
2/19/2018
'''
from GlycanParser import *

def overlap_score(glycan1, glycan2):
    '''
    Finds the largest overlap of the same sugars and bonds between two glycans.
    Overlap scoring is defined as +1 if the same sugars exist in the same position on the tree;
        +0 if the same edges are in the same position on their respective trees; and
        -Inf if edges or nodes are different.
    Note that because each edge must come with a corresponding node, and the mismatch score is
        -Inf, the relative weights of node and edge matching don't matter.
    Inputs: Two Glycan objects
    Outputs:
        i: score of longest overlap
        seeds: all subtrees with an overlap score of i.
    '''
    i=1
    seeds = []
    for n in range(len(glycan1)):
            for m in range(len(glycan2)):
                if glycan1.names[n]==glycan2.names[m]:
                    seeds.append((Glycan([n],[{}]),Glycan([m],[{}])))
    if seeds==[]:
        return 0,[]
    for i in range(2,min(len(glycan1)+1,len(glycan2)+1)):
        new_seeds=[]
        for aln in seeds:
            aln1,aln2 = aln
            # n,m indices in glycans12, names in aln12.
            for n in range(len(glycan1)):
                if not n in aln1.names:
                    for m in range(len(glycan2)):
                       if not m in aln2.names and glycan1.names[n] == glycan2.names[m]:
                           # b,c indices in glycans12, names in aln12
                           for b in glycan1.bonds[n]:
                               if b in aln1.names:
                                   bi = aln1.names.index(b)
                                   c = aln2.names[bi]
                                   if (glycan2.isBond(m,c) and
                                       glycan1.bonds[n][b] == glycan2.bonds[m][c]):
                                       new_aln1_names = aln1.names[:]
                                       new_aln1_names.append(n)
                                       new_aln2_names = aln2.names[:]
                                       new_aln2_names.append(m)
                                       new_aln1_bonds = aln1.bonds[:]
                                       new_aln1_bonds[bi][len(aln1)]=glycan1.bonds[n][b]
                                       new_aln1_bonds.append({bi:glycan1.bonds[n][b]})
                                       new_seeds.append((Glycan(new_aln1_names,new_aln1_bonds),
                                                       Glycan(new_aln2_names,new_aln1_bonds)))
        if new_seeds !=[]:
            seeds=new_seeds
        else:
            i-=1
            break
    return(i,seeds)


