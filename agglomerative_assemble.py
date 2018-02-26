# -*- coding: utf-8 -*-
'''
function agglomerative_assemble
author Owen Chapman auberginekenobi
2/25/2018
'''

from overlap_combine import *
from math import inf

def find_max(adjacency_dict):
    '''
    Finds the maximum element in an adjacency dictionary
    adjacency_dict must be of the form
        obj: [(int, obj, obj)]
    '''
    max_val = -inf
    (key, index) = (None,None)
    for k in adjacency_dict:
        for i in range(len(adjacency_dict[k])):
            val = adjacency_dict[k][i][0]
            if val > max_val:
                max_val = val
                (key,index) = (k,i)
    return (max_val, key, index)


def agglomerative_assemble(glycan_list):
    '''
    Assembles contigs agglomeratively from a list of glycans. Objective function
    is to minimize the number of nodes in the resulting contig set, tiebroken
    by minimizing the number of glycans in the resulting contig set.
    Input: 
        glycan_list - list of Glycan objects
    Output:
        contigs - list of Glycan objects
    '''
    # find pairwise overlap scores across all glycan pairs
    distance_matrix = {} # Glycan : [(d, Glycan,(Glycan,Glycan))]
    for i in range(len(glycan_list)):
        gi = glycan_list[i]
        for j in range(i+1,len(glycan_list)):
            gj = glycan_list[j]
            d,o = overlap_score(gi,gj)
            o = o[0] if len(o)>0  else []
            if d in distance_matrix:
                distance_matrix[gi].append((d,gj,o))
            else:
                distance_matrix[gi] = [(d,gj,o)]

    # while there exists a nonzero overlap score:
    max_score, key, index = find_max(distance_matrix)
    while max_score>0:
        # Find the two glycans with the highest overlap
        # (arbitrary selection)
        g1 = key
        g2 = distance_matrix[key][index][1]
        o = distance_matrix[key][index][2]
        # Combine them
        g3 = overlap_combine(g1,g2,o)
        # update datastructure of distances 
        if g1 in distance_matrix:
            distance_matrix.pop(g1)
        if g2 in distance_matrix:
            distance_matrix.pop(g2)
        for key,value in distance_matrix.items():
            if g1 in value or g2 in value:
                distance_matrix[key].remove(value)
        distance_matrix[g3] = []
        for key in distance_matrix.keys():
            if key == g3:
                continue
            d,o = overlap_score(g3,key)
            o = o[0] if len(o)>0  else []
            distance_matrix[g3].append((d,key,o))
        max_score, key, index = find_max(distance_matrix)
    return list(distance_matrix.keys())