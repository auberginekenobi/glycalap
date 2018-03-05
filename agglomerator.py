# -*- coding: utf-8 -*-
"""
Created on Sun Mar 04 15:34:39 2018

@author: jhyun_000
"""

from GlycanParser import Glycan
from overlap_score_fast import overlap_score

def agglomerative_assemble(glycan_list):
    '''
    Assembles contigs agglomeratively from a list of glycans, i.e. looks at
    all pairwise overlaps, selects a pair with largest overlap, fuses the pair,
    recalculate overlaps of fused pair with other fragments, and repeats until
    no more fragments can be fused.
    
    Objective function is to minimize the number of nodes in the resulting 
    contig set, tiebroken    by minimizing the number of glycans in the 
    resulting contig set.
    
    Input: 
        glycan_list - list of Glycan objects
    Output:
        contigs - list of Glycan objects
    '''
    N = len(glycan_list) # number of fragments initially
    
    ''' Compute initial pairwise overlaps '''
    pairwise_overlaps = {} # maps glycan1:glycan2:(score, overlaps)
    for glycan in glycan_list: # initialize all neighbors
        pairwise_overlaps[glycan] = {}
    for i in xrange(N):
        for j in xrange(i):
            # compute overlaps going one way, then invert for reverse direction
            glycan1 = glycan_list[i]
            glycan2 = glycan_list[j]
            score, overlaps = overlap_score(glycan1, glycan2)
            if score > 0: # actual overlap found, add edge
                inverted_overlaps = map(__invert_mapping__, overlaps)
                pairwise_overlaps[glycan1][glycan2] = (score, overlaps)
                pairwise_overlaps[glycan2][glycan1] = (score, inverted_overlaps)
                
    ''' Repeatedly fuse along highest overlap and recompute overlaps to fused
        fragment, until no more overlaps exist among remaining fragments '''
    overlaps_exist = True
    while overlaps_exist:
        ''' Find an edge with the largest overlap in the current network '''
        best_overlap = (0, None, None) # score, glycan1, glycan2
        for glycan1 in pairwise_overlaps:
            for glycan2 in pairwise_overlaps:
                if glycan2 in pairwise_overlaps[glycan1]: # if overlap > 0 for this pair
                    score = pairwise_overlaps[glycan1][glycan2]
                    if score > best_overlap[0]:
                        best_overlap = (score, glycan1, glycan2)
        
        ''' If overlaps still present, pick a best one to fuse '''
        if best_overlap[0] > 0: # if there are overlaps still present
            score, glycan1, glycan2 = best_overlap
            overlap_mapping = pairwise_overlaps[glycan1][glycan2][0] # pick 1st optimal overlap
            fused_glycan = __overlap_combine__(glycan1, glycan2, overlap_mapping)
            
            ''' Remove all mentions of the two glycans that were fused '''
            pairwise_overlaps.pop(glycan1) # remove entry for glycan1
            del pairwise_overlaps(glycan2) # remove entry for glycan2
            for glycan in pairwise_overlaps: # remove mentions of glycan1/2 from other glycans
                if glycan1 in pairwise_overlaps[glycan]:
                    pairwise_overlaps[glycan].pop(glycan1)
                if glycan2 in pairwise_overlaps[glycan]:
                    pairwise_overlaps[glycan].pop(glycan2)
                    
            ''' Add entry for fused glycan and compute new pairwise overlaps '''
            pairwise_overlaps[fused_glycan] = {}
            for glycan in pairwise_overlaps:
                if glycan != fused_glycan: # no self edges
                    score, overlaps = overlap_score(glycan, fused_glycan)
                    if score > 0: # actual overlap found, add edge
                        inverted_overlaps = map(__invert_mapping__, overlaps)
                        pairwise_overlaps[glycan][fused_glycan] = (score, overlaps)
                        pairwise_overlaps[fused_glycan][glycan] = (score, inverted_overlaps)
                    
        else: # no more overlaps exist
            overlaps_exist = False
    
    return pairwise_overlaps.keys() # return glycans not fused away
        
def __overlap_combine__(glycan1, glycan2, overlap_mapping):
    '''
    Combines two glycans using the output from overlap_score.
    Inputs: 
        glycan1, glycan2 - Glycan objects
        overlap_mapping - dict containing {index of glycan1 -> index of glycan2}
    Output: 
        fused_glycan = glycan combining glycan1 and glycan2 through the mapping
    '''
    
    # Start with copy of first glycan
    fused_names = glycan1.names.copy()
    fused_bonds = glycan1.names.copy()
    sugar_count = len(fused_names)
    
    # Assign new indices to glycan2 nodes, based on overlap or not
    new_glycan2_indices = {}
    reverse_mapping = __invert_mapping__(overlap_mapping)
    for g2_index in xrange(len(glycan2)):
        if g2_index in reverse_mapping: # if overlapped node, inherit glycan1 index
            new_glycan2_indices[g2_index] = reverse_mapping[g2_index]
        else: # if node unique to glycan2, shift index based on unique node so far
            fused_names.append(glycan2.names[g2_index])
            fused_bonds.append({})
            new_glycan2_indices[g2_index] = sugar_count                        
            sugar_count += 1
            
    # Add bonds from glycan2 based on new node indices
    for g2_index in xrange(len(glycan2)):
        new_index = new_glycan2_indices[g2_index]
        node_bonds = glycan2.bonds[g2_index]        
        for target in node_bonds:
            bond_info = node_bonds[target]
            new_target_index = new_glycan2_indices[target]
            fused_bonds[new_index][new_target_index] = bond_info
            
    # Create and return fused Glycan
    return Glycan(fused_names, fused_bonds)
        
def __invert_mapping__(mapping):
    ''' Creates a dict where key and value are swapped '''
    inverted = {}
    for k in mapping:
        inverted[ mapping[k] ] = k
    return mapping
    
            
            
    