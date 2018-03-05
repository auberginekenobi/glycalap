# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 11:59:51 2018

@author: jhyun95
"""

def overlap_score(glycan1, glycan2):
    '''
    Polynomial time approach to finding the largest overlap of the same 
    sugars and bonds between two glycans. Overlap scoring is defined as 
        +1 if the same sugars exist in the same position on the tree;
        +0 if the same edges are in the same position on their respective trees; and
        -Inf if edges or nodes are different.
    Note that because each edge must come with a corresponding node, and the mismatch score is
        -Inf, the relative weights of node and edge matching don't matter.
    Inputs: Two Glycan objects
    Outputs:
        i: score of longest overlap
        overlap_mappings: list of dicts, that map indices in glycan1:glycan2
                          that constitute an overlap / shared subtree.
    '''
    
    bonds1 = glycan1.bonds; bonds2 = glycan2.bonds
    names1 = glycan1.names; names2 = glycan2.names
    
    ''' Construct a list of seeds (s1,s2)1, (s1,s2)2, ... in which for all
        sugar[s1] in glycan1 == sugar[s2] in glycan2 by name (sugar matching) '''
    seeds = []
    for s1 in range(len(glycan1)): 
            for s2 in range(len(glycan2)):
                if names1[s1] == names2[s2]:
                    seeds.append( (s1,s2) )
        
    ''' From each seed, extend out through neighbors. If fusing is possible
        remove from seeds newly fused nodes, since they would have generated
        the same fusion. '''
    overlap_mappings = [] # maps indices in glycan 1 to glycan 2
    bad_seeds = [] # using these nodes leads to an alignment with conflicting bonds
    while len(seeds) > 0:
        s1, s2 = seeds.pop()
        if can_nodes_fuse(s1,s2,glycan1,glycan2): # check if the seeds can actually fuse
            current_mapping = {s1:s2}
            unextended_pairs = [(s1,s2)]
            while len(unextended_pairs) > 0: 
                # pairs that are aligned, but not checked for further extension of overlap
                s1,s2 = unextended_pairs.pop()
                for s1_adj in bonds1[s1]: # iterate through s1 immediate neighbors
                    for s2_adj in bonds2[s2]: # iterate through s2 immediate neighbors
                        is_same_sugar = names1[s1_adj] == names2[s2_adj]
                        is_same_bond = bonds1[s1][s1_adj] == bonds2[s2][s2_adj]
                        not_used = (s1_adj,s2_adj) in seeds + bad_seeds
                        if is_same_sugar and is_same_bond and not_used: # candidate fusable nodes
                            is_fusable = can_nodes_fuse(s1_adj,s2_adj,glycan1,glycan2)
                            if is_fusable: # can fuse here, update node mapping
                                current_mapping[s1_adj] = s2_adj
                                unextended_pairs.append( (s1_adj, s2_adj) )
                                seeds.remove( (s1_adj, s2_adj) ) # remove seed, would produce same overlap
                            else: # conflict, bond occupied by different carbons
                                unextended_pairs = [] # stop aligning with this seed
                                current_mapping = {} # alignment no longer valid
                                bad_seeds += current_mapping.items()
            if len(current_mapping) > 0: # if not reset by a conflict
                overlap_mappings.append(current_mapping)
        else: # if the seeds aligning generates a conflicting bond
            bad_seeds.append((s1,s2))
    
    ''' Return overlap mappings of largest size '''
    if len(overlap_mappings) == 0:
        return 0, []
    largest_overlap_size = max(map(len, overlap_mappings))
    maximum_overlap_mappings = []
    for mapping in overlap_mappings: # get only maximum size mappings
        if len(overlap_mappings[mapping] == largest_overlap_size):
            maximum_overlap_mappings.append(mapping)
    return largest_overlap_size, maximum_overlap_mappings
                    
def can_nodes_fuse(v1, v2, glycan1, glycan2):
    ''' Checks immediate neighbors of v1 in glycan1 and v2 in glycan2
        to verify they can fuse. O(1) run time, since maximum degree
        is limited by sugar carbon count. ''' 
    # Special cases where one or more of the nodes is root, i.e. vi = -1
    if v1 * v2 < 0: # if comparing root to non-root, fail
        return False
    elif v1 == -1 and v2 == -1: # comparing root to root
        root_sugar1, root_bond1 = get_root_bond(glycan1)
        root_sugar2, root_bond2 = get_root_bond(glycan2)
        same_root_bond = root_bond1 == root_bond2 
        same_root_sugar = glycan1.names[root_sugar1] == glycan2.names[root_sugar2]
        # check the sugar connected to root is same, and bonded same way
        return same_root_bond and same_root_sugar
    # Construct reverse mappings of glycan1.bonds and glycan2.bonds
    bond_to_node1 = {}; bond_to_node2 = {}
    for v in glycan1.bonds[v1]:
        bond_to_node1[glycan1.bonds[v1][v]] = v
    for v in glycan2.bonds[v2]:
        bond_to_node2[glycan2.bonds[v2][v]] = v
    # Check if immediate neighbors along identical bonds conflict
    for bond1 in bond_to_node1:
        for bond2 in bond_to_node2: 
            same_bond = bond1[:2] == bond2[:2]
            if same_bond: # same starting carbon occupied
                v1_adj = bond_to_node1[bond1]
                v2_adj = bond_to_node2[bond2]
                # name of -1 indexed sugar is always 'root'
                v1_adj_name = glycan1.names[v1_adj] if v1_adj >= 0 else 'root'
                v2_adj_name = glycan2.names[v2_adj] if v2_adj >= 0 else 'root'
                if v1_adj_name != v2_adj_name: # if conflicting sugar along same bond
                    return False
    return True
    
def get_root_bond(glycan):
    ''' Gets the index of the node connected to root and bond info '''
    for i in xrange(len(glycan)):
        if -1 in glycan.bonds[i]: # is root
            return i, glycan.bonds[i][-1]
    return None