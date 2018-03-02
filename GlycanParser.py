# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 01:16:39 2018

@author: jhyun_000
@author: auberginekenobi Owen Chapman
"""
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
                if glycan1.names[n]==glycan2.names[m]: # are the sugars the same
                    # have to check that no bonds conflict
                    carbon_to_bondinfo1 = {v[1]: v for k, v in glycan1.bonds[n].items()}
                    carbon_to_bondinfo2 = {v[1]: v for k, v in glycan2.bonds[m].items()}
                    can_fuse = True
                    for key in carbon_to_bondinfo1:
                        if key in carbon_to_bondinfo2:
                            if carbon_to_bondinfo1[key] != carbon_to_bondinfo2[key]:
                                can_fuse = False
                    if can_fuse:
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

class Glycan:
    '''
    Defines an unrooted glycan as an acyclic undirected graph.
    Implemented as a list of sugar names, and
    a corresponding list of bonds for each sugar.
    Unique indices denote each sugar.
    '''
    bonds = []
    names = []
    def __len__(self):
        return len(self.names)
    def __repr__(self):
        return self.__str__()
    def __str__(self):
        string=''
        for i in range(len(self.names)):
            string+=str(i)+' '+str(self.names[i])+' '+str(self.bonds[i])#+'\n'
        return string
    def equals(self, other):
        '''
        Checks if this Glycan encodes the same molecule as another
        Glycan. This is NOT an __eq__ function, which returns true only
        if they are the same object.
        '''
        return overlap_score(self,other)[0] == max(len(self),len(other))
    def __overload__(*functions):
        return lambda *args, **kwargs: functions[len(args)-1](*args, **kwargs)    
    def __init1__(self):
        pass
    def __init2__(self, iupac):
        self.names, self.bonds = parse_iupac(iupac, output=False)
    def __init3__(self, names, bonds):
        self.bonds=bonds
        self.names=names
    '''
    May initialize with zero, one or two arguments:
    zero: empty glycan
    one: iupac string
    two: names, bonds as list of dictionaries
    '''
    __init__=__overload__(__init1__,__init2__,__init3__)

    def isBond(self,index_i,index_j):
        '''
        Check if bond exists between sugars at indices i and j
        '''
        try:
            return bool(self.bonds[index_i][index_j])
        except:
            return False

def parse_iupac(iupac, output=True):
    ''' Parses an IUPAC extended string to construct a glycan graph, 
        represented as two arrays, names and bonds.
        * names is a list of strings, mapping sugarID to sugar name 
          (i.e. names[i] = name of the ith sugar )
        * bonds is a list of dicts that with bond types and start/end carbons 
          (i.e. if there is a bond from the ith sugar to the jth sugar, 
           then bonds[i][j] = (bond type, start carbon, end carbon)
        * The root sugar will have a bond to index -1. '''
    count = 0 
    stack = [] # previous sugar branches
    active = [] # list of sugars sharing the same parent sugar, to be encountered
    children = [] # store these as children to the next sugar encountered
    
    bonds = []; names = []
    while len(iupac) > 0:
        sugar, remaining = __pop_next_sugar__(iupac)
        sugarName, bondType, bondStart, bondEnd = __parse_sugar_string__(sugar)
        sugarID = count; count += 1
        names.append(sugarName)
        bonds.append( {-1:(bondType, bondStart, bondEnd)} )
        active.append(sugarID)
        
        while len(children) > 0:
            childID = children.pop(0)
            childToParent = bonds[childID][-1] # get temporary parent bond
            bonds[childID][sugarID] = childToParent # assign actual parent ID
            del bonds[childID][-1] # delete parent placeholder
            parentToChild = (childToParent[0], childToParent[2], childToParent[1])
            bonds[sugarID][childID] = parentToChild # store in other direction
        
        if len(remaining) > 0: # not last sugar
            if remaining[0] == '-': # moving to parent sugar
                children = active; active = []
                remaining = remaining[1:]
            elif remaining[0] == '[': # starting new sibling sugar branch
                stack.append(active); active = []
                remaining = remaining[1:]
            elif remaining[0:2] == ']-': # closing a sibling sugar branch
                children = stack.pop() + [sugarID]
                active = []
                remaining = remaining[2:]
            else: # to single sibling sugar
                pass
        iupac = remaining
        
    if output: # print glycan structure
        for i in xrange(len(bonds)):
            print (i, names[i], bonds[i])
        
    return names, bonds

def __parse_sugar_string__(sugar):
    ''' Parses a sugar and bond as in IUPAC extended '''
    bondType = sugar.split("-")[0]
    nonBond, bondPositions = sugar.split("(")
    bondStart = int(bondPositions[0])
    bondEnd = int(bondPositions[3]) if len(bondPositions) > 3 else 0
    sugarName = "-".join(nonBond.split("-")[1:])
    sugarName = sugarName[:-1] if sugarName[-1] == '-' else sugarName
    return sugarName, bondType, bondStart, bondEnd

def __pop_next_sugar__(iupac):
    ''' Pops the next sugar in the glycan IUPAC extended string 
        Returns the string for the sugar and what's left of iupac '''
    closeParen = iupac.find(")")
    if closeParen == -1: # not found
        return iupac, ""
    else:
        return iupac[:closeParen+1], iupac[closeParen+1:]
