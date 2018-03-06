# jtsorren

import pandas as pd
from GlycanFragmentation import *
from GlycanParser import *
from agglomerator import agglomerative_assemble
from overlap_score_fast import overlap_score

df = pd.read_csv('input_glycans.txt',sep='\t')
columns = ['accession number', 'num contigs', 'all contigs subtrees', 'exact match', 'assembled glycan']
newDf = pd.DataFrame(columns=columns)
acEntries = [] # accession number column
agEntries = [] # contigs column
numContigsEntries = [] # number of contigs
subtreeCheckEntries = [] # all contigs subtrees
exactCheckEntries = [] # only 1 contig and is subtree

for index, row in df.iterrows():
    iupacName = row['IUPAC']
    glycan = Glycan(iupacName)
    fragments = GlycanFragmentation(100, glycan)
    acEntries.append(row['accession number'])
    contigs = agglomerative_assemble(fragments)
    agEntries.append(contigs)
    numContigsEntries.append(len(contigs))
    allSubtrees = True
    for contig in contigs:
        score, overlaps = overlap_score(contig, glycan)
        if score != len(contig): # contig is not fully a subtree
            allSubtrees = False; break
    exactMatch = len(contigs) == 1 and overlap_score(contigs[0], glycan)[0] == len(glycan)
    subtreeCheckEntries.append(allSubtrees)
    exactCheckEntries.append(exactMatch)

newDf['accession number'] = acEntries
newDf['num contigs'] = numContigsEntries
newDf['all contigs subtrees'] = subtreeCheckEntries
newDf['exact match'] = exactCheckEntries
newDf['assembled glycan'] = agEntries
newDf.to_csv(r'results2.txt', header=None, index=None, sep='\t', mode='a')
