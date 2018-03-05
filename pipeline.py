# jtsorren

import pandas as pd
from GlycanFragmentation import *
from GlycanParser import *
from agglomerator import agglomerative_assemble

df = pd.read_csv('input_glycans.txt',sep='\t')
columns = ['accession number', 'assembled glycan']
newDf = pd.DataFrame(columns=columns)
acEntries = []
agEntries = []
for index, row in df.iterrows():
    iupacName = row['IUPAC']
    glycan = Glycan(iupacName)
    fragments = GlycanFragmentation(100, glycan)
    acEntries.append(row['accession number'])
    agEntries.append(agglomerative_assemble(fragments))
newDf['accession number'] = acEntries
newDf['assembled glycan'] = agEntries
newDf.to_csv(r'results.txt', header=None, index=None, sep='\t', mode='a')