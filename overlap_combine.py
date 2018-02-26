# -*- coding: utf-8 -*-
'''
function overlap_combine
author Owen Chapman auberginekenobi
2/19/2018
'''
from GlycanParser import *

def overlap_combine(g1,g2,overlap):
	'''
	Combines two glycans using the the output from overlap_score.
	Inputs:
		g1, g2 - Glycan objects
		overlap - an overlap of type (Glycan, Glycan).
			Should be one of the outputs of overlap_score(g1,g2)
	Output:
		a combined glycan g1+g2
	'''
	# create a dictionary g2 index --> g1 index
	g2_g1 = {}
	# map all overlapped indices in g2 to the corresponding g1
	for i in range(len(overlap[1])):
		g2_g1[overlap[1].names[i]] = overlap[0].names[i]
	# map all other indices in g2 to new indices in g1
	# and add nodes in g2 to g1
	i=len(g1)
	for k in range(len(g2)):
		if k in g2_g1:
			continue
		else:
			g2_g1[k] = i
			i+=1
			g1.names.append(g2.names[k])
			g1.bonds.append({})
	# add g2 edges to g1
	for j in range(len(g2)):
		for key in g2.bonds[j].keys():
			if key == -1:
				continue
			elif g2_g1[key] in g1.bonds[g2_g1[j]]:
				assert(g2.bonds[j][key] == g1.bonds[g2_g1[j]][g2_g1[key]])
			else:
				g1.bonds[g2_g1[j]][g2_g1[key]] = g2.bonds[j][key]
	return g1



