# -*- coding: utf-8 -*-
'''
function agglomerative_assemble
author Owen Chapman auberginekenobi
2/25/2018
'''

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
	contigs = glycan_list
	# find pairwise overlap scores across all glycan pairs
	distance_matrix = {} # distance : [(Glycan, Glycan)]
	for i in range(len(glycan_list)):
		gi = glycan_list[i]
		for j in range(i,len(glycan_list)):
			gj = glycan_list[j]
			print(i,j)
			d,_ = overlap_score(gi,gj)
			if d in distance_matrix:
				distance_matrix[d].append((gi,gj))
			else:
				distance_matrix[d] = [(gi,gj)]

	# while there exists a nonzero overlap score:
	max_score = max(distance_matrix.keys())
	while max_score>0:
		# Find the two glycans with the highest overlap
		# (arbitrary selection)
		g1, g2 = distance_matrix[max_score].pop()
		# Combine them
		# update datastructures of distances and contigs
	return contigs