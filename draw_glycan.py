# -*- coding: utf-8 -*-
"""
Created on Sun Mar 04 20:43:44 2018

@author: jhyun_000
"""

from GlycanParser import Glycan
import networkx as nx
import matplotlib.pyplot as plt

def draw_glycan(glycan):
    ''' Uses networkx to draw a glycan as a graph, showing node labels,
        edge labels, and limited node-coloring based on sugar identity.
        Takes as input a Glycan object. '''
    
    ''' Node colormap for common sugars adapted from GlyTouCan '''
    colormap = {'D-GlcpNAc': '#6666FF',
                'L-Fucp': 'red',
                'D-Manp': '#66FF66',
                'D-Galp': 'yellow',
                'D-NeupAc': 'magenta',
                'Root': 'white'}

    ''' Map sugar names to node indices and labels '''
    n = len(glycan)
    #node_labels = {-1:'Root'}
    node_labels = {}
    for i in range(n):
        node_labels[i] = glycan.names[i]

    ''' Map bonds to edges and edge labels '''
    all_edges = []; # all edges, for rendering purposes
    drawable_edges = []; # only draw one way edges
    edge_labels = {}
    for start in range(n):
        for end in glycan.bonds[start]:
            edge = (start,end)
            label = glycan.bonds[start][end]
            all_edges.append(edge)
            if start < end or start*end < 0: # only show edge in one way, or to root
                drawable_edges.append(edge)
                bond_type, start_carbon, end_carbon = label
                bond_type = bond_type.replace("alpha","a")
                bond_type = bond_type.replace("beta","b")
                if start > -1 and end > -1: # non-root edge
                    edge_label = bond_type + "(" + str(start_carbon) + ">" + str(end_carbon) + ")"
                else:
                    edge_label = bond_type + "(" + str(start_carbon) + ">)"
                edge_labels[edge] = edge_label
    G = nx.DiGraph()
    G.add_edges_from(all_edges)
    
    ''' Assign node colors based on sugar identity '''
    node_colors = []
    for node in G.nodes():
        if node_labels[node] in colormap: # common sugar with known color
            node_colors.append(colormap[node_labels[node]])
        else: # unknown sugar, default to gray
            node_colors.append('lightgrey')
            
    ''' Render graph as matplotlib figure'''
    pos=nx.spring_layout(G)
    nx.draw_networkx_nodes(G, pos, with_labels=True, node_color=node_colors, node_size=1200)
    nx.draw_networkx_edges(G, pos, edgelist=drawable_edges)
    nx.draw_networkx_labels(G, pos, labels=node_labels)
    nx.draw_networkx_edge_labels(G, pos, edgelist=drawable_edges, edge_labels=edge_labels)
    ax = plt.gca() # for disabling axes and ticks
    [sp.set_visible(False) for sp in ax.spines.values()]
    ax.set_xticks([]); ax.set_yticks([])