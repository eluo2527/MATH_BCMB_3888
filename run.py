'''
A python file designed to output lists of proteins as per main.ipynb

For fast graphs see original notebook.

'''
################################################################################
# IMPORTS
from logging import warning
import networkx as nx
import markov_clustering as mc
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd
import func
import warnings
from pprint import pprint

################################################################################
# These are the essential proteins that the biochemist have identified 
# https://docs.google.com/document/d/12kaAjgjEsQtCOaRqw6g2ZNeLzN-rlzmLaGApKCdI1uc/edit 
# E3 protein is LPD1

names = ['LPD1', 'PDA1', 'PYC2', 'PDB1', 'PTC1', 'BAT2', 'KGD1', 'AIM22', 'PKP1', 'PTC5', 'LAT1']
important_nodes = func.parser(names)

# Create the network and remove initial nodes 
network_name = "network_info/4932_protein_links_v11_5.txt"
print("Importing Proteins and removing essentials")
G = func.remove_threshold(network_name, 700)

essential_proteins = "network_info/essential_proteins.csv"
G = func.remove_essential(G, essential_proteins)

# Find the clusters
# Here we are going to lose the protein names as the matrix gets assigned to their index. 
# So we recover that with a hash table (dictionary)

with warnings.catch_warnings():
    print("Creating adjacency matrix ")
    warnings.simplefilter("ignore")
    adj_matrix = nx.adjacency_matrix(G) 

#Runs MCL 
# run with default parameters   
print("Performing MCL") 
result = mc.run_mcl(adj_matrix)         
clusters = mc.get_clusters(result) 

# Create a hash table that takes takes a number and returns a protein name
protein_hash = {}
for index, node in enumerate(G.nodes):
    protein_hash[index] = node

# Renaming proteins in the clusters
named_clusters = func.renaming_clusters(clusters, protein_hash)
# Create weighted network of clusters
print("Creating weighted network of clusters")
weighted_network = func.convert_to_weighted(G, named_clusters)
print(weighted_network)


mapping = {node : f"w{node}" for node in weighted_network.nodes}
weighted_network_rename = nx.relabel_nodes(weighted_network, mapping)

#Make a graph with only nodes in a certain cluster 
#cluster number 32 in this case as it has LPD1

print("Creating central clusters subgraph")
cluster_network = func.cluster_graph(G, named_clusters[32], named_clusters[93])

b_centrality = func.between_centrality(G, named_clusters[32], named_clusters[93])

# This is comparing cluster 32 and 93

# sorts the dictionary such that highest centrality is first
sorted_betweeness = dict(sorted(b_centrality.items(), key = lambda x: -x[1]))

sorted_betweeness_names = {}
for key in sorted_betweeness:
    sorted_betweeness_names[func.name_change(key)] = sorted_betweeness[key]
func.json_save(sorted_betweeness_names, "betweeness.json")

# gives first 5 maybe or more to BCMB students
prefered_names = func.parser(sorted_betweeness.keys(), False)
# func.json_save(prefered_names, "prefered_names.json")

# this here is just within cluster 32

print("Performing between-ness centrality on cluster subgraph")
c32 = nx.betweenness_centrality(func.cluster_graph(G, named_clusters[32]))

# sorts the dictionary such that highest centrality is first
c32_sorted = dict(sorted(c32.items(), key = lambda x: -x[1]))
print(c32_sorted)
# gives first 5 maybe or more to BCMB students
c32_sorted = func.parser(c32_sorted.keys(), False)
print(G)
# Iscolated the middle section of the graph as seen in an earlier picture
# Used clsuter 32 as the source as it contains LPD1 

filtered_weight = func.connected_clusters(weighted_network_rename, 'w32')

print("Creating Weighted Centrality")
weighted_centrality = func.weighted_centrality(filtered_weight, 'w32')

#non_zero_weighted = {protein:weight for (protein,weight) in weighted_centrality.items() if weight!=0}

func.json_save(weighted_centrality, "weighted_centrality.json")
