import networkx as nx
import markov_clustering as mc
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd
import func
from pprint import pprint

def importance(confidence):
    network_name = "network_info/4932_protein_links_v11_5.txt"
    G = func.remove_threshold(network_name, confidence)

    essential_proteins = "network_info/essential_proteins.csv"
    G = func.remove_essential(G, essential_proteins)

    protein_hash = {}
    for index, node in enumerate(G.nodes):
        protein_hash[index] = node
        
    result = mc.run_mcl(adj_matrix)         
    clusters = mc.get_clusters(result) 
    
    named_clusters = func.renaming_clusters(clusters, protein_hash)
    
    target = func.find_cluster("4932.YFL018C", named_clusters)
    
    connected = []
    for i in range(len(named_clusters)):
        if i != target:
            cluster_network = func.cluster_graph(G, named_clusters[target], named_clusters[i])
            if nx.is_connected(cluster_network):
                connected.append(i)
    
    b_centrals = []
    for i in connected:
        b_centrality = func.between_centrality(G, named_clusters[target], named_clusters[i])
        list_between = list(zip(b_centrality.keys(), b_centrality.values()))
        for tup in list_between:
            b_centrals.append(tup)
        
        sorted_centrals = sorted(b_centrals, key=lambda tup: tup[1])
        sorted_centrals.reverse()
    
    return sorted_centrals
