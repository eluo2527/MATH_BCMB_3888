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
    
#     essential_proteins = "network_info/essential_proteins.csv"
#     G = func.remove_essential(G, essential_proteins)

    adj_matrix = nx.adjacency_matrix(G) 
    protein_hash = {}
    for index, node in enumerate(G.nodes):
        protein_hash[index] = node
        
    result = mc.run_mcl(adj_matrix)         
    clusters = mc.get_clusters(result) 
    
    named_clusters = func.renaming_clusters(clusters, protein_hash)
    
    target = func.find_cluster("4932.YFL018C", named_clusters)
    
    weighted_network = func.convert_to_weighted(G, named_clusters)
    
    mapping = {node : f"w{node}" for node in weighted_network.nodes}
    weighted_network_rename = nx.relabel_nodes(weighted_network, mapping)
    filtered_weight = func.connected_clusters(weighted_network_rename, mapping[target])
    weighted_centrality = func.weighted_centrality(filtered_weight, mapping[target])
    
    important = list(weighted_centrality.items())
    
    return important

def ranking(graph, groups, weight):
    names = func.parser(graph.nodes(), False)
    
    essential = pd.read_csv("network_info/essential_proteins.csv", header = None, usecols = [1])
    essential = essential[1].tolist()
    for i in range(len(essential)):
        key = "4932." + essential[i]
        if key in names.keys():
            essential[i] = names[key]
    
    cluster = weight[0][0]
    index = int(cluster[1:])
    
    new_graph = graph.subgraph(groups[index])
    
    bet = nx.betweenness_centrality(new_graph)
    bet = (sorted(bet.items(), key=lambda item: -item[1]))

    eig = nx.eigenvector_centrality(new_graph)
    eig = (sorted(eig.items(), key=lambda item: -item[1]))
    
    for i in range(len(eig)):
        eig[i] = (names[eig[i][0]],eig[i][1])
        bet[i] = (names[bet[i][0]],bet[i][1])
        
    combine = {}
    for i in range(len(eig)):
        combine[eig[i][0]] = i + 1

    for i in range(len(bet)):
        combine[bet[i][0]] += (i + 1)
    
    combine = (sorted(combine.items(), key=lambda item: item[1]))
    
    for protein in combine:
        if protein[0] in essential:
            combine.remove(protein)
    
    return combine
