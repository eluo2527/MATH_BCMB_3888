'''
Modified version of run.py where it can run on multiple threhold scores
and saves relevant data in JSON files
'''


# from logging import warning
from cProfile import run
import sys, os

import operator
from collections import defaultdict
import pandas as pd
from tqdm import tqdm
import networkx as nx
import markov_clustering as mc
# import numpy as np
# import scipy as sp
# import matplotlib.pyplot as plt
import func
import warnings
import json
from pprint import pprint


from joblib import Parallel, delayed
import multiprocessing

def main(threshold : int, important_nodes : dict):

    # Create the network and remove initial nodes 
    network_name = "network_info/4932_protein_links_v11_5.txt"
    print("Importing Proteins and removing essentials")
    G = func.remove_threshold(network_name, threshold)

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

    # Renames the clusters in the weighted network as w0, w1, w2, ...
    mapping = {node : f"w{node}" for node in weighted_network.nodes}
    weighted_network_rename = nx.relabel_nodes(weighted_network, mapping)

    # Finds the index of the cluster each node corresponds to

    # key = protein name, ie LPD1
    # value = cluster index ie 32
    important_clusters = {}

    for key, value in important_nodes.items():
        index = func.find_cluster(value, named_clusters)
        important_clusters[key] = index

    print(important_clusters)

    # loop over all important clusters and find betweenness scores
    for name, cluster in important_clusters.items():
        print(name, cluster)
        
        # this here creates a connected weighted network with one of the important
        # clusters as a source
        filtered_weight = func.connected_clusters(weighted_network_rename, 'w' + str(cluster))
        weighted_centrality = func.weighted_centrality(filtered_weight, 'w' + str(cluster))

        # saves the weighted centrality 
        #func.json_save(weighted_centrality, 'results/' + str(threshold) + '/' + name + '/' + name + '_weighted_centrality.json')
        df = pd.DataFrame(weighted_centrality.items())
        df.to_csv(f"results/{threshold}_{name}_weighted_centrality.csv",index=False)
        # # key = cluster w0, w1,
        # # value = centrality value
        # my_ls = list((key, value) for key, value in weighted_centrality.items())

        # print('here')

        # # this loops over the first 5 (can be changed) most important nodes
        # for i in range(5):
            
        #     new_graph = func.cluster_graph(G, named_clusters[cluster], named_clusters[int(my_ls[i][0][1:])])
        #     if nx.is_connected(new_graph):

        #         # finds betweeness score on the unweighted graph between important cluster (BCMB) and what we found
        #         b_centrality = func.between_centrality(G, named_clusters[cluster], named_clusters[int(my_ls[i][0][1:])])

        #         sorted_betweeness_names = {}
        #         for key in b_centrality:
        #             sorted_betweeness_names[func.name_change(key)] = b_centrality[key]
        #         func.json_save(sorted_betweeness_names, 'results/' + str(threshold) + '/' + name + '/' + name + '_' + my_ls[i][0] + '_betweeness.json')

def run_main():
    # These are the essential proteins that the biochemist have identified 
    # https://docs.google.com/document/d/12kaAjgjEsQtCOaRqw6g2ZNeLzN-rlzmLaGApKCdI1uc/edit 
    # E3 protein is LPD1
    names = ['LPD1', 'PDA1', 'PYC2', 'PDB1', 'PTC1', 'BAT2', 'KGD1', 'AIM22', 'PKP1', 'PTC5', 'LAT1']
    important_nodes = func.parser(names)

    # print(important_nodes)

    threshold_scores = [824,864]

    for threshold in threshold_scores:
        main(threshold, important_nodes)
    return True

def measure_nearest_cluster(threshold:int, protein_of_interest:str):
    # Create the network and remove initial nodes 
    network_name = "network_info/4932_protein_links_v11_5.txt"
    print("Importing Proteins and removing essentials")
    G = func.remove_threshold(network_name, threshold)

    # essential_proteins = "network_info/essential_proteins.csv"
    # G = func.remove_essential(G, essential_proteins)

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

    # Renames the clusters in the weighted network as w0, w1, w2, ...
    mapping = {node : f"w{node}" for node in weighted_network.nodes}
    weighted_network_rename = nx.relabel_nodes(weighted_network, mapping)

    # Finds the index of the cluster each node corresponds to

    # key = protein name, ie LPD1
    # value = cluster index ie 32


    protein_cluster = func.find_cluster(protein_of_interest, named_clusters)
    # loop over all important clusters and find betweenness scores
    print(protein_of_interest, protein_cluster)
    protein_score = {}
    score = 1 # A bit arbitrary - needs some thought
    protein_subgraph = G.subgraph(list(named_clusters[protein_cluster]))    
    e = nx.eigenvector_centrality(protein_subgraph)
    e.pop(protein_of_interest)
    in_cluster_norm_factor = sum(e.values())
    for (protein,in_cluster_score) in sorted(list(e.items()),key=lambda i: i[1],reverse=True):
        # print(func.name_change(protein),in_cluster_score*score/in_cluster_norm_factor)
        protein_score[func.name_change(protein)] = in_cluster_score*score/in_cluster_norm_factor
    return dict( sorted(protein_score.items(), key=operator.itemgetter(1),reverse=True))

def unified_list(threshold : int, protein_of_interest : str):

    # Create the network and remove initial nodes 
    network_name = "network_info/4932_protein_links_v11_5.txt"
    print("Importing Proteins and removing essentials")
    G = func.remove_threshold(network_name, threshold)

    # essential_proteins = "network_info/essential_proteins.csv"
    # G = func.remove_essential(G, essential_proteins)

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

    # Renames the clusters in the weighted network as w0, w1, w2, ...
    mapping = {node : f"w{node}" for node in weighted_network.nodes}
    weighted_network_rename = nx.relabel_nodes(weighted_network, mapping)

    # Finds the index of the cluster each node corresponds to

    # key = protein name, ie LPD1
    # value = cluster index ie 32


    protein_cluster = func.find_cluster(protein_of_interest, named_clusters)
    # loop over all important clusters and find betweenness scores
    print(protein_of_interest, protein_cluster)
    
    # this here creates a connected weighted network with one of the important
    # clusters as a source
    filtered_weight = func.connected_clusters(weighted_network_rename, 'w' + str(protein_cluster))
    weighted_centrality = func.weighted_centrality(filtered_weight, 'w' + str(protein_cluster))

    # saves the weighted centrality 
    protein_score = {}
    for (cluster,score) in sorted(list(weighted_centrality.items()),key=lambda i: i[1],reverse=True)[0:10]:
        cluster_index = cluster[1:]
        proteins = named_clusters[int(cluster_index)]
        # print(proteins)

        protein_subgraph = G.subgraph(list(proteins))
        e = nx.eigenvector_centrality(protein_subgraph)
        in_cluster_norm_factor = sum(e.values())
        for (protein,in_cluster_score) in sorted(list(e.items()),key=lambda i: i[1],reverse=True):
            protein_score[func.name_change(protein)] = in_cluster_score*score/in_cluster_norm_factor
    
    # Check cluster that protein is in
    score = 1 # A bit arbitrary - needs some thought
    protein_subgraph = G.subgraph(list(named_clusters[protein_cluster]))    
    e = nx.eigenvector_centrality(protein_subgraph)
    e.pop(protein_of_interest)
    in_cluster_norm_factor = sum(e.values())
    for (protein,in_cluster_score) in sorted(list(e.items()),key=lambda i: i[1],reverse=True):
        # print(func.name_change(protein),in_cluster_score*score/in_cluster_norm_factor)
        protein_score[func.name_change(protein)] = in_cluster_score*score/in_cluster_norm_factor
    return dict( sorted(protein_score.items(), key=operator.itemgetter(1),reverse=True))

def run_sequential(thresholds,protein_of_interest,base_file_name="proteins_by_threshold"):
    results_by_threshold = []
    for threshold in tqdm(thresholds):
        results_by_threshold.append(measure_nearest_cluster(threshold, func.name_change(protein_of_interest)))

    data = defaultdict(lambda: len(list(thresholds))*[float("nan")])
    data['threshold'] = list(thresholds)
    for index, result in enumerate(results_by_threshold): # you can list as many input dicts as you want here
        for key, value in result.items():
            data[key][index] = value

    df = pd.DataFrame(dict(data))
    # normalize df
    threshold=df['threshold']

    df=df.drop('threshold', axis=1)
    df=df.div(df.sum(axis=1), axis=0)
    df.insert(loc=0, column='threshold', value=threshold)

    print(df)
    df.to_csv(f'results/data/{protein_of_interest}_{base_file_name}.csv', index=False)

def get_multiprocess_func(protein_of_interest):
    def blockPrint():
        sys.stdout = open(os.devnull, 'w')

    # Restore
    def enablePrint():
        sys.stdout = sys.__stdout__

    def output_func(threshold):
        print("Starting threshold ", threshold)
        blockPrint()
        output = (threshold,unified_list(threshold, func.name_change(protein_of_interest)))
        enablePrint()
        print("Finishing threshold ", threshold)
        return output


    return output_func

def run_parallel(threshold,protein_of_interest,base_file_name="proteins_by_threshold"):
    
    multiprocess_func = get_multiprocess_func(protein_of_interest)
    num_cores = multiprocessing.cpu_count()-1
    results = Parallel(n_jobs=num_cores)(delayed(multiprocess_func)(threshold) for threshold in thresholds)
    results_by_threshold = [result[1] for result in sorted(results, key=lambda result: result[0])]


    data = defaultdict(lambda: len(list(thresholds))*[float("nan")])
    data['threshold'] = list(thresholds)
    for index, result in enumerate(results_by_threshold): # you can list as many input dicts as you want here
        for key, value in result.items():
            data[key][index] = value

    df = pd.DataFrame(dict(data))
    # normalize df
    threshold=df['threshold']

    df=df.drop('threshold', axis=1)
    df=df.div(df.sum(axis=1), axis=0)
    df.insert(loc=0, column='threshold', value=threshold)

    print(df)
    df.to_csv(f'results/data/{protein_of_interest}_{base_file_name}.csv', index=False)
    return df


if __name__ == '__main__':
    # names = ['LPD1', 'PDA1', 'PYC2', 'PDB1', 'PTC1', 'BAT2', 'KGD1', 'AIM22', 'PKP1', 'PTC5', 'LAT1'] # https://docs.google.com/document/d/12kaAjgjEsQtCOaRqw6g2ZNeLzN-rlzmLaGApKCdI1uc/edit 
    # # base_file_name = "proteins_by_threshold_824-864-1_25_proteins"
    # # for name in names:
    # #     thresholds = range(824,865,1)
    # #     run_parallel(thresholds,name,base_file_name=base_file_name)
    # df = pd.DataFrame(unified_list(824, func.name_change('PDA1')).items())
    # df.to_csv(f'results/data/PDA1_whole_network_at_824', index=False)
    run_main()