import pandas as pd
import networkx as nx

'''
This file here stores all necessary functions that need to be
reused throughout the different notebooks we have
'''

def remove_threshold(file_name : str, num : int) -> nx.Graph:
    '''
    Creates a nx.graph and removes all edges with a score of less than num

    file_name -> the file that stores the network
    num -> the threshold to remove at

    returns a nx.Graph after removal of edges
    '''

    df = pd.read_csv(file_name, sep=" ")
    nodes = set(df["protein1"])

    # keeps edges with a score of greater or equal to num
    df = df[df.combined_score >= num]
    G = nx.from_pandas_edgelist(df, "protein1", "protein2")

    # the above will also remove nodes taht shouldn't be removed
    # so these nodes are added back here
    new = nodes - set(G.nodes)
    for node in new:
        G.add_node(node)

    return G 

def remove_essential(G : nx.Graph, file_name : str) -> nx.Graph:
    '''
    Removes all essential proteins
    
    G -> the nx.Graph of our network
    file_name -> the file that stores all the essential proteins

    returns the graph without these proteins
    '''

    # this file was in the Ed which details all the essential nodes to remove
    essential = pd.read_csv(file_name, header = None, usecols = [1])
    essential = essential[1].tolist()

    # here is the actual removal process
    for protein in essential:
        name = "4932." + protein
        if name in G.nodes():
            G.remove_node(name)

    return G

def parser(names : list) -> dict:
    '''
    Converts all protein names to the names in the network
    
    names -> the names of the proetins needed to change

    returns a dictionary with the corresponding names
    '''

    # proteins -> a full list of all the proteins information in our network
    proteins = pd.read_csv("network_info/4932.protein.info.v11.5.txt", sep = "\t")

    nodes = {}
    for name in names:
        nodes[name] = (proteins.loc[proteins['preferred_name'] == name])['#string_protein_id'].iloc[0]

    return nodes

def renaming_clusters(clusters : list, protein_hash : dict) -> list:
    '''
    Transforms index based names to protein names
    
    clusters -> list of the clusters, where each cluster is a list
                of proteins by index
    protein_hash -> a dictionary with mapping from index to protein name 
     
    returns cluster list with protein names
    '''

    named_clusters = []
    lpd1pos = 0
    for index,cluster in enumerate(clusters):
        named_cluster = tuple([protein_hash[node] for node in cluster])
        if "4932.YFL018C" in named_cluster:
            lpd1pos = index
        named_clusters.append(named_cluster)

    return named_clusters

def get_weight(G : nx.Graph, cluster1, cluster2) -> int:
    '''
    G -> the orginal unweighted network
    cluster1 -> a cluster after MCL
    cluster2 -> a different cluster after MCL 

    returns the number of edges from cluster1 to cluster2
    '''
    weight = 0

    for node1 in cluster1:
        for node2 in cluster2:
            if G.has_edge(node1, node2):
                weight += 1
    
    return weight

def convert_to_weighted(G_orginal, clusters):
    '''
    Transforms an unweighted graph to a weighted graph based on the number of edges it
    has between different clusters

    G_orginal -> the orginal unweighted graph
    clusters -> clusters identified when running MCL

    returns nx.Graph
    '''

    G = nx.Graph()
    
    #Creates a graph with nodes as each cluster
    nodes = [i for i in range(len(clusters))]
    G.add_nodes_from(nodes)

    #loops over each p[air of clusters and determines the corresponding weight
    i = 0
    while i < len(clusters) - 1:

        j = i + 1
        while j < len(clusters):
            weight = get_weight(G_orginal, clusters[i], clusters[j])

            # only add an edge if there orginaly were edges
            if weight != 0:
                G.add_edge(i ,j, weight = weight)
            j += 1

        i += 1

    return G


