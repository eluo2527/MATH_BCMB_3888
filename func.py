from typing import OrderedDict
import pandas as pd
import networkx as nx

'''
This file here stores all necessary functions that need to be
reused throughout the different notebooks we have
'''

def remove_threshold(file_name : str, num : int) -> nx.Graph:
    '''
    Creates a nx.graph and removes all edges with a score greater than or equal to num

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

def parser(names : list, val = True) -> dict:
    '''
    Converts all protein names to the names in the network
    
    names -> the names of the proetins needed to change
    val -> determies the direction of the mapping

    returns a dictionary with the corresponding names
    '''

    # proteins -> a full list of all the proteins information in our network
    proteins = pd.read_csv("network_info/4932.protein.info.v11.5.txt", sep = "\t")

    nodes = {}
    for name in names:

        if val:
            # changes from perferred to 4xxx
            nodes[name] = (proteins.loc[proteins['preferred_name'] == name])['#string_protein_id'].iloc[0]
        else:
            # changes from 4xxx to perferred 
            nodes[name] = (proteins.loc[proteins['#string_protein_id'] == name])['preferred_name'].iloc[0]

    return nodes

def name_change(node):
    '''
    Changes the names of the clusters to the names of the proteins
    '''
    proteins = pd.read_csv("network_info/4932.protein.info.v11.5.txt", sep = "\t")
    if node[0:5] == "4932.":
        # changes from 4xxx to perferred 
        return(proteins.loc[proteins['#string_protein_id'] == node])['preferred_name'].iloc[0]
    else:
        # changes from perferred to 4xxx
        return(proteins.loc[proteins['preferred_name'] == node])['#string_protein_id'].iloc[0]
       

def renaming_clusters(clusters : list, protein_hash : dict) -> list:
    '''
    Transforms index based names to protein names
    
    clusters -> list of the clusters, where each cluster is a tuple
                of proteins by index
    protein_hash -> a dictionary with mapping from index to protein name 
     
    returns cluster list with protein names
    '''

    named_clusters = []
    for cluster in clusters:
        named_cluster = tuple([protein_hash[node] for node in cluster])
        named_clusters.append(named_cluster)

    return named_clusters

def get_weight(G : nx.Graph, cluster1 : tuple, cluster2 : tuple) -> int:
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
    if weight == 0:
        return weight
    return 1/weight

def convert_to_weighted(G_orginal : nx.Graph, clusters : list) -> nx.Graph:
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

def find_cluster(name : str, clusters : list) -> int:
    '''
    Finds the cluster that the protein is associated with

    name -> name of the protein
    clusters -> clusters returned from the MCL algorithm 
    '''
    idx = 0
    while idx < len(clusters):
        if name in clusters[idx]:
            return idx
        idx += 1

    return idx

def cluster_graph(G : nx.Graph, *clusters : tuple) -> nx.Graph:
    '''
    Creates a subgraph with nodes only in the clusters specified

    G -> the orginal unweighted graph
    clusters -> the clusters of proteins that need to be examined
    '''

    # this combines all the clusters into a single cluster 
    all_cluster = [node for cluster in clusters for node in cluster]

    # creates a new graph with only the nodes specified along with edges
    new_graph = G.subgraph(all_cluster)

    if nx.is_connected(new_graph)==False:
        print("Graph is not connected. Pruning to largest connected component")
        largest_component = [x for x in nx.connected_components(new_graph)][0]
        new_graph = new_graph.subgraph(largest_component)
    return new_graph

def between_centrality(G_orginal : nx.Graph, cluster1 : tuple, cluster2 : tuple) -> dict:
    '''
    Calculates the betweenness centrality for each node in the 2 clusters specified.

    G_orginal -> the orginal unweighted graph
    cluster1 -> a cluster after MCL
    cluster2 -> a different cluster after MCL 

    returns a dictionary with nodes as keys and its betweenness cemtrality as its value
    '''

    # calls the function cluster_graph to create a graph with only nodes in the 2 clusters
    G = cluster_graph(G_orginal, cluster1, cluster2)

    # stores all paths identified with 
    # -> key in the form (source, sink)
    # -> value which is a list of the paths 
    all_paths = {}

    # the total amount of paths 
    total = 0
    
    # initerates over all pairs in each cluster
    for node1 in cluster1:
        for node2 in cluster2:

            # finds all shortest paths from source (node1) to sink (node2)
            paths = list(nx.all_shortest_paths(G, node1, node2))
            all_paths[(node1, node2)] = paths
            total += len(paths)
    
    # initialises a dictionary with the keys as nodes and the value as 0
    merged = cluster1 + cluster2
    centrality_value = dict.fromkeys(merged, 0)

    # finds the amount of times a node appears in the shortest paths exlusive of the source and sink
    for key, val in all_paths.items():
        for path in val:
            for node in path[1:-1]:
                centrality_value[node] += 1

    # divide each value by the total to normalise
    for key, val in centrality_value.items():
        centrality_value[key] /= total

    return centrality_value

def shortest_path(G: nx.Graph, source : str, sink : str) -> list:
    '''
    Finds the shortest path from the source to the sink
    '''
    pass

def connected_clusters(G : nx.Graph, source : int) -> nx.Graph:
    '''
    This function takes in the orginal weighted graph and a source cluster, and finds the
    connected component which contains the cluster

    G -> the orginal unweighted graph
    source -> the cluster which determines the connected component (the index of the cluster)

    returns the connected component of that node
    '''

    nodes = nx.node_connected_component(G, source)
    new_graph = G.subgraph(nodes)

    return new_graph

def weighted_centrality(G_weighted : nx.Graph, source : tuple) :
    '''
    Takes in a weighted graph and finds the betweenness centrality between a source cluster
    and all other clusters within the connected network

    G_weighted -> a connected weighted graph
    source -> a cluster in the graph

    returns a list of tuples, where each tuple is (source)
    '''

    # this will hold the total number of paths
    total = 0
    all_paths = {}

    #need to have all paths starting from source
    for node in G_weighted.nodes:
        if node != source:
            paths = list(nx.all_shortest_paths(G_weighted, source, node, weight = 'weight'))
            total += len(paths)
            all_paths[node] = paths
    

    centrality_value = dict.fromkeys(G_weighted.nodes, 0)

    # finds the amount of times a node appears in the shortest paths exlusive of the source and sink
    for key, val in all_paths.items():
        for path in val:
            for node in path[1:-1]:
                centrality_value[node] += 1

    # divide each value by the total to normalise
    for key, val in centrality_value.items():
        centrality_value[key] /= total

    centrality_value_sort = dict(sorted(centrality_value.items(), key=lambda item: -item[1]))

    return centrality_value_sort


def pretty_draw(G: nx.Graph,  *highlight_nodes: str, with_highlight_labels = False, with_other_labels= False,only_connected=False) -> None:
    '''
    Provide the function with a node,
    Add highlight nodes and it will paint these red,
    Add with_labels = True to show the labels
    Add only_connected = True to only show the connected nodes
    '''
    import matplotlib.pyplot as plt

    if only_connected:
        components = list(nx.connected_components(G))
        nodes = set()
        for component in components:
            if set(highlight_nodes) & component != set():
                nodes |= component
        G0 = G.subgraph(nodes)

        if len(highlight_nodes)==0:
            G0 = G.subgraph(max(components, key=len))
    else:
        G0 = G
    node_colourmap = []
    node_sizemap = []
    highlight_labels = {}
    other_labels = {}
    for node in G0:
        if node in highlight_nodes:
            node_colourmap.append((1,0.2,0.2,0.8))
            node_sizemap.append(300)
            if with_highlight_labels: highlight_labels[node] = node
        else:
            node_colourmap.append((0.4,0.4,0.8,0.6))
            node_sizemap.append(200)
            if with_other_labels: other_labels[node] = node
    
    edge_colourmap = []


    for edge in G0.edges:
        if len(highlight_nodes) == 1:
            if (set([highlight_nodes]) & set(edge)) != set():
                edge_colourmap.append((1,0,0,1))
            else:
                edge_colourmap.append((0,0,0,0.2))
        elif len(highlight_nodes) >  1:
            if (set(highlight_nodes) & set(edge)) != set():
                edge_colourmap.append((1,0,0,1))
            else:
                edge_colourmap.append((0,0,0,0.2))
        else:
            edge_colourmap.append((0,0,0,0.2))
    plt.figure(1,figsize=(30,30))
    seed = 1
    pos = nx.spring_layout(G0, seed=seed)
    nodes = nx.draw_networkx_nodes(G0, pos=pos, node_size=node_sizemap, node_color=node_colourmap)
    edges = nx.draw_networkx_edges(
        G0,
        pos=pos,
        node_size=300,
        edge_color=edge_colourmap,
        width=1,
    )
    if with_highlight_labels: nx.draw_networkx_labels(G0,pos,labels=highlight_labels,font_size=10,bbox=dict(color=(1,1,1,0.6)),font_weight='bold',font_family='sans-serif')
    if with_other_labels: nx.draw_networkx_labels(G0,pos,labels=other_labels,font_size=10,font_weight='bold',font_family='sans-serif')

    return plt.show()

def json_save(string: str, filename: str) -> None:
    '''
    Saves a string to a json file
    '''
    import json
    with open(filename, 'w') as f:
        json.dump(string, f, indent=4)
