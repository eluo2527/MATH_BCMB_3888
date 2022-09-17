import pandas as pd
import networkx as nx

'''
This file here stores all necessary functions that need to be
reused throughout the different notebooks we have
'''

def remove_threshold(file_name : str, num : int) -> nx.graph:
    '''
    Creates a nx.graph and removes all edges with a score of less than num

    file_name -> the file that stores the network
    num -> the threshold to remove at

    returns a nx.graph after removal of edges
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

def remove_essential(G : nx.graph, file_name : str) -> nx.graph:
    '''
    Removes all essential proteins
    
    G -> the nx.graph of our network
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

def parser(names : list) -> list:
    '''
    Converts all protein names to the names in the network
    
    names -> the names of the proetins needed to change

    returns a list with the corresponding names
    '''

    # proteins -> a full list of all the proteins information in our network
    proteins = pd.read_csv("network_info/4932.protein.info.v11.5.txt", sep = "\t")

    nodes = [] 
    for name in names:
        nodes.append((proteins.loc[proteins['preferred_name'] == name])['#string_protein_id'].iloc[0])

    return nodes



