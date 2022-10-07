# MATH_BCMB_3888

MATH3888.ipynb is roughly the same as the colab one except I have modified the file neames so that it can run

main.ipynb is the one I'm mainly using 

func.py is a file that stores a lot of reusable functions so write functions there and comment them

## Description
`func.py` is a library of functions used in our algorithm
`run.py` is the main code, which is executed at runtime
`.json` files are the **output** files

## @Bio-chemists
`.json` files have the relevant values, let us know what can give them more clarity.
`betweeness.json` contains the list of proteins that have the strongest effect on LPD1.

## Plan for Mid-Sem Break
1. Scott will mess around with threshold score
    - Identifying proteins that are impo rtant at multiple threshold scores
2. Patrick is implement other cluster algorithm - namely walktrap
    - Using package cdlib and exploring certain
    - ! Move stuff over to normal .py
    - ! Keep juypter for graphs
    - Optional function to save graph
3. Edward refactor code to always find connected code. Shake the code tree.
4. Avon looking at implementing alternative centrality measures
    - Possibly eigenvector centrality, spectral clustering? (Not very useful)
## Plans for October 7-14th
1. Remove subgraph part of alogorithm.

    After clustering we have nodes 1 to n. Now for some target - say LAPD1 - we can assign centrality scores w_j, for each cluster j in relation to LAPD1.
    - Currently we are using betweeness - maybe others are better
    
    The issue we have is we then create a subgraph that contains only the cluster of LAPD1 and the j with the heighest w_j. Then we perform between-ness on just the subgraph. Georg said this was bad and the better way to do it was to
    a. Find the between-ness of proteins in cluster j with the whole network G.
    b. Use the same idea as a. but some other centrality measure, such as bottleneck. Or,
    c. Perform degree centrality, page-rank or something else on the subgraph of just j. i.e which protein is most important to the structure of j.

    - a. doesn't seem to hard to implement so we do that first. b. and c. should come later.
2. We need to start looking at implementing different centrality measures and plotting them against eachother. (Avon is working on this)
    a. Main point is that if there is a monotonic relationship between degree and between-ness then we learn nothing new. We should try to find patterns and graphs. This will be useful for the algorithm and in the report according to Georg.
3. Finally, we need to sort out how we output the data so that it is as simple for the Biochemists to use. (I'll be doing this)
   - My currently plan is to tidy up the results folder into a single json (maybe csv). Then to perform one last step which is to create a total score for each protein. Namely, if we multiply the score of the protein by the score of the cluster the protein is in, then hopefully we can compare proteins from different clusters and create a single list.
4. Relatively minor, but I'll do my best to fix the file naming stuff. Why does info use `4932.protein` and links use `4932_` etc. etc.