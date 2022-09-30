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