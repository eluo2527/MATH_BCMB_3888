# MATH_BCMB_3888

This github repository contains research on community finding algorithms, k shortest paths and eigenvector centrality measures to identify potentially important proteins in known diseases. Our research is based on the Maple Syrup Urine Disease (MSUD), with our network being modelled through the use of protein-protein interactions of Yeast saccharomyces cerevisiae. 

## Navigation
To navigate to relevant information the structure of our code is as follows:
1. All relevant results can be found in the results folder
    - The `deprecated` folder stores our initial results for threshold scores of 600, 700, 800 and 900. Within each of these directories there are subdirectories relating to the important proteins that the BCMB students have given us. 
    - The `data` folder stores all the current information related to our final algorithm with the `general extra results - includes long run results` folder containing all the extra information that might be useful
    - Lastly, the `graphs` folder contains plots related to centrality measures and threshold values of different proteins

2. All relevant code can be found at:
    - `run2.py` stores the main source code for our final results. Note that there are other files relevant to analysis, however this was the main script which ran our final algorithm
    - `func.py` stores all the relevant functions that we used to restructure and apply centrality measures to the network.
