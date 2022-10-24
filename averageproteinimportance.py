import pandas as pd
import numpy as np
import csv
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
def average_importance(protein_of_interest):
    df = pd.read_csv(f'results/data/{protein_of_interest}_proteins_by_threshold.csv')
    average = {}
    for column in df:
        values = np.array(df[column])
        values[np.isnan(values)] = 0
        if column != 'threshold':
            average[column] = np.median(values)
    for key,value in average.items():
        average[key] = value/sum(average.values())

    average = dict(sorted(average.items(), reverse=True, key= lambda x: x[1]))

    df = pd.DataFrame(average.items(), columns=['threshold', 'importance'])
    df.to_csv(f'results/data/{protein_of_interest}_average_importance.csv', index=False)
    return True

def average_importance_plot(protein_of_interest,names):
    df = pd.read_csv(f'results/data/{protein_of_interest}_proteins_by_threshold_824-864-1_25_proteins.csv')
    average = {}
    for column in df:
        values = np.array(df[column])
        values[np.isnan(values)] = 0
        if column != 'threshold':
            average[column] = np.median(values)
    for key,value in average.items():
        average[key] = value/sum(average.values())

    average = dict(sorted(average.items(), reverse=True, key= lambda x: x[1]))

    # plot the keys of average and the values of average as a bar chart
    x = list(average.keys())[0:25]
    y = list(average.values())[0:25]
    colour = lambda xs: ['#00719e' if x in names else '#f2493d' for x in xs]
    legend_elements = [Line2D([0], [0], marker='o', color='w', label='Known connections', markerfacecolor='#00719e', markersize=10),
                          Line2D([0], [0], marker='o', color='w', label='Novel proteins', markerfacecolor='#f2493d', markersize=10)]
    plt.figure(figsize=(16,10))
    plt.legend(handles=legend_elements, loc='upper right')
    
    plt.bar(x,y,color=colour(x))
    plt.xlabel("Protein",fontsize=14)
    plt.ylabel("Median importance",fontsize=14)
    plt.title(f"Top 25 proteins by median importance to {protein_of_interest} in range 824-864",fontsize=18)
    plt.savefig(f'results/graphs/{protein_of_interest}_median_importance.png',dpi=300)
    
    return True

names = ['LPD1', 'PDA1', 'PYC2', 'PDB1', 'PTC1', 'BAT2', 'KGD1', 'AIM22', 'PKP1', 'PTC5', 'LAT1'] # https://docs.google.com/document/d/12kaAjgjEsQtCOaRqw6g2ZNeLzN-rlzmLaGApKCdI1uc/edit 

for protein in names:
    average_importance_plot(protein,names)
