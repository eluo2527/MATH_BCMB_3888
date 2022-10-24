import pandas as pd
import numpy as np
import csv
def average_importance(protein_of_interest):
    df = pd.read_csv(f'results/data/old/{protein_of_interest}_proteins_by_threshold.csv')
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

names = ['LPD1', 'PDA1', 'PYC2', 'PDB1', 'PTC1', 'BAT2', 'KGD1', 'AIM22', 'PKP1', 'PTC5', 'LAT1'] # https://docs.google.com/document/d/12kaAjgjEsQtCOaRqw6g2ZNeLzN-rlzmLaGApKCdI1uc/edit 

for protein in names:
    average_importance(protein)
