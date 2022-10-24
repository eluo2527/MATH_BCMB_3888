import pandas as pd
import numpy as np
df = pd.read_csv('results\data\PDA1_proteins_by_threshold.csv')
average = {}
for column in df:
    values = np.array(df[column])
    values[np.isnan(values)] = 0
    if column != 'threshold':
        average[column] = np.median(values)
for key,value in average.items():
    average[key] = value/sum(average.values())

average = dict(sorted(average.items(), reverse=True, key= lambda x: x[1]))
print(average)