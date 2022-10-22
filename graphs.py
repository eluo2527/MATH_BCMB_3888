import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("results\proteins_by_threshold.csv")

for index, row in df.iterrows():
    print(list(row.index))
    plt.plot(list(row.index),row.to_numpy())

plt.show()