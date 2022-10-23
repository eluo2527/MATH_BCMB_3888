import func 
import pandas as pd
import matplotlib.pyplot as plt

network_name = "network_info/4932_protein_links_v11_5.txt"
df = pd.read_csv("results/proteins_by_threshold800-900_detailed.csv")

mean_importance = df.mean().drop('Unnamed: 0').drop('threshold').sort_values(ascending=False)

print(sum(mean_importance))