import func 
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
import random
from collections import OrderedDict, defaultdict
import numpy as np
from scipy.stats import linregress
network_name = "network_info/4932_protein_links_v11_5.txt"
df = pd.read_csv("results/proteins_by_threshold_detailed.csv")

data = []
colour_map = {}
lin_reg_data = defaultdict(list)

for index, row in df.iterrows():
    if index%10 == 0:
        print(f"{round(index/len(df)*100)}%")
        row = row.drop('Unnamed: 0')
        threshold = row['threshold']
        row = row.drop('threshold')
        row = row.sort_values(ascending=False)
        importance_values = {}
        for index, value in row.iteritems():
            if value > 0:
                importance_values[index] = value
        colour_map[threshold] = "#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])

        G = func.remove_threshold(network_name, threshold)


        for key in importance_values.keys():
            data.append([G.degree(func.name_change(key)), importance_values[key],threshold])
            lin_reg_data[threshold].append((G.degree(func.name_change(key)), importance_values[key]))
for datum in data:
    plt.plot(datum[0], datum[1], 'o', color=colour_map[datum[2]],label=f"Threshold: {datum[2]}")

# plot linear regressions for each threshold
# print(lin_reg_data)
for threshold in lin_reg_data.keys():
    x = np.array([datum[0] for datum in lin_reg_data[threshold]])
    y = np.array([datum[1] for datum in lin_reg_data[threshold]])
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    plt.plot(x, intercept + slope*x, colour_map[threshold ], label=f"r_sqd = {round(r_value**2,3)}",linestyle=(0,(5,10)), linewidth=0.5)

x = np.array([datum[0] for datum in data])
y =  np.array([datum[1] for datum in data])
a, b, r_value, p_value, std_err = linregress(x, y)

plt.plot(x, a*x + b, color="grey", label=f"Cumulated: r_sqd = {round(r_value**2,3)}",linestyle=(0,(3,10,1,10)), linewidth=1.5)
# plt.text(0.9, 0.9, f"R^2 = {round(r_value**2, 3)}", fontsize=12, transform=plt.gca().transAxes)


handles, labels = plt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys())
plt.title("Degree vs Importance")
plt.xlabel("Degree")
plt.ylabel("Importance")
plt.show()
