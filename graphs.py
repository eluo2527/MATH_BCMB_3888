from tabnanny import check
from turtle import color
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import OrderedDict
import func
import random

df = pd.read_csv("results/proteins_by_threshold600-800_detailed.csv")

threshold=df['threshold']

df=df.drop('Unnamed: 0', axis=1)
df=df.drop('threshold', axis=1)
df=df.div(df.sum(axis=1), axis=0)
df.insert(loc=0, column='threshold', value=threshold)

ys = list(df.loc[0].index)
ys.remove("threshold")
# ys = func.remove_essentials_from_list(ys,"network_info\essential_proteins.csv")


colour_map = {}
for y in ys:
    colour_map[y] = "#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])


ax = df.plot(x="threshold", y=ys, kind="line",color=[colour_map.get(the_y, '#333333') for the_y in ys])
ax.get_legend().remove()
# handles, labels = plt.gca().get_legend_handles_labels()
# by_label = OrderedDict(zip(labels, handles))


plt.xlabel("Threshold")
plt.ylabel("Normalized cluster_score*protein_score")
plt.title("Protein importance by threshold")

# make labels




for y in ys:

    checklist = list((df['threshold'][index], val) for index, val in enumerate(df[y]) if val>=0)
    if len(checklist)>1:
        # random choice for xpos, ypos
        choice = random.choice(range(len(checklist)))
        xpos,ypos = checklist[choice]
        random_offset = 1+(random.random()-0.5)/50
        xpos, ypos = random_offset*xpos, random_offset*ypos
        plt.annotate(y, (xpos, ypos),color=colour_map[y])

    # print(index,x_final)
# plt.legend(by_label.values(), by_label.keys(),loc='upper center', bbox_to_anchor=(0.5, 0.95),
#           ncol=15, fancybox=True, shadow=True)
plt.show()



