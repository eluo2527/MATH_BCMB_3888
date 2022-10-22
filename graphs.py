import pandas as pd
import matplotlib.pyplot as plt
import numpy
from collections import OrderedDict
import func

df = pd.read_csv("results/proteins_by_threshold.csv")
ys = list(df.loc[0].index)
ys.remove("threshold")
ys.remove('Unnamed: 0')
ys = func.remove_essentials_from_list(ys,"network_info\essential_proteins.csv")

df.plot(x="threshold", y=ys, kind="line")

handles, labels = plt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))


plt.xlabel("Threshold")
plt.ylabel("Normalized cluster_score*protein_score")
plt.title("Protein importance by threshold")
plt.legend(by_label.values(), by_label.keys(),loc='upper center', bbox_to_anchor=(0.5, 0.95),
          ncol=15, fancybox=True, shadow=True)
plt.show()