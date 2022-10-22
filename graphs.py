import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("results\proteins_by_threshold.csv")
print(df.T)


# for index, row in df.iterrows():
#     print(row["threshold"])
#     plt.plot(row["threshold"], row["proteins"], 'ro')
#     plt.show()

# plt.show()
