from collections import defaultdict
import matplotlib as mlp
mlp.use("Agg")

import pandas as pd
import numpy as np
import argparse
import glob
import matplotlib.pyplot as plt
from matplotlib import cm

from sklearn.decomposition import PCA
from functools import reduce


parser = argparse.ArgumentParser()
parser.add_argument("--group", nargs="+", action='append')
parser.add_argument("--symbols", nargs="+")
parser.add_argument("-o")
args = parser.parse_args()
groups = args.group
symbols = args.symbols

# filenames = glob.glob("results/methylation/A006200150_1473*")
# groups = [filenames]

sample_to_group_idx = defaultdict()
sample_names = []
data = []
print("loading data")
sample_id = 0
for g_idx, filenames in enumerate(groups):
    for filename in filenames:
        sample_names.append(filename.rsplit("/", 1)[-1].rsplit(".",1)[0].split("_")[-2]) 
        sample_to_group_idx[sample_id] = g_idx
        sample_id += 1
        print(filename)
        d = pd.read_csv(
            filename,
            skiprows=1,
            sep="\t",
            names=["chrom", "start", "end", "methylation", "n", "m"],
            dtype={
                "chrom": "<S10",
                "start": np.int32,
                "end": np.int32,
                "methylation": np.int8,
                "n": np.int16,
                "m": np.int32,
            }
        )
        d.set_index(["chrom", "start"], inplace=True)
        valid_cov = (d["m"] + d["n"]) >= 8
        data.append(d[valid_cov])

print("unifying positions")
index = reduce(lambda a,b: a.intersection(b), map(lambda a: a.index, data))
data = map(lambda d: d.loc[index], data)
x = np.vstack([d["methylation"].values for d in data])

print("performing pca")
pca = PCA(n_components=2)
t = pca.fit_transform(x)

cmap = cm.get_cmap("hsv", len(groups))

ax = plt.subplot(1,1,1)

for i, (single, symbol, name) in enumerate(zip(t, symbols, sample_names)):
    s = ax.scatter(single[0], single[1], c=[cmap(sample_to_group_idx[i])], marker=symbol, label=name)
    # ax.legend([s], [f'label{i}'])

box = ax.get_position()
ax.set_position([box.x0 - box.width * 0.05, box.y0,
                 box.width * 0.9, box.height])

# Put a legend below current axis
ax.legend(loc='upper center', bbox_to_anchor=(1.2, 1.1),
          fancybox=True, shadow=True)

plt.axis('equal')
plt.savefig(args.o)

    # X = pca.transform(X)


    # print(x)
    # print("pca")
    # pca.fit(x)
    # print(pca.explained_variance_ratio_)
    # print(pca.components_)