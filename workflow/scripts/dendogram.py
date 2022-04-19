import matplotlib as mpl
mpl.use("agg")

import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
import scipy.spatial as sp, scipy.cluster.hierarchy as hc
import os


dmr_filename = snakemake.input.dmr
# dmr_filename = "/vol/nano/smith/rrbs/dna-seq-methylation/results/dmr/metilene/dmrs/tumor_vs_healthy_tissue.bed"
dmr_df = pd.read_csv(dmr_filename, sep="\t", usecols=(0,1,2), names=("chrom", "start", "stop"))

m_filename = snakemake.input.mvalues
# m_filename = "/vol/nano/smith/rrbs/dna-seq-methylation/results/dmr/metilene/input/tumor_vs_healthy_tissue.bg"
m_df = pd.read_csv(m_filename, sep="\t", na_values=".")
outdir = snakemake.output.outdir
# outdir = "/vol/nano/smith/rrbs/dna-seq-methylation/results/dmr/metilene/plots/tumor_vs_healthy_tissue"

os.makedirs(outdir, exist_ok=True)

for i, row in dmr_df.iterrows():
    print(i)
    chrom, start, stop = row[["chrom", "start", "stop"]]
    m_values = m_df[(m_df["chr"] == str(chrom)) & (m_df["pos"] >= start) & (m_df["pos"] <= stop)]
    m_values = m_values.set_index(["chr", "pos"])
    # print(chrom, start, stop)
    
    m_values[m_values.isnull()] = 0.5
    # print(m_values.isnull())
    # exit()
    # m_values.drop(["chr", "pos"], axis=1, inplace=True)
    m_values = m_values.transpose()
    # row_dism = 1 - m_values.T.corr()
    # row_linkage = hc.linkage(sp.distance.squareform(row_dism), method='complete')
    # col_dism = 1 - m_values.corr()
    # col_linkage = hc.linkage(sp.distance.squareform(col_dism), method='complete')

    # m_values.drop("pos", axis=1)
    g = sns.clustermap(m_values, mask=m_values.isnull(), cmap="RdBu_r", col_cluster=False, row_cluster=False, vmin=0, vmax=1, facecolor = 'black')
    plt.savefig(os.path.join(outdir, f"{chrom}_{start}_{stop}.pdf"))
    plt.clf()
    # exit()

    #na_values