import pandas as pd
from statsmodels.stats.multitest import multipletests

df = pd.read_csv("stats_analysis-v2.tsv", sep = "\t")
df.dropna(subset=['p-value'], inplace=True)

p_values = df['p-value'].tolist()

reject, corrected_p_values, _, _ = multipletests(p_values, method='fdr_bh')

df['adj_p-value'] = corrected_p_values

df.to_csv('adj_stats_analysis.tsv', sep='\t', index=False)