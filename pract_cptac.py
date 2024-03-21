import cptac
import cptac.utils as ut
import csv
import numpy as np
cptac.download("luad")
lu = cptac.Luad()

lu_clin_tran = lu.join_metadata_to_omics(metadata_df_name='clinical',
  omics_df_name='transcriptomics',
  metadata_cols = ["Sample_Tumor_Normal"])
#Some functions require that we drop the Database_ID column.
lu_clin_tran.to_csv('transcriptomics_data.csv')
#luad_no_id = ut.reduce_multiindex(df=lu_clin_`prot, levels_to_drop="Database_ID")

#lu_transcriptomics = lu.get_transcriptomics()

#lu_filtered = luad_no_id.loc[luad_no_id["Sample_Tumor_Normal"] == "Normal"]
lu_filtered = lu_clin_tran.loc[lu_clin_tran["Sample_Tumor_Normal"] == "Normal"]
print(lu_filtered)
gene_dict = {}
for gene in lu_filtered:
  if gene == "Sample_Tumor_Normal":
    continue
  levels = lu_filtered[gene]
  levels = levels.dropna()
  levels = levels.sort_values()
  median = np.median(levels)
  
  gene = gene[:gene.find("_")]
  gene_dict[gene] = median

sort_dict = dict(sorted(gene_dict.items(), key=lambda item: item[1], reverse=True))
genes = list(sort_dict.keys())
print(len(genes))
n = len(genes)//4
high_exp = genes[:n+1]
print(len(high_exp))
with open("high_exp_genes_og.txt", 'w') as out:
  for gene in high_exp:
    out.write(gene + "\n")

# for name in list(lu_filtered.index):
#   row = lu_filtered.loc[name]
#   row = row.drop(labels='Sample_Tumor_Normal')
#   row = row.sort_values(ascending=False)
#   row = row.dropna()
#   n = len(row)//3
#   row = row[:n]
#   row.to_csv('row.csv')

#   with open(f'{name}.tsv', 'w') as out:
#     writer = csv.writer(out, delimiter = '\t')
#     writer.writerow(['gene','assoc_DB_IDs','protein_expression'])
#     for info in list(row.index):
#       exp = row.loc[info]
#       gene = info[0]
#       gene = gene[:gene.find("_")]
#       proteinID = info[1]
#       writer.writerow([gene, proteinID, exp])
  
