import pandas as pd
import cptac
import cptac.utils as ut

cptac.download("luad")
lu = cptac.Luad()

lu_clin_tran = lu.join_metadata_to_omics(metadata_df_name='clinical',
  omics_df_name='transcriptomics',
  metadata_cols=["Sample_Tumor_Normal"])

lu_filtered = lu_clin_tran.loc[lu_clin_tran["Sample_Tumor_Normal"] == "Normal"]
df = pd.read_csv("geneRampAnalysis.tsv", sep='\t', header = [0], index_col = [0])

patientIDs = lu_filtered.index.tolist()
patients = df.index.tolist()

for pat in patientIDs:
  if pat not in patients:
    lu_filtered = lu_filtered.drop(pat)

for pat in patients:
  if pat not in patientIDs:
    df = df.drop(pat)
df = df.add_suffix("_rampStatus")
lu_filtered = lu_filtered.drop("Sample_Tumor_Normal", axis = 1)

result = pd.concat([df, lu_filtered], axis=1)
result = result.sort_index(axis=1)
result.to_csv("final_data_transcriptomics.tsv", sep='\t', header=True)