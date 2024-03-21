import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pingouin as pg
from scipy.stats import mannwhitneyu

df = pd.read_csv("final_data_transcriptomics.tsv", sep='\t', header = [0], index_col = [0])
col_list = df.columns.tolist()

ramp_dict = {"0/0":[], "0/1":[], "1/1":[]}
i = 0
while i < len(col_list):
    column = col_list[i]
    ramp_stat = df[column]
    i += 1
    column = col_list[i]
    exp_level = df[column]

    #Making sure that the ramp_status column is connected to an expression level column
    a = exp_level[0]
    if a == "0/0" or a == "0/1" or a == "1/1" or a == "NA/0" or a == "NA/1" or a == "NA/NA":
        continue
    
    i += 1
    for j in range(len(ramp_stat)):
        if ramp_stat[j] in ramp_dict.keys():
            if not pd.isna(exp_level[j]):  # Check if the value is not NaN
                temp_list = ramp_dict[ramp_stat[j]]
                temp_list.append(float(exp_level[j]))
                ramp_dict[ramp_stat[j]] = temp_list

mean_group1 = np.mean(ramp_dict["0/0"])
mean_group2 = np.mean(ramp_dict["0/1"])
mean_group3 = np.mean(ramp_dict["1/1"])

print(mean_group1)
print(mean_group2)
print(mean_group3)
# Create box plot
plt.boxplot([ramp_dict["0/0"], ramp_dict["0/1"], ramp_dict["1/1"]], labels=['0/0', '0/1', '1/1'])

plt.title('Comparing Transcriptomics Expression Levels')
plt.ylabel('Expression Level')
plt.xlabel('Ramp Status')
plt.grid(True)
plt.savefig('all_points.png')

# Perform Mann-Whitney U test
statistic, p_value = mannwhitneyu(ramp_dict["0/0"], ramp_dict["1/1"])

# Define significance level
alpha = 0.05

# Print the results
print("Mann-Whitney U Test:")
print("Statistic:", statistic)
print("P-value:", p_value)

# Check if the p-value is less than the significance level
if p_value < alpha:
    print("Reject null hypothesis: There is a significant difference between the two groups.")
else:
    print("Fail to reject null hypothesis: There is no significant difference between the two groups.")

from scipy.stats import ranksums

# Perform Wilcoxon rank-sum test
z_score, p_value = ranksums(ramp_dict["0/0"], ramp_dict["1/1"])

# Calculate rank-biserial correlation coefficient
rank_biserial_corr = z_score / np.sqrt(len(ramp_dict["0/0"]) + len(ramp_dict["1/1"]))

print("Rank-Biserial Correlation Coefficient:", rank_biserial_corr)
