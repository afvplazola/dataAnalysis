import os
import csv
import pandas as pd

genes = []
# get the genes for the headers of the output file
with open("allUniqueGenes.txt") as geneFile:
    for line in geneFile:
        gene = line.strip()
        genes.append(gene)

genes.sort()

# get the patient paths 
patientsHaplo = os.listdir("/home/adelynfs/fsl_groups/fslg_Hawkins/compute/Projects/Ramps/ramp-expression-correlation/haplo_results")
patientsHaplo.sort()

# create the data matrix with default values of 0/0
data = [["0/0"] * len(genes) for _ in range(len(patientsHaplo)//2)]


# get the patient names
patients = []
i = 0
for pat in patientsHaplo:
    if pat[:-2] not in patients:
        patients.append(pat[:-2])
        data[i].insert(0,pat[:-2])
        i+=1

# for each patient
for i in range(len(patients)):
    haplo1 = []
    haplo2 = []
    string1 = patients[i] + ".1"
    string2 = patients[i] + ".2"
    # open both the haplo ramp files and get genes in each file
    with open(os.path.join("/home/adelynfs/fsl_groups/fslg_Hawkins/compute/Projects/Ramps/ramp-expression-correlation/haplo_results", string1, "ramps.fa")) as h1, open(os.path.join("/home/adelynfs/fsl_groups/fslg_Hawkins/compute/Projects/Ramps/ramp-expression-correlation/haplo_results", string2, "ramps.fa")) as h2:
        for line in h1:
            if line.startswith(">"):
                gene = line.strip()[line.find("gene=")+5:line.find("product=")-1]
                haplo1.append(gene)
        for line in h2:
            if line.startswith(">"):
                gene = line.strip()[line.find("gene=")+5:line.find("product=")-1]
                haplo2.append(gene)
    # for each gene, get the ramp status
    for j in range(len(genes)):
        if genes[j] in haplo1 and genes[j] in haplo2:
            data[i][j+1] = "1/1"
        elif genes[j] in haplo1 or genes[j] in haplo2:
            data[i][j+1] = "0/1"

# get the paths to the full CDS files
patientsFull = os.listdir("/home/adelynfs/fsl_groups/fslg_Hawkins/compute/Projects/Ramps/data/LUAD/normal/FASTAs")
patientsFull.sort()

# for each patient
index = 0
for i in range(len(patientsFull)):
    # ignore the N.2 files
    if "N.2" in patientsFull[i]:
        index += 1
        continue
    geneList1 = []
    geneList2 = []

    #EDIT Path
    with open(os.path.join("/home/adelynfs/fsl_groups/fslg_Hawkins/compute/Projects/Ramps/data/LUAD/normal/FASTAs", patientsFull[i])) as inF:
        for line in inF:
            if line.startswith(">"):
                gene = line.strip()[line.find("gene=")+5:line.find("product=")-1]
                if gene not in geneList1:
                    geneList1.append(gene)
    
    secondString = patientsFull[i][:-7] + "2.fasta"

    with open(os.path.join("/home/adelynfs/fsl_groups/fslg_Hawkins/compute/Projects/Ramps/data/LUAD/normal/FASTAs", secondString)) as inF:
        for line in inF:
            if line.startswith(">"):
                gene = line.strip()[line.find("gene=")+5:line.find("product=")-1]
                if gene not in geneList2:
                    geneList2.append(gene)
    
    # modify ramp status if genes are not in the respective full CDS haplo files
    for j in range(len(genes)):
        if genes[j] not in geneList1 and genes[j] not in geneList2:
            data[index][j+1] = "NA/NA"
        elif genes[j] not in geneList1 or genes[j] not in geneList2:
            if data[index][j+1] == "1/1" or data[index][j+1] == "0/1":
                data[index][j+1] = "NA/1"
            else:
                data[index][j+1] = "NA/0"


# add patient column
genes.insert(0,"Patient")

df = pd.DataFrame(data, columns = genes)

# open and write to the output file
with open("geneRampAnalysis.tsv", 'w') as out:
    df.to_csv(out, sep="\t", index = False)
