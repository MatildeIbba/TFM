# Valencia, 11-03-2026
# This script was written by Matilde Ibba
# This script is for convert a df obtained in bash from a ResFam analysis
# to a matrix, filtering rare or very common genes, useful for generate a
# heatmap in R

import sys
import pandas as pd
import matplotlib.pyplot as plt

#---------------------------------------------------
# 1. Input file from terminal
#---------------------------------------------------

input_file = sys.argv[1]

#---------------------------------------------------
# 2. Load dataframes
#---------------------------------------------------

df = pd.read_csv(
    input_file, 
    sep="\t", 
    header=None,
    usecols=[0,1]   # Fuerza usar solo las dos primera columnas del df
)
df.columns = ["genome", "gene"]

#---------------------------------------------------
# 3. Clean genes names
#---------------------------------------------------

df["gene"] = df["gene"].str.strip()
df["gene"] = df["gene"].str.replace(" ", "_")
df["gene"] = df["gene"].str.replace("'", "", regex=False)

# Remove duplicated genome-gene pairs
df = df.drop_duplicates()

#---------------------------------------------------
# 4. Create presence/absence matrix
#---------------------------------------------------

matrix = pd.crosstab(df["genome"], df["gene"])
matrix[matrix>0] = 1

print("Matrix size:") 
print(matrix.shape)

#---------------------------------------------------
# 5. Filter rare (≤2 genomes) o very common genes
#---------------------------------------------------

min_genomes = 2
max_fraction = 0.9

gene_counts = matrix.sum(axis=0)

matrix_filtered = matrix.loc[:,
    (gene_counts >= min_genomes) &
    (gene_counts <= max_fraction * len(matrix))
]

print("Genomes:", len(matrix))
print("Genes before filtering:", matrix.shape[1])
print("Genes after filtering:", matrix_filtered.shape[1])

#---------------------------------------------------
# 6. Sort genes by frequency
#---------------------------------------------------

matrix_filtered = matrix_filtered.loc[:, 
    matrix_filtered.sum().sort_values(ascending=False).index
]

#---------------------------------------------------
# 7. Save matrices in cvs format
#---------------------------------------------------

matrix.to_csv("resfmas_matrix_all.tsv", sep="\t")
matrix_filtered.to_csv("resfams_matrix_filtered.tsv", sep="\t")

#---------------------------------------------------
# 8. Gene frequency analysis
#---------------------------------------------------

gene_counts = matrix_filtered.sum(axis=0)
gene_freq = gene_counts / len(matrix_filtered)

summary = pd.DataFrame({
    "genomes_with_gene": gene_counts,
    "frequency": gene_freq
})

summary["frequency"] = summary["frequency"].round(3)

summary = summary.sort_values("genomes_with_gene", ascending=False)

print(summary)

# Convert index to column (avoid "gene" header problem)
summary = summary.reset_index()
summary = summary.rename(columns={"index": "gene"})   # Rename column explicity

# Save frequency table
summary.to_csv("resfams_gene_frequency.tsv", sep="\t", index=False)

#---------------------------------------------------
# 9. Plot gene frequency
#---------------------------------------------------

summary["percent"] = summary["frequency"] * 100

summary.plot(
    x="gene",
    y="genomes_with_gene",
    kind="bar",
    legend=False
)

plt.ylabel("Prevalence (%)")
plt.xlabel("Resistance gene")
plt.title("Distribution of resistance genes")

plt.tight_layout()
plt.savefig("resfams_gene_frequency.png", dpi=300)

print("Analysis complete.")