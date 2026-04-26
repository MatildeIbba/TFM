# Valencia, 19-04-2025
# This script was written by Matilde Ibba
# This script can be used to study the pangenome 

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt


# --------------------------------------------------
# 1. Detect genome columns
# --------------------------------------------------
def detect_genome_columns(df):
 
    metadata_cols = [
        "Gene",
        "Non-unique Gene name",
        "Annotation",
        "No. isolates",
        "No. sequences",
        "Avg sequences per isolate",
        "Genome Fragment",
        "Order within Fragment",
        "Accessory Fragment",
        "Accessory Order with Fragment",
        "QC",
        "Min group size nuc",
        "Max group size nuc",
        "Avg group size nuc"
    ]

    genome_cols = [col for col in df.columns if col not in metadata_cols]
    genome_cols = [col for col in genome_cols if not str(col).startswith("Unnamed")]

    return genome_cols

# --------------------------------------------------
# 2. Classify genes
# --------------------------------------------------
def classify_gene(freq, n_genomes):
    proportion = freq / n_genomes
    if proportion >= 0.99:
        return "core"
    elif proportion >= 0.95:
        return "soft_core"
    elif proportion >= 0.15:
        return "shell"
    else:
        return "cloud"

# --------------------------------------------------
# 3. Main
# --------------------------------------------------
def main():

    # ----------------------------------------------
    # Input file
    # ----------------------------------------------
    if len(sys.argv) < 2:
        print("Uso:")
        print("python panaroo_pangenome_analysis.py gene_presence_absence.csv")
        sys.exit(1)

    input_file = sys.argv[1]

    if not os.path.exists(input_file):
        print(f"Error: no se encuentra el archivo '{input_file}'")
        sys.exit(1)

    # ----------------------------------------------
    # Read file
    # ----------------------------------------------
    print(f"Cargando archivo: {input_file}")
    df = pd.read_csv(input_file)

    # ----------------------------------------------
    # Detect genome columns
    # ----------------------------------------------
    genome_cols = detect_genome_columns(df)
    n_genomes = len(genome_cols)

    if n_genomes == 0:
        print("Error: no se detectaron columnas de genomas.")
        sys.exit(1)

    print(f"Número de genomas detectados: {n_genomes}")

    # ----------------------------------------------
    # Build binary presence/absence matrix
    # ----------------------------------------------
    print("Construyendo matriz binaria de presencia/ausencia...")
    pa_matrix = df[genome_cols].notna().astype(int)

    # ----------------------------------------------
    # Calculate gene frequency
    # ----------------------------------------------
    print("Calculando frecuencia génica...")
    df["gene_frequency"] = pa_matrix.sum(axis=1)
    df["gene_fraction"] = df["gene_frequency"] / n_genomes

    # ----------------------------------------------
    # Classify pangenome categories
    # ----------------------------------------------
    print("Clasificando genes del pangenoma...")
    df["pangenome_category"] = df["gene_frequency"].apply(
        lambda x: classify_gene(x, n_genomes)
    )

    # ----------------------------------------------
    # Summary tables
    # ----------------------------------------------
    total_genes = len(df)
    core_genes = (df["pangenome_category"] == "core").sum()
    soft_core_genes = (df["pangenome_category"] == "soft_core").sum()
    shell_genes = (df["pangenome_category"] == "shell").sum()
    cloud_genes = (df["pangenome_category"] == "cloud").sum()

    summary_stats = pd.DataFrame({
        "metric": [
            "number_of_genomes",
            "total_genes_pangenome",
            "core_genes",
            "soft_core_genes",
            "shell_genes",
            "cloud_genes"
        ],
        "value": [
            n_genomes,
            total_genes,
            core_genes,
            soft_core_genes,
            shell_genes,
            cloud_genes
        ]
    })

    category_summary = (
        df["pangenome_category"]
        .value_counts()
        .rename_axis("category")
        .reset_index(name="n_genes")
    )

    # ----------------------------------------------
    # Save outputs
    # ----------------------------------------------
    print("Guardando resultados...")

    pa_matrix.to_csv("presence_absence_matrix.csv", index=False)
    df.to_csv("gene_presence_absence_annotated.csv", index=False)
    summary_stats.to_csv("pangenome_stats.csv", index=False)
    category_summary.to_csv("pangenome_category_summary.csv", index=False)

    df[df["pangenome_category"] == "core"].to_csv("core_genes.csv", index=False)
    df[df["pangenome_category"] == "soft_core"].to_csv("soft_core_genes.csv", index=False)
    df[df["pangenome_category"] == "shell"].to_csv("shell_genes.csv", index=False)
    df[df["pangenome_category"] == "cloud"].to_csv("cloud_genes.csv", index=False)

    # ----------------------------------------------
    # Plot 1: Histogram of gene frequencies
    # ----------------------------------------------
    print("Generando histograma de frecuencias...")
    plt.figure(figsize=(10, 6))
    plt.hist(df["gene_frequency"], bins=30, edgecolor="black")
    plt.xlabel("Número de genomas en los que aparece cada gen")
    plt.ylabel("Número de genes")
    plt.title("Distribución de la frecuencia génica")
    plt.tight_layout()
    plt.savefig("gene_frequency_histogram.png", dpi=300)
    plt.close()

    # ----------------------------------------------
    # Plot 2: Barplot of pangenome categories
    # ----------------------------------------------
    print("Generando gráfico de categorías del pangenoma...")
    ordered_categories = ["core", "soft_core", "shell", "cloud"]
    ordered_counts = [
        (df["pangenome_category"] == cat).sum()
        for cat in ordered_categories
    ]

    plt.figure(figsize=(8, 6))
    plt.bar(ordered_categories, ordered_counts, edgecolor="black")
    plt.xlabel("Categoría")
    plt.ylabel("Número de genes")
    plt.title("Composición del pangenoma")
    plt.tight_layout()
    plt.savefig("pangenome_category_barplot.png", dpi=300)
    plt.close()

    # ----------------------------------------------
    # Final summary
    # ----------------------------------------------
    print("\nAnálisis completado.\n")
    print("Resumen del pangenoma:")
    print(summary_stats.to_string(index=False))

    print("\nArchivos generados:")
    print("- presence_absence_matrix.csv")
    print("- gene_presence_absence_annotated.csv")
    print("- pangenome_stats.csv")
    print("- pangenome_category_summary.csv")
    print("- core_genes.csv")
    print("- soft_core_genes.csv")
    print("- shell_genes.csv")
    print("- cloud_genes.csv")
    print("- gene_frequency_histogram.png")
    print("- pangenome_category_barplot.png")


if __name__ == "__main__":
    main()