# Valencia, 10-03-2025
# This script was written by Matilde Ibba
# This script can be used to graphic a heatmap of antibiotic resistance 
# and a histogram of gene prevalence distribution

#--------------------------------------------------
# 1. Set working directory
#--------------------------------------------------
setwd("~/Desktop/TFM")

#--------------------------------------------------
# 2. Load resistome matrix
#--------------------------------------------------
matrix <- read.table(
  "results/Results/resfams_matrix_all.tsv",
  header = TRUE,
  row.names = 1,
  sep = "\t"
)

dim(matrix)
head(matrix)

#--------------------------------------------------
# 3. Calculate gene frequency
#--------------------------------------------------
gene_freq <- colSums(matrix)
gene_fraction <- gene_freq / nrow(matrix)

#--------------------------------------------------
# 4. Filter informative genes
# Remove genes present in > 95% genomes and genes present in < 2 genomes
#--------------------------------------------------
matrix_filtered <- matrix[, gene_fraction < 0.95 & gene_freq > 0.05]
dim(matrix_filtered)

#--------------------------------------------------
# 5. Sort genes by frequency
#--------------------------------------------------
gene_freq_filtered <- colSums(matrix_filtered)
matrix_filtered <- matrix_filtered[, order(gene_freq_filtered,
                                           decreasing = TRUE)]

#--------------------------------------------------
# 6. Create metadata table
#--------------------------------------------------
genomes <- row.names(matrix_filtered)

metadata <- data.frame(
  genome = genomes,
  source = ifelse(grepl("^id-", genomes), "PubMLST", "Paper")
)

# Convert to rownames format for pheatmap
rownames(metadata) <- metadata$genome
metadata$genome <- NULL
head(metadata)

#--------------------------------------------------
# 7. Sort genomes by source
#--------------------------------------------------
sorted_genomes <- order(metadata$source)

matrix_filtered <- matrix_filtered[sorted_genomes, ]
metadata <- metadata[sorted_genomes, , drop = FALSE]

#--------------------------------------------------
# 8. Find where genome groups change
#--------------------------------------------------
group_changes <- which(diff(as.numeric(as.factor(metadata$source))) !=0)

#--------------------------------------------------
# 9. Gene annotation (mechanism)
# Only genes that appear in the filtered matrix
#--------------------------------------------------
gene_mechanism <- c(
    L1 = "beta-lactamase",
    'AAC6.II' = "aminoglycoside modification",
    SPM = "metallo-beta-lactamase",
    TetD = "tetracycline efflux",
    ArmA_Rmt = "rRNA methylation",
    'CMY.LAT.MOX.ACT.MIR.FOX' = 'beta-lactamase'
  )

mechanism_vector <- gene_mechanism[colnames(matrix_filtered)]

# Replace possible missing annotation
mechanism_vector[is.na(mechanism_vector)] <- "other"

gene_annotation <- data.frame(
  mechanism = mechanism_vector
)

rownames(gene_annotation) <- colnames(matrix_filtered)

#--------------------------------------------------
# 10 Define annotation colors
#--------------------------------------------------
ann_colors <- list(
  
  source = c(
    Paper = "orange",
    PubMLST = "skyblue"
  ),
  
  mechanism = c(
    "beta-lactamase" = "steelblue",
    "metallo-beta-lactamase" = "navy",
    "aminoglycoside modification" = "darkgreen",
    "rRNA methylation" = "gold",
    "tetracycline efflux" = "purple"
  )
)

#--------------------------------------------------
# 11. Create heatmap
#--------------------------------------------------
library(pheatmap)

pheatmap(matrix_filtered,
         color = colorRampPalette(c("white", "#fcbba1", "#cb181d"))(50),
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         fontsize_col = 12,
         fontsize_row = 6,
         angle_col = 45,
         cellwidth = 25,
         cellheight = 4,
         annotation_row = metadata,
         annotation_col = gene_annotation,
         annotation_colors = ann_colors,
         gaps_row = group_changes,
         border_color = NA
         )

# Saving the heatmap as an image (pdf)
pdf("results/Results/resfams_heatmap.pdf", width = 6, height = 8)

#--------------------------------------------------
# 12. Create gene frequency plot
#--------------------------------------------------
library(ggplot2)

# Creating dataframe
gene_df <- data.frame(
  gene = names(gene_freq),
  genomes = gene_freq,
  frequency = gene_fraction
)

# Sort by frequency
gene_df <- gene_df[order(gene_df$genomes, decreasing = TRUE), ]
head(gene_df)

# Creating histogram (core genes)
ggplot(data.frame(freq = gene_fraction), aes(x = freq)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "black") +
  geom_vline(xintercept = c(0.05, 0.95),
             linetype = "dashed", color = "red") +
  labs(
    title = "Gene Prevalence Distribution",
    x = "Fraction of Genomes",
    y = "Number of Genes"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Saving the histogram as an image (pdf)
pdf("results/Results/histogram_core_genes.pdf", width = 6, height = 8)

# Data frame for the barplot
gene_df_plot <- gene_df[
  gene_df$frequency > 0.05 & gene_df$frequency < 0.95,
]

# Creating the barpolot
ggplot(gene_df_plot, aes(x = reorder(gene, genomes), y = genomes)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top resistance genes by prevalence",
    x = "Gene",
    y = "Number of genomes"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Saving the barplot as an image (pdf)
pdf("results/Results/barplot_top_resistance.pdf", width = 6, height = 8)