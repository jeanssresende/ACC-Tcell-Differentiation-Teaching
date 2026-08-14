#==============================================================================
# Projeto:
# ACC-Tcell-Differentiation
#
# Aula 02
# Exploração inicial da coorte TCGA-ACC e controle de qualidade
#
# Objetivos:
# 1. Carregar os dados processados da Aula 01
# 2. Explorar dimensões da matriz de expressão
# 3. Avaliar metadados clínicos
# 4. Investigar distribuição das contagens
# 5. Avaliar tamanho das bibliotecas
# 6. Filtrar genes de baixa expressão
# 7. Realizar PCA exploratória
# 8. Salvar resultados e figuras
#
# Autor: Jean Resende
#==============================================================================

#==============================================================================
# 1. Carregar pacotes
#==============================================================================

library(tidyverse)
library(ggplot2)

#==============================================================================
# 2. Criar pastas de saída
#==============================================================================

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

#==============================================================================
# 3. Carregar objetos da Aula 01
#==============================================================================

counts <- readRDS("data/TCGA_ACC_Counts_Tumor.rds")
sample_info <- readRDS("data/TCGA_ACC_SampleInfo_Tumor.rds")
gene_info <- readRDS("data/TCGA_ACC_GeneInfo.rds")

#==============================================================================
# 4. Informações gerais da coorte
#==============================================================================

cat("==============================\\n")
cat("Informações gerais da coorte\\n")
cat("==============================\\n")

cat("Número de genes:", nrow(counts))
cat("Número de amostras:", ncol(counts))

# Pacientes únicos
patients <- substr(colnames(counts), 1, 12)

cat("Número de pacientes únicos:",
    length(unique(patients)))

# Amostras duplicadas
cat("Amostras duplicadas:",
    sum(duplicated(colnames(counts))))

#==============================================================================
# 5. Exploração dos metadados clínicos
#==============================================================================

cat("==============================\\n")
cat("Variáveis clínicas disponíveis\\n")
cat("==============================\\n")

print(colnames(sample_info))

# Exemplo de variáveis frequentes
if ("race" %in% colnames(sample_info)) {
  cat("Distribuição por raca:")
  print(table(sample_info$race))
}

if ("vital_status" %in% colnames(sample_info)) {
  cat("Status vital:")
  print(table(sample_info$vital_status))
}

#==============================================================================
# 6. Distribuição das contagens brutas
#==============================================================================

pdf("figures/01_boxplot_raw_counts.pdf", width = 12, height = 6)

boxplot(log2(counts + 1),
        outline = FALSE,
        las = 2,
        main = "Distribuição das contagens brutas",
        ylab = "log2(counts + 1)",
        col = "lightblue")

dev.off()

#==============================================================================
# 7. Tamanho das bibliotecas
#==============================================================================

library_size <- colSums(counts)

library_df <- tibble(
  sample = names(library_size),
  reads = library_size
)

p_library <- ggplot(library_df,
                    aes(x = reorder(sample, reads),
                        y = reads)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Tamanho das bibliotecas",
       x = "Amostra",
       y = "Total de reads") +
  theme_minimal()

ggsave("figures/02_library_size.pdf",
       p_library,
       width = 8,
       height = 10)

# Estatísticas descritivas
summary(library_size)

#==============================================================================
# 8. Filtragem de genes de baixa expressão
#==============================================================================

# Critério exploratório:
# gene expresso (>=10 reads) em pelo menos 5 amostras

keep <- rowSums(counts >= 10) >= 5

counts_filtered <- counts[keep, ]

cat("\\n==============================\\n")
cat("Filtragem de genes\\n")
cat("==============================\\n")

cat("Genes antes:", nrow(counts), "\\n")
cat("Genes após filtragem:", nrow(counts_filtered), "\\n")
cat("Proporção mantida:",
    round(nrow(counts_filtered) / nrow(counts) * 100, 2),
    "%\\n")

#==============================================================================
# 9. PCA exploratória
#==============================================================================

log_counts <- log2(counts_filtered + 1)

pca <- prcomp(t(log_counts), scale. = TRUE)

# Variância explicada
variance <- pca$sdev^2 / sum(pca$sdev^2)

pca_df <- data.frame(
  sample = colnames(counts_filtered),
  PC1 = pca$x[,1],
  PC2 = pca$x[,2]
)

p_pca <- ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(size = 3, color = "darkred") +
  labs(title = "PCA exploratória — TCGA-ACC",
       x = paste0("PC1 (",
                  round(variance[1] * 100, 1),
                  "%)"),
       y = paste0("PC2 (",
                  round(variance[2] * 100, 1),
                  "%)")) +
  theme_minimal()

ggsave("figures/03_pca_exploratoria.pdf",
       p_pca,
       width = 7,
       height = 6)

#==============================================================================
# 10. Identificação preliminar de possíveis outliers
#==============================================================================

distance_from_center <- sqrt(pca_df$PC1^2 + pca_df$PC2^2)

threshold <- quantile(distance_from_center, 0.95)

outliers <- pca_df %>%
  mutate(distance = distance_from_center) %>%
  filter(distance > threshold)

cat("\\n==============================\\n")
cat("Possíveis outliers (PCA)\\n")
cat("==============================\\n")

print(outliers)

#==============================================================================
# 11. Salvar resultados
#==============================================================================

write.csv(pca_df,
          "results/pca_coordinates.csv",
          row.names = FALSE)

write.csv(outliers,
          "results/pca_outliers.csv",
          row.names = FALSE)

saveRDS(counts_filtered,
        "data/TCGA_ACC_Counts_Filtered.rds")

#==============================================================================
# Fim da Aula 02
#==============================================================================