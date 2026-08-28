###############################################################################
#
# Projeto:
# ACC-Tcell-Differentiation
#
# Aula 03
# Exploração dos marcadores de diferenciação de células T
# na matriz de expressão do TCGA-ACC
#
# Objetivos:
#
# 1. Compreender como marcadores gênicos podem representar estados
#    funcionais de células T
#
# 2. Explorar a expressão de marcadores de células T na coorte TCGA-ACC
#
# 3. Trabalhar com identificadores Ensembl e símbolos gênicos
#
# 4. Utilizar a anotação gênica presente no objeto SummarizedExperiment
#
# 5. Comparar programas associados a:
#       - células T naïve
#       - TSCM
#       - TCM
#       - TEM
#       - TEFF
#       - TPEX
#       - TEX
#       - citotoxicidade
#
# 6. Visualizar a expressão desses genes entre as amostras
#
# 7. Compreender as limitações da análise de marcadores em RNA-seq bulk
#
# IMPORTANTE:
#
# Nesta aula, a expressão dos marcadores será utilizada para uma análise
# EXPLORATÓRIA.
#
# A expressão de um ou poucos genes NÃO permite estimar diretamente a
# proporção de uma população celular.
#
# As estimativas de composição celular e os escores de assinaturas serão
# desenvolvidos em etapas posteriores do projeto.
#
# Autor:
# Jean Resende
#
###############################################################################


#==============================================================================
# 1. Instalação dos pacotes
#==============================================================================

# Executar apenas uma vez

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(
  "edgeR"
)

install.packages(
  c(
    "tidyverse",
    "pheatmap"
  )
)


#==============================================================================
# 2. Carregar os pacotes
#==============================================================================

library(edgeR)
library(tidyverse)
library(pheatmap)
library(SummarizedExperiment)


#==============================================================================
# 3. Criar estrutura de diretórios
#==============================================================================

dir.create(
  "results",
  showWarnings = FALSE
)

dir.create(
  "figures",
  showWarnings = FALSE
)


#==============================================================================
# 4. Carregar os dados
#==============================================================================

# Matriz de contagens das amostras tumorais

counts_tumor <- readRDS(
  "data/TCGA_ACC_Counts_Tumor.rds"
)


# Objeto SummarizedExperiment contendo a anotação dos genes

tcgaACC <- readRDS(
  "data/TCGA_ACC_SummarizedExperiment.rds"
)


#==============================================================================
# 5. Explorar o objeto SummarizedExperiment
#==============================================================================

tcgaACC


# Número de genes

nrow(tcgaACC)


# Número de amostras

ncol(tcgaACC)


#==============================================================================
# 6. Extrair a anotação dos genes
#==============================================================================

# A anotação está armazenada em rowRanges()

rowRanges(tcgaACC)


# Extrair as colunas de anotação

gene_info <- as.data.frame(
  mcols(rowRanges(tcgaACC))
)


# Visualizar as primeiras linhas

head(gene_info)


# Visualizar os nomes das colunas

colnames(gene_info)


#==============================================================================
# 7. Construir tabela de anotação
#==============================================================================

gene_annotation <- gene_info %>%
  select(
    gene_id,
    gene_name,
    gene_type
  ) %>%
  mutate(
    ensembl_id = sub(
      "\\..*$",
      "",
      gene_id
    )
  )


# Visualizar

head(gene_annotation)


#==============================================================================
# 8. Verificar os identificadores Ensembl
#==============================================================================

head(
  gene_annotation$gene_id
)

head(
  gene_annotation$ensembl_id
)


#==============================================================================
# 9. Remover a versão dos Ensembl IDs da matriz
#==============================================================================

# Exemplo:
#
# ENSG00000000003.15
#
# torna-se:
#
# ENSG00000000003

rownames(counts_tumor) <- sub(
  "\\..*$",
  "",
  rownames(counts_tumor)
)


# Verificar

head(
  rownames(counts_tumor)
)


#==============================================================================
# 10. Verificar correspondência entre matriz e anotação
#==============================================================================

genes_matrix <- rownames(
  counts_tumor
)

genes_annotation <- gene_annotation$ensembl_id


sum(
  genes_matrix %in% genes_annotation
)


# Percentual de genes da matriz presentes na anotação

mean(
  genes_matrix %in% genes_annotation
) * 100


#==============================================================================
# 11. Preparar a tabela de anotação
#==============================================================================

# Garantir um identificador Ensembl único

gene_annotation <- gene_annotation %>%
  distinct(
    ensembl_id,
    .keep_all = TRUE
  )


#==============================================================================
# 12. Função para converter símbolo gênico em Ensembl ID
#==============================================================================

symbol_to_ensembl <- function(
    symbols,
    annotation
) {
  
  annotation %>%
    filter(
      gene_name %in% symbols
    ) %>%
    select(
      gene_name,
      ensembl_id
    ) %>%
    distinct()
}


#==============================================================================
# 13. Definir marcadores pelos símbolos gênicos
#==============================================================================

# Os marcadores são definidos por seus símbolos HGNC.
#
# Posteriormente, eles serão convertidos para Ensembl IDs utilizando
# a anotação do próprio objeto TCGA-ACC.


#------------------------------------------------------------------------------
# T cells - geral
#------------------------------------------------------------------------------

T_cell_general_symbols <- c(
  "CD3D",
  "CD3E",
  "CD3G",
  "TRBC1",
  "TRBC2"
)


#------------------------------------------------------------------------------
# T cells naïve
#------------------------------------------------------------------------------

T_naive_symbols <- c(
  "CCR7",
  "SELL",
  "IL7R",
  "LTB",
  "MAL",
  "LEF1",
  "TCF7"
)


#------------------------------------------------------------------------------
# TSCM - T memory stem cells
#------------------------------------------------------------------------------

TSCM_symbols <- c(
  "TCF7",
  "LEF1",
  "CCR7",
  "SELL",
  "IL7R",
  "LTB",
  "BCL2",
  "MAL"
)


#------------------------------------------------------------------------------
# TCM - Central memory
#------------------------------------------------------------------------------

TCM_symbols <- c(
  "CCR7",
  "SELL",
  "IL7R",
  "LTB",
  "TCF7",
  "LEF1"
)


#------------------------------------------------------------------------------
# TEM - Effector memory
#------------------------------------------------------------------------------

TEM_symbols <- c(
  "GZMK",
  "CCL5",
  "CXCR3",
  "IL7R",
  "CD27",
  "CD28"
)


#------------------------------------------------------------------------------
# TEFF - Effector T cells
#------------------------------------------------------------------------------

TEFF_symbols <- c(
  "CCL5",
  "NKG7",
  "GZMA",
  "GZMH",
  "GZMB",
  "PRF1",
  "GNLY",
  "IFNG"
)


#------------------------------------------------------------------------------
# TPEX - Progenitor exhausted
#------------------------------------------------------------------------------

TPEX_symbols <- c(
  "TCF7",
  "PDCD1",
  "CXCR5",
  "SLAMF6",
  "TOX",
  "MAL",
  "IL7R"
)


#------------------------------------------------------------------------------
# TEX - Exhausted T cells
#------------------------------------------------------------------------------

TEX_symbols <- c(
  "PDCD1",
  "HAVCR2",
  "LAG3",
  "TIGIT",
  "CTLA4",
  "TOX",
  "ENTPD1",
  "CXCL13"
)


#------------------------------------------------------------------------------
# Programa citotóxico
#------------------------------------------------------------------------------

Cytotoxicity_symbols <- c(
  "NKG7",
  "CCL5",
  "GZMA",
  "GZMB",
  "GZMH",
  "PRF1",
  "GNLY",
  "GZMK"
)


#==============================================================================
# 14. Criar lista dos marcadores
#==============================================================================

marker_symbols <- list(
  
  T_cell = T_cell_general_symbols,
  
  T_naive = T_naive_symbols,
  
  TSCM = TSCM_symbols,
  
  TCM = TCM_symbols,
  
  TEM = TEM_symbols,
  
  TEFF = TEFF_symbols,
  
  TPEX = TPEX_symbols,
  
  TEX = TEX_symbols,
  
  Cytotoxicity = Cytotoxicity_symbols
  
)


#==============================================================================
# 15. Converter símbolos para Ensembl IDs
#==============================================================================

marker_sets <- lapply(
  marker_symbols,
  symbol_to_ensembl,
  annotation = gene_annotation
)


# Visualizar

marker_sets


#==============================================================================
# 16. Criar tabela completa de marcadores
#==============================================================================

marker_table <- bind_rows(
  lapply(
    names(marker_sets),
    function(set_name) {
      
      marker_sets[[set_name]] %>%
        mutate(
          marker_set = set_name
        )
      
    }
  )
) %>%
  select(
    marker_set,
    gene_name,
    ensembl_id
  )


# Visualizar

marker_table


#==============================================================================
# 17. Verificar quais marcadores estão presentes na matriz
#==============================================================================

marker_table <- marker_table %>%
  mutate(
    present_in_matrix =
      ensembl_id %in% rownames(counts_tumor)
  )


marker_table


#==============================================================================
# 18. Verificar marcadores ausentes
#==============================================================================

marker_table %>%
  filter(
    !present_in_matrix
  )


#==============================================================================
# 19. Resumo da presença dos marcadores
#==============================================================================

marker_summary <- marker_table %>%
  group_by(
    marker_set
  ) %>%
  summarise(
    total_markers = n(),
    markers_present = sum(
      present_in_matrix
    ),
    markers_absent = sum(
      !present_in_matrix
    )
  )


marker_summary


#==============================================================================
# 20. Salvar tabela de marcadores
#==============================================================================

write.csv(
  marker_table,
  "results/TCGA_ACC_Tcell_marker_annotation.csv",
  row.names = FALSE
)


write.csv(
  marker_summary,
  "results/TCGA_ACC_Tcell_marker_summary.csv",
  row.names = FALSE
)


#==============================================================================
# 21. Selecionar somente marcadores presentes
#==============================================================================

genes_present <- marker_table %>%
  filter(
    present_in_matrix
  ) %>%
  pull(
    ensembl_id
  ) %>%
  unique()


#==============================================================================
# 22. Criar matriz de expressão dos marcadores
#==============================================================================

marker_expression_counts <- counts_tumor[
  genes_present,
  ,
  drop = FALSE
]


# Verificar dimensão

dim(
  marker_expression_counts
)


#==============================================================================
# 23. Transformação log2-CPM
#==============================================================================

# A transformação será utilizada apenas para visualização e exploração.

dge <- DGEList(
  counts = counts_tumor
)


logCPM <- cpm(
  dge,
  log = TRUE,
  prior.count = 1
)


#==============================================================================
# 24. Criar matriz log2-CPM dos marcadores
#==============================================================================

marker_expression <- logCPM[
  genes_present,
  ,
  drop = FALSE
]


#==============================================================================
# 25. Criar nomes dos genes para visualização
#==============================================================================

marker_labels <- marker_table %>%
  filter(
    present_in_matrix
  ) %>%
  select(
    ensembl_id,
    gene_name
  ) %>%
  distinct(
    ensembl_id,
    .keep_all = TRUE
  )


# Organizar na mesma ordem da matriz

marker_labels <- marker_labels[
  match(
    rownames(marker_expression),
    marker_labels$ensembl_id
  ),
]


#==============================================================================
# 26. Heatmap dos marcadores
#==============================================================================

pdf(
  "figures/01_Tcell_marker_heatmap.pdf",
  width = 12,
  height = 10
)


pheatmap(
  marker_expression,
  scale = "row",
  show_colnames = FALSE,
  labels_row = marker_labels$gene_name,
  fontsize_row = 9,
  main = "T-cell differentiation markers - TCGA-ACC"
)


dev.off()


#==============================================================================
# 27. Genes selecionados para exploração individual
#==============================================================================

genes_key_symbols <- c(
  "CD3D",
  "CD3E",
  "TCF7",
  "CCR7",
  "IL7R",
  "GZMK",
  "CCL5",
  "NKG7",
  "GZMB",
  "PRF1",
  "PDCD1",
  "TOX",
  "LAG3",
  "TIGIT",
  "HAVCR2"
)


#==============================================================================
# 28. Converter genes selecionados para Ensembl
#==============================================================================

genes_key <- symbol_to_ensembl(
  genes_key_symbols,
  gene_annotation
)


# Manter apenas genes presentes na matriz

genes_key <- genes_key %>%
  filter(
    ensembl_id %in% rownames(logCPM)
  )


genes_key


#==============================================================================
# 29. Criar tabela para os gráficos
#==============================================================================

plot_data <- logCPM[
  genes_key$ensembl_id,
  ,
  drop = FALSE
] %>%
  as.data.frame() %>%
  rownames_to_column(
    "ensembl_id"
  ) %>%
  left_join(
    genes_key,
    by = "ensembl_id"
  ) %>%
  pivot_longer(
    cols = -c(
      ensembl_id,
      gene_name
    ),
    names_to = "sample",
    values_to = "expression"
  )


#==============================================================================
# 30. Boxplots dos marcadores selecionados
#==============================================================================

pdf(
  "figures/02_Tcell_marker_boxplots.pdf",
  width = 12,
  height = 8
)


ggplot(
  plot_data,
  aes(
    x = gene_name,
    y = expression
  )
) +
  
  geom_boxplot() +
  
  theme_bw() +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  
  labs(
    title = "Expression of selected T-cell markers",
    x = "Gene",
    y = "log2-CPM"
  )


dev.off()


#==============================================================================
# 31. Correlação entre marcadores
#==============================================================================

cor_genes <- cor(
  t(marker_expression),
  method = "spearman"
)


#==============================================================================
# 32. Heatmap de correlação
#==============================================================================

pdf(
  "figures/03_Tcell_marker_correlation.pdf",
  width = 10,
  height = 9
)


pheatmap(
  cor_genes,
  labels_row = marker_labels$gene_name,
  labels_col = marker_labels$gene_name,
  main = "Correlation among T-cell markers"
)


dev.off()


#==============================================================================
# 33. Programa de citotoxicidade
#==============================================================================

cytotoxicity_annotation <- marker_table %>%
  filter(
    marker_set == "Cytotoxicity",
    present_in_matrix
  )


cytotoxicity_expression <- logCPM[
  cytotoxicity_annotation$ensembl_id,
  ,
  drop = FALSE
]


#==============================================================================
# 34. Escore exploratório de citotoxicidade
#==============================================================================

# IMPORTANTE:
#
# Este é um escore exploratório.
#
# Ele NÃO representa a proporção de células citotóxicas.
#
# O objetivo é resumir a expressão dos genes selecionados em cada amostra.

cytotoxicity_score <- colMeans(
  cytotoxicity_expression
)


cytotoxicity_score <- data.frame(
  sample = names(
    cytotoxicity_score
  ),
  cytotoxicity_score = as.numeric(
    cytotoxicity_score
  )
)


#==============================================================================
# 35. Salvar escore de citotoxicidade
#==============================================================================

write.csv(
  cytotoxicity_score,
  "results/cytotoxicity_score_exploratory.csv",
  row.names = FALSE
)


#==============================================================================
# 36. Distribuição do escore de citotoxicidade
#==============================================================================

pdf(
  "figures/04_cytotoxicity_score.pdf",
  width = 8,
  height = 6
)


ggplot(
  cytotoxicity_score,
  aes(
    x = cytotoxicity_score
  )
) +
  
  geom_histogram(
    bins = 30
  ) +
  
  theme_bw() +
  
  labs(
    title = "Exploratory cytotoxicity score",
    x = "Mean log2-CPM",
    y = "Number of samples"
  )


dev.off()


#==============================================================================
# 37. Salvar matriz de expressão dos marcadores
#==============================================================================

saveRDS(
  marker_expression,
  "results/TCGA_ACC_Tcell_MarkerExpression.rds"
)


#==============================================================================
# 38. Salvar conjuntos de marcadores
#==============================================================================

saveRDS(
  marker_sets,
  "results/TCGA_ACC_Tcell_marker_sets.rds"
)


#==============================================================================
# 39. Salvar anotação gênica utilizada
#==============================================================================

saveRDS(
  gene_annotation,
  "results/TCGA_ACC_GeneAnnotation.rds"
)


###############################################################################
#
# Fim da Aula 03
#
###############################################################################