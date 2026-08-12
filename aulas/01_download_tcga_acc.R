#==============================================================================
# Projeto:
# ACC-Tcell-Differentiation
#
# Aula 01
# Download dos dados TCGA-ACC utilizando o pacote TCGAbiolinks
#
# Objetivos:
# 1. Conhecer a coorte TCGA-ACC
# 2. Aprender a utilizar o TCGAbiolinks
# 3. Baixar dados de RNA-seq do TCGA-ACC
# 4. Explorar um objeto SummarizedExperiment
# 5. Organizar dados para análises de imunoinformática
#
# Autor: Jean Resende
#==============================================================================

#==============================================================================
# 1. Instalação dos pacotes (executar apenas uma vez)
#==============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "TCGAbiolinks",
  "SummarizedExperiment",
  "DESeq2"
))

install.packages("tidyverse")

#==============================================================================
# 2. Carregar os pacotes
#==============================================================================

library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)

#==============================================================================
# 3. Criar estrutura do projeto
#==============================================================================

dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("scripts", showWarnings = FALSE)

#==============================================================================
# 4. Construindo a consulta ao GDC
#==============================================================================

query <- GDCquery(
  project = "TCGA-ACC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# Visualizar o objeto da consulta
query

#==============================================================================
# 5. Download dos dados
#==============================================================================

GDCdownload(
  query,
  method = "api",
  files.per.chunk = 20
)

#==============================================================================
# 6. Preparar os dados
#==============================================================================

acc <- GDCprepare(query)

#==============================================================================
# 7. Explorando o objeto
#==============================================================================

class(acc)

acc

# Dimensão do objeto
dim(acc)

#==============================================================================
# 8. Matriz de expressão
#==============================================================================

counts <- assay(acc)

dim(counts)

head(counts[, 1:5])

#==============================================================================
# 9. Informações dos genes
#==============================================================================

gene_info <- rowData(acc)

head(gene_info)

#==============================================================================
# 10. Informações das amostras
#==============================================================================

sample_info <- colData(acc) %>%
  as.data.frame()

head(sample_info)

#==============================================================================
# 11. Metadados do objeto
#==============================================================================

metadata(acc)

#==============================================================================
# 12. Identificadores das amostras
#==============================================================================

samples <- colnames(counts)

head(samples)

#==============================================================================
# 13. Identificadores dos genes
#==============================================================================

head(rownames(counts))

#==============================================================================
# 14. Removendo a versão do Ensembl
#==============================================================================

rownames(counts) <- sub("\\\\..*", "", rownames(counts))

head(rownames(counts))

#==============================================================================
# 15. Selecionando apenas amostras tumorais primárias
#==============================================================================

tumor_samples <- TCGAquery_SampleTypes(
  barcode = colnames(counts),
  typesample = "TP"
)

counts_tumor <- counts[, tumor_samples]

sample_info_tumor <- sample_info %>%
  filter(barcode %in% tumor_samples)

dim(counts_tumor)

#==============================================================================
# 16. Verificar presença de dados clínicos relevantes
#==============================================================================

colnames(sample_info_tumor)

#==============================================================================
# 17. Salvar objetos
#==============================================================================

saveRDS(acc,
        "data/TCGA_ACC_SummarizedExperiment.rds")

saveRDS(counts_tumor,
        "data/TCGA_ACC_Counts_Tumor.rds")

saveRDS(sample_info_tumor,
        "data/TCGA_ACC_SampleInfo_Tumor.rds")

saveRDS(gene_info,
        "data/TCGA_ACC_GeneInfo.rds")

#==============================================================================
# Fim da Aula 01
#==============================================================================