# ================================================
# CRC scRNA-seq Analysis
# Script 01: Load Data
# Author: Xin Liu
# Date: 2025
# ================================================

library(Seurat)
library(tidyverse)

# 1. 读取元数据
cat("正在读取元数据...\n")
metadata <- read.table(
  "data/GSE146771_CRC.Leukocyte.10x.Metadata.txt.gz",
  header = TRUE,
  sep = "\t",
  row.names = 1
)

cat("元数据维度：", dim(metadata), "\n")
cat("元数据列名：\n")
print(colnames(metadata))

# 2. 预览元数据
head(metadata)

# 看看有哪些组织类型
table(metadata$Tissue)

# 看看有哪些细胞大类
table(metadata$Global_Cluster)

# 看看有几个患者样本
length(unique(metadata$Sample))