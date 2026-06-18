# ================================================
# CRC scRNA-seq Analysis
# Script 04: SPP1 Analysis in TME
# Author: Xin Liu
# Date: 2025
# ================================================

library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)

# 1. 读取数据
cat("读取数据...\n")
seurat_obj <- readRDS("data/seurat_raw.rds")
metadata <- readRDS("data/metadata.rds")

# 2. 获取SPP1表达量
cat("提取SPP1表达数据...\n")
spp1_expr <- FetchData(seurat_obj, vars = "SPP1")
seurat_obj$SPP1_expr <- spp1_expr$SPP1
seurat_obj$SPP1_positive <- ifelse(spp1_expr$SPP1 > 0, "SPP1+", "SPP1-")

# 3. SPP1在不同细胞类型中的表达
cat("\nSPP1在各细胞类型的表达情况：\n")
spp1_by_cluster <- data.frame(
  cell_type = seurat_obj$Global_Cluster,
  SPP1 = seurat_obj$SPP1_expr,
  tissue = seurat_obj$Tissue
)

# 各细胞类型SPP1阳性细胞数
result_cluster <- spp1_by_cluster %>%
  group_by(cell_type) %>%
  summarise(
    total_cells = n(),
    SPP1_positive = sum(SPP1 > 0),
    SPP1_ratio = round(SPP1_positive / total_cells * 100, 2)
  ) %>%
  arrange(desc(SPP1_ratio))

print(result_cluster)

# 4. 肿瘤 vs 正常组织SPP1表达差异
cat("\n肿瘤 vs 正常组织SPP1阳性比例：\n")
result_tissue <- spp1_by_cluster %>%
  filter(tissue %in% c("T", "N")) %>%
  group_by(tissue, cell_type) %>%
  summarise(
    total_cells = n(),
    SPP1_positive = sum(SPP1 > 0),
    SPP1_ratio = round(SPP1_positive / total_cells * 100, 2),
    .groups = "drop"
  ) %>%
  arrange(cell_type, tissue)

print(result_tissue)

# 5. 保存结果表格
write.csv(result_cluster, 
          "results/tables/SPP1_by_celltype.csv", 
          row.names = FALSE)
write.csv(result_tissue, 
          "results/tables/SPP1_tumor_vs_normal.csv", 
          row.names = FALSE)

cat("\n结果已保存！\n")