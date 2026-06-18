# ================================================
# CRC scRNA-seq Analysis
# Script 03: Create Seurat Object
# Author: Xin Liu
# Date: 2025
# ================================================

library(Seurat)
library(tidyverse)

# ── 参数 ──────────────────────────────────────
MIN_CELLS <- 3      # 基因至少在3个细胞里表达
MIN_FEATURES <- 200 # 每个细胞至少检测到200个基因
# ─────────────────────────────────────────────

# 1. 读取保存好的数据
cat("读取数据...\n")
counts <- readRDS("data/counts_raw.rds")
metadata <- readRDS("data/metadata.rds")

# 2. 创建Seurat对象
cat("创建Seurat对象...\n")
seurat_obj <- CreateSeuratObject(
  counts = counts,
  meta.data = metadata,
  min.cells = MIN_CELLS,
  min.features = MIN_FEATURES,
  project = "CRC_scRNA"
)

cat("Seurat对象创建完成！\n")
print(seurat_obj)

# 3. 检查LOXL2
cat("\nLOXL2表达情况：\n")
loxl2_expr <- FetchData(seurat_obj, vars = "LOXL2")
cat("表达LOXL2的细胞数：", sum(loxl2_expr > 0), "\n")
cat("LOXL2表达范围：", range(loxl2_expr$LOXL2), "\n")

# 4. 保存
saveRDS(seurat_obj, "data/seurat_raw.rds")
cat("保存完成！\n")