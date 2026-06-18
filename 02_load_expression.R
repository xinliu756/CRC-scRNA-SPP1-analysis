# ================================================
# CRC scRNA-seq Analysis
# Script 02: Load Expression Matrix
# Author: Xin Liu
# Date: 2025
# ================================================

library(Seurat)
library(tidyverse)

# ── 参数 ──────────────────────────────────────
TPM_PATH <- "data/GSE146771_CRC.Leukocyte.10x.TPM.txt.gz"
META_PATH <- "data/GSE146771_CRC.Leukocyte.10x.Metadata.txt.gz"
# ─────────────────────────────────────────────

# 1. 读取元数据
cat("读取元数据...\n")
metadata <- read.table(META_PATH, header=TRUE, 
                       sep="\t", row.names=1)
cat("元数据：", nrow(metadata), "个细胞\n")

# 2. 读取表达矩阵（较大，需要几分钟）
cat("读取表达矩阵，请稍等...\n")
start_time <- Sys.time()
#先看一下数据类型再设置空格符，这一步数据读取这一步容易出错
counts <- read.table(TPM_PATH, header=TRUE,
                     sep=" ",
                     row.names=1,
                     quote="\"")

end_time <- Sys.time()
cat("读取完成！耗时：", round(end_time - start_time, 1), "秒\n")
cat("表达矩阵维度：", dim(counts), "\n")
cat("基因数：", nrow(counts), "\n")
cat("细胞数：", ncol(counts), "\n")

# 3. 检查LOXL2是否在数据里
cat("\nLOXL2是否存在：", "LOXL2" %in% rownames(counts), "\n")
cat("ETV4是否存在：", "ETV4" %in% rownames(counts), "\n")

# 保存为RDS格式，下次秒读
cat("保存数据中...\n")
saveRDS(counts, "data/counts_raw.rds")
saveRDS(metadata, "data/metadata.rds")
cat("保存完成！\n")