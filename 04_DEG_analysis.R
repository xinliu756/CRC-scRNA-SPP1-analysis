
# ================================================
# CRC scRNA-seq Analysis
# Script 06: DEG Analysis - Tumor vs Normal
# Author: Xin Liu
# Date: 2025
# ================================================

library(Seurat)
library(tidyverse)
library(ggplot2)

# 1. 读取数据
cat("读取数据...\n")
seurat_obj <- readRDS("data/seurat_normalized.rds")#此下载TPM是已标准化的数据，但是是counts里面，data里面是空的，重新标准化了保存了直接再读取进行DEG，也可以把数据同时放到data层就能读了。

# 2. 只保留Myeloid cell，肿瘤和正常组织
cat("提取Myeloid cell...\n")
myeloid <- subset(seurat_obj, 
                  Global_Cluster == "Myeloid cell" & 
                    Tissue %in% c("T", "N"))

cat("肿瘤Myeloid cell数量：", sum(myeloid$Tissue == "T"), "\n")
cat("正常Myeloid cell数量：", sum(myeloid$Tissue == "N"), "\n")

# 3. 设置分组标识
Idents(myeloid) <- myeloid$Tissue

# 4. 做DEG分析
cat("正在做DEG分析，请稍等...\n")
deg_results <- FindMarkers(
  myeloid,
  ident.1 = "T",      # 肿瘤组
  ident.2 = "N",      # 正常组
  min.pct = 0.1,      # 至少在10%的细胞里表达
  logfc.threshold = 0.25  # 最小fold change阈值
)

# 5. 查看结果
cat("\nDEG分析完成！\n")
cat("总差异基因数：", nrow(deg_results), "\n")

# 加上基因名列
deg_results$gene <- rownames(deg_results)

# 看top20上调基因
cat("\n肿瘤组织上调Top20基因：\n")
top20_up <- deg_results %>%
  filter(avg_log2FC > 0, p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  head(20)
print(top20_up[, c("gene", "avg_log2FC", "p_val_adj")])

# SPP1在不在里面？
cat("\nSPP1是否在差异基因里：", "SPP1" %in% deg_results$gene, "\n")
if("SPP1" %in% deg_results$gene){
  cat("SPP1的结果：\n")
  print(deg_results["SPP1", ])
}

# 6. 保存结果
write.csv(deg_results, 
          "results/tables/DEG_Myeloid_Tumor_vs_Normal.csv",
          row.names = TRUE)
cat("\n结果已保存！\n")