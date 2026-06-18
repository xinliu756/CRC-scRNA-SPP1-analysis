# ================================================
# CRC scRNA-seq Analysis
# Script 07: Volcano Plot
# Author: Xin Liu
# Date: 2025
# ================================================

library(tidyverse)
library(ggplot2)

# 1. 读取DEG结果
deg_results <- read.csv(
  "results/tables/DEG_Myeloid_Tumor_vs_Normal.csv",
  row.names = 1
)
deg_results$gene <- rownames(deg_results)

# 2. 添加分类标签
deg_results <- deg_results %>%
  mutate(
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0.25  ~ "Up in Tumor",
      p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Down in Tumor",
      TRUE ~ "Not significant"
    )
  )

# 3. 标注感兴趣的基因
genes_to_label <- c("SPP1", "FN1", "TREM2", "ANGPTL4", 
                    "APOE", "MMP12", "CLEC5A", "NUPR1")

deg_results <- deg_results %>%
  mutate(label = ifelse(gene %in% genes_to_label, gene, ""))

# 4. 画火山图
p <- ggplot(deg_results, 
            aes(x = avg_log2FC, 
                y = -log10(p_val_adj),
                color = significance,
                label = label)) +
  
  geom_point(alpha = 0.6, size = 1.5) +
  
  # 标注基因名
  ggrepel::geom_text_repel(
    data = subset(deg_results, label != ""),
    size = 3.5,
    fontface = "bold",
    max.overlaps = 20,
    box.padding = 0.5
  ) +
  
  # 颜色
  scale_color_manual(values = c(
    "Up in Tumor"     = "#E85555",
    "Down in Tumor"   = "#4A9FE8",
    "Not significant" = "grey70"
  )) +
  
  # 阈值线
  geom_vline(xintercept = c(-0.25, 0.25), 
             linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), 
             linetype = "dashed", color = "grey50") +
  
  # 标题
  labs(
    title = "DEG: Myeloid Cells in Tumor vs Normal",
    subtitle = "CRC Tumor Microenvironment | GSE146771",
    x = "avg log2 Fold Change",
    y = "-log10 (adjusted p-value)",
    color = ""
  ) +
  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle = element_text(hjust = 0.5, color = "grey50"),
    legend.position = "top"
  )

# 5. 保存
ggsave("results/figures/volcano_DEG_myeloid.png",
       p, width = 8, height = 7, dpi = 300)

cat("火山图已保存！\n")
print(p)