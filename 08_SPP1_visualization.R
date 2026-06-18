# ================================================
# CRC scRNA-seq Analysis
# Script 05: Visualization
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

# 重新计算SPP1数据
spp1_expr <- FetchData(seurat_obj, vars = "SPP1")
spp1_by_cluster <- data.frame(
  cell_type = seurat_obj$Global_Cluster,
  SPP1 = spp1_expr$SPP1,
  tissue = seurat_obj$Tissue
)

# ── 图1：各细胞类型SPP1阳性比例（柱状图）──────
result_tissue <- spp1_by_cluster %>%
  filter(tissue %in% c("T", "N")) %>%
  group_by(tissue, cell_type) %>%
  summarise(
    total_cells = n(),
    SPP1_positive = sum(SPP1 > 0),
    SPP1_ratio = round(SPP1_positive / total_cells * 100, 2),
    .groups = "drop"
  )

p1 <- ggplot(result_tissue, 
             aes(x = cell_type, y = SPP1_ratio, fill = tissue)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("N" = "#4A9FE8", "T" = "#E85555"),
                    labels = c("N" = "Normal", "T" = "Tumor")) +
  labs(
    title = "SPP1+ Cell Ratio in Tumor vs Normal Tissue",
    x = "Cell Type",
    y = "SPP1+ Cells (%)",
    fill = "Tissue"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# ── 图2：Myeloid cell中SPP1表达分布（小提琴图）──
myeloid_data <- spp1_by_cluster %>%
  filter(cell_type == "Myeloid cell",
         tissue %in% c("T", "N"))

p2 <- ggplot(myeloid_data, 
             aes(x = tissue, y = SPP1, fill = tissue)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.8) +
  scale_fill_manual(values = c("N" = "#4A9FE8", "T" = "#E85555"),
                    labels = c("N" = "Normal", "T" = "Tumor")) +
  scale_x_discrete(labels = c("N" = "Normal", "T" = "Tumor")) +
  labs(
    title = "SPP1 Expression in Myeloid Cells",
    x = "Tissue",
    y = "SPP1 Expression (TPM)",
    fill = "Tissue"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

# ── 合并两张图并保存 ───────────────────────────
combined_plot <- p1 + p2 +
  plot_annotation(
    title = "SPP1 Expression in CRC Tumor Microenvironment",
    subtitle = "GSE146771 | 10x scRNA-seq | 10 CRC patients",
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray50")
    )
  )

ggsave("results/figures/SPP1_analysis.png", 
       combined_plot,
       width = 12, height = 6, dpi = 300)

cat("图片已保存到 results/figures/SPP1_analysis.png\n")

# 顺便看一下图
print(combined_plot)