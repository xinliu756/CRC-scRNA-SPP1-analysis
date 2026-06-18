#检查需要的包是不是已安装
packages <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot")

for(pkg in packages){
  if(!requireNamespace(pkg, quietly = TRUE)){
    cat(pkg, "未安装\n")
  } else {
    cat(pkg, "已安装 ✓\n")
  }
}

# ================================================
# CRC scRNA-seq Analysis
# Script 08: GSEA Analysis
# Author: Xin Liu
# Date: 2025
# ================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(tidyverse)
library(ggplot2)

# 1. 读取DEG结果
cat("读取DEG结果...\n")
deg_results <- read.csv(
  "results/tables/DEG_Myeloid_Tumor_vs_Normal.csv",
  row.names = 1
)
deg_results$gene <- rownames(deg_results)

# 2. 准备GSEA输入：按log2FC排序的基因列表
cat("准备基因列表...\n")
gene_list <- deg_results$avg_log2FC
names(gene_list) <- deg_results$gene
gene_list <- sort(gene_list, decreasing = TRUE)

cat("基因列表长度：", length(gene_list), "\n")
cat("最高FC基因：", names(gene_list)[1], "\n")
cat("最低FC基因：", names(gene_list)[length(gene_list)], "\n")

# 3. 运行GSEA - GO分析
cat("\n运行GSEA-GO分析，请稍等...\n")
gsea_go <- gseGO(
  geneList = gene_list,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",           # 生物学过程
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

cat("GSEA-GO完成！显著通路数：", nrow(gsea_go@result), "\n")

# 4. 查看top通路
cat("\n上调top10通路：\n")
top_up <- gsea_go@result %>%
  filter(NES > 0) %>%
  arrange(desc(NES)) %>%
  head(10)
print(top_up[, c("Description", "NES", "p.adjust")])

cat("\n下调top10通路：\n")
top_down <- gsea_go@result %>%
  filter(NES < 0) %>%
  arrange(NES) %>%
  head(10)
print(top_down[, c("Description", "NES", "p.adjust")])

# 5. 可视化
cat("\n生成图表...\n")

# 点图
p1 <- dotplot(gsea_go, 
              showCategory = 15,
              split = ".sign") +
  facet_grid(. ~ .sign) +
  labs(title = "GSEA-GO: Myeloid Cells Tumor vs Normal") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("results/figures/GSEA_GO_dotplot.png",
       p1, width = 12, height = 8, dpi = 300)

# 6. 保存结果
write.csv(gsea_go@result,
          "results/tables/GSEA_GO_results.csv",
          row.names = FALSE)

cat("所有结果已保存！\n")


# 7. 单通路富集曲线 - angiogenesis
cat("生成angiogenesis富集曲线...\n")

# 找到angiogenesis的ID
angio_id <- gsea_go@result %>%
  filter(grepl("angiogenesis", Description, ignore.case = TRUE)) %>%
  pull(ID) %>%
  head(1)

cat("angiogenesis通路ID：", angio_id, "\n")

# 画富集曲线
p2 <- gseaplot2(
  gsea_go,
  geneSetID = angio_id,
  title = "GSEA: Angiogenesis",
  color = "#E85555",
  pvalue_table = TRUE
)

ggsave("results/figures/GSEA_angiogenesis_curve.png",
       p2, width = 8, height = 6, dpi = 300)

cat("富集曲线已保存！\n")
print(p2)