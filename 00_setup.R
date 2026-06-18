# ================================================
# CRC scRNA-seq Analysis
# Script 00: Environment Setup
# Author: Xin Liu
# Date: 2025
# ================================================

# 1. 检查并加载必要的包
packages <- c("Seurat", "tidyverse", "ggplot2", 
              "patchwork", "harmony")

for(pkg in packages){
  if(!requireNamespace(pkg, quietly = TRUE)){
    cat(pkg, "未安装，正在安装...\n")
    install.packages(pkg)
  } else {
    cat(pkg, "已安装 ✓\n")
  }
}

# 2. 加载包
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(harmony)

cat("\n所有包加载完成！\n")

# 3. 创建项目目录结构
dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("results/figures", showWarnings = FALSE)
dir.create("results/tables", showWarnings = FALSE)
dir.create("scripts", showWarnings = FALSE)

cat("目录结构创建完成！\n")

# 4. 检查R版本
cat("\nR版本信息：\n")
R.version.string