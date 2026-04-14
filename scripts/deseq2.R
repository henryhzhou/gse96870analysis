# ============================================================
# deseq2.R — 差异表达分析
# GSE96870: 上呼吸道感染对CNS转录组的影响
# ============================================================

library(DESeq2)
library(BiocParallel)
library(ggplot2)
library(pheatmap)

register(MulticoreParam(4))

# ============================================================
# 1. 读取count矩阵
# ============================================================
counts_file <- "/home/users/nus/e1668109/scratch/5004/results/counts/counts.txt"

counts <- read.table(counts_file,
                     header = TRUE,
                     skip = 1,
                     row.names = 1,
                     check.names = FALSE)

# 只保留样本列
counts <- counts[, 6:ncol(counts)]

# 清理列名，只保留SRR编号
colnames(counts) <- gsub(".Aligned.sortedByCoord.out.bam", "", colnames(counts))
colnames(counts) <- basename(colnames(counts))

cat("Count矩阵维度:", dim(counts), "\n")

# ============================================================
# 2. 读取metadata
# ============================================================
meta <- read.csv("/home/users/nus/e1668109/scratch/5004/pipeline/scripts/metadata.csv",
                 row.names = 1)

# 确保顺序一致
meta <- meta[colnames(counts), ]

# 转换为factor
meta$tissue    <- factor(meta$tissue)
meta$time      <- factor(meta$time, levels = c("Day 0", "Day 4", "Day 8"))
meta$infection <- factor(meta$infection, levels = c("Non-Infected", "Influenza A"))

cat("样本分组:\n")
print(table(meta$tissue, meta$time))

# ============================================================
# 3. 分组分析：Cerebellum 和 Spinal cord 分开跑
# ============================================================
results_dir <- "/home/users/nus/e1668109/scratch/5004/results/deseq2"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

run_deseq2 <- function(tissue_name) {
    cat("\n========================================\n")
    cat("分析组织:", tissue_name, "\n")
    
    # 筛选该组织的样本
    idx    <- meta$tissue == tissue_name
    counts_sub <- counts[, idx]
    meta_sub   <- meta[idx, ]
    
    # 创建DESeq2对象
    dds <- DESeqDataSetFromMatrix(
        countData = counts_sub,
        colData   = meta_sub,
        design    = ~ time
    )
    
    # 过滤低表达基因
    keep <- rowSums(counts(dds) > 1) >= 3
    dds  <- dds[keep, ]
    cat("过滤后基因数:", nrow(dds), "\n")
    
    # 运行DESeq2
    dds <- DESeq(dds, parallel = TRUE)
    
    # 提取多个contrast并行
    contrasts <- list(
        Day4_vs_Day0 = c("time", "Day 4", "Day 0"),
        Day8_vs_Day0 = c("time", "Day 8", "Day 0"),
        Day8_vs_Day4 = c("time", "Day 8", "Day 4")
    )
    
    res_list <- bplapply(names(contrasts), function(name) {
        res <- results(dds, contrast = contrasts[[name]], alpha = 0.05)
        res <- as.data.frame(res)
        res$gene <- rownames(res)
        tissue_clean <- gsub(" ", "_", tissue_name)
        write.csv(res,
                  file.path(results_dir, paste0(tissue_clean, "_", name, ".csv")),
                  row.names = FALSE)
        n_sig <- sum(!is.na(res$padj) & res$padj < 0.05, na.rm = TRUE)
        cat(name, "— 显著差异基因:", n_sig, "\n")
        return(res)
    }, BPPARAM = MulticoreParam(3))
    
    names(res_list) <- names(contrasts)
    
    # PCA图
    vsd <- vst(dds, blind = FALSE)
    pca_data <- plotPCA(vsd, intgroup = c("time", "infection"), returnData = TRUE)
    pca_var  <- round(100 * attr(pca_data, "percentVar"))
    
    p <- ggplot(pca_data, aes(PC1, PC2, color = time, shape = infection)) +
        geom_point(size = 4) +
        xlab(paste0("PC1: ", pca_var[1], "%")) +
        ylab(paste0("PC2: ", pca_var[2], "%")) +
        ggtitle(paste("PCA —", tissue_name)) +
        theme_bw()
    
    tissue_clean <- gsub(" ", "_", tissue_name)
    ggsave(file.path(results_dir, paste0(tissue_clean, "_PCA.pdf")), p,
           width = 8, height = 6)
    
    return(res_list)
}

# 两个组织并行分析
tissues <- c("Cerebellum", "Spinal cord")
all_results <- bplapply(tissues, run_deseq2, BPPARAM = MulticoreParam(2))
names(all_results) <- tissues

cat("\nDESeq2分析完成！结果保存在:", results_dir, "\n")
