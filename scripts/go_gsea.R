# ============================================================
# go_gsea.R — 通路富集分析
# GSE96870: 上呼吸道感染对CNS转录组的影响
# ============================================================

library(clusterProfiler)
library(org.Mm.eg.db)      # 小鼠基因数据库
library(BiocParallel)
library(ggplot2)
library(enrichplot)

# 并行设置
register(MulticoreParam(4))

# ============================================================
# 1. 读取DESeq2结果
# ============================================================
results_dir <- "/home/users/nus/e1668109/scratch/5004/results/deseq2"
gsea_dir    <- "/home/users/nus/e1668109/scratch/5004/results/gsea"
dir.create(gsea_dir, showWarnings = FALSE, recursive = TRUE)

# 读取第一个contrast结果
res <- read.csv(file.path(results_dir, "Cerebellum_Day8_vs_Day0.csv"))

# 过滤显著差异基因
sig_genes <- res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]
cat("显著差异基因数:", nrow(sig_genes), "\n")

# ============================================================
# 2. 基因ID转换（ENSEMBL → ENTREZ）
# clusterProfiler需要ENTREZ ID
# ============================================================
gene_ids <- bitr(sig_genes$gene,
                 fromType = "ENSEMBL",
                 toType   = "ENTREZID",
                 OrgDb    = org.Mm.eg.db)

cat("成功转换基因数:", nrow(gene_ids), "\n")

# ============================================================
# 3. GO富集分析 — BP/MF/CC 三个层次并行
# ============================================================
cat("\n开始GO富集分析（并行跑BP/MF/CC）...\n")

go_ontologies <- c("BP", "MF", "CC")

go_results <- bplapply(go_ontologies, function(ont) {
    cat("正在分析:", ont, "\n")
    result <- enrichGO(
        gene          = gene_ids$ENTREZID,
        OrgDb         = org.Mm.eg.db,
        ont           = ont,
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.05,
        readable      = TRUE
    )
    
    # 保存结果
    if (!is.null(result) && nrow(result@result) > 0) {
        write.csv(result@result,
                  file.path(gsea_dir, paste0("GO_", ont, ".csv")),
                  row.names = FALSE)
        
        # 画图
        p <- dotplot(result, showCategory = 20, title = paste("GO", ont))
        ggsave(file.path(gsea_dir, paste0("GO_", ont, "_dotplot.pdf")),
               p, width = 10, height = 8)
        
        cat(ont, "找到", nrow(result@result), "个显著通路\n")
    }
    return(result)
}, BPPARAM = MulticoreParam(3))    # 3个本体同时跑

names(go_results) <- go_ontologies

# ============================================================
# 4. GSEA分析 — 用所有基因的fold change排序
# ============================================================
cat("\n开始GSEA分析...\n")

# 准备ranked gene list（所有基因，按log2FC排序）
all_genes <- res[!is.na(res$log2FoldChange) & !is.na(res$padj), ]

# 转换基因ID
all_gene_ids <- bitr(all_genes$gene,
                     fromType = "ENSEMBL",
                     toType   = "ENTREZID",
                     OrgDb    = org.Mm.eg.db)

all_genes <- merge(all_genes, all_gene_ids,
                   by.x = "gene", by.y = "ENSEMBL")

# 创建ranked list
gene_list <- all_genes$log2FoldChange
names(gene_list) <- all_genes$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)
gene_list <- gene_list[!duplicated(names(gene_list))]

cat("GSEA基因数:", length(gene_list), "\n")

# 跑GSEA
gsea_result <- gseGO(
    geneList      = gene_list,
    OrgDb         = org.Mm.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    verbose       = TRUE
)

if (!is.null(gsea_result) && nrow(gsea_result@result) > 0) {
    write.csv(gsea_result@result,
              file.path(gsea_dir, "GSEA_BP.csv"),
              row.names = FALSE)
    
    # GSEA enrichment plot（前5个通路）
    p <- gseaplot2(gsea_result,
                   geneSetID = 1:min(5, nrow(gsea_result@result)),
                   title = "GSEA — Infected vs Mock (Day4)")
    ggsave(file.path(gsea_dir, "GSEA_plot.pdf"), p, width = 12, height = 8)
    
    cat("GSEA找到", nrow(gsea_result@result), "个显著通路\n")
}

cat("\nGO/GSEA分析完成！结果保存在:", gsea_dir, "\n")
