#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(optparse)
})

option_list <- list(
  make_option(c("--counts"), type="character", help="counts/Rloop_highconf_universe_counts.tsv"),
  make_option(c("--sample-table"), type="character", help="sample_table.tsv"),
  make_option(c("--condition-a"), type="character", help="numerator condition, e.g. KO"),
  make_option(c("--condition-b"), type="character", help="denominator condition, e.g. WT"),
  make_option(c("--out-prefix"), type="character", default="diff_Rloop")
)

opt <- parse_args(OptionParser(option_list=option_list))

counts <- read.table(opt$counts, sep="\t", header=TRUE, check.names=FALSE)
meta <- read.table(opt$sample.table, sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)

peak_id <- paste(counts$chr, counts$start, counts$end, sep=":")
count_mat <- counts[, !(colnames(counts) %in% c("chr", "start", "end")), drop=FALSE]
rownames(count_mat) <- peak_id

meta <- meta[meta$sample %in% colnames(count_mat), ]
meta <- meta[tolower(meta$call_peak) %in% c("yes", "true", "1"), ]
rownames(meta) <- meta$sample
count_mat <- count_mat[, rownames(meta), drop=FALSE]

meta$condition <- factor(meta$condition)
dds <- DESeqDataSetFromMatrix(
  countData=round(as.matrix(count_mat)),
  colData=meta,
  design=~ condition
)
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition", opt$condition.a, opt$condition.b))
res_df <- as.data.frame(res)
res_df$peak_id <- rownames(res_df)
res_df$chr <- sub(":.*", "", res_df$peak_id)
res_df$start <- as.integer(sub(".*:(.*):.*", "\\1", res_df$peak_id))
res_df$end <- as.integer(sub(".*:", "", res_df$peak_id))
res_df <- res_df[, c("chr", "start", "end", "peak_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

write.table(res_df, paste0(opt$out.prefix, ".tsv"), sep="\t", row.names=FALSE, quote=FALSE)

sig <- subset(res_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
write.table(sig, paste0(opt$out.prefix, ".padj0.05_lfc1.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
