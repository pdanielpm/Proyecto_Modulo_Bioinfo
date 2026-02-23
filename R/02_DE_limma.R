load("results/datos_DLBC_filtrados.RData")

library(Matrix)
library(edgeR)
library(limma)
library(ggplot2)
library(matrixStats)
library(pheatmap)


v <- tolower(trimws(as.character(coldata_final[[vital_col]])))

group <- ifelse(
  v %in% c("alive", "living"),
  "alive",
  ifelse(v %in% c("dead", "deceased"), "dead", NA)
)

keep_ad <- !is.na(group)
cts_ad <- cts_final[, keep_ad]
cd_ad <- coldata_final[keep_ad, , drop = FALSE]
group <- factor(group[keep_ad], levels = c("alive", "dead"))

cat("Tabla alive/dead:\n")
print(table(group))


y <- DGEList(counts = cts_ad, group = group)
keep_gene2 <- filterByExpr(y, group = group)
y <- y[keep_gene2, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)


svg("plots/Plot_limma.svg", width = 6, height = 5)
v <- voom(y, design, plot = TRUE)

dev.off()

fit <- lmFit(v, design)


# contraste: dead vs alive (logFC>0 => más alto en dead)
contr <- makeContrasts(DeadVsAlive = dead - alive, levels = design)
fit.contr <- contrasts.fit(fit, contr)

fit.eb <- eBayes(fit.contr)
res <- topTable(fit.eb, coef = "DeadVsAlive", number = Inf, sort.by = "P")

## --- anotar símbolo si existe
gene_col_candidates <- intersect(
  c("gene_name", "gene", "symbol", "external_gene_name"),
  colnames(rowdata_final)
)
gene_col <- if (length(gene_col_candidates) > 0) {
  gene_col_candidates[1]
} else {
  NULL
}

res$gene <- if (!is.null(gene_col)) {
  rowdata_final[[gene_col]][match(rownames(res), rownames(rowdata_final))]
} else {
  NA_character_
}

head(res, 20)

svg("plots/Volcano_limma.svg", width = 6, height = 5)
vc_plot <- volcanoplot(
  fit.eb,
  coef = "DeadVsAlive",
  highlight = 100,
  names = res$gene,
)
dev.off()


A <- res$AveExpr
M <- res$logFC


svg("plots/MA_limma.svg", width = 6, height = 5)
MAplt <- plot(
  A,
  M,
  pch = 16,
  cex = 0.45,
  xlab = "Average expression (logCPM)",
  ylab = "logFC (dead vs alive)",
  ylim = c(-5, 5)
)
abline(h = 0, col = "red")
dev.off()


logcpm <- v$E

top_n <- 20
vars <- matrixStats::rowVars(logcpm)
top_idx <- order(vars, decreasing = TRUE)[seq_len(min(top_n, nrow(logcpm)))]
mat <- logcpm[top_idx, , drop = FALSE]

## renombrar filas con símbolo si existe
if (!is.null(gene_col)) {
  labs <- rowdata_final[[gene_col]][match(
    rownames(mat),
    rownames(rowdata_final)
  )]
  labs[is.na(labs) | labs == ""] <- rownames(mat)[is.na(labs) | labs == ""]
  rownames(mat) <- make.unique(labs)
}

## anotar columnas
ann_col <- data.frame(vital_status = group)
rownames(ann_col) <- colnames(mat)

## ordenar por grupo
ord <- order(group)
mat <- mat[, ord, drop = FALSE]
ann_col <- ann_col[ord, , drop = FALSE]

pheatmap(
  mat,
  scale = "row",
  show_colnames = FALSE,
  annotation_col = ann_col,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  main = paste0(
    "Top ",
    top_n,
    " genes más variables (logCPM; z-score por gen)"
  ),
  filename = "plots/Heatmap_DE.png",
  width = 6,
  height = 8
)

save(
  res,
  gene_col,
  logcpm,
  group,
  file = "results/DE_limma_results.RData"
)
cat("¡Análisis de DE con limma guardado con éxito!")
