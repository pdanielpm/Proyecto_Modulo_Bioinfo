load("results/DE_limma_results.RData")

library(ggplot2)
library(matrixStats)
library(survival)
library(survminer)
library(ggpubr)


## localizar SPN por símbolo
if (is.null(gene_col)) {
  stop("No hay columna de símbolo en rowData para ubicar SPN.")
}

spn_row <- which(rowdata_final[[gene_col]] == "SPN")
if (length(spn_row) == 0) {
  stop("SPN no encontrado en rowData (símbolos).")
}

SPN_expr <- as.numeric(logcpm[spn_row, ])

df_spn <- data.frame(
  status = group,
  expr = SPN_expr
)
df_spn <- df_spn[!is.na(df_spn$status), ]

violinplt <- ggplot(df_spn, aes(status, expr, fill = status)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.10, alpha = 0.6, size = 1.2) +
  stat_compare_means(method = "wilcox.test") +
  theme_classic() +
  labs(
    title = "Expresión de SPN (CD43) vs estado vital",
    y = "logCPM (TMM)",
    x = ""
  )

print(violinplt)
ggsave(
  "plots/SPN_violin.png",
  width = 6,
  height = 5,
  dpi = 300,
  plot = violinplt
)

d_death <- suppressWarnings(as.numeric(as.character(
  coldata_final$tcga.gdc_cases.diagnoses.days_to_death
)))
d_fu <- suppressWarnings(as.numeric(as.character(
  coldata_final$tcga.gdc_cases.diagnoses.days_to_last_follow_up
)))

time <- ifelse(!is.na(d_death), d_death, d_fu)
event <- ifelse(group == "dead", 1, 0)

expr_group <- ifelse(SPN_expr > median(SPN_expr, na.rm = TRUE), "High", "Low")
expr_group <- factor(expr_group, levels = c("Low", "High"))

surv_df <- data.frame(
  time = as.numeric(time),
  event = as.integer(event),
  expr_group = expr_group
)

surv_df <- surv_df[complete.cases(surv_df), ]
surv_df <- surv_df[surv_df$time >= 0, ]

fit <- survfit(Surv(time, event) ~ expr_group, data = surv_df)

surv_plot <- ggsurvplot(
  fit,
  data = surv_df,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  ggtheme = theme_classic(),
  title = "Overall survival por SPN (High vs Low)"
)

ggsave(
  "plots/SPN_survival.png",
  width = 6,
  height = 5,
  dpi = 300,
  plot = surv_plot
)
pca <- prcomp(t(logcpm), scale. = TRUE)

df_pc <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  SPN = SPN_expr,
  status = group
)

ggplot(df_pc, aes(SPN, PC1, color = status)) +
  geom_point(size = 2.7) +
  theme_classic() +
  labs(
    title = "Relación entre SPN y el eje principal de variación transcriptómica (PC1)",
    x = "SPN (logCPM)",
    y = "PC1"
  )
