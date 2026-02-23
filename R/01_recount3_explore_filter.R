library(recount3)
library(SummarizedExperiment)

human_projects <- available_projects()

project_info <- subset(
  human_projects,
  project == "DLBC" &
    project_type == "data_sources" &
    project_home == "data_sources/tcga"
)

rse_DLBC <- create_rse(project_info)

assay(rse_DLBC, "counts") <- compute_read_counts(rse_DLBC)


cts <- assay(rse_DLBC, "counts")
cd <- as.data.frame(colData(rse_DLBC))
rd <- as.data.frame(rowData(rse_DLBC))

cat("Muestras iniciales:", ncol(cts), "\n")
cat("Genes iniciales:", nrow(cts), "\n")


## --- detectar columna de estado vital (robusto)
vital_candidates <- grep(
  "vital|vital_status|deceased|dead|alive",
  names(cd),
  ignore.case = TRUE,
  value = TRUE
)

if (length(vital_candidates) == 0) {
  stop(
    "No se encontraron columnas tipo vital_status en colData(). Revisa names(cd)."
  )
}

# Elegir la columna más informativa (por defecto, la primera)
vital_col <- vital_candidates[1]
cat("Usando columna vital:", vital_col, "\n")
print(table(cd[[vital_col]], useNA = "ifany"))

## --- remover NA / unknown / not reported
bad_vital <- is.na(cd[[vital_col]]) |
  trimws(tolower(as.character(cd[[vital_col]]))) %in%
    c("", "na", "nan", "not reported", "unknown")

keep_vital <- !bad_vital

rse_vital <- rse_DLBC[, keep_vital]
cts_vital <- assay(rse_vital, "counts")
cd_vital <- as.data.frame(colData(rse_vital))

cat("Muestras tras filtrar vital:", ncol(cts_vital), "\n")
print(table(cd_vital[[vital_col]], useNA = "ifany"))

## --- QC por muestra
lib_size <- Matrix::colSums(cts_vital)
genes_det <- Matrix::colSums(cts_vital > 0)
zero_frac <- Matrix::colMeans(cts_vital == 0)

qc <- data.frame(
  sample = colnames(cts_vital),
  lib_size = as.numeric(lib_size),
  genes_detected = as.numeric(genes_det),
  zero_fraction = as.numeric(zero_frac)
)

summary(qc)

min_lib <- 5e6
min_genes <- 10000
max_zeros <- 0.85

keep_qc <- with(
  qc,
  lib_size >= min_lib & genes_detected >= min_genes & zero_fraction <= max_zeros
)
cat("Muestras que pasan QC:", sum(keep_qc), "de", length(keep_qc), "\n")

rse_qc <- rse_vital[, keep_qc]
cts_qc <- assay(rse_qc, "counts")

## --- filtro de genes (>=10 counts en >=10% de muestras)
n <- ncol(cts_qc)
keep_gene <- Matrix::rowSums(cts_qc >= 10) >= ceiling(0.10 * n)

rse_filt <- rse_qc[keep_gene, ]
cts_final <- assay(rse_filt, "counts")
coldata_final <- as.data.frame(colData(rse_filt))
rowdata_final <- as.data.frame(rowData(rse_filt))

cat("Genes finales:", nrow(cts_final), "\n")
cat("Muestras finales:", ncol(cts_final), "\n")

save(
  rse_filt,
  cts_final,
  coldata_final,
  rowdata_final,
  file = "results/datos_DLBC_filtrados.RData"
)
cat("¡Datos guardados con éxito!")
