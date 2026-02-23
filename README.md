# Expresión Diferencial en DLBC (TCGA) y Análisis de SPN/CD43 con limma-voom

## Descripción
En este documento se analiza expresión diferencial a partir de datos públicos de RNA-seq del consorcio TCGA, accesibles vía recount3, comparando pacientes alive vs dead. 
Además, se explora el comportamiento del gen SPN (que codifica CD43/sialoforina), dada la evidencia de que CD43 puede contribuir a fenotipos pro-tumorales en distintos contextos (adhesión, crecimiento y señalización) y, de manera relevante para DLBCL, su expresión inmunohistoquímica se ha asociado con peores desenlaces clínicos.

## Estructura del Repositorio

📁 R/: Scripts de R numerados según el flujo de trabajo (01_recount3_explore_filter.R, 02_DE_limma.R, 03_CD43_dist.R).

📁 plots/: Gráficos generados (Heatmaps, Volcano plots, análisis de supervivencia de SPN, etc.).

📁 results/: Objetos de RData con los resultados crudos y filtrados (datos_DLBC_filtrados.RData, DE_limma_results.RData).

📁 docs/: Archivos de resumen de resultados(Rmd), para su visualizacion en una pagina web. Puede ser consultada en el apartado de descripcion del proyecto o aqui: https://pdanielpm.github.io/Proyecto_Modulo_Bioinfo/

## Requisitos (Dependencias)
  library(recount3)
  library(SummarizedExperiment)
  library(Matrix)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(ggpubr)
  library(matrixStats)
  library(pheatmap)
  library(survival)
  library(survminer)

## Instrucciones de Uso / Reproducibilidad

1. Ejecutar el script 01 para descargar y preprocesar los datos.

2. Correr el script 02 para el modelo lineal y la extracción de genes diferencialmente expresados.

3. Finalizar con el 03 para el análisis específico de CD43.
