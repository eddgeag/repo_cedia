generate_exome_report <- function(
    muestra,
    hist_file,
    stats_file,
    exome_file,
    metrics_dir = "./metrics"
) {
  
  suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggpubr)
    library(officer)
    library(flextable)
  })
  
  dir.create(metrics_dir, showWarnings = FALSE, recursive = TRUE)
  
  out_fig  <- file.path(metrics_dir, paste0(muestra, "_coverage_plots.jpeg"))
  out_docx <- file.path(metrics_dir, paste0(muestra, "_reporte_cobertura_exoma.docx"))
  
  # ==================================================
  # 1. Métricas clínicas de cobertura
  # ==================================================
  stats <- fread(stats_file, header = FALSE)
  setnames(stats, c("metric", "value"))
  
  mean_cov <- as.numeric(stats[metric == "Mean_depth", value])
  pct20    <- as.numeric(stats[metric == "Pct_ge_20X", value])
  pct30    <- as.numeric(stats[metric == "Pct_ge_30X", value])
  
  # ==================================================
  # 2. Histograma de cobertura
  # ==================================================
  df <- fread(hist_file, col.names = c("depth", "bases"))
  df <- df[order(depth)]
  
  total_bases <- sum(df$bases)
  
  df[, pct_ge := 100 * (1 - cumsum(bases) / total_bases)]
  df[, pct_at := 100 * bases / total_bases]
  
  df_plot <- df[depth <= 200]
  
  # ==================================================
  # 3. Gráficas
  # ==================================================
  title_global <- paste0(
    muestra,
    " | Media: ", round(mean_cov,2), "X",
    " | ≥20X: ", round(pct20,2), "%",
    " | ≥30X: ", round(pct30,2), "%"
  )
  
  p1 <- ggplot(df_plot, aes(depth, pct_ge)) +
    geom_line(color = "steelblue", linewidth = 1.2) +
    geom_vline(xintercept = c(20, 30), linetype = "dashed", color = "grey50") +
    labs(
      x = "Profundidad de cobertura (X)",
      y = "% de bases ≥ X"
    ) +
    theme_classic(base_size = 10)
  
  p2 <- ggplot(df_plot, aes(depth, pct_at)) +
    geom_col(fill = "steelblue") +
    labs(
      x = "Profundidad de cobertura (X)",
      y = "% de bases"
    ) +
    theme_classic(base_size = 13)
  
  p_final <- ggarrange(p1, p2, ncol = 2)
  p_final <- annotate_figure(
    p_final,
    top = text_grob(title_global, face = "bold", size = 12)
  )
  
  ggsave(
    filename = out_fig,
    plot     = p_final,
    width    = 12,
    height   = 5,
    dpi      = 300
  )
  
  # ==================================================
  # 4. Estadísticas del exoma curado
  # ==================================================
  exome <- fread(exome_file, na.strings = c("", ".", "NA"))
  
  stopifnot("DP_INFO" %in% names(exome))
  stopifnot("GENES_AFFECTED" %in% names(exome))
  
  exome <- exome[GENES_AFFECTED != "" & !is.na(GENES_AFFECTED)]
  
  mean_dp   <- mean(exome$DP_INFO, na.rm = TRUE)
  median_dp <- median(exome$DP_INFO, na.rm = TRUE)
  
  genes_dt <- exome[, .(mean_dp_gene = mean(DP_INFO, na.rm = TRUE)), by = GENES_AFFECTED]
  
  genes_totales <- nrow(genes_dt)
  genes_20      <- genes_dt[mean_dp_gene >= 20, .N]
  genes_por_20  <- 100 * genes_20 / genes_totales
  
  genes_mayor_media <- genes_dt[mean_dp_gene >= mean_dp, .N]
  genes_por_media   <- 100 * genes_mayor_media / genes_totales
  
  variantes_totales <- nrow(exome)
  variantes_20      <- exome[DP_INFO >= 20, .N]
  variantes_20_x    <- 100 * variantes_20 / variantes_totales
  
  variantes_media   <- exome[DP_INFO >= mean_dp, .N]
  variantes_media_x <- 100 * variantes_media / variantes_totales
  
  exome_stats <- data.table(
    Descripción = c(
      "Cobertura media (variantes)",
      "Cobertura mediana (variantes)",
      "Genes totales",
      "Genes ≥ 20X",
      "Genes % ≥ 20X",
      "Genes ≥ media",
      "Genes % ≥ media",
      "Variantes totales",
      "Variantes ≥ 20X",
      "Variantes % ≥ 20X",
      "Variantes ≥ media",
      "Variantes % ≥ media"
    ),
    Valor = round(c(
      mean_dp,
      median_dp,
      genes_totales,
      genes_20,
      genes_por_20,
      genes_mayor_media,
      genes_por_media,
      variantes_totales,
      variantes_20,
      variantes_20_x,
      variantes_media,
      variantes_media_x
    ), 2)
  )
  
  ft_stats <- autofit(flextable(exome_stats))
  
  # ==================================================
  # 5. Textos
  # ==================================================
  texto_clinico <- paste0(
    "Esta muestra se ha estudiado mediante secuenciación masiva en paralelo del exoma completo. ",
    "Se analizaron ", genes_totales, " genes. ",
    "La sensibilidad y especificidad del método son superiores al 98% para variantes SNV e INDELs menores a 20 pb. ",
    "El ", round(genes_por_20,2), "% de los genes presentó una cobertura media ≥20X. ",
    "Del total de ", variantes_totales, " variantes identificadas, ",
    variantes_20, " presentaron una profundidad ≥20X, lo que corresponde al ",
    round(variantes_20_x,2), "% del total, garantizando una adecuada confiabilidad diagnóstica."
  )
  
  texto_metodologia <- paste(
    "La muestra fue procesada mediante un pipeline bioinformático estandarizado para análisis de exoma completo. ",
    "Se realizó control de calidad de lecturas (FastQC), alineamiento al genoma de referencia hg38 mediante BWA-MEM, ",
    "marcado de duplicados y recalibración de calidad de bases (BQSR) utilizando GATK. ",
    "El llamado de variantes se realizó con HaplotypeCaller y GenotypeGVCFs, seguido de filtrado duro. ",
    "Las variantes fueron anotadas con SnpEff, ClinVar, dbNSFP y recursos poblacionales. ",
    "La curación bioinformática incluyó análisis de profundidad, consistencia alélica, zigosidad, ",
    "herencia compuesta y clasificación ACMG/AMP automatizada y revisable clínicamente."
  )
  
  texto_graficos <- paste(
    "El gráfico de la izquierda muestra la proporción acumulada de bases con una profundidad de cobertura ",
    "igual o superior a un umbral dado, mientras que el gráfico de la derecha representa la distribución ",
    "exacta de profundidades de cobertura. Estas visualizaciones permiten evaluar la calidad global del ",
    "exoma y la robustez del análisis clínico."
  )
  
  
  # ==================================================
  # 6. Word
  # ==================================================
  doc <- read_docx() |>
    body_add_par("INFORME DE COBERTURA – EXOMA COMPLETO", style = "heading 1") |>
    body_add_par(paste("Muestra:", muestra)) |>
    body_add_par("Metodología", style = "heading 2") |>
    body_add_par(texto_metodologia) |>
    body_add_par("Resultados de cobertura", style = "heading 2") |>
    body_add_par(texto_graficos) |>
    body_add_img(src = out_fig, width = 6.5, height = 3) |>
    body_add_par("Resumen estadístico del exoma", style = "heading 2") |>
    body_add_flextable(ft_stats) |>
    body_add_par("Interpretación clínica", style = "heading 2") |>
    body_add_par(texto_clinico)
  
  print(doc, target = out_docx)
  
  message("✔ Informe Word generado: ", out_docx)
}
