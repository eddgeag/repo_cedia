suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggpubr)
})

# ========= PARÁMETROS =========
coverage_file <- "DX042-25_coverage.hist.txt"
output_dir    <- "coverage_and_stats"
muestra       <- "DX042-25"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ========= 1) LEER ARCHIVO =========
cov_raw <- fread(coverage_file, header = FALSE, sep = "\t", fill = TRUE)

cov_all <- cov_raw[V1 == "all", .(
  depth          = as.numeric(V2),
  bases_at_depth = as.numeric(V3),
  total_bases    = as.numeric(V4),
  fraction       = as.numeric(V5)
)]

if (nrow(cov_all) == 0)
  stop("No hay filas 'all' en el archivo de cobertura")

cov_all <- cov_all[order(depth)]

total_bases <- cov_all$total_bases[1]

# ========= 2) MÉTRICAS GLOBALES (incluye 0X) =========
mean_cov_global <- sum(cov_all$depth * cov_all$bases_at_depth) / total_bases

cum_bases_global <- cumsum(cov_all$bases_at_depth) / total_bases
median_cov_global <- cov_all$depth[which(cum_bases_global >= 0.5)[1]]

sd_cov_global <- sqrt(
  sum(cov_all$bases_at_depth * (cov_all$depth - mean_cov_global)^2) / total_bases
)

# ========= 3) MÉTRICAS CLÍNICAS (>0X) =========
cov_pos <- cov_all[cov_all$depth > 0, ]

total_bases_pos <- sum(cov_pos$bases_at_depth)

mean_cov_clin <- sum(cov_pos$depth * cov_pos$bases_at_depth) / total_bases_pos

cum_bases_pos <- cumsum(cov_pos$bases_at_depth) / total_bases_pos
median_cov_clin <- cov_pos$depth[which(cum_bases_pos >= 0.5)[1]]

sd_cov_clin <- sqrt(
  sum(cov_pos$bases_at_depth * (cov_pos$depth - mean_cov_clin)^2) / total_bases_pos
)

# ========= 4) MÉTRICAS ≥20X / ≥30X =========
pct_20x <- 100 * sum(cov_all$bases_at_depth[cov_all$depth >= 20]) / total_bases
pct_30x <- 100 * sum(cov_all$bases_at_depth[cov_all$depth >= 30]) / total_bases

# ========= 5) REDONDEO =========
mean_cov_global  <- round(mean_cov_global, 2)
median_cov_global<- round(median_cov_global, 2)
sd_cov_global    <- round(sd_cov_global, 2)

mean_cov_clin    <- round(mean_cov_clin, 2)
median_cov_clin  <- round(median_cov_clin, 2)
sd_cov_clin      <- round(sd_cov_clin, 2)

pct_20x <- round(pct_20x, 2)
pct_30x <- round(pct_30x, 2)

# ========= 6) DATOS PARA GRÁFICOS =========
plot_df <- data.frame(
  depth = cov_all$depth,
  pct_ge_depth = 100 * (1 - cumsum(cov_all$bases_at_depth) / total_bases),
  pct_at_depth = 100 * cov_all$bases_at_depth / total_bases
)

plot_df <- plot_df[plot_df$depth <= 200, ]

# ========= 7) GRÁFICOS =========
p1 <- ggplot(plot_df, aes(depth, pct_ge_depth)) +
  geom_line(color = "steelblue", linewidth = 1.2) +
  labs(
    x = "Profundidad de cobertura (X)",
    y = "% de la región ≥ profundidad",
    title = paste0(
      muestra,
      " | Media clínica: ", mean_cov_clin, "X",
      " | Mediana: ", median_cov_clin, "X",
      " | SD: ", sd_cov_clin, "X",
      "\n≥20X: ", pct_20x, "%  ≥30X: ", pct_30x, "%"
    )
  ) +
  theme_classic()

p2 <- ggplot(plot_df, aes(depth, pct_at_depth)) +
  geom_col(fill = "steelblue") +
  labs(
    x = "Profundidad de cobertura (X)",
    y = "% de la región"
  ) +
  theme_classic()

p_final <- ggarrange(p1, p2, ncol = 2)

# ========= 8) GUARDAR =========
ggsave(
  filename = file.path(output_dir, "coverage_depth_distribution.jpeg"),
  plot     = p_final,
  width    = 12,
  height   = 5
)

stats <- data.frame(
  metric = c(
    "mean_global", "median_global", "sd_global",
    "mean_clinical", "median_clinical", "sd_clinical",
    "pct_20x", "pct_30x"
  ),
  value = c(
    mean_cov_global, median_cov_global, sd_cov_global,
    mean_cov_clin, median_cov_clin, sd_cov_clin,
    pct_20x, pct_30x
  )
)

write.csv(
  stats,
  file.path(output_dir, "coverage_basic_stats.csv"),
  row.names = FALSE
)

# ========= 9) LOG =========
message("✔ Cobertura computada correctamente")
message("  Media global:    ", mean_cov_global, "X")
message("  Media clínica:   ", mean_cov_clin, "X")
message("  ≥20X:            ", pct_20x, "%")
message("  ≥30X:            ", pct_30x, "%")
