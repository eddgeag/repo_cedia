#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tools)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Uso: Rscript main.R <MUESTRA>\nEj: Rscript main.R DX046-25")
}
muestra <- args[1]
t_start <- Sys.time()
## ===============================
## 1. Paths base
## ===============================
NAS_ROOT <- path.expand("~/NAS_NGS")

yy <- sub(".*-(\\d{2})$", "\\1", muestra)
year_full <- paste0("20", yy)

muestra_dir <- file.path(NAS_ROOT, year_full, muestra)
output_dir  <- file.path(muestra_dir, "output_dir")

if (!dir.exists(output_dir))
  stop("output_dir no existe: ", output_dir)

## ===============================
## 2. Whitelist (reglas duras)
## ===============================
keep_files <- c(
  file.path(output_dir, "mapping_output", paste0(muestra, ".sorted.rg.mark_dup_bqsr.bam")),
  file.path(output_dir, "mapping_output", paste0(muestra, ".sorted.rg.mark_dup_bqsr.bai")),
  file.path(output_dir, "anotation", paste0(muestra, "_05_dbnsfp_reduced.vcf.gz")),
  file.path(output_dir, "anotation", paste0(muestra, "_05_dbnsfp_reduced.vcf.gz.tbi")),
  file.path(output_dir, "variantCalling", paste0(muestra, ".hardfiltered.all.vcf.gz"))
)

keep_dirs <- file.path(
  output_dir,
  c("QC", "post_process_results", "metrics", "bam_metrics")
)
keep_dirs_structure <- file.path(
  output_dir,
  c("mapping_output", "anotation", "variantCalling")
)
dir.create(keep_dirs_structure, showWarnings = FALSE, recursive = TRUE)

## ===============================
## 3. Función de limpieza
## ===============================
cleanup_output_dir <- function(output_dir, keep_files, keep_dirs, keep_dirs_structure) {
  
  norm <- function(x) normalizePath(x, winslash = "/", mustWork = FALSE)
  
  keep_files_n <- norm(keep_files)
  keep_dirs_n  <- norm(keep_dirs)
  keep_dirs_structure_n <- norm(keep_dirs_structure)
  
  all_paths <- list.files(
    output_dir,
    recursive   = TRUE,
    full.names = TRUE,
    all.files  = TRUE,
    include.dirs = TRUE,
    no.. = TRUE
  )
  
  all_n <- norm(all_paths)
  
  is_keep_file <- all_n %in% keep_files_n
  
  is_keep_dir <- vapply(
    all_n,
    function(p) {
      # conservar carpetas completas
      any(startsWith(p, paste0(keep_dirs_n, "/")) | p %in% keep_dirs_n) ||
        # conservar SOLO la estructura (no el contenido)
        p %in% keep_dirs_structure_n
    },
    logical(1)
  )
  
  to_delete <- all_paths[!(is_keep_file | is_keep_dir)]
  
  unlink(to_delete[file.exists(to_delete) & !dir.exists(to_delete)], force = TRUE)
  unlink(to_delete[dir.exists(to_delete)], recursive = TRUE, force = TRUE)
}

## ===============================
## 4. Comprobación de artefactos críticos
## ===============================
missing_core <- keep_files[!file.exists(keep_files)]

if (length(missing_core) > 0L) {
  
  message(">>> Faltan artefactos críticos del exoma:")
  message(paste0(" - ", missing_core, collapse = "\n"))
  
  message(">>> Limpieza previa (dejando solo QC / metrics / bam_metrics si existen)")
  cleanup_output_dir(output_dir, keep_files, keep_dirs)
  
  message(">>> Ejecutando pipeline_lab.R (pipeline completo GATK)")
  
  pipeline_script <- "pipeline_lab.R"
  if (!file.exists(pipeline_script))
    stop("No existe pipeline_lab.R")
  
  res <- system2(
    Sys.which("Rscript"),
    c(pipeline_script, muestra),
    stdout = TRUE,
    stderr = TRUE
  )
  message(paste(res, collapse = "\n"))
  
  ## Revalidación dura
  missing_core <- keep_files[!file.exists(keep_files)]
  if (length(missing_core) > 0L) {
    stop(
      "Pipeline GATK terminó pero siguen faltando artefactos:\n",
      paste0(" - ", missing_core, collapse = "\n")
    )
  }
}

## ===============================
## 5. Curación bioinformática
## ===============================
vcf_file <- keep_files[grepl("_05_dbnsfp_reduced.vcf.gz$", keep_files)]
bam_file <- keep_files[grepl("mark_dup_bqsr.bam$", keep_files)]

genes_horizon <- "./genes_ACMG_HORIZON.csv"
base_lab      <- "./bd.rds"
hpo_file      <- "./genes_to_phenotype.txt"

post_dir <- file.path(output_dir, "post_process_results")

expected_outputs <- file.path(
  post_dir,
  c(
    "file_ready_analysis_optimized.csv",
    "file_ready_analysis_optimized_UNICAS.csv",
    "file_ready_analysis_optimized_UNICAS_ACMG_HORIZON.csv"
  )
)

run_vcf_process <- !all(file.exists(expected_outputs))

source("process_vcf_function.R")

if (run_vcf_process) {
  message(">>> Ejecutando vcf_process")
  vcf_process(
    vcf_file      = vcf_file,
    genes_horizon = genes_horizon,
    base_lab      = base_lab,
    hpo_file      = hpo_file
  )
} else {
  message(">>> Curación ya existe, se omite")
}

## ===============================
## 6. Cobertura
## ===============================
bed_kit <- "./MGI_Exome_Capture_V5.hg38.sorted.merged.bed"
fasta   <- file.path(NAS_ROOT, "datos_exomas/datos_gatk/hg38/hg38.fasta")

coverage_out <- file.path(
  output_dir, "metrics",
  paste0(muestra, "_coverage.hist.txt")
)

dir.create(dirname(coverage_out), showWarnings = FALSE, recursive = TRUE)

coverage_script <- "./computar_coverage.sh"
if (!file.exists(coverage_script))
  stop("No existe computar_coverage.sh")

system2(
  coverage_script,
  c(bam_file, bed_kit, fasta, coverage_out),
  stdout = TRUE,
  stderr = TRUE
)

## ===============================
## 7. Reporte
## ===============================
metrics_dir <- file.path(output_dir, "metrics")

hist_file  <- file.path(metrics_dir, paste0(muestra, "_coverage.depth.hist.tsv"))
stats_file <- file.path(metrics_dir, paste0(muestra, "_coverage.depth.stats.txt"))
exome_file <- file.path(post_dir, "file_ready_analysis_optimized.csv")

if (!file.exists(hist_file))  stop("No existe histograma: ", hist_file)
if (!file.exists(stats_file)) stop("No existe stats: ", stats_file)
if (!file.exists(exome_file)) stop("No existe exoma curado")

source("generate_exome_report.R")

generate_exome_report(
  muestra     = muestra,
  hist_file   = hist_file,
  stats_file  = stats_file,
  exome_file  = exome_file,
  metrics_dir = metrics_dir
)

## ===============================
## 8. Limpieza final
## ===============================
message(">>> Limpieza final del output_dir")
cleanup_output_dir(
  output_dir,
  keep_files,
  keep_dirs,
  keep_dirs_structure
)

message(">>> main.R finalizado correctamente")
t_end <- Sys.time()
elapsed <- difftime(t_end, t_start, units = "mins")

message(sprintf(
  ">>> Tiempo total de ejecución: %.2f minutos",
  as.numeric(elapsed)
))
