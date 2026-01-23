#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tools)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Uso: Rscript main.R <MUESTRA>\nEj: Rscript main.R DX046-25")
}

muestra <- args[1]

## ===============================
## 1. Derivar paths base (MISMA lógica que main.R)
## ===============================
NAS_ROOT <- path.expand("~/NAS_NGS")

yy <- sub(".*-(\\d{2})$", "\\1", muestra)
year_full <- paste0("20", yy)

muestra_dir <- file.path(NAS_ROOT, year_full, muestra)
output_dir  <- file.path(muestra_dir, "output_dir")

if (!dir.exists(output_dir))
  stop("output_dir no existe: ", output_dir)

## ===============================
## 2. Derivar inputs automáticamente
## ===============================

# --- VCF final (ajusta patrón si cambia)
vcf_file <- file.path(
  output_dir,
  "anotation",
  paste0(muestra, "_05_dbnsfp_reduced.vcf.gz")
)

if (!file.exists(vcf_file))
  stop("VCF final no existe: ", vcf_file)

# --- BAM final post-BQSR
bam_file <- file.path(
  output_dir,
  "mapping_output",
  paste0(muestra, ".sorted.rg.mark_dup_bqsr.bam")
)

if (!file.exists(bam_file))
  stop("BAM final no existe: ", bam_file)

# --- Recursos fijos
genes_horizon <- "./genes_ACMG_HORIZON.csv"
base_lab      <- "./bd.rds"
hpo_file      <- "./genes_to_phenotype.txt"

bed_kit <- "./MGI_Exome_Capture_V5.hg38.sorted.merged.bed"
fasta   <- file.path(NAS_ROOT, "datos_exomas/datos_gatk/hg38/hg38.fasta")

coverage_out <- file.path(
  output_dir,
  "metrics",
  paste0(muestra, "_coverage.hist.txt")
)

dir.create(dirname(coverage_out), showWarnings = FALSE, recursive = TRUE)

## ===============================
## 3. Curación bioinformática (R)
## ===============================
message(">>> Ejecutando curación bioinformática")
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
  message(">>> Ejecutando curación bioinformática (vcf_process)")
  
  vcf_process(
    vcf_file      = vcf_file,
    genes_horizon = genes_horizon,
    base_lab      = base_lab,
    hpo_file      = hpo_file
  )
  
} else {
  message(">>> Curación bioinformática ya existe, se omite vcf_process")
}

## ===============================
## 4. Cobertura on-target (shell)
## ===============================
message(">>> Ejecutando cómputo de cobertura")

coverage_script <- "./computar_coverage.sh"

if (!file.exists(coverage_script))
  stop("No existe computar_coverage.sh")

cmd <- coverage_script
args <- c(
  bam_file,
  bed_kit,
  fasta,
  coverage_out
)

res <- system2(
  command = cmd,
  args    = args,
  stdout  = TRUE,
  stderr  = TRUE
)

message(res)

message(">>> Pipeline post-proceso FINALIZADO")
metrics_dir <- file.path(output_dir, "metrics")
post_dir    <- file.path(output_dir, "post_process_results")

## ===============================
## 2. Derivar inputs automáticamente
## ===============================
hist_file <- file.path(
  metrics_dir,
  paste0(muestra, "_coverage.depth.hist.tsv")
)

stats_file <- file.path(
  metrics_dir,
  paste0(muestra, "_coverage.depth.stats.txt")
)

exome_file <- file.path(
  post_dir,
  "file_ready_analysis_optimized.csv"
)

## Chequeos mínimos
if (!file.exists(hist_file))
  stop("No existe histograma de cobertura: ", hist_file)

if (!file.exists(stats_file))
  stop("No existe stats de cobertura: ", stats_file)

if (!file.exists(exome_file))
  stop("No existe exoma curado: ", exome_file)

## ===============================
## 3. Generar reporte
## ===============================
message(">>> Generando informe de cobertura y exoma")

source("generate_exome_report.R")

generate_exome_report(
  muestra     = muestra,
  hist_file   = hist_file,
  stats_file  = stats_file,
  exome_file  = exome_file,
  metrics_dir = metrics_dir
)

message(">>> Informe generado correctamente")