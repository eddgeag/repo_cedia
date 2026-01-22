#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tools)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Uso: Rscript post_pipeline.R <MUESTRA>\nEj: Rscript post_pipeline.R DX046-25")
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
  "bqsr",
  paste0(muestra, ".recal.bam")
)

if (!file.exists(bam_file))
  stop("BAM final no existe: ", bam_file)

# --- Recursos fijos
genes_horizon <- "./genes_ACMG_HORIZON.csv"
base_lab      <- "./bd.rds"
hpo_file      <- "./hpo_file.txt"

bed_kit <- "./MGI_Exome_Capture_V5.hg38.sorted.merged.bed"
fasta   <- file.path(NAS_ROOT, "datos_exomas/datos_gatk/hg38/hg38.fa")

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

source("process_vcf_function.R")

vcf_process(
  vcf_file      = vcf_file,
  genes_horizon = genes_horizon,
  base_lab      = base_lab,
  hpo_file      = hpo_file
)

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
