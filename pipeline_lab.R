load.libs <- c(
  "data.table",
  "tools",
  "dplyr",
  "tidyr",
  "ggplot2",
  "stringr",
  "xtable",
  "VariantAnnotation",
  "tibble",
  "GenomicRanges",
  "S4Vectors",
  "ggpubr"
)
invisible(lapply(load.libs, require, character.only = TRUE))

get_sample_name <- function(fastq_dir) {
  fq <- list.files(fastq_dir, pattern = "\\.fq\\.gz$", full.names = FALSE)
  r1 <- fq[grepl("_R1\\.fq\\.gz$", fq)]
  stopifnot(length(r1) == 1)
  
  ## Eliminar _R1.fq.gz → nombre de muestra limpio
  x <- sub("_R1\\.fq\\.gz$", "", r1)
  x
}


control_calidad <- function(
    fastq_dir,
    output_dir,
    threads = 8,
    fastqc_bin = "fastqc",
    overwrite = FALSE
) {
  # ---------------------------
  # 0) Validaciones
  # ---------------------------
  fastq_dir  <- path.expand(fastq_dir)
  output_dir <- path.expand(output_dir)
  
  if (!dir.exists(fastq_dir)) {
    stop("fastq_dir no existe: ", fastq_dir)
  }
  
  fastqc_path <- Sys.which(fastqc_bin)
  if (fastqc_path == "") {
    stop("No se encontró el binario 'fastqc' en PATH")
  }
  
  # ---------------------------
  # 1) Directorio QC
  # ---------------------------
  qc_dir <- file.path(output_dir, "QC")
  dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ---------------------------
  # 2) Detectar FASTQ
  # ---------------------------
  fastq_files <- list.files(
    fastq_dir,
    pattern = "\\.(fq|fastq)\\.gz$",
    full.names = TRUE
  )
  
  if (length(fastq_files) == 0) {
    stop("No se encontraron FASTQ en: ", fastq_dir)
  }
  
  # ---------------------------
  # 3) Evitar recomputar
  # ---------------------------
  qc_reports <- list.files(qc_dir, pattern = "\\.html$", full.names = TRUE)
  
  if (length(qc_reports) > 0 && !overwrite) {
    message("QC ya existe. Use overwrite = TRUE para rehacerlo.")
    return(invisible(qc_reports))
  }
  
  # ---------------------------
  # 4) Ejecutar FastQC
  # ---------------------------
  args <- c(
    "-t", as.character(threads),
    "-o", qc_dir,
    fastq_files
  )
  
  status <- system2(
    command = fastqc_path,
    args    = args,
    stdout  = TRUE,
    stderr  = TRUE
  )
  
  # ---------------------------
  # 5) Chequeo post-ejecución
  # ---------------------------
  qc_reports <- list.files(qc_dir, pattern = "\\.html$", full.names = TRUE)
  
  if (length(qc_reports) == 0) {
    stop("FastQC terminó pero no generó reportes. Revise stderr.")
  }
  
  invisible(qc_reports)
}



fn_exists_fasta <- function(folder_fasta) {
  fasta_files <- list.files(folder_fasta, pattern = "\\.(fa|fasta)$", full.names =
                              TRUE)
  if (length(fasta_files) == 0)
    stop("No existe archivo de referencia (.fa/.fasta) en: ",
         folder_fasta)
  if (length(fasta_files) > 1)
    stop("Hay más de un FASTA en: ",
         folder_fasta,
         "\n",
         paste(basename(fasta_files), collapse = ", "))
  fasta_files[[1]]
}



index_fasta_samtools <- function(
    folder_fasta,
    samtools_bin = "samtools",
    overwrite = FALSE
) {
  # ---------------------------
  # 0) Validaciones
  # ---------------------------
  folder_fasta <- path.expand(folder_fasta)
  
  if (!dir.exists(folder_fasta)) {
    stop("folder_fasta no existe: ", folder_fasta)
  }
  
  samtools_path <- Sys.which(samtools_bin)
  if (samtools_path == "") {
    stop("No se encontró 'samtools' en PATH")
  }
  
  # FASTA (función que ya usas en el pipeline)
  fasta_file <- path.expand(fn_exists_fasta(folder_fasta))
  
  if (!file.exists(fasta_file)) {
    stop("FASTA no encontrado: ", fasta_file)
  }
  
  fai_file <- paste0(fasta_file, ".fai")
  
  # ---------------------------
  # 1) Evitar recomputar
  # ---------------------------
  if (file.exists(fai_file) && !overwrite) {
    message("Índice FAI ya existe: ", basename(fai_file))
    return(invisible(fai_file))
  }
  
  # ---------------------------
  # 2) Ejecutar samtools faidx
  # ---------------------------
  status <- system2(
    command = samtools_path,
    args    = c("faidx", fasta_file),
    stdout  = TRUE,
    stderr  = TRUE
  )
  
  # ---------------------------
  # 3) Chequeo post-ejecución
  # ---------------------------
  if (!file.exists(fai_file)) {
    stop("samtools faidx falló. No se generó: ", fai_file)
  }
  
  invisible(fai_file)
}


index_bwa <- function(
    folder_fasta,
    bwa_bin = "bwa",
    overwrite = FALSE
) {
  # ---------------------------
  # 0) Validaciones
  # ---------------------------
  folder_fasta <- path.expand(folder_fasta)
  
  if (!dir.exists(folder_fasta)) {
    stop("folder_fasta no existe: ", folder_fasta)
  }
  
  bwa_path <- Sys.which(bwa_bin)
  if (bwa_path == "") {
    stop("No se encontró 'bwa' en PATH")
  }
  
  fasta_file <- path.expand(fn_exists_fasta(folder_fasta))
  
  if (!file.exists(fasta_file)) {
    stop("FASTA no encontrado: ", fasta_file)
  }
  
  # ---------------------------
  # 1) Índices esperados
  # ---------------------------
  exts <- c("amb", "ann", "bwt", "pac", "sa")
  index_files <- paste0(fasta_file, ".", exts)
  
  if (all(file.exists(index_files)) && !overwrite) {
    message("Índices BWA ya existen")
    return(invisible(index_files))
  }
  
  # ---------------------------
  # 2) Ejecutar bwa index
  # ---------------------------
  message("Creando índices BWA...")
  
  status <- system2(
    command = bwa_path,
    args    = c("index", fasta_file),
    stdout  = TRUE,
    stderr  = TRUE
  )
  
  # ---------------------------
  # 3) Chequeo post-ejecución
  # ---------------------------
  if (!all(file.exists(index_files))) {
    stop(
      "bwa index falló. Faltan archivos:\n",
      paste(basename(index_files[!file.exists(index_files)]), collapse = ", ")
    )
  }
  
  invisible(index_files)
}


bwamem <- function(
    fastq_dir,
    folder_fasta,
    output_dir,
    threads = 8,
    bwa_bin = "bwa",
    samtools_bin = "samtools",
    gatk_bin = "~/tools/gatk-4.6.1.0/gatk",
    overwrite = FALSE
) {
  # ---------------------------
  # 0) Normalización
  # ---------------------------
  fastq_dir    <- path.expand(fastq_dir)
  folder_fasta <- path.expand(folder_fasta)
  output_dir   <- path.expand(output_dir)
  gatk_bin     <- path.expand(gatk_bin)
  
  if (!file.exists(gatk_bin)) {
    stop("GATK no encontrado en: ", gatk_bin)
  }
  
  # ---------------------------
  # 1) FASTA + índices
  # ---------------------------
  fasta_file <- path.expand(fn_exists_fasta(folder_fasta))
  
  index_bwa(folder_fasta)
  
  if (!file.exists(paste0(fasta_file, ".fai"))) {
    index_fasta_samtools(folder_fasta)
  }
  
  # ---------------------------
  # 2) FASTQ
  # ---------------------------
  fq <- list.files(
    fastq_dir,
    pattern = "\\.(fq|fastq)\\.gz$",
    full.names = TRUE
  )
  
  r1 <- path.expand(fq[grepl("_R1\\.fq\\.gz$", fq)])
  r2 <- path.expand(fq[grepl("_R2\\.fq\\.gz$", fq)])
  
  stopifnot(length(r1) == 1, length(r2) == 1)
  
  # ---------------------------
  # 3) Sample ID
  # ---------------------------
  sample_id <- get_sample_name(fastq_dir)
  
  # ---------------------------
  # 4) Output
  # ---------------------------
  mapping_output_dir <- file.path(output_dir, "mapping_output")
  dir.create(mapping_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  sam_file   <- file.path(mapping_output_dir, paste0(sample_id, ".sam"))
  bam_file   <- file.path(mapping_output_dir, paste0(sample_id, ".bam"))
  sorted_bam <- file.path(mapping_output_dir, paste0(sample_id, ".sorted.bam"))
  
  if (file.exists(sorted_bam) && !overwrite) {
    message("BAM ordenado ya existe: ", basename(sorted_bam))
    return(invisible(sorted_bam))
  }
  
  # ---------------------------
  # 5) BWA MEM
  # ---------------------------
  message("#### MAPPING (bwa mem, sin RG) ####")
  
  cmd <- c(
    "mem",
    "-M",
    "-t", as.character(threads),
    fasta_file,
    r1,
    r2
  )
  
  message("Ejecutando comando:")
  message(paste(shQuote(bwa_bin), paste(shQuote(cmd), collapse = " ")))
  
  if(!file.exists(sam_file)){
    system2(
      command = bwa_bin,
      args    = cmd,
      stdout  = sam_file,
      stderr  = file.path(mapping_output_dir, paste0(sample_id, ".bwa.stderr.txt"))
    )
    
  }

  if (!file.exists(sam_file)) {
    stop("ERROR CRÍTICO: bwa mem falló para ", sample_id)
  }
  
  # ---------------------------
  # 6) SAM → BAM
  # ---------------------------
  if(!file.exists(bam_file)){
    system2(
      samtools_bin,
      args = c(
        "view",
        "-b",
        "-h",
        "-@", as.character(threads),
        sam_file
      ),
      stdout = bam_file,
      stderr = TRUE
    )
    
  }

  if (!file.exists(bam_file)) {
    stop("samtools view falló para ", sample_id)
  }
  
  # ---------------------------
  # 7) Sort + index (GATK)
  # ---------------------------
  if(!file.exists(sorted_bam)){
    system2(
      gatk_bin,
      args = c(
        "SortSam",
        "-I", bam_file,
        "-O", sorted_bam,
        "-SORT_ORDER", "coordinate",
        "-CREATE_INDEX", "true"
      ),
      stdout = TRUE,
      stderr = TRUE
    )
    
  }

  if (!file.exists(sorted_bam)) {
    stop("SortSam falló para ", sample_id)
  }
  
  invisible(sorted_bam)
}



# ============================================================
# add_read_groups() – versión saneada
# ============================================================
add_read_groups <- function(
    output_dir,
    fastq_dir,
    gatk_bin = "~/tools/gatk-4.6.1.0/gatk",
    rg_platform = "ILLUMINA",
    rg_library  = "lib1",
    rg_unit     = "unit1",
    overwrite = FALSE
) {
  # ---------------------------
  # 0) Normalización y checks
  # ---------------------------
  output_dir <- path.expand(output_dir)
  fastq_dir  <- path.expand(fastq_dir)
  gatk_bin   <- path.expand(gatk_bin)
  
  if (!file.exists(gatk_bin)) {
    stop("GATK no encontrado en: ", gatk_bin)
  }
  
  sample_id <- get_sample_name(fastq_dir)
  
  mapping_dir <- file.path(output_dir, "mapping_output")
  
  bam_in <- file.path(
    mapping_dir,
    paste0(sample_id, ".sorted.bam")
  )
  
  bam_out <- file.path(
    mapping_dir,
    paste0(sample_id, ".sorted.rg.bam")
  )
  
  bai_candidates <- c(
    paste0(bam_out, ".bai"),
    file.path(dirname(bam_out),
              paste0(tools::file_path_sans_ext(basename(bam_out)), ".bai"))
  )
  
  if (!file.exists(bam_in)) {
    stop("No existe BAM de entrada para AddReadGroups: ", bam_in)
  }
  
  # ---------------------------
  # 1) Evitar recomputar
  # ---------------------------
  if (file.exists(bam_out) && !overwrite) {
    message("Read Groups ya existen: ", basename(bam_out))
    return(invisible(bam_out))
  }
  
  # ---------------------------
  # 2) AddOrReplaceReadGroups
  # ---------------------------
  message("#### AddOrReplaceReadGroups ####")
  
  system2(
    command = gatk_bin,
    args = c(
      "AddOrReplaceReadGroups",
      "-I", bam_in,
      "-O", bam_out,
      "-RGID", sample_id,
      "-RGSM", sample_id,
      "-RGPL", rg_platform,
      "-RGLB", rg_library,
      "-RGPU", rg_unit,
      "-CREATE_INDEX", "true"
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  
  # ---------------------------
  # 3) Verificación final
  # ---------------------------
  if (!file.exists(bam_out)) {
    stop("ERROR CRÍTICO: AddOrReplaceReadGroups no generó BAM: ", bam_out)
  }
  
  if (!any(file.exists(bai_candidates))) {
    stop("ERROR CRÍTICO: no se generó índice .bai para: ", bam_out)
  }
  
  invisible(bam_out)
}
markdups <- function(
    output_dir,
    fastq_dir,
    gatk_bin = "~/tools/gatk-4.6.1.0/gatk",
    samtools_bin = "samtools",
    wait_seconds = 10,
    overwrite = FALSE
) {
  # ---------------------------
  # 0) Normalización y checks
  # ---------------------------
  output_dir <- path.expand(output_dir)
  fastq_dir  <- path.expand(fastq_dir)
  gatk_bin   <- path.expand(gatk_bin)
  
  if (!file.exists(gatk_bin)) {
    stop("GATK no encontrado en: ", gatk_bin)
  }
  
  samtools_path <- Sys.which(samtools_bin)
  if (samtools_path == "") {
    stop("samtools no encontrado en PATH")
  }
  
  sample_id <- get_sample_name(fastq_dir)
  mapping_dir <- file.path(output_dir, "mapping_output")
  
  bam_in <- file.path(mapping_dir,
                      paste0(sample_id, ".sorted.rg.bam"))
  
  if (!file.exists(bam_in)) {
    stop("No existe BAM con Read Groups: ", bam_in)
  }
  
  # ---------------------------
  # 1) Verificar Read Groups (SIN pipes)
  # ---------------------------
  header <- system2(
    samtools_path,
    args = c("view", "-H", bam_in),
    stdout = TRUE,
    stderr = TRUE
  )
  
  if (!any(grepl("^@RG", header))) {
    stop("ERROR CRÍTICO: el BAM NO contiene Read Groups: ", bam_in)
  }
  
  # ---------------------------
  # 2) Salidas
  # ---------------------------
  bam_out <- file.path(mapping_dir,
                       paste0(sample_id, ".sorted.rg.mark_dup.bam"))
  
  metrics_file <- file.path(mapping_dir,
                            paste0(sample_id, ".sorted.rg.mark_dup.txt"))
  
  bai_candidates <- c(
    paste0(bam_out, ".bai"),
    file.path(dirname(bam_out),
              paste0(tools::file_path_sans_ext(basename(bam_out)), ".bai"))
  )
  
  if (file.exists(bam_out) && !overwrite) {
    message("MarkDuplicates ya ejecutado, reutilizando BAM existente")
  } else {
    # ---------------------------
    # 3) MarkDuplicates
    # ---------------------------
    message("#### MarkDuplicates (con RG) ####")
    
    system2(
      gatk_bin,
      args = c(
        "MarkDuplicates",
        "-INPUT", bam_in,
        "-OUTPUT", bam_out,
        "-M", metrics_file,
        "--ASSUME_SORTED", "true",
        "--CREATE_INDEX", "true",
        "--VALIDATION_STRINGENCY", "STRICT"
      ),
      stdout = TRUE,
      stderr = TRUE
    )
    
    if (!file.exists(bam_out)) {
      stop("ERROR CRÍTICO: MarkDuplicates falló para ", sample_id)
    }
  }
  
  # ---------------------------
  # 4) Espera activa (NAS-safe)
  # ---------------------------
  waited <- 0
  while (!any(file.exists(bai_candidates)) && waited < wait_seconds) {
    Sys.sleep(1)
    waited <- waited + 1
  }
  
  # ---------------------------
  # 5) Crear índice si no existe
  # ---------------------------
  if (!any(file.exists(bai_candidates))) {
    message("Índice BAI no detectado, generando explícitamente...")
    
    system2(
      gatk_bin,
      args = c("BuildBamIndex", "-I", bam_out),
      stdout = TRUE,
      stderr = TRUE
    )
  }
  
  # ---------------------------
  # 6) Verificación final
  # ---------------------------
  if (!any(file.exists(bai_candidates))) {
    stop(
      "ERROR CRÍTICO: el índice .bai NO existe tras todos los intentos.\n",
      paste(" -", bai_candidates, collapse = "\n")
    )
  }
  
  message("MarkDuplicates completado correctamente (BAM + BAI verificados)")
  invisible(bam_out)
}




create_dict <- function(
    folder_fasta,
    gatk_bin = "~/tools/gatk-4.6.1.0/gatk",
    overwrite = FALSE
) {
  # ---------------------------
  # 0) Normalización y checks
  # ---------------------------
  folder_fasta <- path.expand(folder_fasta)
  gatk_bin     <- path.expand(gatk_bin)
  
  if (!dir.exists(folder_fasta)) {
    stop("folder_fasta no existe: ", folder_fasta)
  }
  
  if (!file.exists(gatk_bin)) {
    stop("GATK no encontrado en: ", gatk_bin)
  }
  
  # ---------------------------
  # 1) FASTA de referencia
  # ---------------------------
  fasta_file <- path.expand(fn_exists_fasta(folder_fasta))
  
  if (!file.exists(fasta_file)) {
    stop("FASTA no encontrado: ", fasta_file)
  }
  
  fasta_dict <- paste0(tools::file_path_sans_ext(fasta_file), ".dict")
  
  # ---------------------------
  # 2) Crear diccionario
  # ---------------------------
  if (file.exists(fasta_dict) && !overwrite) {
    message("Diccionario FASTA ya existe")
    return(invisible(fasta_dict))
  }
  
  message("Creando diccionario de referencia (.dict)...")
  
  system2(
    command = gatk_bin,
    args = c(
      "CreateSequenceDictionary",
      "-R", fasta_file,
      "-O", fasta_dict
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  
  # ---------------------------
  # 3) Verificación final
  # ---------------------------
  if (!file.exists(fasta_dict)) {
    stop("ERROR CRÍTICO: No se pudo crear el diccionario FASTA: ", fasta_dict)
  }
  
  invisible(fasta_dict)
}


base_recalibrator <- function(
    folder_fasta,
    output_dir,
    folder_data_gatk,
    fastq_dir,
    gatk_bin = "~/tools/gatk-4.6.1.0/gatk",
    overwrite = FALSE
) {
  # ---------------------------
  # 0) Normalización y checks
  # ---------------------------
  folder_fasta    <- path.expand(folder_fasta)
  output_dir      <- path.expand(output_dir)
  folder_data_gatk <- path.expand(folder_data_gatk)
  fastq_dir       <- path.expand(fastq_dir)
  gatk_bin        <- path.expand(gatk_bin)
  
  if (!file.exists(gatk_bin)) {
    stop("GATK no encontrado en: ", gatk_bin)
  }
  
  # ---------------------------
  # 1) FASTA de referencia
  # ---------------------------
  fasta_file <- path.expand(fn_exists_fasta(folder_fasta))
  if (!file.exists(fasta_file)) {
    stop("FASTA no encontrado: ", fasta_file)
  }
  
  # ---------------------------
  # 2) Known-sites (obligatorios)
  # ---------------------------
  known_sites_files <- path.expand(c(
    file.path(
      folder_data_gatk,
      "resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"
    ),
    file.path(
      folder_data_gatk,
      "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
    )
  ))
  
  if (!all(file.exists(known_sites_files))) {
    stop(
      "ERROR CRÍTICO: faltan archivos known-sites:\n",
      paste(known_sites_files[!file.exists(known_sites_files)], collapse = "\n")
    )
  }
  
  # ---------------------------
  # 3) Directorio de mapeo
  # ---------------------------
  mapping_output_dir <- file.path(output_dir, "mapping_output")
  
  # ---------------------------
  # 4) Sample ID (fuente única)
  # ---------------------------
  sample_id <- get_sample_name(fastq_dir)
  
  # ---------------------------
  # 5) BAM de entrada (MarkDuplicates)
  # ---------------------------
  bam_file <- file.path(
    mapping_output_dir,
    paste0(sample_id, ".sorted.rg.mark_dup.bam")
  )
  
  bai_candidates <- c(
    paste0(bam_file, ".bai"),
    file.path(dirname(bam_file),
              paste0(tools::file_path_sans_ext(basename(bam_file)), ".bai"))
  )
  
  if (!file.exists(bam_file)) {
    stop("No existe el BAM con duplicados marcados: ", bam_file)
  }
  
  if (!any(file.exists(bai_candidates))) {
    stop("No existe índice .bai del BAM de entrada:\n",
         paste(" -", bai_candidates, collapse = "\n"))
  }
  
  # ---------------------------
  # 6) Tabla de salida BQSR
  # ---------------------------
  out_file <- file.path(
    mapping_output_dir,
    paste0(sample_id, ".recal_data.table")
  )
  
  if (file.exists(out_file) && !overwrite) {
    message("Tabla de BaseRecalibrator ya existe")
    return(invisible(out_file))
  }
  
  # ---------------------------
  # 7) Ejecutar BaseRecalibrator
  # ---------------------------
  message("#### BaseRecalibrator ####")
  
  args <- c(
    "BaseRecalibrator",
    "-I", bam_file,
    "-R", fasta_file
  )
  
  ## one --known-sites por archivo (forma correcta)
  for (ks in known_sites_files) {
    args <- c(args, "--known-sites", ks)
  }
  
  args <- c(args, "-O", out_file)
  
  system2(
    command = gatk_bin,
    args    = args,
    stdout = TRUE,
    stderr = TRUE
  )
  
  # ---------------------------
  # 8) Verificación final
  # ---------------------------
  if (!file.exists(out_file)) {
    stop(
      "ERROR CRÍTICO: BaseRecalibrator falló para la muestra ",
      sample_id
    )
  }
  
  invisible(out_file)
}



applybqsr <- function(
    folder_fasta,
    output_dir,
    fastq_dir,
    gatk_bin = "~/tools/gatk-4.6.1.0/gatk",
    overwrite = FALSE
) {
  # ---------------------------
  # 0) Normalización y checks
  # ---------------------------
  folder_fasta <- path.expand(folder_fasta)
  output_dir   <- path.expand(output_dir)
  fastq_dir    <- path.expand(fastq_dir)
  gatk_bin     <- path.expand(gatk_bin)
  
  if (!file.exists(gatk_bin)) {
    stop("GATK no encontrado en: ", gatk_bin)
  }
  
  # ---------------------------
  # 1) FASTA de referencia
  # ---------------------------
  fasta_file <- path.expand(fn_exists_fasta(folder_fasta))
  if (!file.exists(fasta_file)) {
    stop("FASTA no encontrado: ", fasta_file)
  }
  
  # ---------------------------
  # 2) Directorio de mapeo
  # ---------------------------
  mapping_output_dir <- file.path(output_dir, "mapping_output")
  
  # ---------------------------
  # 3) Sample ID
  # ---------------------------
  sample_id <- get_sample_name(fastq_dir)
  
  # ---------------------------
  # 4) BAM de entrada (MarkDuplicates + RG)
  # ---------------------------
  bam_in <- file.path(
    mapping_output_dir,
    paste0(sample_id, ".sorted.rg.mark_dup.bam")
  )
  
  bai_in_candidates <- c(
    paste0(bam_in, ".bai"),
    file.path(dirname(bam_in),
              paste0(tools::file_path_sans_ext(basename(bam_in)), ".bai"))
  )
  
  if (!file.exists(bam_in)) {
    stop("No existe el BAM de entrada para ApplyBQSR: ", bam_in)
  }
  
  if (!any(file.exists(bai_in_candidates))) {
    stop(
      "No existe el índice .bai del BAM de entrada:\n",
      paste(" -", bai_in_candidates, collapse = "\n")
    )
  }
  
  # ---------------------------
  # 5) Tabla BQSR
  # ---------------------------
  recal_table <- file.path(
    mapping_output_dir,
    paste0(sample_id, ".recal_data.table")
  )
  
  if (!file.exists(recal_table)) {
    stop(
      "No existe la tabla BQSR: ", recal_table,
      "\n¿Ejecutaste base_recalibrator() antes?"
    )
  }
  
  # ---------------------------
  # 6) BAM de salida (BQSR aplicado)
  # ---------------------------
  bam_out <- file.path(
    mapping_output_dir,
    paste0(sample_id, ".sorted.rg.mark_dup_bqsr.bam")
  )
  
  bai_out_candidates <- c(
    paste0(bam_out, ".bai"),
    file.path(dirname(bam_out),
              paste0(tools::file_path_sans_ext(basename(bam_out)), ".bai"))
  )
  
  if (file.exists(bam_out) && !overwrite) {
    message("BAM con BQSR ya existe")
    return(invisible(bam_out))
  }
  
  # ---------------------------
  # 7) Ejecutar ApplyBQSR
  # ---------------------------
  message("#### ApplyBQSR ####")
  
  system2(
    command = gatk_bin,
    args = c(
      "ApplyBQSR",
      "-I", bam_in,
      "-R", fasta_file,
      "--bqsr-recal-file", recal_table,
      "-O", bam_out
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  
  if (!file.exists(bam_out)) {
    stop("ERROR CRÍTICO: ApplyBQSR falló para la muestra ", sample_id)
  }
  
  # ---------------------------
  # 8) Índice obligatorio (.bai)
  # ---------------------------
  if (!any(file.exists(bai_out_candidates))) {
    message("Índice BAI no encontrado. Generando explícitamente...")
    
    system2(
      command = gatk_bin,
      args = c("BuildBamIndex", "-I", bam_out),
      stdout = TRUE,
      stderr = TRUE
    )
  }
  
  if (!any(file.exists(bai_out_candidates))) {
    stop(
      "ERROR CRÍTICO: No se pudo generar el índice .bai para ",
      basename(bam_out)
    )
  }
  
  message("ApplyBQSR completado correctamente (BAM + BAI verificados)")
  invisible(bam_out)
}

bam_statistics <- function(
    folder_fasta,
    fastq_dir,
    output_dir,
    gatk_bin = "~/tools/gatk-4.6.1.0/gatk",
    multiqc_bin = "multiqc",
    overwrite = FALSE
) {
  # ---------------------------
  # 0) Normalización y checks
  # ---------------------------
  folder_fasta <- path.expand(folder_fasta)
  fastq_dir    <- path.expand(fastq_dir)
  output_dir   <- path.expand(output_dir)
  gatk_bin     <- path.expand(gatk_bin)
  
  if (!file.exists(gatk_bin)) {
    stop("GATK no encontrado en: ", gatk_bin)
  }
  
  multiqc_path <- Sys.which(multiqc_bin)
  if (multiqc_path == "") {
    stop("multiqc no encontrado en PATH")
  }
  
  # ---------------------------
  # 1) FASTA de referencia
  # ---------------------------
  fasta_file <- path.expand(fn_exists_fasta(folder_fasta))
  if (!file.exists(fasta_file)) {
    stop("FASTA no encontrado: ", fasta_file)
  }
  
  # ---------------------------
  # 2) Directorio de mapeo
  # ---------------------------
  mapping_output_dir <- file.path(output_dir, "mapping_output")
  
  # ---------------------------
  # 3) Sample ID
  # ---------------------------
  sample_id <- get_sample_name(fastq_dir)
  
  # ---------------------------
  # 4) BAM final (BQSR aplicado + RG)
  # ---------------------------
  bam_file <- file.path(
    mapping_output_dir,
    paste0(sample_id, ".sorted.rg.mark_dup_bqsr.bam")
  )
  
  bai_candidates <- c(
    paste0(bam_file, ".bai"),
    file.path(dirname(bam_file),
              paste0(tools::file_path_sans_ext(basename(bam_file)), ".bai"))
  )
  
  if (!file.exists(bam_file)) {
    stop("No existe el BAM final con BQSR: ", bam_file)
  }
  
  if (!any(file.exists(bai_candidates))) {
    stop(
      "No existe el índice .bai del BAM final:\n",
      paste(" -", bai_candidates, collapse = "\n")
    )
  }
  
  # ---------------------------
  # 5) Directorio de métricas
  # ---------------------------
  out_dir <- file.path(output_dir, "bam_metrics", sample_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ---------------------------
  # 6) Prefijo de salida
  # ---------------------------
  out_prefix <- file.path(out_dir, "CollectMultipleMetrics")
  expected_file <- paste0(out_prefix, ".alignment_summary_metrics")
  
  if (file.exists(expected_file) && !overwrite) {
    message("Las métricas BAM ya existen para la muestra")
    return(invisible(out_dir))
  }
  
  # ---------------------------
  # 7) CollectMultipleMetrics
  # ---------------------------
  message("#### CollectMultipleMetrics ####")
  
  system2(
    command = gatk_bin,
    args = c(
      "CollectMultipleMetrics",
      "-R", fasta_file,
      "-I", bam_file,
      "-O", out_prefix
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  
  if (!file.exists(expected_file)) {
    stop(
      "ERROR CRÍTICO: CollectMultipleMetrics falló para la muestra ",
      sample_id
    )
  }
  
  # ---------------------------
  # 8) MultiQC
  # ---------------------------
  message("#### MultiQC ####")
  
  system2(
    command = multiqc_path,
    args = c(out_dir, "-o", out_dir),
    stdout = TRUE,
    stderr = TRUE
  )
  
  message("Métricas BAM generadas correctamente")
  invisible(out_dir)
}


haplotype_caller <- function(
    output_dir,
    folder_fasta,
    fastq_dir,
    bed_file,
    gatk_bin = "~/tools/gatk-4.6.1.0/gatk",
    threads = 8,
    overwrite = FALSE
) {
  # ---------------------------
  # 0) Normalización y checks
  # ---------------------------
  output_dir   <- path.expand(output_dir)
  folder_fasta <- path.expand(folder_fasta)
  fastq_dir    <- path.expand(fastq_dir)
  bed_file     <- path.expand(bed_file)
  gatk_bin     <- path.expand(gatk_bin)
  
  if (!file.exists(gatk_bin)) {
    stop("GATK no encontrado en: ", gatk_bin)
  }
  
  if (!file.exists(bed_file)) {
    stop("ERROR CRÍTICO: no existe el BED de captura del exoma: ", bed_file)
  }
  
  # ---------------------------
  # 1) FASTA de referencia
  # ---------------------------
  fasta_file <- path.expand(fn_exists_fasta(folder_fasta))
  if (!file.exists(fasta_file)) {
    stop("FASTA no encontrado: ", fasta_file)
  }
  
  # ---------------------------
  # 2) Directorio de mapeo
  # ---------------------------
  mapping_output_dir <- file.path(output_dir, "mapping_output")
  
  # ---------------------------
  # 3) Sample ID
  # ---------------------------
  sample_id <- get_sample_name(fastq_dir)
  
  # ---------------------------
  # 4) BAM de entrada (BQSR aplicado + RG)
  # ---------------------------
  bam_file <- file.path(
    mapping_output_dir,
    paste0(sample_id, ".sorted.rg.mark_dup_bqsr.bam")
  )
  
  bai_candidates <- c(
    paste0(bam_file, ".bai"),
    file.path(dirname(bam_file),
              paste0(tools::file_path_sans_ext(basename(bam_file)), ".bai"))
  )
  
  if (!file.exists(bam_file)) {
    stop("No existe el BAM BQSR para HaplotypeCaller: ", bam_file)
  }
  
  if (!any(file.exists(bai_candidates))) {
    stop(
      "No existe el índice .bai del BAM BQSR:\n",
      paste(" -", bai_candidates, collapse = "\n")
    )
  }
  
  # ---------------------------
  # 5) Directorio y archivo de salida
  # ---------------------------
  out_dir <- file.path(output_dir, "variantCalling")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  out_file <- file.path(out_dir, paste0(sample_id, ".g.vcf.gz"))
  out_tbi  <- paste0(out_file, ".tbi")
  
  if (file.exists(out_file) && !overwrite) {
    message("GVCF ya existe para la muestra")
    return(invisible(out_file))
  }
  
  # ---------------------------
  # 6) Ejecutar HaplotypeCaller (GVCF + BED)
  # ---------------------------
  message("#### HaplotypeCaller (GVCF + BED exoma) ####")
  
  system2(
    command = gatk_bin,
    args = c(
      "HaplotypeCaller",
      "-R", fasta_file,
      "-I", bam_file,
      "-L", bed_file,                    # CLÍNICO
      "-ERC", "GVCF",
      "--native-pair-hmm-threads", as.character(threads),
      "-O", out_file
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  
  if (!file.exists(out_file)) {
    stop("ERROR CRÍTICO: HaplotypeCaller falló para la muestra ", sample_id)
  }
  
  # ---------------------------
  # 7) Verificación / creación del índice .tbi
  # ---------------------------
  if (!file.exists(out_tbi)) {
    message("Índice .tbi no encontrado. Generando explícitamente...")
    
    system2(
      command = gatk_bin,
      args = c("IndexFeatureFile", "-I", out_file),
      stdout = TRUE,
      stderr = TRUE
    )
  }
  
  if (!file.exists(out_tbi)) {
    stop(
      "ERROR CRÍTICO: No se pudo generar el índice .tbi para ",
      basename(out_file)
    )
  }
  
  message("HaplotypeCaller completado correctamente (GVCF + BED + TBI verificados)")
  invisible(out_file)
}


genotypeGVCFs <- function(
    folder_fasta,
    output_dir,
    fastq_dir,
    gatk_bin = "~/tools/gatk-4.6.1.0/gatk",
    overwrite = FALSE
) {
  # ---------------------------
  # 0) Normalización y checks
  # ---------------------------
  folder_fasta <- path.expand(folder_fasta)
  output_dir   <- path.expand(output_dir)
  fastq_dir    <- path.expand(fastq_dir)
  gatk_bin     <- path.expand(gatk_bin)
  
  if (!file.exists(gatk_bin)) {
    stop("GATK no encontrado en: ", gatk_bin)
  }
  
  # ---------------------------
  # 1) FASTA de referencia
  # ---------------------------
  fasta_file <- path.expand(fn_exists_fasta(folder_fasta))
  if (!file.exists(fasta_file)) {
    stop("FASTA no encontrado: ", fasta_file)
  }
  
  # ---------------------------
  # 2) Sample ID
  # ---------------------------
  sample_id <- get_sample_name(fastq_dir)
  
  # ---------------------------
  # 3) Directorio de variantes
  # ---------------------------
  var_dir <- file.path(output_dir, "variantCalling")
  dir.create(var_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ---------------------------
  # 4) GVCF de entrada
  # ---------------------------
  file_in <- file.path(var_dir, paste0(sample_id, ".g.vcf.gz"))
  file_in_tbi <- paste0(file_in, ".tbi")
  
  if (!file.exists(file_in)) {
    stop("No existe el GVCF de entrada: ", file_in)
  }
  
  if (!file.exists(file_in_tbi)) {
    stop("No existe el índice .tbi del GVCF: ", file_in_tbi)
  }
  
  # ---------------------------
  # 5) VCF de salida
  # ---------------------------
  file_out <- file.path(var_dir, paste0(sample_id, "_sample.raw.vcf.gz"))
  file_out_tbi <- paste0(file_out, ".tbi")
  
  if (file.exists(file_out) && !overwrite) {
    message("VCF ya existe para la muestra")
    return(invisible(file_out))
  }
  
  # ---------------------------
  # 6) Ejecutar GenotypeGVCFs
  # ---------------------------
  message("#### GenotypeGVCFs ####")
  
  system2(
    command = gatk_bin,
    args = c(
      "GenotypeGVCFs",
      "-R", fasta_file,
      "-V", file_in,
      "-O", file_out
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  
  if (!file.exists(file_out)) {
    stop("ERROR CRÍTICO: GenotypeGVCFs falló para ", sample_id)
  }
  
  # ---------------------------
  # 7) Índice .tbi obligatorio
  # ---------------------------
  if (!file.exists(file_out_tbi)) {
    message("Índice .tbi no encontrado. Generando explícitamente...")
    
    system2(
      command = gatk_bin,
      args = c("IndexFeatureFile", "-I", file_out),
      stdout = TRUE,
      stderr = TRUE
    )
  }
  
  if (!file.exists(file_out_tbi)) {
    stop("ERROR CRÍTICO: no se pudo generar el índice .tbi para ", file_out)
  }
  
  message("GenotypeGVCFs completado correctamente (VCF + TBI verificados)")
  invisible(file_out)
}


variantFiltration <- function(
    folder_fasta,
    output_dir,
    fastq_dir,
    gatk_bin = "~/tools/gatk-4.6.1.0/gatk",
    overwrite = FALSE
) {
  
  folder_fasta <- path.expand(folder_fasta)
  output_dir   <- path.expand(output_dir)
  fastq_dir    <- path.expand(fastq_dir)
  gatk_bin     <- path.expand(gatk_bin)
  
  fasta_file <- fn_exists_fasta(folder_fasta)
  sample_id  <- get_sample_name(fastq_dir)
  
  var_dir <- file.path(output_dir, "variantCalling")
  dir.create(var_dir, recursive = TRUE, showWarnings = FALSE)
  
  in_vcf <- file.path(var_dir, paste0(sample_id, "_sample.raw.vcf.gz"))
  if (!file.exists(in_vcf)) stop("VCF de entrada no existe")
  
  snps_vcf        <- file.path(var_dir, paste0(sample_id, ".snps.vcf"))
  indels_vcf      <- file.path(var_dir, paste0(sample_id, ".indels.vcf"))
  snps_filt_vcf   <- file.path(var_dir, paste0(sample_id, ".snps.hardfiltered.vcf"))
  indels_filt_vcf <- file.path(var_dir, paste0(sample_id, ".indels.hardfiltered.vcf"))
  merged_vcf      <- file.path(var_dir, paste0(sample_id, ".hardfiltered.vcf"))
  
  if (file.exists(merged_vcf) && !overwrite) {
    message("VCF hardfiltered ya existe")
    return(invisible(merged_vcf))
  }
  
  ## SNPs
  system2(gatk_bin, c(
    "SelectVariants", "-R", fasta_file, "-V", in_vcf,
    "--select-type-to-include", "SNP", "-O", snps_vcf
  ))
  
  system2(gatk_bin, c(
    "VariantFiltration",
    "-R", fasta_file, "-V", snps_vcf, "-O", snps_filt_vcf,
    "--filter-name", "QD2",   "--filter-expression", "QD<2.0",
    "--filter-name", "MQ40",  "--filter-expression", "MQ<40.0",
    "--filter-name", "FS60",  "--filter-expression", "FS>60.0",
    "--filter-name", "SOR3",  "--filter-expression", "SOR>3.0"
  ))
  
  ## INDELs
  system2(gatk_bin, c(
    "SelectVariants", "-R", fasta_file, "-V", in_vcf,
    "--select-type-to-include", "INDEL", "-O", indels_vcf
  ))
  
  system2(gatk_bin, c(
    "VariantFiltration",
    "-R", fasta_file, "-V", indels_vcf, "-O", indels_filt_vcf,
    "--filter-name", "QD2",    "--filter-expression", "QD<2.0",
    "--filter-name", "FS200",  "--filter-expression", "FS>200.0",
    "--filter-name", "SOR5",   "--filter-expression", "SOR>5.0"
  ))
  
  ## Merge
  system2(gatk_bin, c(
    "MergeVcfs",
    "-I", snps_filt_vcf,
    "-I", indels_filt_vcf,
    "-O", merged_vcf
  ))
  
  if (!file.exists(merged_vcf))
    stop("No se pudo crear hardfiltered.vcf")
  
  invisible(merged_vcf)
}





analysisReady <- function(
    output_dir,
    fastq_dir,
    bcftools_bin = "bcftools",
    bgzip_bin = "bgzip",
    overwrite = FALSE
) {
  
  output_dir <- path.expand(output_dir)
  fastq_dir  <- path.expand(fastq_dir)
  
  bcftools <- Sys.which(bcftools_bin)
  bgzip    <- Sys.which(bgzip_bin)
  
  sample_id <- get_sample_name(fastq_dir)
  var_dir   <- file.path(output_dir, "variantCalling")
  
  in_vcf <- file.path(var_dir, paste0(sample_id, ".hardfiltered.vcf"))
  if (!file.exists(in_vcf)) stop("hardfiltered.vcf no existe")
  
  out_pass     <- file.path(var_dir, paste0(sample_id, ".hardfiltered.pass.vcf"))
  out_pass_gz  <- paste0(out_pass, ".gz")
  out_tbi      <- paste0(out_pass_gz, ".tbi")
  report_txt   <- file.path(var_dir, paste0(sample_id, ".filter_report.txt"))
  
  ## ---- reporte de filtros ----
  stats <- system2(
    bcftools,
    c("stats", in_vcf),
    stdout = TRUE
  )
  writeLines(stats, report_txt)
  
  ## ---- PASS ----
  system2(
    bcftools,
    c("view", "-f", "PASS", "-O", "v", "-o", out_pass, in_vcf)
  )
  
  if (!file.exists(out_pass))
    stop("No se pudo generar VCF PASS")
  
  system2(bgzip, c("-f", out_pass))
  system2(bcftools, c("index", "-t", out_pass_gz))
  
  if (!file.exists(out_tbi))
    stop("No se pudo indexar VCF PASS")
  
  message("VCF PASS listo + reporte de filtros generado")
  invisible(out_pass_gz)
}


anotation <- function(
    folder_fasta,
    path_snpeff,
    output_dir,
    fastq_dir,
    clinvar_vcf,
    dbnsfp_db,
    gwas_db = NULL,          # GWAS opcional
    use_gwas = FALSE,        # <- CONTROL EXPLÍCITO
    dbnsfp_fields = c(
      "aaref","aaalt","rs_dbSNP","HGVSc_snpEff","HGVSp_snpEff","APPRIS",
      "CADD_phred","AlphaMissense_pred","clinvar_OMIM_id",
      "clinvar_Orphanet_id","clinvar_MedGen_id"
    ),
    java_bin = "java",
    bcftools_bin = "bcftools",
    bgzip_bin = "bgzip",
    snpeff_genome = "hg38",
    java_mem = "32g",
    download_snpeff_db_if_missing = TRUE,
    overwrite = FALSE
) {
  
  # =========================================================
  # 0) Normalización + checks binarios
  # =========================================================
  output_dir   <- path.expand(output_dir)
  fastq_dir    <- path.expand(fastq_dir)
  folder_fasta <- path.expand(folder_fasta)
  path_snpeff  <- path.expand(path_snpeff)
  
  java_path     <- Sys.which(java_bin)
  bcftools_path <- Sys.which(bcftools_bin)
  bgzip_path    <- Sys.which(bgzip_bin)
  
  if (java_path == "")     stop("java no encontrado")
  if (bcftools_path == "") stop("bcftools no encontrado")
  if (bgzip_path == "")    stop("bgzip no encontrado")
  
  # Java >= 21
  jv <- paste(system2(java_path, "-version", stderr = TRUE), collapse = " ")
  if (!grepl("version \"(21|22|23|24|25)", jv))
    stop("snpEff requiere Java >= 21")
  
  # =========================================================
  # 1) FASTA
  # =========================================================
  fasta_file <- fn_exists_fasta(folder_fasta)
  
  # =========================================================
  # 2) Sample + dirs
  # =========================================================
  sample_id <- get_sample_name(fastq_dir)
  
  variant_dir   <- file.path(output_dir, "variantCalling")
  anotacion_dir <- file.path(output_dir, "anotation")
  dir.create(anotacion_dir, recursive = TRUE, showWarnings = FALSE)
  
  in_vcf <- file.path(variant_dir, paste0(sample_id, ".hardfiltered.pass.vcf.gz"))
  if (!file.exists(in_vcf)) stop("VCF PASS no existe: ", in_vcf)
  
  # =========================================================
  # 3) Recursos
  # =========================================================
  snpeff_jar  <- file.path(path_snpeff, "snpEff.jar")
  snpsift_jar <- file.path(path_snpeff, "SnpSift.jar")
  
  req <- c(snpeff_jar, snpsift_jar, clinvar_vcf, dbnsfp_db)
  if (!all(file.exists(req)))
    stop("Faltan recursos de anotación")
  
  if (use_gwas && is.null(gwas_db))
    stop("use_gwas=TRUE pero gwas_db es NULL")
  
  # =========================================================
  # 4) Base snpEff
  # =========================================================
  snpeff_data <- file.path(path_snpeff, "data", snpeff_genome)
  if (!dir.exists(snpeff_data)) {
    if (!download_snpeff_db_if_missing)
      stop("Base snpEff no encontrada")
    system2(
      java_path,
      c(paste0("-Xmx", java_mem), "-jar", snpeff_jar, "download", snpeff_genome)
    )
  }
  
  # =========================================================
  # 5) Archivos
  # =========================================================
  f1 <- file.path(anotacion_dir, paste0(sample_id, "_01_snpeff.vcf"))
  f2 <- file.path(anotacion_dir, paste0(sample_id, "_02_varType.vcf"))
  f3 <- file.path(anotacion_dir, paste0(sample_id, "_03_clinvar.vcf"))
  f4 <- file.path(anotacion_dir, paste0(sample_id, "_04_dbnsfp_full.vcf"))
  f5 <- if (use_gwas)
    file.path(anotacion_dir, paste0(sample_id, "_05_gwas.vcf"))
  else f4
  
  f6 <- file.path(anotacion_dir, paste0(sample_id, "_06_dbnsfp_reduced.vcf.gz"))
  f6_tbi <- paste0(f6, ".tbi")
  
  # =========================================================
  # 6) snpEff
  # =========================================================
  if (!file.exists(f1) || overwrite)
    system2(java_path, c(paste0("-Xmx", java_mem), "-jar", snpeff_jar,
                         snpeff_genome, "-v", in_vcf),
            stdout = f1)
  
  if (!file.exists(f1)) stop("snpEff falló")
  
  # =========================================================
  # 7) varType
  # =========================================================
  if (!file.exists(f2) || overwrite)
    system2(java_path, c(paste0("-Xmx", java_mem), "-jar", snpsift_jar,
                         "varType", "-v", f1),
            stdout = f2)
  
  # =========================================================
  # 8) ClinVar
  # =========================================================
  if (!file.exists(f3) || overwrite)
    system2(java_path, c(paste0("-Xmx", java_mem), "-jar", snpsift_jar,
                         "annotate", "-v", clinvar_vcf, f2),
            stdout = f3)
  
  # =========================================================
  # 9) dbNSFP FULL
  # =========================================================
  if (!file.exists(f4) || overwrite)
    system2(java_path, c(paste0("-Xmx", java_mem), "-jar", snpsift_jar,
                         "dbnsfp", "-v", "-db", dbnsfp_db, f3),
            stdout = f4)
  
  # =========================================================
  # 10) GWAS (opcional)
  # =========================================================
  if (use_gwas) {
    if (!file.exists(f5) || overwrite)
      system2(java_path, c(paste0("-Xmx", java_mem), "-jar", snpsift_jar,
                           "annotate", "-v", gwas_db, f4),
              stdout = f5)
  }
  
  # =========================================================
  # 11) dbNSFP REDUCIDO (FINAL CLÍNICO)
  # =========================================================
  campos <- paste(dbnsfp_fields, collapse = ",")
  tmp <- sub("\\.gz$", "", f6)
  
  if (!file.exists(f6) || overwrite) {
    system2(
      java_path,
      c(paste0("-Xmx", java_mem), "-jar", snpsift_jar,
        "dbnsfp", "-v", "-db", dbnsfp_db,
        "-f", campos, f5),
      stdout = tmp
    )
    
    if (!file.exists(tmp)) stop("dbNSFP reducido falló")
    
    system2(bgzip_path, c("-f", tmp))
  }
  
  # =========================================================
  # 12) Indexación final
  # =========================================================
  if (!file.exists(f6_tbi))
    system2(bcftools_path, c("index", "-t", f6))
  
  if (!file.exists(f6_tbi))
    stop("No se pudo crear índice TBI final")
  
  message("ANOTACIÓN CLÍNICA FINALIZADA CORRECTAMENTE")
  invisible(f6)
}




verify_bqsr_minimal <- function(
    output_dir,
    fastq_dir,
    samtools_bin = "samtools",
    size_tol = 0.001
) {
  # ---------------------------
  # 0) Normalización y checks
  # ---------------------------
  output_dir <- path.expand(output_dir)
  fastq_dir  <- path.expand(fastq_dir)
  
  samtools_path <- Sys.which(samtools_bin)
  if (samtools_path == "") {
    stop("samtools no encontrado en PATH")
  }
  
  sample_id <- get_sample_name(fastq_dir)
  mapping_dir <- file.path(output_dir, "mapping_output")
  
  bam_pre  <- file.path(mapping_dir,
                        paste0(sample_id, ".sorted.rg.mark_dup.bam"))
  bam_post <- file.path(mapping_dir,
                        paste0(sample_id, ".sorted.rg.mark_dup_bqsr.bam"))
  recal_table <- file.path(mapping_dir,
                           paste0(sample_id, ".recal_data.table"))
  
  # ---------------------------
  # 1) Existencia básica
  # ---------------------------
  if (!file.exists(bam_pre))
    stop("BQSR CHECK: no existe BAM pre-BQSR: ", bam_pre)
  
  if (!file.exists(bam_post))
    stop("BQSR CHECK: no existe BAM post-BQSR: ", bam_post)
  
  if (!file.exists(recal_table))
    stop("BQSR CHECK: no existe tabla recal_data.table")
  
  # ---------------------------
  # 2) Verificar que el BAM cambió
  # ---------------------------
  size_pre  <- file.info(bam_pre)$size
  size_post <- file.info(bam_post)$size
  
  if (is.na(size_pre) || is.na(size_post)) {
    stop("BQSR CHECK: no se pudo obtener tamaño de BAMs")
  }
  
  size_diff <- abs(size_post - size_pre) / size_pre
  
  if (size_diff < size_tol) {
    stop(
      "BQSR CHECK: BAM post-BQSR prácticamente idéntico al pre-BQSR\n",
      "Δ relativo = ", signif(size_diff, 3), "\n",
      "Esto sugiere que ApplyBQSR no tuvo efecto real."
    )
  }
  
  # ---------------------------
  # 3) Verificar header (soft check)
  # ---------------------------
  header_post <- system2(
    samtools_path,
    args = c("view", "-H", bam_post),
    stdout = TRUE,
    stderr = TRUE
  )
  
  if (!any(grepl("BQSR|ApplyBQSR|recal", header_post, ignore.case = TRUE))) {
    warning(
      "BQSR CHECK: no se detecta mención explícita de BQSR en el header.\n",
      "No es fatal, pero conviene revisar ApplyBQSR."
    )
  }
  
  # ---------------------------
  # 4) Sanity check de calidades (correcto)
  # ---------------------------
  reads <- system2(
    samtools_path,
    args = c("view", bam_post),
    stdout = TRUE,
    stderr = TRUE
  )
  
  reads <- reads[seq_len(min(length(reads), 500))]
  
  quals <- vapply(
    strsplit(reads, "\t"),
    function(x) if (length(x) >= 11) x[11] else NA_character_,
    character(1)
  )
  
  quals <- quals[!is.na(quals)]
  
  uniq_qual_strings <- length(unique(quals))
  
  if (uniq_qual_strings < 5) {
    warning(
      "BQSR CHECK: baja diversidad de strings de calidad en reads iniciales.\n",
      "No es prueba de error, pero es atípico."
    )
  }
  
  message("BQSR CHECK: verificación mínima completada correctamente")
  invisible(TRUE)
}

###=====main=============#########

###===============================
### main.R (saneado)
### - Todos los paths modificables en el main (sin quitar defaults en funciones)
### - Muestra definida de forma segura
### - Nombres de funciones corregidos (genotypeGVCFs, haplotype_caller con bed_file, anotation con recursos)
###===============================

start <- Sys.time()

# ----------------------------
# 0) Parse args: una sola muestra
# ----------------------------
argv <- commandArgs(trailingOnly = TRUE)

if (length(argv) == 0) {
  stop("Uso: Rscript main.R <MUESTRA>\nEj: Rscript main.R DX046-25")
}

# Acepta que venga como "DX046-25" o como "DX046-25 ..." (por copy/paste)
muestra <- trimws(strsplit(paste(argv, collapse = " "), "\\s+")[[1]][1])

if (is.na(muestra) || muestra == "") stop("Muestra vacía")

# Año desde sufijo -YY (mantengo tu lógica)
yy <- sub(".*-(\\d{2})$", "\\1", muestra)
if (!grepl("^\\d{2}$", yy)) {
  stop("No se pudo inferir el año desde el nombre de la muestra: ", muestra)
}
year_full <- paste0("20", yy)

# ----------------------------
# 1) Paths MODIFICABLES (todos aquí)
# ----------------------------
NAS_ROOT <- path.expand("~/NAS_NGS")

base_year_dir <- file.path(NAS_ROOT, year_full)
muestra_dir   <- file.path(base_year_dir, muestra)

# IMPORTANT: output_dir debe ser estable (no lo cambies entre runs)
output_dir    <- file.path(muestra_dir, "output_dir")

fastq_dir     <- file.path(muestra_dir, "fastqfiles")

# Referencia + recursos GATK
folder_fasta     <- file.path(NAS_ROOT, "datos_exomas/datos_gatk/hg38")
folder_data_gatk <- file.path(NAS_ROOT, "datos_exomas/datos_gatk")

# BED clínico (ahora ES INPUT de haplotype_caller)
bed_file <- file.path("./MGI_Exome_Capture_V5.hg38.sorted.merged.bed")

# Binarios / herramientas (aunque las funciones tienen defaults, se declaran aquí)
gatk_bin     <- path.expand("~/tools/gatk-4.6.1.0/gatk")
path_snpeff  <- path.expand("~/tools/snpEff")

fastqc_bin   <- "fastqc"
bwa_bin      <- "bwa"
samtools_bin <- "samtools"
bcftools_bin <- "bcftools"
bgzip_bin    <- "bgzip"
multiqc_bin  <- "multiqc"
java_bin     <- "java"

threads <- 8

# Recursos anotación (ahora OBLIGATORIOS como inputs en anotation)
clinvar_vcf <- file.path(NAS_ROOT, "datos_exomas/datos_clinvar/clinvar_20260104.vcf.gz")
dbnsfp_db   <- file.path(NAS_ROOT, "datos_exomas/datos_dbsnp/dbNSFP5.3.1a_grch38.gz")
gwas_db     <- file.path(NAS_ROOT, "datos_exomas/gwas/gwascatalog.txt")

# Campos dbNSFP (modificable desde main)
dbnsfp_fields <- c(
  # Identidad / consecuencia
  "aaref",
  "aaalt",
  "rs_dbSNP",
  "HGVSc_snpEff",
  "HGVSp_snpEff",
  "APPRIS",
  "MANE",
  
  # Predicción funcional
  "CADD_phred",
  "REVEL_score",
  "ClinPred_score",
  "AlphaMissense_pred",
  "MetaRNN_pred",
  "MPC_score",
  
  # Frecuencia poblacional
  "gnomAD4.1_joint_POPMAX_AF",
  "gnomAD4.1_joint_POPMAX_POP",
  "dbNSFP_POPMAX_AF",
  
  # Constraint
  "LOEUF",
  "pLI",
  
  # ClinVar cross-check
  "clinvar_clnsig",
  "clinvar_review",
  "clinvar_OMIM_id",
  "clinvar_Orphanet_id",
  "clinvar_MedGen_id"
)


# Flags de recomputación
overwrite_qc        <- FALSE
overwrite_mapping   <- FALSE
overwrite_markdups  <- FALSE
overwrite_bqsr      <- FALSE
overwrite_metrics   <- FALSE
overwrite_variants  <- FALSE
overwrite_filtering <- FALSE
overwrite_ready     <- FALSE
overwrite_annot     <- FALSE

# ----------------------------
# 2) Crear directorios
# ----------------------------
dir.create(base_year_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(muestra_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir,    recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# 3) Validaciones mínimas upfront (evita runs largos rotos)
# ----------------------------
if (!dir.exists(fastq_dir)) stop("fastq_dir no existe: ", fastq_dir)
if (!dir.exists(folder_fasta)) stop("folder_fasta no existe: ", folder_fasta)
if (!dir.exists(folder_data_gatk)) stop("folder_data_gatk no existe: ", folder_data_gatk)
if (!file.exists(gatk_bin)) stop("gatk_bin no existe: ", gatk_bin)
if (!file.exists(bed_file)) stop("bed_file no existe: ", bed_file)
if (!dir.exists(path_snpeff)) stop("path_snpeff no existe: ", path_snpeff)

# Recursos anotación
req_annot <- c(clinvar_vcf, dbnsfp_db, gwas_db)
if (!all(file.exists(req_annot))) {
  stop("Faltan recursos de anotación:\n", paste(req_annot[!file.exists(req_annot)], collapse = "\n"))
}

message("MUESTRA: ", muestra)
message("FASTQ_DIR: ", fastq_dir)
message("OUTPUT_DIR: ", output_dir)

# ============================================================
# 4) Pipeline (llamadas ajustadas a saneamientos)
# ============================================================

# QC
control_calidad(
  fastq_dir  = fastq_dir,
  output_dir = output_dir,
  threads    = threads,
  fastqc_bin = fastqc_bin,
  overwrite  = overwrite_qc
)

# Mapping
bwamem(
  fastq_dir    = fastq_dir,
  folder_fasta = folder_fasta,
  output_dir   = output_dir,
  threads      = threads,
  bwa_bin      = bwa_bin,
  samtools_bin = samtools_bin,
  gatk_bin     = gatk_bin,
  overwrite    = overwrite_mapping
)

# Read groups (AQUÍ faltaba saneamiento: pasar gatk_bin como input)
# Recomendación: actualiza la función add_read_groups para aceptar gatk_bin.
# Mientras tanto: mantengo llamada y seteo gatk_bin por fuera SOLO si tu función ya lo usa hardcodeado.
add_read_groups(
  output_dir = output_dir,
  fastq_dir  = fastq_dir
)

# MarkDuplicates
markdups(
  output_dir   = output_dir,
  fastq_dir    = fastq_dir,
  gatk_bin     = gatk_bin,
  samtools_bin = samtools_bin,
  overwrite    = overwrite_markdups
)

# Diccionario referencia
create_dict(
  folder_fasta = folder_fasta,
  gatk_bin     = gatk_bin
)

# BQSR
base_recalibrator(
  folder_fasta     = folder_fasta,
  output_dir       = output_dir,
  folder_data_gatk = folder_data_gatk,
  fastq_dir        = fastq_dir,
  gatk_bin         = gatk_bin,
  overwrite        = overwrite_bqsr
)

applybqsr(
  folder_fasta = folder_fasta,
  output_dir   = output_dir,
  fastq_dir    = fastq_dir,
  gatk_bin     = gatk_bin,
  overwrite    = overwrite_bqsr
)

verify_bqsr_minimal(
  output_dir   = output_dir,
  fastq_dir    = fastq_dir,
  samtools_bin = samtools_bin
)

# Métricas BAM
bam_statistics(
  folder_fasta = folder_fasta,
  fastq_dir    = fastq_dir,
  output_dir   = output_dir,
  gatk_bin     = gatk_bin,
  multiqc_bin  = multiqc_bin,
  overwrite    = overwrite_metrics
)

# Variant calling (ajuste: haplotype_caller requiere bed_file)
haplotype_caller(
  output_dir   = output_dir,
  folder_fasta = folder_fasta,
  fastq_dir    = fastq_dir,
  bed_file     = bed_file,
  gatk_bin     = gatk_bin,
  threads      = threads,
  overwrite    = overwrite_variants
)

# Genotype (ajuste: nombre correcto genotypeGVCFs)
genotypeGVCFs(
  folder_fasta = folder_fasta,
  output_dir   = output_dir,
  fastq_dir    = fastq_dir,
  gatk_bin     = gatk_bin,
  overwrite    = overwrite_variants
)

# Hard-filter
variantFiltration(
  folder_fasta = folder_fasta,
  output_dir   = output_dir,
  fastq_dir    = fastq_dir,
  gatk_bin     = gatk_bin,
  overwrite    = overwrite_filtering
)

# PASS bgzip + index
analysisReady(
  output_dir    = output_dir,
  fastq_dir     = fastq_dir,
  bcftools_bin  = bcftools_bin,
  bgzip_bin     = bgzip_bin,
  overwrite     = overwrite_ready
)

# Anotación (ajuste: recursos + campos como inputs)
anotation(
  folder_fasta  = folder_fasta,
  path_snpeff   = path_snpeff,
  output_dir    = output_dir,
  fastq_dir     = fastq_dir,
  clinvar_vcf   = clinvar_vcf,
  dbnsfp_db     = dbnsfp_db,
  gwas_db       = gwas_db,
  dbnsfp_fields = dbnsfp_fields,
  java_bin      = java_bin,
  bcftools_bin  = bcftools_bin,
  bgzip_bin     = bgzip_bin,
  overwrite     = overwrite_annot
)

message("Tiempo total: ", round(as.numeric(difftime(Sys.time(), start, units = "mins")), 2), " min")
