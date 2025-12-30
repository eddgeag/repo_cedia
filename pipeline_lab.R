load.libs <- c("data.table","tools","dplyr","tidyr","ggplot2","stringr","xtable",
               "VariantAnnotation","tibble","GenomicRanges","S4Vectors","ggpubr")
invisible(lapply(load.libs, require, character.only = TRUE))

get_sample_name <- function(fastq_dir){
  fq <- list.files(fastq_dir, pattern="\\.fq\\.gz$", full.names = FALSE)
  r1 <- fq[grepl("_R1\\.fq\\.gz$", fq)]
  stopifnot(length(r1) == 1)
  
  ## Eliminar _R1.fq.gz → nombre de muestra limpio
  x <- sub("_R1\\.fq\\.gz$", "", r1)
  x
}


control_calidad <- function(fastq_dir, output_dir) {
  
  qc_dir <- file.path(output_dir, "QC")
  dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)
  
  if (length(list.files(qc_dir)) == 0) {
    
    fastq_files <- list.files(
      fastq_dir,
      pattern = "\\.(fq|fastq)\\.gz$",
      full.names = TRUE
    )
    
    if (length(fastq_files) == 0) {
      stop("No se encontraron FASTQ en: ", fastq_dir)
    }
    
    command <- paste(
      "fastqc -t 4",
      paste(fastq_files, collapse = " "),
      "-o", qc_dir
    )
    
    system(command)
    
  } else {
    message("Ya se ha hecho el control de calidad")
  }
}


fn_exists_fasta <- function(folder_fasta) {
  fasta_files <- list.files(folder_fasta, pattern="\\.(fa|fasta)$", full.names=TRUE)
  if (length(fasta_files) == 0) stop("No existe archivo de referencia (.fa/.fasta) en: ", folder_fasta)
  if (length(fasta_files) > 1)  stop("Hay más de un FASTA en: ", folder_fasta, "\n", paste(basename(fasta_files), collapse=", "))
  fasta_files[[1]]
}



index_fasta_samtools <- function(folder_fasta) {
  fasta_file <- fn_exists_fasta(folder_fasta)
  if (length(list.files(folder_fasta, pattern="\\.fai$", full.names=TRUE)) == 0) {
    cmd <- paste("samtools faidx", shQuote(fasta_file))
    print(cmd); system(cmd)
  } else message("Ya esta el index fai")
}


index_bwa <- function(folder_fasta) {
  fasta_file <- fn_exists_fasta(folder_fasta)
  exts <- c("amb","ann","bwt","pac","sa")
  expected <- paste0(fasta_file, ".", exts)
  if (!all(file.exists(expected))) {
    message("Creando índices BWA...")
    system(paste("bwa index", shQuote(fasta_file)))
  } else message("Ya existen los índices BWA")
}


bwamem <- function(fastq_dir,
                   folder_fasta,
                   output_dir) {
  
  ## =========================
  ## 1) FASTA de referencia + índices
  ## =========================
  fasta_file <- fn_exists_fasta(folder_fasta)
  index_bwa(folder_fasta)
  
  ## Asegurar índice .fai
  if (length(list.files(dirname(fasta_file), pattern = "\\.fai$")) == 0) {
    index_fasta_samtools(folder_fasta)
  }
  
  ## =========================
  ## 2) FASTQ R1 / R2
  ## =========================
  fq <- list.files(fastq_dir, pattern = "\\.fq\\.gz$", full.names = TRUE)
  r1 <- fq[grepl("_R1\\.fq\\.gz$", basename(fq))]
  r2 <- fq[grepl("_R2\\.fq\\.gz$", basename(fq))]
  stopifnot(length(r1) == 1, length(r2) == 1)
  
  ## =========================
  ## 3) Nombre de muestra (CANÓNICO)
  ## =========================
  output_file_name <- get_sample_name(fastq_dir)
  
  ## =========================
  ## 4) Directorio de salida
  ## =========================
  mapping_output_dir <- file.path(output_dir, "mapping_output")
  dir.create(mapping_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  ## =========================
  ## 5) Archivos de salida
  ## =========================
  output_file_sam <- file.path(mapping_output_dir, paste0(output_file_name, ".sam"))
  output_file_bam <- file.path(mapping_output_dir, paste0(output_file_name, ".bam"))
  output_file_sorted_bam <- file.path(
    mapping_output_dir,
    paste0(output_file_name, ".sorted.bam")
  )
  output_file_sorted_bai <- paste0(output_file_sorted_bam, ".bai")
  
  ## =========================
  ## 6) Ejecutar mapeo si no existe BAM ordenado
  ## =========================
  if (!file.exists(output_file_sorted_bam)) {
    
    message("#### MAPPING (bwa mem) ####")
    
    comando_map <- paste(
      "bwa mem -M -t 8",
      "-R '@RG\\tID:", output_file_name,
      "\\tSM:", output_file_name,
      "\\tPL:ILLUMINA\\tLB:lib1\\tPU:unit1'",
      shQuote(fasta_file),
      shQuote(r1),
      shQuote(r2),
      ">",
      shQuote(output_file_sam)
    )
    
    print(comando_map)
    ret <- system(comando_map)
    if (ret != 0 || !file.exists(output_file_sam)) {
      stop("ERROR CRÍTICO: bwa mem falló para la muestra ", output_file_name)
    }
    
    message("#### SAM -> BAM ####")
    
    comando_sam2bam <- paste(
      "samtools view -b -h -@ 8",
      shQuote(output_file_sam),
      "-o",
      shQuote(output_file_bam)
    )
    
    print(comando_sam2bam)
    ret <- system(comando_sam2bam)
    if (ret != 0 || !file.exists(output_file_bam)) {
      stop("ERROR CRÍTICO: samtools view falló para ", output_file_name)
    }
    
    message("#### BAM -> SORTED BAM ####")
    
    comando_sort <- paste(
      "~/tools/gatk-4.6.1.0/gatk SortSam",
      "-CREATE_INDEX true",
      "-INPUT",  shQuote(output_file_bam),
      "-OUTPUT", shQuote(output_file_sorted_bam),
      "-SORT_ORDER coordinate",
      "-VALIDATION_STRINGENCY STRICT"
    )
    
    print(comando_sort)
    ret <- system(comando_sort)
    if (ret != 0 || !file.exists(output_file_sorted_bam)) {
      stop("ERROR CRÍTICO: SortSam falló para ", output_file_name)
    }
    
  } else {
    message("BAM ordenado ya existe")
  }
  
  ## =========================
  ## 7) Verificación OBLIGATORIA del índice .bai
  ## =========================
  if (!file.exists(output_file_sorted_bai)) {
    
    message("Índice .bai no encontrado para BAM ordenado. Generando índice...")
    
    cmd_index <- paste(
      "~/tools/gatk-4.6.1.0/gatk BuildBamIndex",
      "-I", shQuote(output_file_sorted_bam)
    )
    
    print(cmd_index)
    ret <- system(cmd_index)
    
    if (ret != 0 || !file.exists(output_file_sorted_bai)) {
      stop(
        "ERROR CRÍTICO: No se pudo generar el índice .bai para ",
        basename(output_file_sorted_bam)
      )
    }
  }
  
  message("BWA-MEM + SortSam completados correctamente (BAM + BAI verificados)")
}




markdups <- function(output_dir,
                     fastq_dir) {
  
  ## =========================
  ## 1) Nombre base de la muestra (FUENTE ÚNICA)
  ## =========================
  output_file_name <- get_sample_name(fastq_dir)
  
  ## =========================
  ## 2) Directorio de mapeo
  ## =========================
  mapping_output_dir <- file.path(output_dir, "mapping_output")
  
  ## =========================
  ## 3) BAM ordenado de entrada
  ## =========================
  bam_file <- file.path(
    mapping_output_dir,
    paste0(output_file_name, ".sorted.bam")
  )
  
  if (!file.exists(bam_file)) {
    stop("No existe BAM ordenado para MarkDuplicates: ", bam_file)
  }
  
  ## =========================
  ## 4) BAM de salida con duplicados marcados
  ## =========================
  mark_file <- file.path(
    mapping_output_dir,
    paste0(output_file_name, ".sorted.mark_dup.bam")
  )
  
  mark_bai <- paste0(mark_file, ".bai")
  
  ## =========================
  ## 5) Archivo de métricas
  ## =========================
  metrics_file <- file.path(
    mapping_output_dir,
    paste0(output_file_name, ".sorted.mark_dup.txt")
  )
  
  ## =========================
  ## 6) Ejecutar MarkDuplicates
  ## =========================
  if (!file.exists(mark_file)) {
    
    command <- paste(
      "~/tools/gatk-4.6.1.0/gatk MarkDuplicates",
      "-CREATE_INDEX true",
      "--ASSUME_SORTED true",
      "-INPUT",  shQuote(bam_file),
      "-OUTPUT", shQuote(mark_file),
      "-M",      shQuote(metrics_file),
      "--VALIDATION_STRINGENCY STRICT"
    )
    
    print(command)
    
    ret <- system(command)
    if (ret != 0) {
      stop("MarkDuplicates falló para la muestra: ", output_file_name)
    }
  } else {
    message("BAM con duplicados ya existe")
  }
  
  ## =========================
  ## 7) Verificación OBLIGATORIA del índice .bai
  ## =========================
  if (!file.exists(mark_bai)) {
    
    message("Índice .bai no encontrado para BAM con duplicados. Generando índice...")
    
    cmd_index <- paste(
      "~/tools/gatk-4.6.1.0/gatk BuildBamIndex",
      "-I", shQuote(mark_file)
    )
    
    print(cmd_index)
    
    ret <- system(cmd_index)
    if (ret != 0 || !file.exists(mark_bai)) {
      stop(
        "ERROR CRÍTICO: No se pudo generar el índice .bai para ",
        basename(mark_file)
      )
    }
  }
  
  message("MarkDuplicates completado correctamente (BAM + BAI verificados)")
}




create_dict <- function(folder_fasta) {
  fasta_file <- fn_exists_fasta(folder_fasta)
  fasta_dict <- paste0(tools::file_path_sans_ext(fasta_file), ".dict")
  if (!file.exists(fasta_dict)) {
    cmd <- paste("~/tools/gatk-4.6.1.0/gatk CreateSequenceDictionary",
                 "-R", shQuote(fasta_file),
                 "-O", shQuote(fasta_dict))
    print(cmd); system(cmd)
  } else message("Ya se ha creado el diccionario fasta")
}

# 
# creacion_readgroup <-
#   function(output_dir = output_dir,
#            fastq_dir = fastq_dir) {
#     ## el directorio del mapeo esta creado pero lo tenemos que guardar de nuevo
#     mapping_output_dir <- file.path(output_dir, "mapping_output")
#     fastq_files <- list.files(fastq_dir, full.names = F)
#     output_file_name <-
#       unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
#     ## quitamos la extension del archivo
#     output_file_name <- file_path_sans_ext(output_file_name)
#     ## obtenemos el archivo marcado con duplicados
#     mark_file <-
#       file.path(mapping_output_dir,
#                 paste0(output_file_name, ".sorted.mark_dup.bam"))
#     ## quitamos la extension y renombramos al archivo de salida
#     out_file <- paste0(file_path_sans_ext(mark_file), "_RG.bam")
#     
#     ## si ya se ha creado el archivo con grupo
#     
#     if (!file.exists(out_file)) {
#       command <-
#         paste(
#           "~/tools/gatk-4.6.1.0/gatk AddOrReplaceReadGroups I=",
#           mark_file,
#           "O=",
#           out_file,
#           "RGID=1 RGLB=lib2 RGPL=ILLUMINA RGPU=unit1 RGSM=1"
#         )
#       system(command = command, intern = T)
#       
#     } else{
#       message("Ya estan los grupos")
#     }
#     
#     
#     
#     
#   }



base_recalibrator <- function(folder_fasta,
                              output_dir,
                              folder_data_gatk,
                              fastq_dir) {
  
  ## =========================
  ## 1) Known-sites (uno por flag)
  ## =========================
  known_sites_file <- c(
    file.path(folder_data_gatk,
              "resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"),
    file.path(folder_data_gatk,
              "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz")
  )
  
  known_sites_args <- paste(
    paste0("--known-sites ", shQuote(known_sites_file)),
    collapse = " "
  )
  
  ## =========================
  ## 2) FASTA de referencia
  ## =========================
  fasta_file <- fn_exists_fasta(folder_fasta)
  
  ## =========================
  ## 3) Directorio de mapeo
  ## =========================
  mapping_output_dir <- file.path(output_dir, "mapping_output")
  
  ## =========================
  ## 4) Nombre base de la muestra (FUENTE ÚNICA)
  ## =========================
  output_file_name <- get_sample_name(fastq_dir)
  
  ## =========================
  ## 5) BAM de entrada (mark duplicates)
  ## =========================
  RG_file <- file.path(
    mapping_output_dir,
    paste0(output_file_name, ".sorted.mark_dup.bam")
  )
  
  if (!file.exists(RG_file)) {
    stop("No existe el BAM con duplicados marcados: ", RG_file)
  }
  
  ## =========================
  ## 6) Tabla de salida BQSR
  ## =========================
  out_file <- file.path(
    mapping_output_dir,
    paste0(output_file_name, ".recal_data.table")
  )
  
  ## =========================
  ## 7) Ejecutar BaseRecalibrator
  ## =========================
  if (!file.exists(out_file)) {
    
    command <- paste(
      "~/tools/gatk-4.6.1.0/gatk BaseRecalibrator",
      "-I", shQuote(RG_file),
      "-R", shQuote(fasta_file),
      known_sites_args,
      "-O", shQuote(out_file)
    )
    
    print(command)
    system(command)
    
  } else {
    message("Ya existe la tabla de BaseRecalibrator")
  }
}



applybqsr <- function(folder_fasta,
                      output_dir,
                      fastq_dir) {
  
  ## =========================
  ## 1) FASTA de referencia
  ## =========================
  fasta_file <- fn_exists_fasta(folder_fasta)
  
  ## =========================
  ## 2) Directorio de mapeo
  ## =========================
  mapping_output_dir <- file.path(output_dir, "mapping_output")
  
  ## =========================
  ## 3) Nombre base de la muestra (FUENTE ÚNICA)
  ## =========================
  output_file_name <- get_sample_name(fastq_dir)
  
  ## =========================
  ## 4) BAM de entrada (mark duplicates)
  ## =========================
  bam_file <- file.path(
    mapping_output_dir,
    paste0(output_file_name, ".sorted.mark_dup.bam")
  )
  
  if (!file.exists(bam_file)) {
    stop("No existe el BAM de entrada para ApplyBQSR: ", bam_file)
  }
  
  ## =========================
  ## 5) Tabla BQSR
  ## =========================
  recal_data_table <- file.path(
    mapping_output_dir,
    paste0(output_file_name, ".recal_data.table")
  )
  
  if (!file.exists(recal_data_table)) {
    stop(
      "No existe la tabla BQSR: ", recal_data_table,
      "\n¿Ejecutaste base_recalibrator() antes?"
    )
  }
  
  ## =========================
  ## 6) BAM de salida (BQSR aplicado)
  ## =========================
  out_file <- file.path(
    mapping_output_dir,
    paste0(output_file_name, ".sorted.mark_dup_bqsr.bam")
  )
  
  out_bai <- paste0(out_file, ".bai")
  
  ## =========================
  ## 7) Ejecutar ApplyBQSR
  ## =========================
  if (!file.exists(out_file)) {
    
    command <- paste(
      "~/tools/gatk-4.6.1.0/gatk ApplyBQSR",
      "-I", shQuote(bam_file),
      "-R", shQuote(fasta_file),
      "--bqsr-recal-file", shQuote(recal_data_table),
      "-O", shQuote(out_file)
    )
    
    print(command)
    
    ret <- system(command)
    if (ret != 0) {
      stop("ApplyBQSR falló para la muestra: ", output_file_name)
    }
  } else {
    message("BAM con BQSR ya existe")
  }
  
  ## =========================
  ## 8) Verificación OBLIGATORIA del índice .bai
  ## =========================
  if (!file.exists(out_bai)) {
    
    message("Índice .bai no encontrado. Generando índice...")
    
    cmd_index <- paste(
      "~/tools/gatk-4.6.1.0/gatk BuildBamIndex",
      "-I", shQuote(out_file)
    )
    
    print(cmd_index)
    
    ret <- system(cmd_index)
    if (ret != 0 || !file.exists(out_bai)) {
      stop(
        "ERROR CRÍTICO: No se pudo generar el índice .bai para ",
        basename(out_file)
      )
    }
  }
  
  message("ApplyBQSR completado correctamente (BAM + BAI verificados)")
}



bam_statistics <- function(folder_fasta,
                           fastq_dir,
                           output_dir) {
  
  ## FASTA de referencia
  fasta_file <- fn_exists_fasta(folder_fasta)
  
  ## Directorio de mapeo
  mapping_output_dir <- file.path(output_dir, "mapping_output")
  
  ## Nombre base de la muestra (CANÓNICO)
  output_file_name <- get_sample_name(fastq_dir)
  
  ## BAM final (BQSR aplicado)
  bam_file <- file.path(
    mapping_output_dir,
    paste0(output_file_name, ".sorted.mark_dup_bqsr.bam")
  )
  
  if (!file.exists(bam_file)) {
    stop("No existe el BAM final con BQSR: ", bam_file)
  }
  
  ## Directorio de salida de métricas
  out_dir <- file.path(
    output_dir,
    "bam_metrics",
    output_file_name
  )
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  ## Prefijo de salida
  out_prefix <- file.path(out_dir, "CollectMultipleMetrics")
  
  ## Archivo clave esperado
  expected_file <- paste0(out_prefix, ".alignment_summary_metrics")
  
  ## Ejecutar solo si no existen métricas
  if (!file.exists(expected_file)) {
    
    command <- paste(
      "~/tools/gatk-4.6.1.0/gatk CollectMultipleMetrics",
      "-R", fasta_file,
      "-I", bam_file,
      "-O", out_prefix
    )
    
    print(command)
    system(command)
    
    ## MultiQC
    command_mqc <- paste(
      "multiqc",
      out_dir,
      "-o",
      out_dir
    )
    print(command_mqc)
    system(command_mqc)
    
  } else {
    message("Las métricas BAM ya existen para la muestra")
  }
}



haplotype_caller <- function(output_dir,
                             folder_fasta,
                             fastq_dir) {
  
  ## =========================
  ## 1) FASTA de referencia
  ## =========================
  fasta_file <- fn_exists_fasta(folder_fasta)
  
  ## =========================
  ## 2) Directorio de mapeo
  ## =========================
  mapping_output_dir <- file.path(output_dir, "mapping_output")
  
  ## =========================
  ## 3) Nombre base de la muestra (FUENTE ÚNICA)
  ## =========================
  output_file_name <- get_sample_name(fastq_dir)
  
  ## =========================
  ## 4) BAM de entrada (BQSR aplicado)
  ## =========================
  bam_file <- file.path(
    mapping_output_dir,
    paste0(output_file_name, ".sorted.mark_dup_bqsr.bam")
  )
  
  bam_bai <- paste0(bam_file, ".bai")
  
  if (!file.exists(bam_file)) {
    stop("No existe el BAM BQSR para HaplotypeCaller: ", bam_file)
  }
  
  if (!file.exists(bam_bai)) {
    stop("No existe el índice .bai del BAM BQSR: ", bam_bai)
  }
  
  ## =========================
  ## 5) Directorio y archivo de salida
  ## =========================
  out_dir <- file.path(output_dir, "variantCalling")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  out_file <- file.path(
    out_dir,
    paste0(output_file_name, ".g.vcf.gz")
  )
  
  out_tbi <- paste0(out_file, ".tbi")
  
  ## =========================
  ## 6) Ejecutar HaplotypeCaller (GVCF)
  ## =========================
  if (!file.exists(out_file)) {
    
    command <- paste(
      "~/tools/gatk-4.6.1.0/gatk HaplotypeCaller",
      "-I", shQuote(bam_file),
      "-R", shQuote(fasta_file),
      "-ERC GVCF",
      "-O", shQuote(out_file),
      "--native-pair-hmm-threads 8"
    )
    
    print(command)
    
    ret <- system(command)
    if (ret != 0 || !file.exists(out_file)) {
      stop("ERROR CRÍTICO: HaplotypeCaller falló para la muestra ", output_file_name)
    }
  } else {
    message("GVCF ya existe para la muestra")
  }
  
  ## =========================
  ## 7) Verificación OBLIGATORIA del índice .tbi
  ## =========================
  if (!file.exists(out_tbi)) {
    
    message("Índice .tbi no encontrado para GVCF. Generando índice...")
    
    cmd_index <- paste(
      "~/tools/gatk-4.6.1.0/gatk IndexFeatureFile",
      "-I", shQuote(out_file)
    )
    
    print(cmd_index)
    
    ret <- system(cmd_index)
    if (ret != 0 || !file.exists(out_tbi)) {
      stop(
        "ERROR CRÍTICO: No se pudo generar el índice .tbi para ",
        basename(out_file)
      )
    }
  }
  
  message("HaplotypeCaller completado correctamente (GVCF + TBI verificados)")
}

genotypeGVCF <- function(folder_fasta,
                         output_dir,
                         fastq_dir) {
  
  ## =========================
  ## 1) FASTA de referencia
  ## =========================
  fasta_file <- fn_exists_fasta(folder_fasta)
  
  ## =========================
  ## 2) Nombre base de la muestra (FUENTE ÚNICA)
  ## =========================
  output_file_name <- get_sample_name(fastq_dir)
  
  ## =========================
  ## 3) Directorio de variantes
  ## =========================
  var_dir <- file.path(output_dir, "variantCalling")
  
  ## =========================
  ## 4) GVCF de entrada (HaplotypeCaller)
  ## =========================
  file_in <- file.path(
    var_dir,
    paste0(output_file_name, ".g.vcf.gz")
  )
  
  file_in_tbi <- paste0(file_in, ".tbi")
  
  if (!file.exists(file_in)) {
    stop("No existe el GVCF de entrada para GenotypeGVCFs: ", file_in)
  }
  
  if (!file.exists(file_in_tbi)) {
    stop("No existe el índice .tbi del GVCF de entrada: ", file_in_tbi)
  }
  
  ## =========================
  ## 5) VCF genotipado de salida
  ## =========================
  file_out <- file.path(
    var_dir,
    paste0(output_file_name, "_sample.raw.vcf.gz")
  )
  
  file_out_tbi <- paste0(file_out, ".tbi")
  
  ## =========================
  ## 6) Ejecutar GenotypeGVCFs
  ## =========================
  if (!file.exists(file_out)) {
    
    command <- paste(
      "~/tools/gatk-4.6.1.0/gatk GenotypeGVCFs",
      "-R", shQuote(fasta_file),
      "-V", shQuote(file_in),
      "-O", shQuote(file_out)
    )
    
    print(command)
    
    ret <- system(command)
    if (ret != 0 || !file.exists(file_out)) {
      stop(
        "ERROR CRÍTICO: GenotypeGVCFs falló para la muestra ",
        output_file_name
      )
    }
    
  } else {
    message("VCF genotipado ya existe para la muestra")
  }
  
  ## =========================
  ## 7) Verificación OBLIGATORIA del índice .tbi
  ## =========================
  if (!file.exists(file_out_tbi)) {
    
    message("Índice .tbi no encontrado para VCF genotipado. Generando índice...")
    
    cmd_index <- paste(
      "~/tools/gatk-4.6.1.0/gatk IndexFeatureFile",
      "-I", shQuote(file_out)
    )
    
    print(cmd_index)
    
    ret <- system(cmd_index)
    if (ret != 0 || !file.exists(file_out_tbi)) {
      stop(
        "ERROR CRÍTICO: No se pudo generar el índice .tbi para ",
        basename(file_out)
      )
    }
  }
  
  message("GenotypeGVCFs completado correctamente (VCF + TBI verificados)")
}

# variantRecalibrator <- function(fastq_dir,
#                                 folder_fasta,
#                                 folder_data_gatk,
#                                 output_dir) {
#   
#   ## FASTA de referencia
#   fasta_file <- fn_exists_fasta(folder_fasta)
#   
#   ## Nombre base de la muestra
#   fastq_files <- list.files(fastq_dir, full.names = FALSE)
#   output_file_name <-
#     unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
#   output_file_name <- file_path_sans_ext(output_file_name)
#   
#   ## VCF de entrada (GenotypeGVCFs)
#   in_file <- file.path(
#     output_dir,
#     "variantCalling",
#     paste0(output_file_name, "_apply_genotypeGVCF.g.vcf.gz")
#   )
#   
#   ## Archivos de salida (SNPs)
#   snps_recal_file <- file.path(
#     output_dir,
#     "variantCalling",
#     paste0(output_file_name, ".snps.recal")
#   )
#   
#   tranches_file <- file.path(
#     output_dir,
#     "variantCalling",
#     paste0(output_file_name, ".snps.tranches")
#   )
#   
#   if (!file.exists(snps_recal_file) || !file.exists(tranches_file)) {
#     
#     command <- paste(
#       "~/tools/gatk-4.6.1.0/gatk VariantRecalibrator",
#       "-R", fasta_file,
#       "-V", in_file,
#       
#       "--resource:hapmap,known=false,training=true,truth=true,prior=15.0",
#       file.path(folder_data_gatk,
#                 "resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz"),
#       
#       "--resource:omni,known=false,training=true,truth=true,prior=12.0",
#       file.path(folder_data_gatk,
#                 "resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz"),
#       
#       "--resource:1000G,known=false,training=true,truth=false,prior=10.0",
#       file.path(folder_data_gatk,
#                 "resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz"),
#       
#       "--resource:dbsnp,known=true,training=false,truth=false,prior=2.0",
#       file.path(folder_data_gatk,
#                 "resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"),
#       
#       "-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP",
#       "-mode SNP",
#       "-O", snps_recal_file,
#       "--tranches-file", tranches_file
#     )
#     
#     print(command)
#     system(command)
#     
#   } else {
#     message("Ya se ha hecho el recalibrado VQSR de SNPs")
#   }
# }

# 
# applyVQSR <- function(folder_fasta,
#                       fastq_dir,
#                       output_dir) {
#   
#   ## FASTA de referencia
#   fasta_file <- fn_exists_fasta(folder_fasta)
#   
#   ## Nombre base de la muestra
#   fastq_files <- list.files(fastq_dir, full.names = FALSE)
#   output_file_name <-
#     unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
#   output_file_name <- file_path_sans_ext(output_file_name)
#   
#   ## VCF de entrada (GenotypeGVCFs)
#   in_file <- file.path(
#     output_dir,
#     "variantCalling",
#     paste0(output_file_name, "_apply_genotypeGVCF.g.vcf.gz")
#   )
#   
#   ## Archivos VQSR (generados por VariantRecalibrator)
#   snps_recal_file <- file.path(
#     output_dir,
#     "variantCalling",
#     paste0(output_file_name, ".snps.recal")
#   )
#   
#   tranches_file <- file.path(
#     output_dir,
#     "variantCalling",
#     paste0(output_file_name, ".snps.tranches")
#   )
#   
#   ## Archivo de salida
#   out_file <- file.path(
#     output_dir,
#     "variantCalling",
#     paste0(output_file_name, ".snps.vqsr.vcf")
#   )
#   
#   ## Aplicar VQSR solo si no existe
#   if (!file.exists(out_file)) {
#     
#     command <- paste(
#       "~/tools/gatk-4.6.1.0/gatk ApplyVQSR",
#       "-R", fasta_file,
#       "-V", in_file,
#       "-O", out_file,
#       "--truth-sensitivity-filter-level 99.0",
#       "--tranches-file", tranches_file,
#       "--recal-file", snps_recal_file,
#       "-mode SNP"
#     )
#     
#     print(command)
#     system(command)
#     
#   } else {
#     message("Ya se ha aplicado VQSR (SNPs)")
#   }
#   
# }

variantFiltration <- function(folder_fasta,
                              output_dir,
                              fastq_dir) {
  
  ## =========================
  ## 1) FASTA de referencia
  ## =========================
  fasta_file <- fn_exists_fasta(folder_fasta)
  
  ## =========================
  ## 2) Nombre base de la muestra (FUENTE ÚNICA)
  ## =========================
  output_file_name <- get_sample_name(fastq_dir)
  
  ## =========================
  ## 3) Directorio de variantes
  ## =========================
  var_dir <- file.path(output_dir, "variantCalling")
  
  ## =========================
  ## 4) VCF de entrada (GenotypeGVCFs)
  ## =========================
  in_vcf <- file.path(
    var_dir,
    paste0(output_file_name, "_sample.raw.vcf.gz")
  )
  in_vcf_tbi <- paste0(in_vcf, ".tbi")
  
  if (!file.exists(in_vcf)) {
    stop("No existe el VCF de entrada para hard-filter: ", in_vcf)
  }
  if (!file.exists(in_vcf_tbi)) {
    stop("No existe el índice .tbi del VCF de entrada: ", in_vcf_tbi)
  }
  
  ## =========================
  ## 5) Comprobación VCF NO vacío
  ## =========================
  n_variants <- as.integer(
    system(paste("bcftools view -H", shQuote(in_vcf), "| head -n 1 | wc -l"),
           intern = TRUE)
  )
  if (n_variants == 0) {
    stop("ERROR CRÍTICO: el VCF de entrada está vacío: ", in_vcf)
  }
  
  ## =========================
  ## 6) Archivos intermedios
  ## =========================
  snps_vcf        <- file.path(var_dir, paste0(output_file_name, ".snps.vcf"))
  indels_vcf      <- file.path(var_dir, paste0(output_file_name, ".indels.vcf"))
  
  snps_filt_vcf   <- file.path(var_dir, paste0(output_file_name, ".snps.hardfiltered.vcf"))
  indels_filt_vcf <- file.path(var_dir, paste0(output_file_name, ".indels.hardfiltered.vcf"))
  
  merged_vcf      <- file.path(var_dir, paste0(output_file_name, ".hardfiltered.vcf"))
  merged_vcf_tbi  <- paste0(merged_vcf, ".tbi")
  
  ## =========================
  ## 7) Seleccionar SNPs
  ## =========================
  if (!file.exists(snps_vcf)) {
    cmd_snps <- paste(
      "~/tools/gatk-4.6.1.0/gatk SelectVariants",
      "-R", shQuote(fasta_file),
      "-V", shQuote(in_vcf),
      "--select-type-to-include SNP",
      "-O", shQuote(snps_vcf)
    )
    print(cmd_snps)
    if (system(cmd_snps) != 0 || !file.exists(snps_vcf)) {
      stop("Fallo al seleccionar SNPs")
    }
  }
  
  ## =========================
  ## 8) Seleccionar INDELs
  ## =========================
  if (!file.exists(indels_vcf)) {
    cmd_indels <- paste(
      "~/tools/gatk-4.6.1.0/gatk SelectVariants",
      "-R", shQuote(fasta_file),
      "-V", shQuote(in_vcf),
      "--select-type-to-include INDEL",
      "-O", shQuote(indels_vcf)
    )
    print(cmd_indels)
    if (system(cmd_indels) != 0 || !file.exists(indels_vcf)) {
      stop("Fallo al seleccionar INDELs")
    }
  }
  
  ## =========================
  ## 9) Hard-filter SNPs
  ## =========================
  if (!file.exists(snps_filt_vcf)) {
    cmd_snps_filt <- paste(
      "~/tools/gatk-4.6.1.0/gatk VariantFiltration",
      "-R", shQuote(fasta_file),
      "-V", shQuote(snps_vcf),
      "-O", shQuote(snps_filt_vcf),
      "--filter-name QD2 --filter-expression 'QD < 2.0'",
      "--filter-name FS60 --filter-expression 'FS > 60.0'",
      "--filter-name MQ40 --filter-expression 'MQ < 40.0'",
      "--filter-name MQRS-12.5 --filter-expression 'MQRankSum < -12.5'",
      "--filter-name RPRS-8 --filter-expression 'ReadPosRankSum < -8.0'",
      "--filter-name SOR3 --filter-expression 'SOR > 3.0'"
    )
    print(cmd_snps_filt)
    if (system(cmd_snps_filt) != 0 || !file.exists(snps_filt_vcf)) {
      stop("Fallo en hard-filter de SNPs")
    }
  }
  
  ## =========================
  ## 10) Hard-filter INDELs
  ## =========================
  if (!file.exists(indels_filt_vcf)) {
    cmd_indels_filt <- paste(
      "~/tools/gatk-4.6.1.0/gatk VariantFiltration",
      "-R", shQuote(fasta_file),
      "-V", shQuote(indels_vcf),
      "-O", shQuote(indels_filt_vcf),
      "--filter-name QD2 --filter-expression 'QD < 2.0'",
      "--filter-name FS200 --filter-expression 'FS > 200.0'",
      "--filter-name RPRS-20 --filter-expression 'ReadPosRankSum < -20.0'",
      "--filter-name SOR10 --filter-expression 'SOR > 10.0'"
    )
    print(cmd_indels_filt)
    if (system(cmd_indels_filt) != 0 || !file.exists(indels_filt_vcf)) {
      stop("Fallo en hard-filter de INDELs")
    }
  }
  
  ## =========================
  ## 11) Merge SNPs + INDELs
  ## =========================
  if (!file.exists(merged_vcf)) {
    cmd_merge <- paste(
      "~/tools/gatk-4.6.1.0/gatk MergeVcfs",
      "-I", shQuote(snps_filt_vcf),
      "-I", shQuote(indels_filt_vcf),
      "-O", shQuote(merged_vcf)
    )
    print(cmd_merge)
    if (system(cmd_merge) != 0 || !file.exists(merged_vcf)) {
      stop("Fallo al fusionar SNPs + INDELs")
    }
  }
  
  ## =========================
  ## 12) Índice .tbi obligatorio
  ## =========================
  if (!file.exists(merged_vcf_tbi)) {
    cmd_index <- paste(
      "~/tools/gatk-4.6.1.0/gatk IndexFeatureFile",
      "-I", shQuote(merged_vcf)
    )
    print(cmd_index)
    if (system(cmd_index) != 0 || !file.exists(merged_vcf_tbi)) {
      stop("No se pudo generar el índice .tbi del VCF final")
    }
  }
  
  ## =========================
  ## 13) Comprobación variantes PASS
  ## =========================
  n_pass <- as.integer(
    system(
      paste("bcftools view -f PASS -H", shQuote(merged_vcf),
            "| head -n 1 | wc -l"),
      intern = TRUE
    )
  )
  
  if (n_pass == 0) {
    stop(
      "ERROR CRÍTICO: el VCF final no contiene variantes PASS.\n",
      "Revisa calidad, cobertura o filtros."
    )
  }
  
  message("Hard-filter SNPs + INDELs completado correctamente (VCF válido y con PASS)")
}





analysisReady <- function(folder_fasta,
                          output_dir,
                          fastq_dir) {
  
  ## =========================
  ## 1) Nombre base de la muestra (FUENTE ÚNICA)
  ## =========================
  output_file_name <- get_sample_name(fastq_dir)
  
  ## =========================
  ## 2) Directorio de variantes
  ## =========================
  var_dir <- file.path(output_dir, "variantCalling")
  
  ## =========================
  ## 3) VCF de entrada (hard-filter SNPs + INDELs, mergeado)
  ## =========================
  in_file <- file.path(
    var_dir,
    paste0(output_file_name, ".hardfiltered.vcf")
  )
  
  if (!file.exists(in_file)) {
    stop("No existe el VCF hardfiltered de entrada: ", in_file)
  }
  
  ## =========================
  ## 4) VCF PASS final
  ## =========================
  out_file <- file.path(
    var_dir,
    paste0(output_file_name, ".hardfiltered.pass.vcf")
  )
  
  ## =========================
  ## 5) Filtrar solo variantes PASS
  ## =========================
  if (!file.exists(out_file)) {
    
    command <- paste(
      "bcftools view -f PASS -O v -o",
      shQuote(out_file),
      shQuote(in_file)
    )
    
    print(command)
    system(command)
    
  } else {
    message("El VCF PASS final ya existe")
  }
}

anotation <- function(folder_fasta,
                      path_snpeff,
                      output_dir,
                      fastq_dir) {
  
  ## =========================
  ## 1) Nombre base de la muestra (ÚNICO)
  ## =========================
  output_file_name <- get_sample_name(fastq_dir)
  
  ## =========================
  ## 2) Directorio de anotación
  ## =========================
  anotacion_dir <- file.path(output_dir, "anotation")
  if (!dir.exists(anotacion_dir)) {
    dir.create(anotacion_dir, recursive = TRUE)
  }
  
  ## =========================
  ## 3) VCF de entrada (PASS final)
  ## =========================
  in_file <- file.path(
    output_dir,
    "variantCalling",
    paste0(output_file_name, ".hardfiltered.pass.vcf")
  )
  
  if (!file.exists(in_file)) {
    stop("No existe el VCF PASS para anotación: ", in_file)
  }
  
  ## =========================
  ## 4) Archivos de salida
  ## =========================
  output_file_anno1 <- file.path(anotacion_dir,
                                 paste0(output_file_name, "_anno_snpeff.vcf"))
  
  output_file_anno1_1 <- file.path(anotacion_dir,
                                   paste0(output_file_name, "_anno_snpeff_varType.vcf"))
  
  output_file_anno2 <- file.path(anotacion_dir,
                                 paste0(output_file_name, "_anno_snpeff_clinvar.vcf"))
  
  output_file_anno3 <- file.path(anotacion_dir,
                                 paste0(output_file_name, "_anno_snpeff_clinvar_dbnsfp.vcf"))
  
  output_file_anno4 <- file.path(anotacion_dir,
                                 paste0(output_file_name, "_anno_snpeff_clinvar_dbnsfp_gwas.vcf"))
  
  output_file_anno5 <- file.path(anotacion_dir,
                                 paste0(output_file_name, "_anno_final.vcf"))
  
  ## =========================
  ## 5) Comandos
  ## =========================
  
  ## snpEff
  cmd1 <- paste(
    "java -Xmx32g -jar",
    shQuote(file.path(path_snpeff, "snpEff.jar")),
    "hg38 -v",
    shQuote(in_file),
    ">",
    shQuote(output_file_anno1)
  )
  
  ## SnpSift varType
  cmd2 <- paste(
    "java -Xmx32g -jar",
    shQuote(file.path(path_snpeff, "SnpSift.jar")),
    "varType -v",
    shQuote(output_file_anno1),
    ">",
    shQuote(output_file_anno1_1)
  )
  
  ## ClinVar
  cmd3 <- paste(
    "java -Xmx32g -jar",
    shQuote(file.path(path_snpeff, "SnpSift.jar")),
    "annotate -v",
    shQuote("~/NAS_NGS/datos_exomas/datos_clinvar/clinvar.vcf.gz"),
    shQuote(output_file_anno1_1),
    ">",
    shQuote(output_file_anno2)
  )
  
  ## dbNSFP (completo)
  cmd4 <- paste(
    "java -Xmx32g -jar",
    shQuote(file.path(path_snpeff, "SnpSift.jar")),
    "dbnsfp -v -db",
    shQuote("~/NAS_NGS/datos_exomas/datos_dbsnp/dbNSFP5.1a_grch38.gz"),
    shQuote(output_file_anno2),
    ">",
    shQuote(output_file_anno3)
  )
  
  ## GWAS Catalog
  cmd5 <- paste(
    "java -Xmx32g -jar",
    shQuote(file.path(path_snpeff, "SnpSift.jar")),
    "gwasCat -db",
    shQuote("~/NAS_NGS/datos_exomas/gwas/gwascatalog.txt"),
    shQuote(output_file_anno3),
    ">",
    shQuote(output_file_anno4)
  )
  
  ## dbNSFP reducido (campos clave)
  campos <- paste(
    "aaref,aaalt,rs_dbSNP,HGVSc_snpEff,HGVSp_snpEff,APPRIS,M-CAP_pred,",
    "CADD_phred,clinvar_OMIM_id,clinvar_Orphanet_id,clinvar_MedGen_id,",
    "AlphaMissense_pred,Reliability_index",
    sep = ""
  )
  
  cmd6 <- paste(
    "java -Xmx32g -jar",
    shQuote(file.path(path_snpeff, "SnpSift.jar")),
    "dbnsfp -v -db",
    shQuote("~/NAS_NGS/datos_exomas/datos_dbsnp/dbNSFP5.1a_grch38.gz"),
    "-f", shQuote(campos),
    shQuote(output_file_anno4),
    ">",
    shQuote(output_file_anno5)
  )
  
  ## =========================
  ## 6) Ejecución secuencial
  ## =========================
  if (!file.exists(output_file_anno1))  { print(cmd1); system(cmd1) }
  if (!file.exists(output_file_anno1_1)){ print(cmd2); system(cmd2) }
  if (!file.exists(output_file_anno2))  { print(cmd3); system(cmd3) }
  if (!file.exists(output_file_anno3))  { print(cmd4); system(cmd4) }
  if (!file.exists(output_file_anno4))  { print(cmd5); system(cmd5) }
  if (!file.exists(output_file_anno5))  { print(cmd6); system(cmd6) }
  
  message("ANOTACIÓN COMPLETA FINALIZADA")
}



obtener_exoma_overlap <- function(bd_list, exoma_df, cromosomas, muestra) {
  # 1) Unir bd_list en un solo data.frame, convirtiendo todo a caracteres y filtrando cromosomas
  cat(
    "\nCheckpoint 1: Uniendo bd_list. Cols bd_list[[1]]:",
    paste(names(bd_list[[1]]), collapse = ", "),
    "\n"
  )
  
  
  
  bd_combined <-
    lapply(bd_list, function(x) {
      df <- as.data.frame(lapply(x, as.character), stringsAsFactors = FALSE)
      df[df$Chr %in% cromosomas, ]
    }) %>%
    bind_rows()
  cat("Checkpoint 2: bd_combined columns:",
      paste(names(bd_combined), collapse = ", "),
      "\n")
  
  # *** Aquí verifica la existencia de Gene.refGene ***
  if (!"Gene.refGene" %in% names(bd_combined)) {
    stop(
      "La columna 'Gene.refGene' NO existe. Las columnas disponibles son:\n",
      paste(names(bd_combined), collapse = ", ")
    )
  }
  
  # 2) Agrupar por Chr, Start, End, Gene.refGene; concatenar "codigo" y contar N
  
  bd_grouped <- bd_combined %>%
    group_by(Chr, Start, End, Gene.refGene) %>%
    summarise(
      n = n(),
      gene_name = dplyr::first(Gene.refGene),
      # Aquí
      paste_m = toString(codigo),
      .groups = "drop"
    ) %>%
    mutate(paste_m = sapply(strsplit(paste_m, ",\\s*"), function(v)
      paste(unique(v), collapse = ",")),
      N = sapply(strsplit(paste_m, ","), length))
  # Ya no necesitas rename(gene_name = Gene.refGene)
  cat("Checkpoint 3: bd_grouped columns:",
      paste(names(bd_grouped), collapse = ", "),
      "\n")
  
  
  
  # 3) Filtrar únicamente las filas que contienen la muestra de interés
  bd_filtrado_muestra <- bd_grouped[grep(muestra, bd_grouped$paste_m), ]
  cat("Checkpoint 4: Filtrado muestra. Rows:",
      nrow(bd_filtrado_muestra),
      "\n")
  
  # 4) Asegurar que Start y End sean numéricos
  bd_filtrado_muestra <- bd_filtrado_muestra %>%
    mutate(Start = as.numeric(Start), End   = as.numeric(End))
  cat("Checkpoint 5: bd_filtrado_muestra columns:",
      paste(names(bd_filtrado_muestra), collapse = ", "),
      "\n")
  
  
  # 5) Crear objeto GRanges a partir de bd_filtrado_muestra
  gr_bd <- makeGRangesFromDataFrame(
    bd_filtrado_muestra,
    seqnames.field     = "Chr",
    start.field        = "Start",
    end.field          = "End",
    keep.extra.columns = TRUE
  )
  
  # 6) Convertir exoma_df a GRanges (asegurar POS y END como numéricos)
  exoma_numeric <- exoma_df %>%
    mutate(POS = as.numeric(START), END = as.numeric(END))
  
  gr_exoma <- makeGRangesFromDataFrame(
    exoma_numeric,
    seqnames.field     = "CHROM",
    start.field        = "POS",
    end.field          = "END",
    keep.extra.columns = TRUE
  )
  
  cat("Checkpoint 6: GRanges creados\n")
  
  # 7) Encontrar solapamientos
  hits <- findOverlaps(gr_bd, gr_exoma)
  
  # 8) Crear data.frame de hits con identificadores únicos para detectar duplicados
  query_ids <- paste0(seqnames(gr_bd)[queryHits(hits)], ":", start(gr_bd)[queryHits(hits)], "-", end(gr_bd)[queryHits(hits)])
  subject_ids <- paste0(seqnames(gr_exoma)[subjectHits(hits)],
                        ":",
                        start(gr_exoma)[subjectHits(hits)],
                        "-",
                        end(gr_exoma)[subjectHits(hits)])
  
  hits_df <- data.frame(
    query = queryHits(hits),
    subject = subjectHits(hits),
    query_id = query_ids,
    subject_id = subject_ids,
    stringsAsFactors = FALSE
  )
  
  # 9) Eliminar duplicados en los pares query-subject para mantener correspondencia 1 a 1
  # Puedes decidir qué criterio usar para duplicados, aquí eliminamos duplicados en query_id + subject_id
  hits_df_unique <- hits_df[!duplicated(paste0(hits_df$query_id, "_", hits_df$subject_id)), ]
  
  cat("Duplicados eliminados en hits:",
      nrow(hits_df) - nrow(hits_df_unique),
      "\n")
  
  # 10) Extraer granges filtrados según hits únicos
  gr_bd_hits    <- gr_bd[hits_df_unique$query]
  gr_exoma_hits <- gr_exoma[hits_df_unique$subject]
  
  # 11) Extraer metadatos para gr_bd_hits y gr_exoma_hits
  bd_meta <- as_tibble(mcols(gr_bd_hits)[, c("gene_name", "N", "paste_m")])
  exoma_meta <- as_tibble(mcols(gr_exoma_hits))
  
  # 12) Construir tibble resultados
  resultados <- tibble(
    CHROM     = as.character(seqnames(gr_bd_hits)),
    Start     = start(gr_bd_hits),
    End       = end(gr_bd_hits),
    N         = bd_meta$N,
    samples   = bd_meta$paste_m,
    !!!exoma_meta
  )
  
  resultados <- as.data.frame(resultados)
  
  cat("Checkpoint 7: resultados columns:",
      paste(names(resultados), collapse = ", "),
      "\n")
  resultados_unicos <- resultados[!duplicated(resultados[, c("gene_name", "Start", "End", "CHROM")]), ]
  resultados_unicos$N <- sapply(strsplit(resultados_unicos$samples, ","), length)
  colnames(resultados_unicos)[colnames(resultados_unicos) == "Start"]   <- "POS"
  colnames(resultados_unicos)[colnames(resultados_unicos) == "End"]     <- "END"
  
  primeros <- c("CHROM", "POS", "END", "gene_name", "N", "samples")
  restantes <- setdiff(names(resultados_unicos), primeros)
  bd_exoma_overlap <- resultados_unicos[, c(primeros, restantes)]
  cat("Checkpoint 8: bd_exoma_overlap columns:",
      paste(names(bd_exoma_overlap), collapse = ", "),
      "\n")
  
  return(bd_exoma_overlap)
}

buscar_herencia <- function(df) {
  vector_hpo <- df$hpo  # Vector de la columna "hpo"
  vector_valores <-
    strsplit(vector_hpo, ";")  # Separar los valores por el delimitador ";"
  vector_resultado1 <-
    sapply(vector_valores, function(x)
      paste(unique(x[grep("inheritance", ignore.case = T, x = x)]), sep = ",", collapse = ","))
  X <- bind_cols(herencia = vector_resultado1, df)
  return(X)
}


process_vcf_to_table <- function(folder_fasta,
                                 output_dir,
                                 fastq_dir,
                                 muestra,
                                 db,
                                 hpo_file_) {
  print("=== Iniciando process_vcf_to_table ===")
  fastq_files <- list.files(fastq_dir, full.names = F)
  print(paste("FASTQ files encontrados:", paste(fastq_files, collapse = ", ")))
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  output_file_name <- file_path_sans_ext(output_file_name)
  anotacion_dir <- file.path(output_dir, "anotation")
  output_file_anno5 <- file.path(anotacion_dir, paste0(output_file_name, "_anno_final.vcf"))
  
  post_process_dir <- file.path(output_dir, "post_process_results")
  if (!dir.exists(post_process_dir)) {
    dir.create(post_process_dir, recursive = T)
    print(paste("Directorio creado:", post_process_dir))
  } else {
    print(paste("Directorio ya existe:", post_process_dir))
  }
  archivo_final <- file.path(post_process_dir, "file_ready_analysis_optimized.csv")
  
  if (!file.exists(archivo_final)) {
    vcf_file <- output_file_anno5
    print(paste("Archivo VCF a procesar:", vcf_file))
    
    # Leer el VCF
    print("Leyendo archivo VCF...")
    vcf <- readVcf(vcf_file, "hg38")
    print("Archivo VCF leído exitosamente.")
    
    # Variantes, info y genotipos
    print("Extrayendo variantes, info y genotipos...")
    variantes_df <- as.data.frame(rowRanges(vcf))
    info_df     <- as.data.frame(info(vcf))
    geno_list   <- geno(vcf)
    geno_df     <- do.call(cbind, lapply(geno_list, as.data.frame))
    print("Dataframes extraídos.")
    print("Primeras filas de info_df:")
    print(head(info_df))
    
    # --- Extraer ANN
    print("Procesando anotación ANN...")
    info_ann_1 <- info_df %>%
      dplyr::mutate(
        ANN_single = if_else(is.na(ANN), NA_character_, sub(",.*", "", ANN)),
        ANN_parts  = stringr::str_split(ANN_single, "\\|")
      ) %>%
      tidyr::unnest_wider(col = ANN_parts, names_sep = "_") %>%
      dplyr::rename_with(
        .cols = starts_with("ANN_parts_"),
        .fn = ~ stringr::str_replace(.x, "ANN_parts_", "ANN_")
      ) %>%   dplyr::select(-ANN_single, -ANN)
    
    info_con_ann_df <- info_ann_1
    print("ANN procesado.")
    print(head(info_con_ann_df))
    
    # --- Unir variantes + info + ANN + geno
    print("Uniendo variantes, info, ANN y genotipos...")
    final_vcf_df <- bind_cols(variantes_df, info_df, info_con_ann_df, geno_df)
    final_vcf_df$ALT <- lapply(final_vcf_df$ALT, function(x)
      as.character(x[[1]]))
    final_vcf_df <- as.data.frame(lapply(final_vcf_df, as.character))
    print("Unión de dataframes realizada.")
    
    # --- Eliminar columnas específicas
    print("Eliminando columnas no deseadas (1ra ronda)...")
    remove_columns <- c(
      "strand",
      "FILTER",
      "BaseQRankSum...14",
      "ExcessHet...17",
      "FS...18",
      "MLEAC...20",
      "MLEAF...21",
      "MQ...22",
      "MQRankSum...23",
      "NEGATIVE_TRAIN_SITE...24",
      "POSITIVE_TRAIN_SITE...25",
      "RAW_MQandDP...27",
      "ReadPosRankSum...28",
      "SOR...29",
      "VQSLOD...30",
      "culprit...31",
      "ANN",
      "HOM"
    )
    busca_exoma <- function(terms, cols)
      unique(unlist(sapply(terms, function(x)
        grep(x, cols, ignore.case = TRUE))))
    final_vcf_df <- final_vcf_df[, -busca_exoma(remove_columns, colnames(final_vcf_df)), drop = FALSE]
    print("Columnas eliminadas (1ra ronda).")
    
    # --- Columnas dbNSFP y frecuencias
    print("Extrayendo y limpiando columnas dbNSFP y de frecuencias...")
    cols_dbnsfp <- c(
      "dbNSFP_rs_dbSNP",
      "dbNSFP_clinvar_OMIM_id",
      "dbNSFP_clinvar_MedGen_id",
      "dbNSFP_HGVSc_snpEff",
      "dbNSFP_HGVSp_snpEff",
      "dbNSFP_clinvar_Orphanet_id",
      "dbNSFP_Reliability_index",
      "dbNSFP_AlphaMissense_pred",
      "dbNSFP_SIFT_pred",
      "dbNSFP_Polyphen2_HVAR_pred",
      "dbNSFP_MutationTaster_pred",
      "dbNSFP_Polyphen2_HDIV_pred",
      "dbNSFP_CADD_phred",
      "dbNSFP_Uniprot_acc",
      "dbNSFP_Interpro_domain"
    )
    cols_freqs  <- c(
      "dbNSFP_1000Gp3_SAS_AF",
      "dbNSFP_1000Gp3_AFR_AF",
      "dbNSFP_1000Gp3_EUR_AF",
      "dbNSFP_1000Gp3_EAS_AF",
      "dbNSFP_1000Gp3_AF",
      "dbNSFP_1000Gp3_AMR_AF"
    )
    clean_dbnsfp_cols <- function(df, cols_target, prefix = "dbNSFP_") {
      colnames_base <- sub("\\.\\.\\.[0-9]+$", "", names(df))
      df_unique <- df[, !duplicated(colnames_base)]
      names(df_unique) <- sub("\\.\\.\\.[0-9]+$", "", names(df_unique))
      df_final <- df_unique %>% dplyr::select(any_of(cols_target))
      names(df_final) <- sub(paste0("^", prefix), "", names(df_final))
      return(df_final)
    }
    df_dbnsfp <- clean_dbnsfp_cols(final_vcf_df[, busca_exoma(cols_dbnsfp, colnames(final_vcf_df)), drop = FALSE], cols_dbnsfp)
    
    # Extraer y limpiar frecuencias
    replace_dot_all_cases <- function(df) {
      df[] <- lapply(df, as.character)
      process_cell <- function(cell) {
        if (is.na(cell))
          return(NA_character_)
        cell_trim <- trimws(cell)
        if (grepl("^c\\(.*\\)$", cell_trim)) {
          content <- gsub("^c\\((.*)\\)$", "\\1", cell_trim)
          elems <- strsplit(content, ",")[[1]]
          elems <- trimws(gsub("^['\"]?|['\"]?$", "", elems))
          elems_clean <- elems[!grepl("^\\.*\\s*\\.\\s*\\.*$", elems)]
          if (length(elems_clean) == 0)
            return(NA_character_)
          return(paste(elems_clean, collapse = ","))
        }
        if (grepl("^\\.*\\s*\\.\\s*\\.*$", cell_trim))
          return(NA_character_)
        return(cell_trim)
      }
      df[] <- lapply(df, function(col)
        vapply(col, process_cell, FUN.VALUE = character(1)))
      df
    }
    df_freqs <- final_vcf_df[, busca_exoma(cols_freqs, colnames(final_vcf_df)), drop = FALSE]
    df_freqs <- df_freqs[, -grep("_AC", colnames(df_freqs)), drop = FALSE]
    df_freqs <- replace_dot_all_cases(df_freqs)
    
    mean_freqs_by_row <- function(df) {
      parse_cell <- function(cell) {
        if (is.na(cell))
          return(NA_real_)
        cell <- as.character(cell)
        if (grepl(",", cell)) {
          nums <- as.numeric(unlist(strsplit(cell, ",")))
          if (all(is.na(nums)))
            return(NA_real_)
          return(mean(nums, na.rm = TRUE))
        }
        as.numeric(cell)
      }
      df_num <- as.data.frame(lapply(df, function(col)
        vapply(col, parse_cell, numeric(1))))
      rowMeans(df_num, na.rm = TRUE)
    }
    df_freqs <- mean_freqs_by_row(df_freqs)
    print("Columnas dbNSFP y frecuencias listas.")
    
    # --- Quitar columnas intermedias
    print("Eliminando columnas intermedias y duplicadas (2da ronda)...")
    temp_Df <- final_vcf_df[, -c(busca_exoma(cols_freqs, colnames(final_vcf_df)),
                                 busca_exoma(cols_dbnsfp, colnames(final_vcf_df)))]
    other_columns <- c(
      "END...16",
      "InbreedingCoeff...19",
      "SNP...36",
      "INS...38",
      "DEL...39",
      "MC...62",
      "AF_EXAC...68",
      "AF_ESP...70",
      "AF_TGP...75",
      "GWASCAT_TRAIT...107"
    )
    temp_Df <- temp_Df[, !colnames(temp_Df) %in% other_columns, drop = FALSE]
    temp_Df <- temp_Df[, -grep("^dbNSFP", colnames(temp_Df)), drop = FALSE]
    temp_Df <- temp_Df[, -grep("GWASCAT", colnames(temp_Df)), drop = FALSE]
    print("Columnas intermedias y duplicadas eliminadas.")
    
    # --- Renombrar y eliminar duplicados con sufijos
    print("Renombrando columnas y eliminando sufijos repetidos...")
    remove_col_suffix_duplicates_and_rename <- function(df) {
      original_names <- names(df)
      is_x1 <- grepl("^X1", original_names)
      df_x1 <- df[, is_x1, drop = FALSE]
      df_rest <- df[, !is_x1, drop = FALSE]
      base_names <- sub("\\.\\.\\.[0-9]+$", "", names(df_rest))
      duplicated_names <- duplicated(base_names) |
        duplicated(base_names, fromLast = TRUE)
      new_names <- base_names
      for (name in unique(base_names[duplicated_names])) {
        indices <- which(base_names == name)
        for (i in seq_along(indices))
          new_names[indices[i]] <- paste0(name, "_", i)
      }
      names(df_rest) <- new_names
      df_final <- cbind(df_x1, df_rest)
      df_final
    }
    temp_Df.tmp <- remove_col_suffix_duplicates_and_rename(temp_Df)
    print("Columnas renombradas.")
    
    # --- Pasar vectores tipo c(1,2) a "1|2"
    print("Procesando vectores a formato pipe (|)...")
    process_pipe_vector <- function(df) {
      df[] <- lapply(df, as.character)
      process_cell <- function(cell) {
        if (is.na(cell) || cell == "character(0)")
          return(NA_character_)
        if (grepl("^c\\(", cell)) {
          content <- gsub("^c\\((.*)\\)$", "\\1", cell)
          items <- strsplit(content, ",")[[1]]
          items <- trimws(gsub("^['\"]?|['\"]?$", "", items))
          return(paste(items, collapse = "|"))
        }
        cell
      }
      df[] <- lapply(df, function(col)
        vapply(col, process_cell, character(1)))
      df
    }
    df_out <- process_pipe_vector(temp_Df.tmp)
    print("Vectores convertidos.")
    
    # --- Quitar columnas finales
    print("Eliminando columnas finales (3ra ronda)...")
    colnames_to_remove <- c(
      "BaseQRankSum",
      "ExcessHet",
      "END",
      "FS",
      "InbreedingCoeff",
      "MLEAC",
      "MLEAF",
      "MQ",
      "MQRankSum",
      "NEGATIVE_TRAIN_SITE",
      "POSITIVE_TRAIN_SITE",
      "RAW_MQandDP",
      "ReadPosRankSum",
      "SOR",
      "VQSLOD",
      "culprit",
      "SNP",
      "INS",
      "DEL",
      "MC",
      "AF_EXAC",
      "AF_ESP",
      "AF_TGP"
    )
    df_out <- df_out[, !colnames(df_out) %in% colnames_to_remove, drop = FALSE]
    equis <- df_out[, grep("^X1", colnames(df_out)), drop = FALSE]
    df_to_modify <- df_out[, -grep("^X1", colnames(df_out)), drop = FALSE]
    print("Columnas finales eliminadas.")
    
    # --- Conservar columna con más info entre duplicadas
    print("Conservando columnas con mayor información (si hay duplicadas)...")
    conservar_mas_informacion <- function(df) {
      nombres_base <- sub("_[0-9]+$", "", names(df))
      cols_a_conservar <- character(0)
      for (nombre in unique(nombres_base)) {
        idx <- which(nombres_base == nombre)
        if (length(idx) == 1) {
          cols_a_conservar <- c(cols_a_conservar, names(df)[idx])
        } else {
          n_info <- sapply(idx, function(i)
            sum(!is.na(df[[i]]) & df[[i]] != ""))
          col_elegida <- names(df)[idx[which.max(n_info)]]
          cols_a_conservar <- c(cols_a_conservar, col_elegida)
        }
      }
      df[, cols_a_conservar, drop = FALSE]
    }
    tmp <- conservar_mas_informacion(df_to_modify)
    geno_df_ <- geno_df[, -grep("SB", colnames(geno_df)), drop = FALSE]
    colnames(geno_df_) <- c("GT",
                            "AD",
                            "DP",
                            "GQ",
                            "MIN_DP",
                            "PGT",
                            "PID",
                            "PL",
                            "PS",
                            "RGQ")
    final <- bind_cols(tmp, geno_df_)
    final <- final[, !grepl("SB\\.1\\.", names(final)), drop = FALSE]
    final$freqs <- df_freqs
    final <- bind_cols(final, df_dbnsfp)
    final <- final[, -grep("^RS_1$", colnames(final)), drop = FALSE]
    print("Columnas con mayor información conservadas.")
    
    # --- ANN columnas
    print("Procesando columnas ANN...1")
    ann_cols <- info_con_ann_df[, grep("^ANN_\\d+$", names(info_con_ann_df), value = TRUE), drop = FALSE]
    print("Procesando columnas ANN...2")
    colnames(ann_cols) <- c(
      "alterno_quitar",
      "effect",
      "impact",
      "gene_name",
      "gene_name_quitar",
      "effect_quitar",
      "annotation_id",
      "gene_biotype",
      "exon_intron_rank",
      "nt_change",
      "aa_change",
      "cDNA_position.cDNA_len",
      "nt_position",
      "aa_position",
      "distance_to_feature",
      "errors"
    )[1:ncol(ann_cols)]
    
    print("Procesando columnas ANN... 3")
    ann_cols <- ann_cols[, !grepl("quitar", colnames(ann_cols)), drop = FALSE]
    final <- bind_cols(final, ann_cols)
    print("Columnas ANN procesadas.")
    
    # --- Nombres y orden final
    print("Reordenando y renombrando columnas finales...")
    colnames(final)[1:3] <- c("CHROM", "START", "END")
    if ("HET_1" %in% colnames(final))
      final$HET_1 <- ifelse(final$HET_1 == T, "Het", "Hom")
    quitar_col_repetidas_y_sufijo <- function(df) {
      base_names <- sub("_[0-9]+$", "", names(df))
      cols_a_conservar <- c()
      ya_vistos <- list()
      for (i in seq_along(df)) {
        nombre_base <- base_names[i]
        contenido   <- df[[i]]
        es_repetida <- FALSE
        if (!is.null(ya_vistos[[nombre_base]])) {
          for (col_existente in ya_vistos[[nombre_base]]) {
            if (isTRUE(all.equal(contenido, col_existente))) {
              es_repetida <- TRUE
              break
            }
          }
        }
        if (!es_repetida) {
          cols_a_conservar <- c(cols_a_conservar, i)
          ya_vistos[[nombre_base]] <- c(ya_vistos[[nombre_base]], list(contenido))
        }
      }
      df_nuevo <- df[, cols_a_conservar, drop = FALSE]
      names(df_nuevo) <- base_names[cols_a_conservar]
      df_nuevo
    }
    df_limpio <- quitar_col_repetidas_y_sufijo(final)
    
    orden_columnas <- c(
      "CHROM",
      "START",
      "END",
      "width",
      "paramRangeID",
      "REF",
      "ALT",
      "QUAL",
      "GT",
      "AD",
      "DP",
      "GQ",
      "MIN_DP",
      "PGT",
      "PID",
      "PL",
      "PS",
      "RGQ",
      "freqs",
      "AC_1",
      "AF_1",
      "AN_1",
      "DP_1",
      "QD_1",
      "LOF_1",
      "NMD_1",
      "VARTYPE_1",
      "MNP_1",
      "MIXED_1",
      "HET_1",
      "DBVARID_1",
      "SCISCV_1",
      "ALLELEID_1",
      "rs_dbSNP",
      "HGVSc_snpEff",
      "HGVSp_snpEff",
      "gene_name",
      "annotation_id",
      "gene_biotype",
      "exon_intron_rank",
      "nt_change",
      "aa_change",
      "cDNA_position.cDNA_len",
      "nt_position",
      "aa_position",
      "distance_to_feature",
      "CLNSIG_1",
      "CLNVCSO_1",
      "SCIDNINCL_1",
      "CLNREVSTAT_1",
      "ONCREVSTAT_1",
      "CLNDNINCL_1",
      "ONC_1",
      "CLNSIGSCV_1",
      "ORIGIN_1",
      "ONCINCL_1",
      "ONCDNINCL_1",
      "ONCDISDB_1",
      "SCIREVSTAT_1",
      "ONCDISDBINCL_1",
      "ONCSCV_1",
      "CLNDN_1",
      "ONCCONF_1",
      "CLNVC_1",
      "SCIDISDB_1",
      "CLNVI_1",
      "ONCDN_1",
      "CLNSIGINCL_1",
      "CLNDISDB_1",
      "GENEINFO_1",
      "CLNDISDBINCL_1",
      "CLNSIGCONF_1",
      "SCIDISDBINCL_1",
      "CLNHGVS_1",
      "SCIINCL_1",
      "SCIDN_1",
      "SCI_1",
      "clinvar_OMIM_id",
      "clinvar_MedGen_id",
      "clinvar_Orphanet_id",
      "Reliability_index",
      "AlphaMissense_pred",
      "SIFT_pred",
      "Polyphen2_HVAR_pred",
      "MutationTaster_pred",
      "Polyphen2_HDIV_pred",
      "CADD_phred",
      "Uniprot_acc",
      "Interpro_domain",
      "effect",
      "impact",
      "errors"
    )
    
    
    #df_limpio <- df_limpio[, orden_columnas[orden_columnas %in% names(df_limpio)], drop = FALSE]
    print("#### columnas df_limpio ########")
    print(colnames(df_limpio))
    print("#### ahora vemos cuales de orden_columnas estan en df_limpio ###")
    print(length(orden_columnas))
    print(orden_columnas)
    print(length(colnames(df_limpio)))
    w <- which(colnames(df_limpio) %in% orden_columnas)
    print(w)
    length(w)
    
    #df_limpio <- df_limpio[,w]
    
    print(colnames(df_limpio))
    print(orden_columnas)
    names(df_limpio) <- sub("_1$", "", names(df_limpio))
    df_limpio <- df_limpio[, names(df_limpio) != "DP.1", drop = FALSE]
    print("Columnas finales reordenadas y renombradas.")
    
    # --- Limpieza y formateo final
    print("Limpieza y filtrado final de errores y formatos...")
    df_limpio <- df_limpio %>%
      mutate(errors = gsub('["\'\\\\]', '', errors),
             errors = trimws(errors)) %>%
      filter(is.na(errors) |
               errors == "" | errors == "INFO_REALIGN_3_PRIME")
    df_limpio$AD <- sapply(df_limpio$AD, function(x)
      paste(x, collapse = "|"))
    df_limpio$PL <- sapply(df_limpio$PL, function(x)
      if (all(is.na(x)))
        NA_character_
      else
        paste(x, collapse = "|"))
    df_limpio <- bind_cols(codigo = muestra, df_limpio)
    
    print("Primeras filas del dataframe final antes de overlap:")
    print(head(df_limpio))
    
    # --- Overlap con exoma
    print("Aplicando overlap con regiones de exoma...")
    bd_list_ <- readRDS(db)
    
    exoma_pos <- df_limpio[, c("CHROM", "START", "END", "gene_name", "codigo")]
    colnames(exoma_pos)[5] <- "codigo"
    colnames(exoma_pos) <- c("Chr", "Start", "End", "Gene.refGene", "codigo")
    
    if (!any(names(bd_list_) == muestra)) {
      bd_list_[[muestra]] <- exoma_pos
      saveRDS(bd_list_, db)
    }
    
    cromosomas_ <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
    
    df_limpio <- obtener_exoma_overlap(
      bd_list = bd_list_,
      exoma_df = df_limpio,
      cromosomas = cromosomas_,
      muestra = muestra
    )
    
    print("Overlap con exoma realizado.")
    print("=== Proceso COMPLETO: process_vcf_to_table ===")
    
    hpo <- read.delim(hpo_file_, skip = 1, header = F)[, c(2, 4)]
    
    
    hpo_ <- aggregate(V4 ~ V2, hpo, FUN = paste, collapse = ";")
    colnames(hpo_) <- c("gene_name", "hpo")
    
    df_limpio <- left_join(df_limpio, hpo_, by = "gene_name")
    df_limpio <- buscar_herencia(df_limpio)
    
    archivo_final <- file.path(post_process_dir, "file_ready_analysis_optimized.csv")
    
    
    df_limpio  <- df_limpio  %>%
      mutate(
        REF_reads = as.numeric(sub(",.*", "", AD)),
        ALT_reads = as.numeric(sub(".*,", "", AD)),
        balance_alelico = ALT_reads / (REF_reads + ALT_reads)
      )
    
    
    orden_columnas <- c(
      "codigo",
      "CHROM",
      "POS",
      "END",
      "gene_name",
      "N",
      "samples",
      "REF",
      "ALT",
      "QUAL",
      "AC",
      "AF",
      "AN",
      "DP",
      "QD",
      "LOF",
      "NMD",
      "VARTYPE",
      "MNP",
      "MIXED",
      "HET",
      "balance_alelico",
      "freqs",
      "rs_dbSNP",
      "herencia",
      "GT",
      "AD",
      "GQ",
      "MIN_DP",
      "PGT",
      "PID",
      "PL",
      "PS",
      "RGQ",
      "clinvar_OMIM_id",
      "clinvar_MedGen_id",
      "HGVSc_snpEff",
      "Interpro_domain",
      "effect",
      "impact",
      "annotation_id",
      "gene_biotype",
      "exon_intron_rank",
      "nt_change",
      "aa_change",
      "cDNA_position.cDNA_len",
      "nt_position",
      "aa_position",
      "distance_to_feature",
      "errors",
      "paramRangeID",
      "DBVARID",
      "SCISCV",
      "ALLELEID",
      "CLNSIG",
      "CLNVCSO",
      "SCIDNINCL",
      "CLNREVSTAT",
      "ONCREVSTAT",
      "CLNDNINCL",
      "ONC",
      "CLNSIGSCV",
      "ORIGIN",
      "ONCINCL",
      "ONCDNINCL",
      "ONCDISDB",
      "SCIREVSTAT",
      "ONCDISDBINCL",
      "ONCSCV",
      "CLNDN",
      "ONCCONF",
      "CLNVC",
      "SCIDISDB",
      "CLNVI",
      "ONCDN",
      "CLNSIGINCL",
      "CLNDISDB",
      "GENEINFO",
      "CLNDISDBINCL",
      "CLNSIGCONF",
      "SCIDISDBINCL",
      "CLNHGVS",
      "SCIINCL",
      "SCIDN",
      "SCI",
      "HGVSp_snpEff",
      "clinvar_Orphanet_id",
      "Reliability_index",
      "AlphaMissense_pred",
      "SIFT_pred",
      "Polyphen2_HVAR_pred",
      "MutationTaster_pred",
      "Polyphen2_HDIV_pred",
      "CADD_phred",
      "Uniprot_acc",
      "hpo"
    )
    
    df_limpio <- df_limpio[, orden_columnas]
    colnames(df_limpio)[grep("^HET$", colnames(df_limpio))] <- "CIGOSIDAD"
    
    
    write.csv(
      df_limpio,
      file.path(
        output_dir,
        "post_process_results",
        "file_ready_analysis_optimized.csv"
      )
    )
    
    unique_variants <- df_limpio[df_limpio$N == 1, ]
    acmg <- read.delim("./acmg.txt", header = T)[, 1]
    
    acmg_v3_3_genes <- c(
      # Cáncer y predisposición tumoral
      "APC",
      "RET",
      "BRCA1",
      "BRCA2",
      "PALB2",
      "SDHD",
      "SDHAF2",
      "SDHC",
      "SDHB",
      "MAX",
      "TMEM127",
      "BMPR1A",
      "SMAD4",
      "TP53",
      "MLH1",
      "MSH2",
      "MSH6",
      "PMS2",
      "MEN1",
      "MUTYH",
      "NF2",
      "STK11",
      "PTEN",
      "RB1",
      "TSC1",
      "TSC2",
      "VHL",
      "WT1",
      
      # Cardiovasculares
      "FBN1",
      "TGFBR1",
      "TGFBR2",
      "SMAD3",
      "ACTA2",
      "MYH11",
      "PKP2",
      "DSP",
      "DSC2",
      "TMEM43",
      "DSG2",
      "RYR2",
      "CASQ2",
      "TRDN",
      "TNNT2",
      "LMNA",
      "FLNC",
      "TTN",
      "BAG3",
      "DES",
      "RBM20",
      "TNNC1",
      "PLN",
      "COL3A1",
      "LDLR",
      "APOB",
      "PCSK9",
      "MYH7",
      "MYBPC3",
      "TNNI3",
      "TPM1",
      "MYL3",
      "ACTC1",
      "PRKAG2",
      "MYL2",
      "KCNQ1",
      "KCNH2",
      "SCN5A",
      "CALM1",
      "CALM2",
      "CALM3",
      
      # Errores innatos del metabolismo
      "BTD",
      "CYP27A1",
      "GLA",
      "HFE",
      "OTC",
      "GAA",
      "ABCD1",
      "ACVRL1",
      "ENG",
      "TTR",
      "RYR1",
      "CACNA1S",
      "HNF1A",
      "RPE65",
      "ATP7B"
    )
    acmg_FINAL <- unique(c(acmg_v3_3_genes, acmg))
    
    
    unicas <- df_limpio[which(df_limpio$N == 1), ]
    
    acmg <- read.delim("./acmg.txt")[, 1]
    
    secondary <- df_limpio[which(df_limpio$gene_name %in% acmg_FINAL), ]
    
    secondary_unicas <- secondary[which(secondary$N == 1), ]
    
    write.csv(df_limpio, archivo_final)
    archivo_final_unicas <- file.path(post_process_dir,
                                      "file_ready_analysis_optimized_UNICAS.csv")
    archivo_final_secundarias <- file.path(post_process_dir, "ACMG_ALL.csv")
    archivo_final_sec_unicas <- file.path(post_process_dir, "ACMG_UNICAS.csv")
    
    write.csv(unicas, archivo_final_unicas)
    write.csv(secondary, archivo_final_secundarias)
    write.csv(secondary_unicas, archivo_final_sec_unicas)
    
  } else{
    print("Ya esta procesado")
  }
  
  #return(df_limpio)
}




compute_depth <- function(folder_fasta,
                          fastq_dir,
                          output_dir) {
  
  ## =========================
  ## 1) Nombre base de la muestra (CANÓNICO)
  ## =========================
  output_file_name <- get_sample_name(fastq_dir)
  
  ## =========================
  ## 2) Directorio de mapeo
  ## =========================
  mapping_output_dir <- file.path(output_dir, "mapping_output")
  
  ## =========================
  ## 3) BAM final (BQSR aplicado)
  ## =========================
  bam_file <- file.path(
    mapping_output_dir,
    paste0(output_file_name, ".sorted.mark_dup_bqsr.bam")
  )
  
  if (!file.exists(bam_file)) {
    stop("No existe el BAM para compute_depth: ", bam_file)
  }
  
  ## =========================
  ## 4) Directorio de salida
  ## =========================
  dir_coverage <- file.path(output_dir, "coverage_and_stats")
  dir.create(dir_coverage, recursive = TRUE, showWarnings = FALSE)
  
  outfile_coverage <- file.path(dir_coverage, "coverage.txt")
  
  ## =========================
  ## 5) BED y FASTA (NO hardcodeado)
  ## =========================
  bed_file  <- "./MGI_Exome_Capture_V5.bed"
  fasta_ref <- fn_exists_fasta(folder_fasta)
  
  if (!file.exists(bed_file)) {
    stop("No existe el BED de captura: ", bed_file)
  }
  
  ## =========================
  ## 6) Ejecutar cobertura
  ## =========================
  comando <- paste(
    "./compute_depthV2.sh",
    shQuote(bam_file),
    shQuote(bed_file),
    shQuote(fasta_ref),
    shQuote(outfile_coverage)
  )
  
  if (!file.exists(outfile_coverage)) {
    message("Computando cobertura...")
    print(comando)
    system(comando)
  } else {
    message("La cobertura ya está computada")
  }
}


compute_stats <- function(fastq_dir, output_dir, muestra) {
  dir_coverage <- file.path(output_dir, "coverage_and_stats")
  if (!file.exists(file.path(dir_coverage, "stats.csv"))) {
    dir_coverage <- file.path(output_dir, "coverage_and_stats")
    cov_file <- file.path(dir_coverage, "coverage.txt")
    process_dir <- file.path(output_dir, "post_process_results")
    exoma_file <- file.path(process_dir, "file_ready_analysis_optimized.csv")
    
    cov_data <- read.delim(cov_file, header = F)
    exoma <- read.csv(exoma_file, na.strings = ".")
    
    cromosomas <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
    exoma <- exoma[exoma$CHROM %in% cromosomas, ]
    df_hist <- cov_data[cov_data$V1 != "all", ]
    
    colnames(df_hist) <- c("chr",
                           "start",
                           "end",
                           "depth",
                           "bases_at_depth",
                           "region_length",
                           "fraction")
    df_hist <- df_hist[df_hist$depth > 0, ]
    # Calcular cobertura media por región
    
    # Filtrar sólo las filas "all"
    
    
    mean_coverage <- mean(mean(df_hist$depth))
    
    mean_coverage20 <- round(mean(df_hist[df_hist$depth > 20, "depth"], na.rm =
                                    T), 2)
    
    exoma_genes_variants <- unique(exoma[, c("CHROM", "POS", "END", "DP", "gene_name")])
    
    
    exoma_ <- exoma[-grep("-", exoma$gene_name), ]
    
    exoma_2  <- exoma_ %>% group_by(gene_name) %>% summarise(mean_dp = mean(DP))
    
    genes_totales <- length(unique(exoma_2$gene_name))
    
    genes_20 <- length(exoma_2[which(exoma_2$mean_dp > 20), ]$gene_name)
    
    genes_por_20 <- round(100 * genes_20 / genes_totales, 2)
    
    
    genes_mayor_media <- length(exoma_2[which(exoma_2$mean_dp > mean_coverage), "gene_name"]$gene_name)
    
    genes_por_media <- round(100 * genes_mayor_media / genes_totales, 2)
    
    
    variantes_totales <- dim(exoma_)[1]
    
    variantes_20 <- dim((exoma_[which(exoma_$DP > 20), ]))[1]
    variantes_20_x <- round(100 * variantes_20 / variantes_totales, 2)
    
    variantes_media <- dim((exoma_[which(exoma_$DP > mean_coverage), ]))[1]
    variantes_media_x <- round(100 * variantes_media / variantes_totales, 2)
    
    res <- as.data.frame(
      c(
        mean_coverage = mean_coverage,
        mean_coverage20 = mean_coverage20,
        genes_totales = genes_totales,
        genes_20 = genes_20,
        genes_por_20 = genes_por_20,
        genes_mayor_media = genes_mayor_media,
        genes_por_media = genes_por_media,
        variantes_totales = variantes_totales,
        variantes_20 = variantes_20,
        variantes_20_x = variantes_20_x,
        variantes_media = variantes_media,
        variantes_media_x
      )
    )
    res$Descripcion <- c(
      "Cobertura media",
      "Cobertura media > 20X",
      "Genes totales",
      "Genes > 20 X",
      "Genes % > 20X",
      "Genes > media",
      "Genes % > media",
      "Variantes totales",
      "Variantes > 20X",
      "Variantes % >20X",
      "Variantes > media X",
      "Variantes % > media X"
    )
    res <- res[, c(2, 1)]
    colnames(res) <- c("Descripcion", "Stats")
    
    
    gcov <- cov_data[cov_data[, 1] == 'all', ]
    
    ###
    longitud <- 20:200
    datos.pre <-
      data.frame(
        X = gcov[longitud, 2],
        Y = 100 * (1 -  cumsum(gcov[longitud, 5])),
        Z =  100 * gcov[longitud, 5],
        relleno = gcov[longitud, 1]
      )
    datos.pre$relleno <- as.factor(datos.pre$relleno)
    p1 <-
      ggplot(data = datos.pre, aes(X, Y, fill = relleno)) + geom_line(color =
                                                                        "steelblue", linewidth = 2) + xlab("Profundidad de Cobertura") + ylab("Porcentaje de la region >= Profunidad") + theme(legend.position = "none") +
      theme_classic()
    p2 <-
      ggplot(data = datos.pre, aes(X, Z, fill = relleno)) + geom_col() + scale_fill_discrete(type =
                                                                                               "steelblue") + xlab("Profundidad de Cobertura") + ylab("Porcentaje de la region") + theme(legend.position = "none") + theme_classic()
    
    p3 <- ggpubr::ggarrange(p1, p2, ncol = 2)
    p3
    figura_file_name <- file.path(dir_coverage, "cobertura.jpeg")
    
    p4 <-
      ggpubr::annotate_figure(p3,
                              top = paste(
                                muestra,
                                "Profundidad media de cobertura:",
                                mean_coverage,
                                "X"
                              ))
    ggsave(filename = figura_file_name, plot = p4)
    res_file <- file.path(dir_coverage, "stats.csv")
    write.csv(res, res_file)
    
    blah <- paste(
      "<p>Esta muestra se ha estudiado por el metodo de secuenciacion masiva en paralelo del exoma completo. Se analizaron",
      res[3, 2],
      "genes. La sensibilidad y la especificidad del metodo son superiores al 98% (SNV< 20 bp INDELS). El porcentaje de genes con una cobertura mayor a 20X es de",
      round(as.numeric(res[5, 2]), 2),
      "%. De todas las variantes identificadas, que son un total de",
      res[8, 2],
      ",",
      res[9, 2],
      "tienen una cobertura mayor a 20X, esto significa un",
      round(as.numeric(res[10, 2]), 2),
      "%.</p>"
    )
    
    res <- print(xtable(res), type = "html")
    
    fileconn <- "./aux1.html"
    writeLines(blah, fileconn)
    fileconn <- "./aux2.html"
    writeLines(res, fileconn)
    
    command <- paste("cat aux1.html aux2.html > ",
                     file.path(dir_coverage, "doc.html"))
    
    system(command)
    
    command <- paste(
      "pandoc --output",
      file.path(dir_coverage, "reporte.docx"),
      file.path(dir_coverage, "doc.html")
    )
    
    system(command)
    
    
    
    
    
  } else{
    print("ya se computaron las estadisticas")
  }
}






start <- Sys.time()
args <- commandArgs(trailingOnly = T)
args <- unlist(strsplit(args, " "))
muestra <- args
# for (muestra in muestras) {
yy <- sub(".*-(\\d{2})$", "\\1", muestra)

if (!grepl("^\\d{2}$", yy)) {
  stop("No se pudo inferir el año desde el nombre de la muestra: ", muestra)
}

year_full <- paste0("20", yy)

base_year_dir <- file.path("~/NAS_NGS", year_full)
muestra_dir   <- file.path(base_year_dir, muestra)
output_dir    <- file.path(muestra_dir, "output_dir")

dir.create(base_year_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(muestra_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir,    recursive = TRUE, showWarnings = FALSE)

hpo_file <- "~/NAS_NGS/datos_exomas/data_pipeline/genes_to_phenotype.txt"
fastq_dir <- file.path(muestra_dir, "fastqfiles")
print("aquiiiiiiiiiiii")
print(fastq_dir)
print(output_dir)

folder_fasta <-
  file.path("~/NAS_NGS/datos_exomas/datos_gatk/hg38")
folder_data_gatk <- file.path("~/NAS_NGS/datos_exomas/datos_gatk")
path_snpeff <- "~/tools/snpEff/"
bd_data <- "./bd.rds"

control_calidad(fastq_dir, output_dir)
bwamem(fastq_dir = fastq_dir, folder_fasta = folder_fasta, output_dir = output_dir)
markdups(output_dir = output_dir, fastq_dir = fastq_dir)
## creamos diccionario
create_dict(folder_fasta)
## anadimos reaad group
# creacion_readgroup(output_dir, fastq_dir)
## Recalibramos
base_recalibrator(folder_fasta, output_dir, folder_data_gatk, fastq_dir)
### aplicamos el recalibrado
applybqsr(folder_fasta, output_dir, fastq_dir)
## estadisticas del pieline bam
bam_statistics(folder_fasta, fastq_dir, output_dir)
## llamamos a las variantes
haplotype_caller(output_dir, folder_fasta, fastq_dir)
## Calculamos la probabilidad posterior del alelo referente
genotypeGVCF(folder_fasta, output_dir, fastq_dir)
## calculamos variant Recalibrator
# variantRecalibrator(fastq_dir, folder_fasta, folder_data_gatk, output_dir)
## apply VQSR
# applyVQSR(folder_fasta, fastq_dir, output_dir)
## primer filtraje
variantFiltration(folder_fasta, output_dir, fastq_dir)
## preparamos el archivo listo para elanalisis
analysisReady(folder_fasta, output_dir, fastq_dir)
##
anotation(folder_fasta, path_snpeff, output_dir, fastq_dir)

# process_vcf_to_table(
#   folder_fasta = folder_fasta,
#   output_dir = output_dir,
#   fastq_dir = fastq_dir,
#   muestra = muestra,
#   db = bd_data,
#   hpo_file = hpo_file
# )
# 
# 
# compute_depth(
#   folder_fasta = folder_fasta,
#   fastq_dir = fastq_dir,
#   output_dir = output_dir
# )
# 
# compute_stats(fastq_dir, output_dir, muestra)
# 

# }
print(Sys.time() - start)
