

compute_depth <- function(folder_fasta, fastq_dir, output_dir) {
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
  bam_file <- file.path(mapping_output_dir,
                        paste0(output_file_name, ".sorted.mark_dup_bqsr.bam"))
  
  if (!file.exists(bam_file)) {
    stop("No existe el BAM para compute_depth: ", bam_file)
  }
  
  ## =========================
  ## 4) Directorio de salida
  ## =========================
  dir_coverage <- file.path(output_dir, "coverage_and_stats")
  dir.create(dir_coverage,
             recursive = TRUE,
             showWarnings = FALSE)
  
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
    samples   = bd_meta$paste_m,!!!exoma_meta
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

