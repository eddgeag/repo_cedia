

vcf_process <- function(vcf_file,
                        # ruta completa al VCF .vcf.gz
                        genes_horizon,
                        # ruta ./archivo_horizon.csv
                        base_lab,
                        # ruta ./bd.rds
                        hpo_file)
  # ruta hpo_file.txt)
{
  suppressPackageStartupMessages({
    library(VariantAnnotation)
    library(GenomicRanges)
    library(GenomeInfoDb)
    library(data.table)
  })
  ## ============================================================
  ## DERIVAR RUTAS DESDE input
  ## ============================================================

  
  vcf_file      <- normalizePath(vcf_file, mustWork = TRUE)
  genes_horizon <- normalizePath(genes_horizon, mustWork = TRUE)
  base_lab      <- normalizePath(base_lab, mustWork = TRUE)
  hpo_file      <- normalizePath(hpo_file, mustWork = TRUE)
  ## ============================================================
  ## DERIVAR RUTAS DESDE vcf_file
  ## ============================================================
  
  exome_code <- sub("_.*$", "", basename(vcf_file))
  
  path_parts <- strsplit(dirname(vcf_file), .Platform$file.sep)[[1]]
  year <- path_parts[which(grepl("^20[0-9]{2}$", path_parts))[1]]
  
  if (is.na(year))
    stop("No se pudo inferir el año desde la ruta del VCF")
  
  base_dir <- file.path(
    .Platform$file.sep,
    paste(path_parts[1:which(path_parts == exome_code)],
          collapse = .Platform$file.sep),
    "output_dir"
  )
  
  out_dir <- file.path(base_dir, "post_process_results")
  if (!dir.exists(out_dir))
    dir.create(out_dir, recursive = TRUE)
  
  message("Output dir: ", out_dir)
  
  ### falta por hacer
  ### limpiar archivos q no se necesitan
  ### computar stats
  ### computar cobertura
  
  
  
  ## ============================================================
  ## 1) CARGA VCF
  ## ============================================================
  ## ============================================================
  ## CARGA DE RECURSOS EXTERNOS
  ## ============================================================
  
  hpo_dt <- read.delim(hpo_file)
  genes_amcg_horizon <- read.csv(genes_horizon)[, 2]
  aux <- readRDS(base_lab)
  
  # aux$`DX042-25` <- NULL
  vcf <- readVcf(vcf_file, genome = "hg38")
  
  ## ============================================================
  ## 2) VARIANTES (1 fila = 1 variante)
  ## ============================================================
  rr <- rowRanges(vcf)
  
  dt_variants <- data.table(
    VAR_IDX = seq_len(length(rr)),
    # clave estable
    CHROM   = as.character(seqnames(rr)),
    START   = start(rr),
    END     = end(rr),
    WIDTH   = width(rr),
    REF     = as.character(ref(vcf)),
    ALT     = vapply(alt(vcf), function(x)
      if (length(x))
        as.character(x[[1]])
      else
        NA_character_, character(1)),
    QUAL    = as.numeric(qual(vcf)),
    FILTER = as.character(rr@elementMetadata$FILTER)
  )
  setkey(dt_variants, VAR_IDX)
  
  rm(rr)
  gc()
  
  message("Bloque VARIANTES listo")
  print(head(dt_variants))
  
  ## ============================================================
  ## 3) GENOTIPOS (single-sample; 1 fila = 1 variante)
  ## ============================================================
  geno_raw <- geno(vcf)
  
  dt_geno <- as.data.table(lapply(geno_raw, function(x) {
    if (is.matrix(x))
      as.character(x[, 1])
    else
      as.character(x)
  }))
  dt_geno[, VAR_IDX := .I]
  setcolorder(dt_geno, c("VAR_IDX", setdiff(names(dt_geno), "VAR_IDX")))
  setkey(dt_geno, VAR_IDX)
  
  message("Bloque GENOTIPOS listo")
  print(names(dt_geno))
  
  ## ============================================================
  ## 4) INFO (SIN ANN; 1 fila = 1 variante)
  ## ============================================================
  info_raw <- info(vcf)
  dt_info  <- as.data.table(info_raw)
  
  if ("ANN" %in% names(dt_info))
    dt_info[, ANN := NULL]
  
  dt_info[, VAR_IDX := .I]
  setcolorder(dt_info, c("VAR_IDX", setdiff(names(dt_info), "VAR_IDX")))
  setkey(dt_info, VAR_IDX)
  
  message("Bloque INFO listo (ANN excluido)")
  print(names(dt_info))
  
  ## ============================================================
  ## 5) ANN LONG (sin pérdida; 1 fila = 1 anotación)
  ## ============================================================
  ann_raw <- info(vcf)$ANN
  if (is.null(ann_raw))
    stop("ERROR CRÍTICO: INFO/ANN no existe en el VCF")
  
  dt_ann <- data.table(
    VAR_IDX = rep(seq_along(ann_raw), lengths(ann_raw)),
    ANN_RAW = unlist(ann_raw, use.names = FALSE)
  )
  dt_ann[, ANN_RAW := as.character(ANN_RAW)]
  setkey(dt_ann, VAR_IDX)
  
  message("Bloque ANN_RAW listo")
  print(head(dt_ann))
  
  ## ============================================================
  ## 6) ANN -> COLUMNAS DINÁMICAS (sin esquema fijo)
  ## ============================================================
  # OJO: usar fixed=TRUE porque el separador es literal '|'
  ann_split <- strsplit(dt_ann$ANN_RAW, "|", fixed = TRUE)
  max_fields <- if (length(ann_split) > 0L) {
    max(lengths(ann_split))
  } else {
    0L
  }
  
  if (max_fields == 0L) {
    stop("ERROR: ANN existe pero no tiene contenido parseable")
  }
  
  
  
  ann_mat <- t(vapply(ann_split, function(x) {
    length(x) <- max_fields
    x
  }, character(max_fields)))
  
  dt_ann_fields <- as.data.table(ann_mat)
  setnames(dt_ann_fields, paste0("ANN_F", seq_len(max_fields)))
  
  dt_ann_full <- cbind(dt_ann[, .(VAR_IDX)], dt_ann_fields)
  setkey(dt_ann_full, VAR_IDX)
  
  message("Bloque ANN expandido listo")
  print(head(dt_ann_full))
  
  ## ============================================================
  ## 7) RESUMEN POR VARIANTE (COLAPSO DE TRANSCRITOS)
  ##     - 1 fila = 1 variante
  ##     - transcritos delimitados por comas
  ## ============================================================
  # ANN_F7 suele ser el "feature_id"/transcript ID en snpEff ANN
  # (en tu ejemplo: NR_027055.1, NM_152486.3, etc.)
  
  tx_col <- "ANN_F7"
  if (!tx_col %in% names(dt_ann_full)){
    stop("ERROR: no existe ANN_F7; revisa max_fields/parseo")}
  
  vcf_ann_by_variant <- dt_ann_full[, .(
    TRANSCRIPTS_AFFECTED = paste(
      unique(get(tx_col)[get(tx_col) != "" & !is.na(get(tx_col))]),
      collapse = ","
    ),
    N_TRANSCRIPTS = uniqueN(
      get(tx_col)[get(tx_col) != "" & !is.na(get(tx_col))]
    ),
    
    EFFECTS = paste(
      unique(ANN_F2[ANN_F2 != "" & !is.na(ANN_F2)]),
      collapse = ","
    ),
    IMPACTS = paste(
      unique(ANN_F3[ANN_F3 != "" & !is.na(ANN_F3)]),
      collapse = ","
    ),
    GENES_AFFECTED = paste(
      unique(ANN_F4[ANN_F4 != "" & !is.na(ANN_F4)]),
      collapse = ","
    ),
    BIOTYPES = paste(
      unique(ANN_F8[ANN_F8 != "" & !is.na(ANN_F8)]),
      collapse = ","
    ),
    NT_CHANGE = paste(
      unique(ANN_F10[ANN_F10 != "" & !is.na(ANN_F10)]),
      collapse = ","
    ),
    AA_CHANGE = paste(
      unique(ANN_F11[ANN_F11 != "" & !is.na(ANN_F11)]),
      collapse = ","
    ),
    DISTANCE_TO_FEATURE = paste(
      unique(ANN_F15[ANN_F15 != "" & !is.na(ANN_F15)]),
      collapse = ","
    ),
    ALIGNMENT = paste(
      unique(ANN_F16[ANN_F16 != "" & !is.na(ANN_F16)]),
      collapse = ","
    )
  ), by = VAR_IDX]
  
  setkey(vcf_ann_by_variant, VAR_IDX)
  
  message("Resumen por variante listo")
  message("Bloque ANN colapsado por variante listo")
  print(head(vcf_ann_by_variant))
  
  ## ============================================================
  ## 8) ENSAMBLE FINAL (VARIANTE-CÉNTRICO; 1 fila = 1 variante)
  ## ============================================================
  vcf_dt <- merge(dt_variants,
                  vcf_ann_by_variant,
                  by = "VAR_IDX",
                  all.x = TRUE)
  vcf_dt <- merge(vcf_dt, dt_geno, by = "VAR_IDX", all.x = TRUE)
  vcf_dt <- merge(vcf_dt, dt_info, by = "VAR_IDX", all.x = TRUE)
  
  message("VCF FINAL (VARIANTE-CÉNTRICO) LISTO")
  message("Dimensiones: ", paste(dim(vcf_dt), collapse = " x "))
  print(head(vcf_dt))
  ## ============================================================
  ## 9) CHEQUEOS CRÍTICOS
  ## ============================================================
  stopifnot(
    nrow(vcf_dt) == nrow(dt_variants),
    "N_TRANSCRIPTS" %in% names(vcf_dt),
    "TRANSCRIPTS_AFFECTED" %in% names(vcf_dt)
  )
  
  
  message("Importación completada correctamente")
  
  ## ============================================================
  ## 10) LIMPIEZA (opcional)
  ## ============================================================
  rm(vcf, geno_raw, info_raw, ann_raw, dt_ann, dt_ann_full)
  gc()
  
  ## ============================================================
  ## ELIMINAR COLUMNAS CONSTANTES
  ## ============================================================
  
  # Detectar columnas constantes
  const_cols <- names(vcf_dt)[vapply(vcf_dt, function(col) {
    # convertir listas a carácter para evaluar unicidad
    if (is.list(col))
      col <- as.character(col)
    uniqueN(col, na.rm = FALSE) <= 1
  }, logical(1))]
  
  message("Columnas constantes detectadas: ", length(const_cols))
  print(const_cols)
  
  # Eliminar columnas constantes
  vcf_dt[, (const_cols) := NULL]
  
  message("Columnas constantes eliminadas")
  message("Dimensiones finales: ", paste(dim(vcf_dt), collapse = " x "))
  
  ## ============================================================
  ## CONVERTIR SOLO "." -> NA (EXACTO)
  ## ============================================================
  
  vcf_dt[, names(vcf_dt) := lapply(.SD, function(col) {
    if (is.list(col))
      col <- as.character(col)
    if (is.character(col))
      col[col == "."] <- NA_character_
    col
  })]
  
  ## ============================================================
  ## NORMALIZACIÓN DEFINITIVA A CHARACTER + LIMPIEZA
  ## ============================================================
  vcf_dt <- vcf_dt[, lapply(.SD, function(col) {
    if (is.list(col)) {
      col <- vapply(col, function(x)
        if (length(x) == 0)
          NA_character_
        else
          paste(x, collapse = ","), character(1))
    }
    col <- as.character(col)
    col[col == "."] <- NA_character_
    col[nchar(col) == 0] <- NA_character_
    col
  })]
  
  vcf_dt[vcf_dt == "character(0)"] <- NA_character_
  
  
  
  ## ============================================================
  ## LIMPIAR STRINGS TIPO c("a","b"),c("c") -> A,B | C
  ## ============================================================
  
  
  vcf_dt[, names(vcf_dt) := lapply(.SD, function(col) {
    if (!is.character(col))
      return(col)
    
    col <- vapply(col, function(x) {
      if (is.na(x))
        return(NA_character_)
      
      # Si no hay c( , devolver tal cual
      if (!grepl("c\\(", x))
        return(x)
      
      # Extraer TODOS los c(...)
      matches <- gregexpr("c\\(([^\\)]*)\\)", x, perl = TRUE)
      parts <- regmatches(x, matches)[[1]]
      
      if (length(parts) == 0)
        return(x)
      
      # Limpiar cada vector
      # Limpiar cada vector (CORRECTO)
      cleaned <- vapply(parts, function(p) {
        p <- sub("^c\\(", "", p)   # quitar c(
        p <- sub("\\)$", "", p)    # quitar )
        p <- gsub("\"", "", p)     # quitar comillas
        p <- gsub("\\s+", "", p)   # quitar espacios
        p
      }, character(1))
      
      
      # Unir vectores con " | "
      paste(cleaned, collapse = " | ")
      
    }, character(1))
    
    col
  })]   # <- paréntesis/corchete final correcto
  
  setnames(
    vcf_dt,
    old = c("DP.x", "DP.y"),
    new = c("DP_FORMAT", "DP_INFO"),
    skip_absent = TRUE
  )
  
  # conservar END.x como END
  setnames(vcf_dt, "END.x", "END", skip_absent = TRUE)
  
  
  vcf_dt[vcf_dt == ".,."] <- NA_character_
  
  
  ## ============================================================
  ## ELIMINAR SNP / INS / DEL (REDUNDANTES)
  ## ============================================================
  
  vcf_dt[, c("SNP", "INS", "DEL") := NULL]
  
  message("Columnas SNP, INS, DEL eliminadas")
  
  ## ============================================================
  ## CREAR COLUMNA ZYGOSITY Y ELIMINAR HOM/HET
  ## ============================================================
  
  vcf_dt[, ZYGOSITY := fifelse(HOM == "TRUE",
                               "HOM",
                               fifelse(HET == "TRUE", "HET", NA_character_))]
  
  vcf_dt[, c("HOM", "HET") := NULL]
  
  message("ZYGOSITY creada y HOM/HET eliminadas")
  
  
  
  ## ============================================================
  ## UNIFICAR RS IDs (RS numérico + dbNSFP_rs_dbSNP con rs)
  ## ============================================================
  
  vcf_dt[, RS_ID := {
    rs1 <- RS
    rs2 <- dbNSFP_rs_dbSNP
    
    ## Normalizar vacíos
    rs1[rs1 == "" |
          is.na(rs1)] <- NA_character_
    rs2[rs2 == "" |
          is.na(rs2)] <- NA_character_
    
    ## RS numérico -> prefijar "rs"
    rs1 <- ifelse(!is.na(rs1), paste0("rs", rs1), NA_character_)
    
    ## Unificación lógica
    fifelse(is.na(rs1) & is.na(rs2),
            NA_character_,
            fifelse(
              !is.na(rs1) & is.na(rs2),
              rs1,
              fifelse(
                is.na(rs1) &
                  !is.na(rs2),
                rs2,
                fifelse(rs1 == rs2, rs1, paste(rs1, rs2, sep = ","))
              )
            ))
  }]
  
  ## Eliminar columnas originales
  vcf_dt[, c("RS", "dbNSFP_rs_dbSNP") := NULL]
  
  message("RS_ID unificado (RS normalizado con 'rs' y dbNSFP integrado)")
  
  
  
  ## 1) Parseo seguro de AD
  ## AD esperado tipo "ref,alt" o "ref,alt1,alt2,..."
  vcf_dt[, c("AD_REF", "AD_ALT1", "DP_AD", "AB_ALT") := {
    ad <- AD
    
    # defaults
    ad_ref <- rep(NA_integer_, .N)
    ad_alt1 <- rep(NA_integer_, .N)
    dp_ad <- rep(NA_integer_, .N)
    ab <- rep(NA_real_, .N)
    
    ok <- !is.na(ad) &
      nzchar(ad) & grepl("^[0-9]+(,[0-9]+)+$", ad)
    if (any(ok)) {
      parts <- strsplit(ad[ok], ",", fixed = TRUE)
      
      refv  <- vapply(parts, function(x) as.integer(x[1]), integer(1))
      alt1v <- vapply(parts, function(x) as.integer(x[2]), integer(1))
      sumv  <- vapply(parts, function(x) sum(as.integer(x)), integer(1))
      
      ad_ref[ok]  <- refv
      ad_alt1[ok] <- alt1v
      dp_ad[ok]   <- sumv
      
      den <- refv + alt1v
      idx <- which(ok)[den > 0]
      
      ab[idx] <- alt1v[den > 0] / den[den > 0]
    }
    list(ad_ref, ad_alt1, dp_ad, ab)
  }]
  
  ## 2) DP final para filtros (preferir DP_AD; fallback DP_FORMAT)
  vcf_dt[, DP_FOR_AB := fifelse(!is.na(DP_AD), DP_AD, as.integer(DP_FORMAT))]
  
  ## 3) Banderas útiles
  vcf_dt[, `:=`(
    FLAG_NO_AD      = is.na(DP_AD),
    # no se pudo usar AD
    FLAG_LOW_DP     = is.na(DP_FOR_AB) |
      DP_FOR_AB < 10,
    FLAG_AB_NA      = is.na(AB_ALT),
    FLAG_DP_MISMATCH = {
      dpF <- suppressWarnings(as.integer(DP_FORMAT))
      ! is.na(DP_AD) &
        !is.na(dpF) & abs(dpF - DP_AD) >= 10
    }
  )]
  
  ## 4) Interpretación simple por genotipo (umbral configurable)
  ## Nota: en exoma recomiendo rangos amplios; luego los ajustas por dataset.
  vcf_dt[, AB_EXPECTATION := fifelse(
    grepl("^0[\\/|]1$|^1[\\/|]0$", GT),
    "HET_expected_mid",
    fifelse(
      grepl("^1[\\/|]1$", GT),
      "HOM_ALT_expected_high",
      fifelse(grepl("^0[\\/|]0$", GT), "HOM_REF_expected_low", "OTHER")
    )
  )]
  
  vcf_dt[, AB_SUSPECT := fifelse(
    FLAG_LOW_DP | FLAG_AB_NA,
    NA,
    # no juzgar si no es confiable
    fifelse(
      AB_EXPECTATION == "HET_expected_mid" &
        (AB_ALT < 0.20 | AB_ALT > 0.80),
      TRUE,
      fifelse(
        AB_EXPECTATION == "HOM_ALT_expected_high" & AB_ALT < 0.90,
        TRUE,
        fifelse(
          AB_EXPECTATION == "HOM_REF_expected_low" & AB_ALT > 0.10,
          TRUE,
          FALSE
        )
      )
    )
  )]
  
  ## ================================
  ## 1) HET real (usando AB)
  ## ================================
  
  vcf_dt[, IS_HET_REAL := {
    is_het_gt <- grepl("^0[\\/|]1$|^1[\\/|]0$", GT)
    is_het_ab <- !is.na(AB_ALT) &
      AB_ALT >= 0.20 &
      AB_ALT <= 0.80 & DP_FOR_AB >= 10
    is_het_gt & is_het_ab
  }]
  
  
  ## ================================
  ## 2) Expandir GENES_AFFECTED
  ## ================================
  
  genes_long <- vcf_dt[IS_HET_REAL == TRUE &
                         !is.na(GENES_AFFECTED), .(VAR_IDX, GENE = unlist(strsplit(GENES_AFFECTED, ",", fixed = TRUE)), PID, PGT)]
  
  genes_long[, GENE := trimws(GENE)]
  
  
  ## ================================
  ## 3) Contar HET por gen
  ## ================================
  
  gene_counts <- genes_long[, .N, by = GENE]
  setnames(gene_counts, "N", "N_HET")
  
  
  ## ================================
  ## 4) Unir a variantes
  ## ================================
  
  vcf_dt <- merge(
    vcf_dt,
    gene_counts,
    by.x = "GENES_AFFECTED",
    by.y = "GENE",
    all.x = TRUE
  )
  
  vcf_dt[is.na(N_HET), N_HET := 0]
  
  
  ## ================================
  ## 5) Evaluar fase (si existe)
  ## ================================
  
  genes_long[, PHASE_BLOCK := paste(PID, PGT, sep = "|")]
  
  phase_info <- genes_long[, .(N_PHASES = uniqueN(PHASE_BLOCK)), by = GENE]
  
  
  vcf_dt <- merge(vcf_dt,
                  phase_info,
                  by.x = "GENES_AFFECTED",
                  by.y = "GENE",
                  all.x = TRUE)
  
  vcf_dt[is.na(N_PHASES), N_PHASES := 0]
  
  
  ## ================================
  ## 6) Clasificación final
  ## ================================
  
  vcf_dt[, COMPOUND_HET := fifelse(IS_HET_REAL == FALSE |
                                     N_HET < 2,
                                   "no",
                                   fifelse(
                                     N_PHASES >= 2,
                                     "confirmed",
                                     # dos haplotipos → verdadero compound het
                                     "possible"                    # sin fase → posible
                                   ))]
  
  rm(
    dt_variants,
    dt_geno,
    dt_info,
    dt_ann,
    dt_ann_fields,
    dt_ann_full,
    vcf_ann_by_variant,
    genes_long,
    gene_counts,
    phase_info
  )
  
  gc(verbose = TRUE)
  ## ============================================================
  ## CONSTRUIR FRECUENCIAS DE COHORTE DESDE bd.rds (aux)
  ##   + registrar DX042-25 desde vcf_dt
  ## ============================================================
  exome_code <- sub("_.*$", "", basename(vcf_file))   # DX042-25
  
  ## ============================================================
  ## 0) REGISTRAR DX042-25 EN aux USANDO vcf_dt
  ## ============================================================
  vcf_dt_ <- copy(vcf_dt)
  
  if (!exome_code %in% names(aux)) {
    message("Registrando exoma nuevo en aux: ", exome_code)
    
    dx <- vcf_dt_[, .(
      CHROM = CHROM,
      START = as.integer(START),
      END   = as.integer(END),
      Gene.refGene = GENES_AFFECTED,
      REF = REF,
      ALT = ALT
    )]
    
    dx$codigo <- exome_code
    aux[[exome_code]] <- as.data.frame(dx)
    
    saveRDS(aux, "../repo_cedia/bd.rds")
    message(exome_code, " registrado en bd.rds")
    
  } else {
    message(exome_code, " ya existe en aux")
  }
  
  
  extract_loci <- function(dt, code) {
    if (is.null(dt) || nrow(dt) == 0)
      return(NULL)
    
    dt <- as.data.table(dt)
    
    setnames(dt, intersect(names(dt), c("Chr", "chr", "CHROM")), "CHROM", skip_absent =
               TRUE)
    setnames(dt, intersect(names(dt), c("Start", "START")), "START", skip_absent =
               TRUE)
    setnames(dt, intersect(names(dt), c("End", "END")), "END", skip_absent =
               TRUE)
    
    data.table(LOCUS  = paste(dt$CHROM, dt$START, dt$END, sep = ":"),
               codigo = code)
  }
  
  all_loci <- rbindlist(lapply(names(aux), function(code) {
    extract_loci(aux[[code]], code)
  }), fill = TRUE)
  
  locus_freq <- all_loci[, .(
    N_EXOMES_LOCUS = uniqueN(codigo),
    EXOMES_LOCUS  = paste(sort(unique(codigo)), collapse = ",")
  ), by = LOCUS]
  
  
  
  vcf_dt_[, LOCUS := paste(CHROM, START, END, sep = ":")]
  
  vcf_dt_ <- merge(vcf_dt_, locus_freq, by = "LOCUS", all.x = TRUE)
  
  
  ## ============================================================
  ## FRECUENCIAS POBLACIONALES (ACMG)
  ## ============================================================
  
  N_TOTAL_EXOMES <- length(names(aux))
  toNum <- function(x)
    suppressWarnings(as.numeric(x))
  
  
  vcf_dt_[, `:=`(
    AF_gnomad = toNum(dbNSFP_gnomAD4.1_joint_POPMAX_AF),
    AF_exac   = toNum(AF_EXAC),
    AF_tgp    = toNum(AF_TGP),
    AF_esp    = toNum(AF_ESP)
  )]
  
  vcf_dt_[, POP_AF := pmax(AF_gnomad, AF_exac, AF_tgp, AF_esp, na.rm = TRUE)]
  vcf_dt_[is.infinite(POP_AF), POP_AF := NA_real_]
  
  vcf_dt_[, COHORT_AF := N_EXOMES_LOCUS / N_TOTAL_EXOMES]
  
  
  vcf_dt_[, `:=`(
    ACMG_PM2 = COHORT_AF < 0.01 & (is.na(POP_AF) | POP_AF < 0.001),
    ACMG_BS1 = POP_AF >= 0.05,
    ACMG_BA1 = POP_AF >= 0.1,
    ACMG_BS2 = COHORT_AF >= 0.05
  )]
  
  
  summary(vcf_dt_$COHORT_AF)
  summary(vcf_dt_$POP_AF)
  table(vcf_dt_$ACMG_PM2, useNA = "ifany")
  table(vcf_dt_$ACMG_BS1, useNA = "ifany")
  
  ## ============================================================
  ## CLNSIG_CLINVAR = unión limpia de CLNSIG + dbNSFP_clinvar_clnsig
  ## ============================================================
  
  ## ============================================================
  ## CLNSIG_CLINVAR = unión limpia de CLNSIG + dbNSFP_clinvar_clnsig
  ## Separadores válidos: |  y  /
  ## ============================================================
  
  clean_merge_clnsig <- function(x, y) {
    v <- c(x, y)
    
    ## eliminar NA, ".", "", "not_provided"
    v <- v[!is.na(v)]
    v <- v[v != ""]
    v <- v[v != "."]
    v <- v[!tolower(v) %in% c("not_provided")]
    
    if (length(v) == 0)
      return(NA_character_)
    
    ## convertir SOLO | y / a coma
    v <- gsub("[|/]", ",", v)
    
    ## separar
    v <- unlist(strsplit(v, ",", fixed = TRUE))
    
    ## limpiar
    v <- trimws(v)
    v <- v[v != ""]
    
    if (length(v) == 0)
      return(NA_character_)
    
    ## únicos por variante, respetando términos completos
    paste(sort(unique(v)), collapse = ",")
  }
  
  vcf_dt_[, CLNSIG_CLINVAR := mapply(clean_merge_clnsig, CLNSIG, dbNSFP_clinvar_clnsig)]
  
  ###=== clasificador acmg========
  
  
  
  ## ============================================================
  ## ACMG PVS1 / PP3 / BP4 — versión extendida con dbNSFP
  ## ============================================================
  ## ============================================================
  ## ACMG PVS1 / PP3 / BP4 — versión extendida con dbNSFP
  ## ============================================================
  
  toNum <- function(x)
    suppressWarnings(as.numeric(x))
  
  split_effects <- function(x) {
    unlist(strsplit(x, "[,&]"))
  }
  
  has_effect <- function(x, pattern) {
    sapply(x, function(z)
      any(grepl(pattern, split_effects(z))))
  }
  
  ## ---------- PVS1 (LoF real según snpEff)
  pvs1_terms <- c(
    "stop_gained",
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "start_lost",
    "stop_lost"
  )
  
  vcf_dt_[, ACMG_PVS1 :=
            has_effect(EFFECTS, paste(pvs1_terms, collapse = "|"))]
  
  ## ---------- Missense puro
  vcf_dt_[, IS_MISSENSE :=
            has_effect(EFFECTS, "^missense_variant$")]
  
  ## ============================================================
  ## PP3 — evidencia computacional de daño
  ## ============================================================
  
  vcf_dt_[, ACMG_PP3 :=
            IS_MISSENSE &
            (
              ## Predictores principales
              toNum(dbNSFP_REVEL_score)       >= 0.50 |
                toNum(dbNSFP_CADD_phred)        >= 20   |
                toNum(dbNSFP_MPC_score)         >= 2    |
                toNum(dbNSFP_ClinPred_score)    >= 0.5  |
                
                ## Clasificadores
                dbNSFP_MetaRNN_pred            == "D"  |
                dbNSFP_AlphaMissense_pred %chin% c("P", "LP")
            )]
  
  ## ============================================================
  ## BP4 — evidencia computacional de benignidad
  ## ============================================================
  
  vcf_dt_[, ACMG_BP4 :=
            IS_MISSENSE &
            (
              toNum(dbNSFP_REVEL_score)    < 0.15 |
                toNum(dbNSFP_CADD_phred)     < 10   |
                toNum(dbNSFP_MPC_score)      < 1    |
                toNum(dbNSFP_ClinPred_score) < 0.20 |
                
                dbNSFP_MetaRNN_pred          == "T" |
                dbNSFP_AlphaMissense_pred    == "B"
            )]
  
  ## ---------- resolver contradicciones
  vcf_dt_[ACMG_PP3 &
            ACMG_BP4, ACMG_BP4 := FALSE]
  
  
  ## Homocigoto raro → compatible con AR
  vcf_dt_[, ACMG_PM3 :=
            ZYGOSITY == "HOM" &
            COHORT_AF < 0.01]
  
  ## Het compuesta posible
  vcf_dt_[, ACMG_PM3_compound :=
            COMPOUND_HET %chin% c("possible", "confirmed")]
  
  ## PM3 efectivo: compound tiene prioridad
  vcf_dt_[, ACMG_PM3_EFFECTIVE :=
            fifelse(ACMG_PM3_compound, TRUE, ACMG_PM3)]
  
  ## PM3 simple NO puede coexistir
  vcf_dt_[ACMG_PM3_compound == TRUE, ACMG_PM3 := FALSE]
  
  vcf_dt_[, `:=`(
    ACMG_CLINVAR_P =
      grepl("Pathogenic", CLNSIG_CLINVAR, ignore.case = TRUE) &
      !grepl("Benign", CLNSIG_CLINVAR, ignore.case = TRUE),
    
    ACMG_CLINVAR_B =
      grepl("Benign", CLNSIG_CLINVAR, ignore.case = TRUE) &
      !grepl("Pathogenic", CLNSIG_CLINVAR, ignore.case = TRUE)
  )]
  
  
  vcf_dt_[, ACMG_PATH_SCORE :=
            2 * ACMG_PVS1 +
            1 * ACMG_PM2 +
            1 * ACMG_PM3_EFFECTIVE +
            1 * ACMG_PP3 +
            2 * ACMG_CLINVAR_P]
  
  vcf_dt_[, ACMG_BENIGN_SCORE :=
            2 * ACMG_BA1 +
            1 * ACMG_BS1 +
            1 * ACMG_BS2 +
            1 * ACMG_BP4 +
            2 * ACMG_CLINVAR_B]
  
  vcf_dt_[, ACMG_CLASS := fifelse(
    ACMG_BENIGN_SCORE >= 3,
    "Benign",
    fifelse(
      ACMG_PATH_SCORE >= 4,
      "Pathogenic",
      fifelse(
        ACMG_PATH_SCORE >= 2,
        "Likely_pathogenic",
        fifelse(ACMG_BENIGN_SCORE >= 2, "Likely_benign", "VUS")
      )
    )
  )]
  
  
  vcf_dt_[, ACMG_EVIDENCE := paste(
    c(
      ifelse(ACMG_PVS1, "PVS1", ""),
      ifelse(ACMG_PM2, "PM2", ""),
      ifelse(ACMG_PM3_EFFECTIVE & ACMG_PM3_compound, "PM3_compound", ""),
      ifelse(ACMG_PM3_EFFECTIVE & !ACMG_PM3_compound, "PM3", ""),
      ifelse(ACMG_PP3, "PP3", ""),
      ifelse(ACMG_BA1, "BA1", ""),
      ifelse(ACMG_BS1, "BS1", ""),
      ifelse(ACMG_BS2, "BS2", ""),
      ifelse(ACMG_BP4, "BP4", ""),
      ifelse(ACMG_CLINVAR_P, "ClinVar_P", ""),
      ifelse(ACMG_CLINVAR_B, "ClinVar_B", "")
    ),
    collapse = ";"
  ), by = seq_len(nrow(vcf_dt_))]
  
  vcf_dt_[, ACMG_EVIDENCE := gsub(";;+", ";", gsub("^;|;$", "", ACMG_EVIDENCE))]
  
  vcf_dt_[, ACMG_CLASS_FINAL := fifelse(
    ## Benigno fuerte
    ACMG_BENIGN_SCORE >= 3 &
      ACMG_PATH_SCORE == 0,
    "Benign",
    
    fifelse(
      ## Probablemente benigno
      ACMG_BENIGN_SCORE >= 2 &
        ACMG_PATH_SCORE <= 1,
      "Likely_benign",
      
      fifelse(
        ## Patogénico fuerte
        ACMG_PATH_SCORE >= 4 &
          ACMG_BENIGN_SCORE == 0,
        "Pathogenic",
        
        fifelse(
          ## Probablemente patogénico
          ACMG_PATH_SCORE >= 2 &
            ACMG_BENIGN_SCORE <= 1,
          "Likely_pathogenic",
          
          "VUS"
        )
      )
    )
  )]
  
  
  vcf_dt_[, ACMG_SUMMARY :=
            paste0(
              ACMG_CLASS_FINAL,
              " | P=",
              ACMG_PATH_SCORE,
              " B=",
              ACMG_BENIGN_SCORE,
              " | ",
              ACMG_EVIDENCE
            )]
  
  
  suppressPackageStartupMessages({
    library(data.table)
    library(brms)
    library(rstan)
    
  })
  ## helpers
  ## ============================================================
  toNum1 <- function(x)
    suppressWarnings(as.numeric(x))
  
  pick_col <- function(dt, candidates) {
    hit <- intersect(candidates, names(dt))
    if (length(hit) == 0)
      return(NA_character_)
    hit[1]
  }
  
  clean_dbnsfp_cat <- function(x, pathogenic_set, benign_set) {
    if (is.na(x) || x == "" || x == ".")
      return(NA_character_)
    v <- unlist(strsplit(x, ",", fixed = TRUE))
    v <- trimws(v)
    v <- v[v != "" & v != "."]
    if (length(v) == 0)
      return(NA_character_)
    if (any(v %in% pathogenic_set))
      return("P")
    if (any(v %in% benign_set))
      return("B")
    NA_character_
  }
  
  clean_numeric_max <- function(x) {
    if (is.na(x) || x == "" || x == ".")
      return(NA_real_)
    v <- unlist(strsplit(x, ",", fixed = TRUE))
    v <- trimws(v)
    v <- v[v != "" & v != "."]
    v <- suppressWarnings(as.numeric(v))
    v <- v[!is.na(v)]
    if (length(v) == 0)
      return(NA_real_)
    max(v)
  }
  
  ## ============================================================
  ## 1) GOLD BINARIO (P vs B)
  ## ============================================================
  build_gold_binary <- function(dt) {
    out <- copy(dt)
    out[, y := fifelse(
      grepl("pathogenic", CLNSIG_CLINVAR, ignore.case = TRUE) &
        !grepl("benign", CLNSIG_CLINVAR, ignore.case = TRUE),
      1L,
      fifelse(
        grepl("benign", CLNSIG_CLINVAR, ignore.case = TRUE) &
          !grepl("pathogenic", CLNSIG_CLINVAR, ignore.case = TRUE),
        0L,
        NA_integer_
      )
    )]
    out
  }
  
  ## ============================================================
  ## 2) NORMALIZAR PREDICTORES dbNSFP (robusto a nombres)
  ## ============================================================
  normalize_dbnsfp_predictors <- function(dt) {
    x <- copy(dt)
    
    ## --- columnas candidatas (porque dbNSFP a veces cambia naming)
    c_am   <- pick_col(x,
                       c(
                         "dbNSFP_AlphaMissense_pred",
                         "AlphaMissense_pred",
                         "AlphaMissense"
                       ))
    c_mrn  <- pick_col(x, c("dbNSFP_MetaRNN_pred", "MetaRNN_pred", "MetaRNN"))
    c_sift <- pick_col(x, c("dbNSFP_SIFT_pred", "SIFT_pred", "SIFT"))
    c_phd  <- pick_col(
      x,
      c(
        "dbNSFP_PolyPhen2_HDIV_pred",
        "dbNSFP_Polyphen2_HDIV_pred",
        "PolyPhen2_HDIV_pred"
      )
    )
    c_phv  <- pick_col(
      x,
      c(
        "dbNSFP_PolyPhen2_HVAR_pred",
        "dbNSFP_Polyphen2_HVAR_pred",
        "PolyPhen2_HVAR_pred"
      )
    )
    
    c_revel <- pick_col(x, c("dbNSFP_REVEL_score", "REVEL_score", "REVEL"))
    c_cadd  <- pick_col(x, c("dbNSFP_CADD_phred", "CADD_phred", "CADD"))
    c_mpc   <- pick_col(x, c("dbNSFP_MPC_score", "MPC_score", "dbNSFP_MPC", "MPC"))
    c_cp    <- pick_col(x,
                        c("dbNSFP_ClinPred_score", "ClinPred_score", "ClinPred"))
    
    ## --- categóricos
    x[, AlphaMissense_clean := if (!is.na(c_am))
      vapply(get(c_am),
             clean_dbnsfp_cat,
             c("P", "LP"),
             c("B", "LB"),
             FUN.VALUE = character(1))
      else
        NA_character_]
    
    x[, MetaRNN_clean := if (!is.na(c_mrn))
      vapply(get(c_mrn),
             clean_dbnsfp_cat,
             c("D"),
             c("T"),
             FUN.VALUE = character(1))
      else
        NA_character_]
    
    x[, SIFT_clean := if (!is.na(c_sift))
      vapply(
        get(c_sift),
        clean_dbnsfp_cat,
        c("D", "deleterious"),
        c("T", "tolerated"),
        FUN.VALUE = character(1)
      )
      else
        NA_character_]
    
    x[, PolyPhen2_HDIV_clean := if (!is.na(c_phd))
      vapply(get(c_phd),
             clean_dbnsfp_cat,
             c("D", "P"),
             c("B"),
             FUN.VALUE = character(1))
      else
        NA_character_]
    
    x[, PolyPhen2_HVAR_clean := if (!is.na(c_phv))
      vapply(get(c_phv),
             clean_dbnsfp_cat,
             c("D", "P"),
             c("B"),
             FUN.VALUE = character(1))
      else
        NA_character_]
    
    ## --- numéricos (max si multivalor)
    x[, REVEL_num    := if (!is.na(c_revel))
      vapply(get(c_revel), clean_numeric_max, numeric(1))
      else
        NA_real_]
    x[, CADD_num     := if (!is.na(c_cadd))
      vapply(get(c_cadd), clean_numeric_max, numeric(1))
      else
        NA_real_]
    x[, MPC_num      := if (!is.na(c_mpc))
      vapply(get(c_mpc), clean_numeric_max, numeric(1))
      else
        NA_real_]
    x[, ClinPred_num := if (!is.na(c_cp))
      vapply(get(c_cp), clean_numeric_max, numeric(1))
      else
        NA_real_]
    
    x
  }
  
  ## ============================================================
  ## 3) MODEL MATRIX (sin drop)
  ## ============================================================
  make_model_dt <- function(dt0) {
    d <- copy(dt0)
    
    ## ----------------
    ## numéricas
    ## ----------------
    num_cols <- c("REVEL_num",
                  "CADD_num",
                  "MPC_num",
                  "ClinPred_num",
                  "POP_AF",
                  "COHORT_AF")
    for (v in num_cols) {
      if (!v %in% names(d))
        d[, (v) := NA_real_]
      d[, (v) := as.numeric(get(v))]
    }
    
    d[, log10_POP_AF :=
        fifelse(!is.na(POP_AF) &
                  POP_AF > 0, log10(POP_AF), NA_real_)]
    d[, log10_COHORT_AF :=
        fifelse(!is.na(COHORT_AF) &
                  COHORT_AF > 0,
                log10(COHORT_AF),
                NA_real_)]
    
    ## ----------------
    ## categóricas ACMG-style (SIN _NA)
    ## ----------------
    cat_cols <- c(
      "AlphaMissense_clean",
      "MetaRNN_clean",
      "SIFT_clean",
      "PolyPhen2_HDIV_clean",
      "PolyPhen2_HVAR_clean"
    )
    
    for (v in cat_cols) {
      if (!v %in% names(d))
        d[, (v) := NA_character_]
      
      d[, paste0(v, "_P") := as.integer(get(v) == "P")]
      d[, paste0(v, "_B") := as.integer(get(v) == "B")]
    }
    
    ## ----------------
    ## missense
    ## ----------------
    if (!"IS_MISSENSE" %in% names(d))
      d[, IS_MISSENSE := NA]
    d[, IS_MISSENSE := as.integer(IS_MISSENSE)]
    
    ## ----------------
    ## columnas finales
    ## ----------------
    keep <- c(
      "y",
      "REVEL_num",
      "CADD_num",
      "MPC_num",
      "ClinPred_num",
      "log10_POP_AF",
      "log10_COHORT_AF",
      paste0(rep(cat_cols, each = 2), c("_P", "_B")),
      "IS_MISSENSE"
    )
    
    d[, ..keep]
  }
  
  ## ============================================================
  ## 4) split + scale
  ## ============================================================
  stratified_split <- function(d,
                               test_frac = 0.2,
                               seed = 1) {
    set.seed(seed)
    d0 <- d[!is.na(y)]
    idxP <- which(d0$y == 1L)
    idxB <- which(d0$y == 0L)
    test_idx <- c(sample(idxP, max(1, floor(
      length(idxP) * test_frac
    ))), sample(idxB, max(1, floor(
      length(idxB) * test_frac
    ))))
    list(train = d0[-test_idx], test = d0[test_idx])
  }
  
  scale_numeric <- function(train, test, vars) {
    mu <- train[, lapply(.SD, mean, na.rm = TRUE), .SDcols = vars]
    sd <- train[, lapply(.SD, sd, na.rm = TRUE), .SDcols = vars]
    scale_one <- function(dt) {
      out <- copy(dt)
      for (v in vars) {
        m <- mu[[v]]
        s <- sd[[v]]
        out[, (v) := (get(v) - m) / ifelse(is.na(s) |
                                             s == 0, 1, s)]
      }
      out
    }
    list(train = scale_one(train), test = scale_one(test))
  }
  
  ## ============================================================
  ## 5) fit horseshoe (usa argumentos, no globals)
  ## ============================================================
  fit_horseshoe <- function(train_dt,
                            iter = 2000,
                            adapt_delta = 0.97,
                            seed = 1) {
    fit <- brm(
      y ~ .,
      # <-- intercepto
      data = sc$train,
      family = bernoulli("logit"),
      prior = c(
        set_prior("normal(0, 2)", class = "Intercept"),
        set_prior("horseshoe(df = 0.5, scale_global = 0.5)", class = "b")
      ),
      chains = 4,
      cores  = 4,
      iter   = 2000,
      control = list(adapt_delta = 0.99, max_treedepth = 15)
    )
    return(fit)
  }
  make_stratified_folds <- function(y, K = 5, seed = 1) {
    set.seed(seed)
    
    idx_pos <- which(y == 1L)
    idx_neg <- which(y == 0L)
    
    if (length(idx_pos) < K) {
      stop("No hay suficientes positivos para K folds")
    }
    
    folds_pos <- split(sample(idx_pos), rep(1:K, length.out = length(idx_pos)))
    folds_neg <- split(sample(idx_neg), rep(1:K, length.out = length(idx_neg)))
    
    folds <- vector("list", K)
    for (k in 1:K) {
      folds[[k]] <- c(folds_pos[[k]], folds_neg[[k]])
    }
    
    folds
  }
  
  compute_metrics <- function(p, y, threshold) {
    ok <- !is.na(p) & !is.na(y)
    p <- p[ok]
    y <- y[ok]
    
    pred <- as.integer(p >= threshold)
    
    TP <- sum(pred == 1 & y == 1)
    FP <- sum(pred == 1 & y == 0)
    FN <- sum(pred == 0 & y == 1)
    TN <- sum(pred == 0 & y == 0)
    
    sensitivity <- if ((TP + FN) == 0)
      NA_real_
    else
      TP / (TP + FN)
    specificity <- if ((TN + FP) == 0)
      NA_real_
    else
      TN / (TN + FP)
    PPV <- if ((TP + FP) == 0)
      NA_real_
    else
      TP / (TP + FP)
    NPV <- if ((TN + FN) == 0)
      NA_real_
    else
      TN / (TN + FN)
    F1  <- if (is.na(PPV) ||
               is.na(sensitivity) ||
               (PPV + sensitivity) == 0)
      NA_real_
    else
      2 * PPV * sensitivity / (PPV + sensitivity)
    
    data.table(
      threshold = threshold,
      TP = TP,
      FP = FP,
      FN = FN,
      TN = TN,
      sensitivity = sensitivity,
      specificity = specificity,
      PPV = PPV,
      NPV = NPV,
      F1 = F1,
      n_called_P = TP + FP
    )
  }
  
  safe_pr_curve <- function(p, y, label = "DATA") {
    ## limpieza NA
    ok <- !is.na(p) & !is.na(y)
    p <- p[ok]
    y <- y[ok]
    
    n_pos <- sum(y == 1L)
    n_neg <- sum(y == 0L)
    
    ## chequeos mínimos
    if (n_pos < 2 || n_neg < 2) {
      message(
        sprintf(
          "⚠️ PR no definida en %s (positivos=%d, negativos=%d)",
          label,
          n_pos,
          n_neg
        )
      )
      return(list(
        auc = NA_real_,
        n_pos = n_pos,
        n_neg = n_neg,
        curve = NULL
      ))
    }
    
    pr <- pr.curve(
      scores.class0 = p[y == 1L],
      scores.class1 = p[y == 0L],
      curve = TRUE
    )
    
    list(
      auc = pr$auc.integral,
      n_pos = n_pos,
      n_neg = n_neg,
      curve = pr$curve
    )
  }
  
  
  
  ## ============================================================
  ## LIMPIEZA INICIAL DE MEMORIA
  ## (conservar solo vcf_dt_)
  ## ============================================================
  gc(verbose = TRUE)
  
  rstan::rstan_options(auto_write = TRUE)
  options(brms.backend = "rstan")
  
  ## ============================================================
  ## 0) DATOS BASE
  ## ============================================================
  
  ## ============================================================
  ## MAIN ROBUSTO (CV + FIT FINAL + SCORE PARA VUS/CONFLICT)
  ## ============================================================
  ## ============================================================
  ## MAIN FINAL – SCORE DE EVIDENCIA (NO CLASIFICADOR)
  ## ============================================================
  
  suppressPackageStartupMessages({
    library(data.table)
    library(brms)
    library(rstan)
  })
  
  rstan::rstan_options(auto_write = TRUE)
  options(brms.backend = "rstan")
  
  ## ============================================================
  ## 0) DATOS BASE
  ## ============================================================
  ## ============================================================
  ## MAIN FINAL — SCORE DE EVIDENCIA COMPUTACIONAL (PP3)
  ## ============================================================
  
  d <- copy(vcf_dt_)
  d <- normalize_dbnsfp_predictors(d)
  d <- build_gold_binary(d)
  d_model <- make_model_dt(d)
  
  ## GOLD para entrenar
  d_gold <- d_model[!is.na(y)]
  n_pos <- sum(d_gold$y == 1L)
  n_neg <- sum(d_gold$y == 0L)
  
  message("GOLD disponibles: P=", n_pos, " | B=", n_neg)
  
  ## columnas que SIEMPRE deben existir
  vcf_dt_[, `:=`(p_hat = NA_real_, PP3_ZONE = NA_character_)]
  
  ## ============================================================
  ## Entrenar SOLO si hay base razonable
  ## ============================================================
  
  if (n_pos >= 10 && n_neg >= 50) {
    ## split simple (solo para estabilidad)
    sp <- stratified_split(d_gold, test_frac = 0.2, seed = 1)
    
    num_vars <- c(
      "REVEL_num",
      "CADD_num",
      "MPC_num",
      "ClinPred_num",
      "log10_POP_AF",
      "log10_COHORT_AF"
    )
    
    sc <- scale_numeric(sp$train, sp$test, vars = num_vars)
    assign("sc", sc, envir = .GlobalEnv)
    
    ## ---- fit congelado
    fit <- fit_horseshoe(
      train_dt = sc$train,
      iter = 3000,
      adapt_delta = 0.99,
      seed = 1
    )
    
    ## ---- score para TODO el dataset
    d_all_model <- make_model_dt(d)
    
    mu <- sc$train[, lapply(.SD, mean, na.rm = TRUE), .SDcols = num_vars]
    sdv <- sc$train[, lapply(.SD, sd, na.rm = TRUE), .SDcols = num_vars]
    
    for (v in num_vars) {
      d_all_model[, (v) :=
                    (get(v) - mu[[v]]) / ifelse(is.na(sdv[[v]]) |
                                                  sdv[[v]] == 0, 1, sdv[[v]])]
    }
    
    vcf_dt_[, p_hat :=
              fitted(fit, newdata = d_all_model, scale = "response")[, "Estimate"]]
    
    ## ============================================================
    ## Calibración empírica SOLO con benignas
    ## ============================================================
    
    p_benign <- vcf_dt_[CLNSIG_CLINVAR %like% "Benign" &
                          !CLNSIG_CLINVAR %like% "Pathogenic", p_hat]
    
    q_ref <- quantile(p_benign,
                      probs = c(0.90, 0.95, 0.99),
                      na.rm = TRUE)
    
    vcf_dt_[, PP3_ZONE := fifelse(p_hat >= q_ref["99%"],
                                  "PP3_strong",
                                  fifelse(
                                    p_hat >= q_ref["95%"],
                                    "PP3_moderate",
                                    fifelse(p_hat >= q_ref["90%"], "PP3_supporting", NA_character_)
                                  ))]
    
  } else {
    message("⚠️ GOLD insuficiente → PP3 computacional desactivado")
  }
  
  ## ============================================================
  ## INTEGRACIÓN ACMG — PP3 desde el modelo
  ## ============================================================
  
  vcf_dt_[, ACMG_PP3 :=
            IS_MISSENSE &
            PP3_ZONE %in% c("PP3_supporting", "PP3_moderate", "PP3_strong")]
  
  vcf_dt_[, ACMG_PP3_strength := PP3_ZONE]
  
  objs <- ls(envir = .GlobalEnv)
  
  is_fun <- vapply(objs, function(nm)
    is.function(get(
      nm, envir = .GlobalEnv, inherits = TRUE
    )), logical(1))
  
  keep <- c("vcf_dt_", objs[is_fun])
  # rm(list = setdiff(objs, keep), envir = .GlobalEnv)
  gc(verbose = TRUE)
  
  exome_dt <- vcf_dt_
  
  setcolorder(
    exome_dt,
    c(
      "LOCUS",
      "CHROM",
      "START",
      "END",
      "WIDTH",
      "REF",
      "ALT",
      "AB_SUSPECT",
      "FILTER",
      "N_EXOMES_LOCUS",
      "EXOMES_LOCUS",
      "VAR_IDX",
      "VARTYPE",
      "ALLELEID",
      "ZYGOSITY",
      "COMPOUND_HET",
      "GENES_AFFECTED",
      "RS_ID",
      "GENEINFO",
      "TRANSCRIPTS_AFFECTED",
      "N_TRANSCRIPTS",
      "BIOTYPES",
      "EFFECTS",
      "IMPACTS",
      "LOF",
      "NMD",
      "MC",
      "CLNHGVS",
      "CLNSIG",
      "CLNSIG_CLINVAR",
      "ACMG_CLASS",
      "ACMG_EVIDENCE",
      "ACMG_CLASS_FINAL",
      "ACMG_SUMMARY",
      "CLNDN",
      "CLNSIGCONF",
      "CLNREVSTAT",
      "CLNSIGSCV",
      "CLNSIGINCL",
      "CLNVCSO",
      "CLNVC",
      "CLNVI",
      "CLNDNINCL",
      "CLNDISDB",
      "CLNDISDBINCL",
      "ONCDN",
      "ONC",
      "ONCDISDB",
      "ONCREVSTAT",
      "ONCSCV",
      "ORIGIN",
      "dbNSFP_HGVSc_snpEff",
      "dbNSFP_HGVSp_snpEff",
      "dbNSFP_aaref",
      "dbNSFP_aaalt",
      "COHORT_AF",
      "AF_EXAC",
      "AF_ESP",
      "AF_TGP",
      "AF_gnomad",
      "POP_AF",
      "QUAL",
      "GT",
      "IS_HET_REAL",
      "N_HET",
      "N_PHASES",
      "PGT",
      "PID",
      "ACMG_PM3",
      "ACMG_PM3_compound",
      "DP_FORMAT",
      "DP_INFO",
      "GQ",
      "PL",
      "PS",
      "BaseQRankSum",
      "FS",
      "MQ",
      "MQRankSum",
      "QD",
      "ReadPosRankSum",
      "SOR",
      "AD",
      "AD_REF",
      "AD_ALT1",
      "DP_AD",
      "AB_ALT",
      "DP_FOR_AB",
      "AB_EXPECTATION",
      "FLAG_NO_AD",
      "FLAG_LOW_DP",
      "FLAG_AB_NA",
      "FLAG_DP_MISMATCH",
      "ACMG_PM2",
      "ACMG_BS1",
      "ACMG_BA1",
      "ACMG_BS2",
      "SCIREVSTAT",
      "ACMG_CLINVAR_P",
      "ACMG_CLINVAR_B",
      "dbNSFP_REVEL_score",
      "dbNSFP_CADD_phred",
      "dbNSFP_MPC_score",
      "dbNSFP_ClinPred_score",
      "dbNSFP_gnomAD4.1_joint_POPMAX_AF",
      "dbNSFP_dbNSFP_POPMAX_AF",
      "dbNSFP_MetaRNN_pred",
      "dbNSFP_AlphaMissense_pred",
      "dbNSFP_SIFT_pred",
      "dbNSFP_Polyphen2_HDIV_pred",
      "dbNSFP_Polyphen2_HVAR_pred",
      "dbNSFP_SIFT_score",
      "dbNSFP_Polyphen2_HDIV_score",
      "dbNSFP_Polyphen2_HVAR_score",
      "dbNSFP_APPRIS",
      "dbNSFP_MANE",
      "IS_MISSENSE",
      "ACMG_PVS1",
      "ACMG_PP3",
      "ACMG_BP4",
      "p_hat",
      "PP3_ZONE",
      "ACMG_PP3_strength",
      "ACMG_PATH_SCORE",
      "ACMG_BENIGN_SCORE"
    )
  )
  
  
  ## ============================================================
  ## NORMALIZACIÓN DEFINITIVA DE p_hat Y PP3 (PRE-HERENCIA)
  ## ============================================================
  
  setDT(exome_dt)
  
  ## 1) Normalizar p_hat inválidos
  exome_dt[is.nan(p_hat) |
             is.infinite(p_hat), p_hat := NA_real_]
  
  ## 2) Forzar PP3_ZONE a NA si no hay p_hat
  exome_dt[is.na(p_hat), PP3_ZONE := NA_character_]
  
  ## 3) Recalcular ACMG_PP3 SOLO desde el modelo bayesiano
  ##    (no desde reglas antiguas)
  exome_dt[, ACMG_PP3 :=
             IS_MISSENSE &
             !is.na(p_hat) &
             PP3_ZONE %in% c("PP3_supporting", "PP3_moderate", "PP3_strong")]
  
  ## 4) Guardar fuerza explícita (auditable)
  exome_dt[, ACMG_PP3_strength :=
             fifelse(ACMG_PP3, PP3_ZONE, NA_character_)]
  
  message("✔ p_hat normalizado y PP3 recalculado (modelo bayesiano activo)")
  ## ============================================================
  ## FIX DEFINITIVO: USO CORRECTO DE p_hat EN ACMG
  ## ============================================================
  
  setDT(exome_dt)
  
  ## 1) Normalizar p_hat inválido
  exome_dt[is.nan(p_hat) |
             is.infinite(p_hat), p_hat := NA_real_]
  
  ## 2) Definir PP3 SOLO si:
  ##    - missense
  ##    - p_hat definido
  ##    - zona válida
  exome_dt[, ACMG_PP3 :=
             IS_MISSENSE &
             !is.na(p_hat) &
             PP3_ZONE %in% c("PP3_supporting", "PP3_moderate", "PP3_strong")]
  
  ## 3) Registrar fuerza de evidencia computacional
  exome_dt[, ACMG_PP3_strength :=
             fifelse(ACMG_PP3, PP3_ZONE, NA_character_)]
  
  ## 4) Recalcular PATH_SCORE (único cambio: PP3 ahora es real)
  exome_dt[, ACMG_PATH_SCORE :=
             2 * ACMG_PVS1 +
             1 * ACMG_PM2 +
             1 * ACMG_PM3 +
             1 * ACMG_PM3_compound +
             1 * ACMG_PP3 +
             2 * ACMG_CLINVAR_P]
  
  
  ## 5) Clasificación final (idéntica lógica)
  exome_dt[, ACMG_CLASS_FINAL := fifelse(
    ACMG_BENIGN_SCORE >= 3 & ACMG_PATH_SCORE == 0,
    "Benign",
    fifelse(
      ACMG_BENIGN_SCORE >= 2 & ACMG_PATH_SCORE <= 1,
      "Likely_benign",
      fifelse(
        ACMG_PATH_SCORE >= 4 & ACMG_BENIGN_SCORE == 0,
        "Pathogenic",
        fifelse(
          ACMG_PATH_SCORE >= 2 & ACMG_BENIGN_SCORE <= 1,
          "Likely_pathogenic",
          "VUS"
        )
      )
    )
  )]
  ## ===============================
  ## 1. HPO: pasar a data.table y colapsar por gen
  ## ===============================
  hpo_df <- read.delim(hpo_file)
  hpo_dt <- as.data.table(hpo_df)
  
  hpo_collapsed <- hpo_dt[
    ,
    .(
      HPO_IDS     = paste(unique(hpo_id), collapse = "|"),
      HPO_NAMES   = paste(unique(hpo_name), collapse = "|"),
      DISEASE_IDS = paste(unique(disease_id), collapse = "|")
    ),
    by = gene_symbol
  ]
  
  ## ===============================
  ## 2. Preparar exome_dt
  ## ===============================
  setDT(exome_dt)
  exome_dt[, row_id := .I]
  
  ## Expandir genes (coma → filas)
  exome_long <- exome_dt[
    ,
    .(gene_symbol = trimws(unlist(strsplit(GENES_AFFECTED, ",")))),
    by = row_id
  ]
  
  ## ===============================
  ## 3. Merge gen–HPO
  ## ===============================
  exome_long <- merge(
    exome_long,
    hpo_collapsed,
    by = "gene_symbol",
    all.x = TRUE
  )
  
  ## ===============================
  ## 4. Recolapsar por variante original
  ## ===============================
  hpo_per_variant <- exome_long[
    ,
    .(
      HPO_IDS     = paste(HPO_IDS, collapse = ","),
      HPO_NAMES   = paste(HPO_NAMES, collapse = ","),
      DISEASE_IDS = paste(DISEASE_IDS, collapse = ",")
    ),
    by = row_id
  ]
  
  ## ===============================
  ## 5. Añadir columnas finales a exome_dt
  ## ===============================
  exome_dt <- merge(
    exome_dt,
    hpo_per_variant,
    by = "row_id",
    all.x = TRUE
  )
  
  ## Opcional: limpiar
  exome_dt[, row_id := NULL]
  ## exoma a escribir
  
  library(data.table)
  
  inferir_herencia_hpo <- function(exome_dt,
                                   hpo_col = "HPO_NAMES",
                                   out_col = "INHERITANCE_HPO") {
    
    stopifnot(hpo_col %in% names(exome_dt))
    
    library(data.table)
    setDT(exome_dt)
    
    extract_inheritance <- function(x) {
      if (is.na(x) || x == "")
        return(NA_character_)
      
      # separar por genes (coma = genes distintos)
      genes <- unlist(strsplit(x, ",", fixed = TRUE))
      
      inh_all <- unlist(lapply(genes, function(g) {
        # separar términos HPO del mismo gen
        terms <- unlist(strsplit(g, "\\|"))
        
        # capturar CUALQUIER término que contenga "inheritance"
        hits <- terms[grepl("inheritance", terms, ignore.case = TRUE)]
        
        hits
      }))
      
      inh_all <- unique(trimws(inh_all))
      
      if (length(inh_all) == 0)
        return(NA_character_)
      
      paste(inh_all, collapse = "|")
    }
    
    exome_dt[, (out_col) := vapply(
      get(hpo_col),
      extract_inheritance,
      character(1)
    )]
    
    exome_dt
  }
  
  # exome_dt <- fread("./file_ready_analysis_optimized.csv")
  
  exome_dt <- inferir_herencia_hpo(exome_dt)
  
  cols <- names(exome_dt)
  i <- match("ALLELEID", cols)
  
  if (is.na(i)) {
    stop("ALLELEID no existe en el data.table")
  }
  
  # quitar INHERITANCE_HPO del vector original
  cols_wo_inh <- cols[cols != "INHERITANCE_HPO"]
  
  setcolorder(
    exome_dt,
    c(
      cols_wo_inh[1:i],                # hasta ALLELEID
      "INHERITANCE_HPO",               # insertar aquí
      cols_wo_inh[(i + 1):length(cols_wo_inh)]
    )
  )
  
  sanitize_inheritance <- function(exome_dt,
                                   in_col  = "INHERITANCE_HPO",
                                   out_col = "INHERITANCE_CANONICAL") {
    
    library(data.table)
    setDT(exome_dt)
    
    map_inheritance <- function(x) {
      if (is.na(x) || x == "")
        return("UNK")
      
      x <- tolower(x)
      
      tags <- character(0)
      
      if (grepl("autosomal recessive", x)) tags <- c(tags, "AR")
      if (grepl("autosomal dominant",  x)) tags <- c(tags, "AD")
      
      if (grepl("x-linked", x)) tags <- c(tags, "XL")
      if (grepl("y-linked", x)) tags <- c(tags, "YL")
      
      if (grepl("mitochondrial", x)) tags <- c(tags, "MT")
      
      if (grepl("digenic|oligogenic", x)) tags <- c(tags, "DG")
      if (grepl("polygenic", x)) tags <- c(tags, "POLY")
      
      if (grepl("non-mendelian", x)) tags <- c(tags, "NM")
      
      tags <- unique(tags)
      
      if (length(tags) == 0)
        return("UNK")
      
      paste(sort(tags), collapse = "|")
    }
    
    exome_dt[, (out_col) := vapply(
      get(in_col),
      map_inheritance,
      character(1)
    )]
    
    exome_dt
  }
  
  exome_dt <- sanitize_inheritance(exome_dt)
  
  integrate_inheritance_into_acmg <- function(
    exome_dt,
    inherit_col = "INHERITANCE_CANONICAL",
    zyg_col = "ZYGOSITY",
    comphet_col = "COMPOUND_HET",
    af_col = "COHORT_AF",
    pm3_af_max = 0.01,
    allow_pm3_if_unk = FALSE,   # clínico: por defecto NO
    update_evidence_string = TRUE,
    evidence_col = "ACMG_EVIDENCE",
    update_scores = TRUE,
    path_score_col = "ACMG_PATH_SCORE",
    benign_score_col = "ACMG_BENIGN_SCORE",
    update_class = TRUE,
    class_final_col = "ACMG_CLASS_FINAL",
    ## --- recomendado: mantener pesos originales
    keep_original_weights = TRUE,
    ## --- opcional: PM3 strength (supporting/moderate/strong)
    use_pm3_strength = FALSE,
    pm3_supporting_weight = 1,   # si activas supporting
    pm3_moderate_weight   = 1,
    pm3_strong_weight     = 2
  ) {
    
    suppressPackageStartupMessages(library(data.table))
    stopifnot(is.data.table(exome_dt) || is.data.frame(exome_dt))
    setDT(exome_dt)
    
    ## -------------------------
    ## Checks mínimos
    ## -------------------------
    need_cols <- c(inherit_col, zyg_col, comphet_col, "ACMG_CLINVAR_B", "ACMG_CLINVAR_P")
    miss <- setdiff(need_cols, names(exome_dt))
    if (length(miss) > 0) stop("Faltan columnas: ", paste(miss, collapse = ", "))
    
    if (!af_col %in% names(exome_dt)) exome_dt[, (af_col) := NA_real_]
    
    ## crear columnas ACMG si no existen (robustez)
    for (cc in c("ACMG_PVS1","ACMG_PM2","ACMG_PP3","ACMG_BA1","ACMG_BS1","ACMG_BS2","ACMG_BP4")) {
      if (!cc %in% names(exome_dt)) exome_dt[, (cc) := FALSE]
    }
    if (!evidence_col %in% names(exome_dt)) exome_dt[, (evidence_col) := NA_character_]
    if (!path_score_col %in% names(exome_dt)) exome_dt[, (path_score_col) := NA_real_]
    if (!benign_score_col %in% names(exome_dt)) exome_dt[, (benign_score_col) := NA_real_]
    if (!class_final_col %in% names(exome_dt)) exome_dt[, (class_final_col) := NA_character_]
    
    ## -------------------------
    ## Helpers
    ## -------------------------
    trim_upper <- function(x) toupper(trimws(as.character(x)))
    has_token  <- function(x, tok) grepl(paste0("(^|\\|)", tok, "(\\||$)"), x)
    as_int_T   <- function(x) as.integer(x %in% TRUE)
    
    ## -------------------------
    ## Normalización base
    ## -------------------------
    exome_dt[, (inherit_col) := fifelse(is.na(get(inherit_col)) | get(inherit_col) == "", "UNK", get(inherit_col))]
    
    exome_dt[, `:=`(
      .inh = as.character(get(inherit_col)),
      .zyg = trim_upper(get(zyg_col)),
      .ch  = tolower(trimws(as.character(get(comphet_col)))),
      .af  = suppressWarnings(as.numeric(get(af_col)))
    )]
    
    exome_dt[, .is_AR  := has_token(.inh, "AR")]
    exome_dt[, .is_UNK := has_token(.inh, "UNK")]
    
    ## -------------------------
    ## AR completo
    ## -------------------------
    exome_dt[, .is_comphet_confirmed := .ch %chin% c("confirmed", "yes", "si", "sí", "true")]
    exome_dt[, .is_comphet_possible  := .ch %chin% c("possible", "posible")]
    exome_dt[, .ar_complete := (.zyg == "HOM") | .is_comphet_confirmed]
    
    ## -------------------------
    ## PM3 / PM3_compound (bloqueado si ClinVar_B)
    ## -------------------------
    exome_dt[, ACMG_PM3 := {
      gate  <- (.is_AR | (allow_pm3_if_unk & .is_UNK))
      af_ok <- is.na(.af) | (.af <= pm3_af_max)
      gate & .ar_complete & af_ok & !ACMG_CLINVAR_B
    }]
    
    exome_dt[, ACMG_PM3_compound := {
      gate <- (.is_AR | (allow_pm3_if_unk & .is_UNK))
      gate & .is_comphet_confirmed & !ACMG_CLINVAR_B
    }]
    
    ## -------------------------
    ## PM3 efectivo (compound > simple)
    ## -------------------------
    exome_dt[, ACMG_PM3_EFFECTIVE := fifelse(ACMG_PM3_compound, TRUE, ACMG_PM3)]
    exome_dt[ACMG_PM3_compound == TRUE, ACMG_PM3 := FALSE]
    
    ## -------------------------
    ## (Opcional) PM3 strength real
    ##  - strong: AR & (HOM o comphet confirmed)
    ##  - moderate: AR & comphet possible
    ##  - supporting: (si decides) AR & PM3_EFFECTIVE pero sin strong/moderate (se deja NA por defecto)
    ## -------------------------
    exome_dt[, `:=`(ACMG_PM3_strength = NA_character_, ACMG_PM3_WEIGHT = 0)]
    if (isTRUE(use_pm3_strength)) {
      
      exome_dt[, ACMG_PM3_strength := {
        gate  <- (.is_AR | (allow_pm3_if_unk & .is_UNK))
        af_ok <- is.na(.af) | (.af <= pm3_af_max)
        
        ok <- gate & af_ok & !ACMG_CLINVAR_B
        
        fifelse(!ok, NA_character_,
                fifelse(.zyg == "HOM" | .is_comphet_confirmed, "PM3_strong",
                        fifelse(.is_comphet_possible, "PM3_moderate",
                                fifelse(ACMG_PM3_EFFECTIVE, "PM3_supporting", NA_character_)
                        )
                )
        )
      }]
      
      exome_dt[, ACMG_PM3_WEIGHT := fifelse(
        ACMG_PM3_strength == "PM3_strong", pm3_strong_weight,
        fifelse(ACMG_PM3_strength == "PM3_moderate", pm3_moderate_weight,
                fifelse(ACMG_PM3_strength == "PM3_supporting", pm3_supporting_weight, 0))
      )]
    }
    
    ## -------------------------
    ## Evidence string (limpieza robusta + reinsert)
    ## -------------------------
    if (isTRUE(update_evidence_string) && evidence_col %in% names(exome_dt)) {
      
      exome_dt[, (evidence_col) := {
        x <- as.character(get(evidence_col))
        x[is.na(x)] <- ""
        
        ## limpiar tokens PM3 previos (incluye variantes de naming)
        x <- gsub("(^|;)(PM3_strong|PM3_moderate|PM3_supporting|PM3_compound|PM3_EFFECTIVE|PM3)(;|$)", ";", x, perl = TRUE)
        x <- gsub(";;+", ";", x)
        x <- gsub("^;|;$", "", x)
        
        add_tok <- function(s, tok) ifelse(is.na(s) | s == "", tok, paste(s, tok, sep = ";"))
        out <- x
        
        if (isTRUE(use_pm3_strength)) {
          idx <- !is.na(ACMG_PM3_strength) & nzchar(ACMG_PM3_strength)
          if (any(idx)) out[idx] <- add_tok(out[idx], ACMG_PM3_strength[idx])
        } else {
          idx1 <- ACMG_PM3_EFFECTIVE & ACMG_PM3_compound
          idx2 <- ACMG_PM3_EFFECTIVE & !ACMG_PM3_compound
          
          if (any(idx1)) out[idx1] <- add_tok(out[idx1], "PM3_compound")
          if (any(idx2)) out[idx2] <- add_tok(out[idx2], "PM3")
        }
        
        gsub("^;|;$", "", gsub(";;+", ";", out))
      }]
    }
    
    ## -------------------------
    ## Recalcular SCORES
    ##  - FIX crítico: mantener pesos 2/1 y ClinVar_* 2
    ## -------------------------
    if (isTRUE(update_scores)) {
      
      if (isTRUE(keep_original_weights)) {
        
        if (isTRUE(use_pm3_strength)) {
          exome_dt[, (path_score_col) :=
                     2L*as_int_T(ACMG_PVS1) +
                     1L*as_int_T(ACMG_PM2) +
                     as.numeric(ACMG_PM3_WEIGHT) +
                     1L*as_int_T(ACMG_PP3) +
                     2L*as_int_T(ACMG_CLINVAR_P)]
        } else {
          exome_dt[, (path_score_col) :=
                     2L*as_int_T(ACMG_PVS1) +
                     1L*as_int_T(ACMG_PM2) +
                     1L*as_int_T(ACMG_PM3_EFFECTIVE) +
                     1L*as_int_T(ACMG_PP3) +
                     2L*as_int_T(ACMG_CLINVAR_P)]
        }
        
        exome_dt[, (benign_score_col) :=
                   2L*as_int_T(ACMG_BA1) +
                   1L*as_int_T(ACMG_BS1) +
                   1L*as_int_T(ACMG_BS2) +
                   1L*as_int_T(ACMG_BP4) +
                   2L*as_int_T(ACMG_CLINVAR_B)]
        
      } else {
        ## fallback: suma binaria (NO recomendado clínicamente; solo por compatibilidad)
        path_cols <- intersect(c("ACMG_PVS1","ACMG_PM2","ACMG_PM3_EFFECTIVE","ACMG_PP3","ACMG_CLINVAR_P"), names(exome_dt))
        benign_cols <- intersect(c("ACMG_BA1","ACMG_BS1","ACMG_BS2","ACMG_BP4","ACMG_CLINVAR_B"), names(exome_dt))
        
        exome_dt[, (path_score_col) :=
                   if (length(path_cols) == 0) NA_real_
                 else Reduce(`+`, lapply(.SD, function(x) as_int_T(x))), .SDcols = path_cols]
        
        exome_dt[, (benign_score_col) :=
                   if (length(benign_cols) == 0) NA_real_
                 else Reduce(`+`, lapply(.SD, function(x) as_int_T(x))), .SDcols = benign_cols]
      }
    }
    
    ## -------------------------
    ## Clasificación FINAL (misma lógica que tienes, pero usando scores ya corregidos)
    ## -------------------------
    if (isTRUE(update_class) &&
        class_final_col %in% names(exome_dt) &&
        path_score_col %in% names(exome_dt) &&
        benign_score_col %in% names(exome_dt)) {
      
      exome_dt[, (class_final_col) := {
        P <- suppressWarnings(as.numeric(get(path_score_col)))
        B <- suppressWarnings(as.numeric(get(benign_score_col)))
        
        out <- fifelse(
          !is.na(B) & B >= 5, "Benign",
          fifelse(
            !is.na(B) & B >= 2 & (is.na(P) | P <= 1), "Likely_benign",
            fifelse(
              !is.na(P) & P >= 5, "Pathogenic",
              fifelse(!is.na(P) & P >= 3, "Likely_pathogenic", "VUS")
            )
          )
        )
        
        ## techo ClinVar Benign
        out[ACMG_CLINVAR_B & out %chin% c("VUS","Likely_pathogenic","Pathogenic")] <- "Likely_benign"
        ## suelo ClinVar Pathogenic
        out[ACMG_CLINVAR_P & out %chin% c("Benign","Likely_benign","VUS")] <- "Likely_pathogenic"
        
        out
      }]
    }
    
    ## -------------------------
    ## Limpieza
    ## -------------------------
    exome_dt[, c(".inh",".zyg",".ch",".af",".is_AR",".is_UNK",
                 ".is_comphet_confirmed",".is_comphet_possible",".ar_complete") := NULL]
    
    exome_dt
  }
  
  
  # ============================
  
  exome_dt[, `:=`(
    ACMG_CLASS_BEFORE_INH   = ACMG_CLASS_FINAL,
    ACMG_PATH_SCORE_BEFORE  = ACMG_PATH_SCORE,
    ACMG_BENIGN_SCORE_BEFORE = ACMG_BENIGN_SCORE,
    ACMG_EVIDENCE_BEFORE_INH = ACMG_EVIDENCE
  )]
  
  
  exome_dt <- integrate_inheritance_into_acmg(exome_dt)
  order_evidence <- function(x) {
    x <- as.character(x)
    
    out <- vapply(x, function(xx) {
      if (is.na(xx) || trimws(xx) == "") return("")
      
      toks <- unlist(strsplit(xx, ";", fixed = TRUE))
      toks <- unique(trimws(toks))
      toks <- toks[toks != ""]
      
      priority <- c(
        "PVS1",
        "PM3_strong","PM3_compound","PM3","PM2",
        "PP3",
        "BA1","BS1","BS2","BP4",
        "ClinVar_P","ClinVar_B"
      )
      
      toks <- c(intersect(priority, toks),
                setdiff(toks, priority))
      
      paste(toks, collapse = ";")
    }, character(1))
    
    out
  }
  exome_dt[, ACMG_SUMMARY := {
    
    cls <- fifelse(
      is.na(ACMG_CLASS_FINAL) | trimws(ACMG_CLASS_FINAL) == "",
      "NA",
      as.character(ACMG_CLASS_FINAL)
    )
    
    P <- fifelse(
      is.na(ACMG_PATH_SCORE),
      "NA",
      as.character(ACMG_PATH_SCORE)
    )
    
    B <- fifelse(
      is.na(ACMG_BENIGN_SCORE),
      "NA",
      as.character(ACMG_BENIGN_SCORE)
    )
    
    ev <- order_evidence(ACMG_EVIDENCE)
    
    paste0(
      cls,
      " | P=", P,
      " B=", B,
      " | ",
      ev
    )
  }]
  
  library(data.table)
  setDT(exome_dt)
  # efectos que, si están solos, suelen ser "baja fuerza" (sin tocar ACMG)
  low_only_terms <- c(
    "synonymous_variant",
    "intron_variant",
    "intergenic_region",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "3_prime_UTR_variant",
    "5_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "splice_region_variant"
  )
  
  # efectos "fuertes" (si aparece cualquiera, sube fuerza)
  strong_terms <- c(
    "missense_variant",
    "frameshift_variant",
    "stop_gained",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "start_lost",
    "stop_lost",
    "inframe_insertion",
    "inframe_deletion",
    "disruptive_inframe_insertion",
    "disruptive_inframe_deletion",
    "conservative_inframe_insertion",
    "conservative_inframe_deletion"
  )
  
  exome_dt[, LAB_FUNCTIONAL_CONFIDENCE := {
    eff <- tolower(fifelse(is.na(EFFECTS), "", EFFECTS))
    imp <- toupper(fifelse(is.na(IMPACTS), "", IMPACTS))
    
    # vector de tokens de efectos por fila (split por coma, & y espacios)
    eff_tokens <- strsplit(gsub("[[:space:]]+", "", eff), "[,&]")
    # para cada fila, evaluar condiciones
    vapply(seq_along(eff_tokens), function(i){
      toks <- eff_tokens[[i]]
      toks <- toks[toks != ""]
      
      has_strong <- any(toks %chin% tolower(strong_terms))
      has_any_effect <- length(toks) > 0
      
      # "solo baja" si TODOS los tokens están en low_only_terms
      low_only <- has_any_effect && all(toks %chin% tolower(low_only_terms))
      
      # IMPACTS: HIGH/MODERATE/LOW/MODIFIER (puede venir como "MODERATE,LOW")
      has_high <- grepl("HIGH", imp[i], fixed = TRUE)
      has_mod  <- grepl("MODERATE", imp[i], fixed = TRUE)
      has_low  <- grepl("LOW", imp[i], fixed = TRUE)
      has_only_low_or_modifier <- has_any_effect && !has_high && !has_mod && (has_low || grepl("MODIFIER", imp[i], fixed=TRUE))
      
      if (has_high) return("HIGH_FUNCTIONAL")
      if (has_strong || has_mod) return("MODERATE_FUNCTIONAL")
      if (low_only || has_only_low_or_modifier) return("LOW_FUNCTIONAL")
      if (!has_any_effect) return("UNCERTAIN_FUNCTIONAL")
      return("LOW_FUNCTIONAL")
    }, character(1))
  }]
  exome_dt[, LAB_PATHOGENIC_CONFIDENCE := fifelse(
    
    # 1️⃣ Solo Pathogenic / Likely_pathogenic
    !ACMG_CLASS_FINAL %chin% c("Pathogenic", "Likely_pathogenic"),
    
    "NOT_APPLICABLE",
    
    # 2️⃣ Penalización dura por exoma de referencia
    fifelse(
      !is.na(EXOMES_LOCUS) &
        grepl("exoma_referencia", tolower(EXOMES_LOCUS)),
      
      "LOW_CONFIDENCE_INTERNAL",
      
      # 3️⃣ Frecuencia en cohorte clínica real
      fifelse(
        is.na(N_EXOMES_LOCUS),
        "UNCERTAIN_CONFIDENCE",
        
        fifelse(
          N_EXOMES_LOCUS <= 2,
          "HIGH_CONFIDENCE",
          
          fifelse(
            N_EXOMES_LOCUS <= 5,
            "MODERATE_CONFIDENCE",
            
            "LOW_CONFIDENCE_INTERNAL"
          )
        )
      )
    )
  )]
  exome_dt[, LAB_FINAL_PRIORITY := fifelse(
    LAB_PATHOGENIC_CONFIDENCE %chin% c("HIGH_CONFIDENCE","MODERATE_CONFIDENCE") &
      LAB_FUNCTIONAL_CONFIDENCE %chin% c("HIGH_FUNCTIONAL","MODERATE_FUNCTIONAL"),
    "PRIORITY",
    "DEPRIORITIZE"
  )]
  
  
  ## ============================================================
  ## EXPORTES FINALES
  ## ============================================================
  
  unicas <- exome_dt[N_EXOMES_LOCUS == 1]
  unicas_second <- unicas[GENES_AFFECTED %in% genes_amcg_horizon]
  
  f_all  <- file.path(out_dir, "file_ready_analysis_optimized.csv")
  f_uni  <- file.path(out_dir, "file_ready_analysis_optimized_UNICAS.csv")
  f_hor  <- file.path(out_dir,
                      "file_ready_analysis_optimized_UNICAS_ACMG_HORIZON.csv")
  
  fwrite(exome_dt, f_all)
  fwrite(unicas, f_uni)
  fwrite(unicas_second, f_hor)
  
  message("Archivos generados:")
  message(" - ", f_all)
  message(" - ", f_uni)
  message(" - ", f_hor)
  
  
  invisible(
    list(
      exome_dt = exome_dt,
      unicas = unicas,
      unicas_acmg_horizon = unicas_second,
      output_dir = out_dir
    )
  )
}
if (sys.nframe() == 0L) {
  argv <- commandArgs(trailingOnly = TRUE)
  if (length(argv) < 4) {
    stop("Uso: vcf_process.R <vcf_file> <genes_horizon.csv> <base_lab.rds> <hpo_file.txt>")
  }
  
  vcf_process(
    vcf_file      = argv[1],
    genes_horizon = argv[2],
    base_lab      = argv[3],
    hpo_file      = argv[4]
  )
}