load.libs <- c(
  "data.table",
  "limma",
  "grid",
  "egg",
  "Rbwa",
  "ggpubr",
  "tools",
  "dplyr",
  "tidyr",
  "ggplot2",
  "stringr",
  "ggplot2",
  "xtable",
  "VariantAnnotation",
  "tibble",
  "stringr",
  "GenomicRanges",
  "S4Vectors"
)

lapply(load.libs, require, character.only = TRUE)


control_calidad <- function(fastq_dir, output_dir) {
  ### Creamos directorio de QC si no existe
  if (!dir.exists(file.path(output_dir, "QC"))) {
    dir.create(file.path(output_dir, "QC"))
  }
  ### Si no existen archivos de salida
  if (length(list.files(file.path(output_dir, "QC"))) == 0) {
    command <-
      paste(
        "fastqc -t 4 ",
        paste0(fastq_dir, "/*.", unique(file_ext(
          list.files(fastq_dir)
        ))),
        "-o",
        file.path(output_dir, "QC")
      )
    system(command, intern = T)
  } else{
    message("Ya se ha hecho el control de Calidad")
  }
  
  
  
}

fn_exists_fasta <- function(folder_fasta) {
  extension = unlist(lapply(list.files(folder_fasta, pattern = "fa"), function(x)
    file_ext(x)))
  extension_fa <- extension[grep("^fa$", extension)]
  extension_fasta <- extension[grep("^fasta$", extension)]
  ## Es fa o fasta, y existe ?
  if (length(extension_fa) == 0 && length(extension_fasta) == 0) {
    stop("No existe archivo de referencia")
    
  } else if (length(extension_fa) != 0 ||
             length(extension_fasta) != 0) {
    if (length(extension_fa) != 0) {
      extension <- extension_fa
    } else if (length(extension_fasta) != 0) {
      extension <- extension_fasta
    }
    
    
  }
  ## el archivo fasta ?
  fasta_file <-
    list.files(folder_fasta,
               pattern = paste0(".", extension, "$"),
               full.names = T)
  return(fasta_file)
  
  
}



index_fasta_samtools <- function(folder_fasta = folder_fasta) {
  ## Vemos si existe el archivo con patron fa
  
  fasta_file <- fn_exists_fasta(folder_fasta)
  if (length(list.files(folder_fasta, pattern = "fai$")) == 0) {
    command <- paste("samtools faidx", fasta_file)
    print(command)
    system(command, intern = T)
  } else{
    message("Ya esta el index fai")
  }
  
  
  
  
}

index_bwa <- function(folder_fasta) {
  fasta_file <- fn_exists_fasta(folder_fasta)
  
  extension <- c("amb", "ann", "bwt", "pac", "sa")
  
  if (length(file.exists(list.files(folder_fasta, extension[1], full.names = T))) ==
      0 |
      length(file.exists(list.files(folder_fasta, extension[2], full.names = T))) ==
      0 |
      length(file.exists(list.files(folder_fasta, extension[3], full.names = T))) ==
      0 |
      length(file.exists(list.files(folder_fasta, extension[4], full.names = T))) ==
      0 |
      length(file.exists(list.files(folder_fasta, extension[5], full.names = T))) ==
      0) {
    ## no existen bwa index
    print("Creando ficheros índices para bwa mem...")
    
    comando <- paste("bwa index", fasta_file)
    system(command = comando, intern = T)
    
  } else if (length(file.exists(list.files(folder_fasta, extension[1], full.names = T))) !=
             0 &
             length(file.exists(list.files(folder_fasta, extension[2], full.names = T))) !=
             0 &
             length(file.exists(list.files(folder_fasta, extension[3], full.names = T))) !=
             0 &
             length(file.exists(list.files(folder_fasta, extension[4], full.names = T))) !=
             0 &
             length(file.exists(list.files(folder_fasta, extension[5], full.names = T))) !=
             0) {
    message("Ya se han creado los ficheros para el alineamiento")
    
  }
  
  
}



bwamem <- function(fastq_dir = fastq_dir ,
                   folder_fasta = folder_fasta) {
  ### conseguimos el archivo fasta
  
  fasta_file <- fn_exists_fasta(folder_fasta)
  
  ####fai index exist ?######
  
  if (!length(file.exists(list.files(dirname(fasta_file), "fai")))) {
    print("### generating fai index...")
    
    index_fasta_samtools(folder_fasta)
    
  }
  
  ## buscamos los archivos fastq
  fastq_files <- list.files(fastq_dir, full.names = T)
  ## Archivos fastq concatenados
  fastq_full_path_files <-
    c(fastq_1 = fastq_files[1], fastq_2 = fastq_files[2])
  ## Creamos directorio de mapeo
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  
  output_file_name <-
    file_path_sans_ext(output_file_name[length(output_file_name)])
  ### Creando el directorio de mapeo
  
  mapping_output_dir <- file.path(output_dir, "mapping_output")
  if (!dir.exists(mapping_output_dir)) {
    dir.create(mapping_output_dir)
  }
  ### out file name sam file
  output_file_sam <-
    file.path(mapping_output_dir, paste0(output_file_name, ".sam"))
  ### out file name bam file
  
  
  output_file_bam <-
    file.path(mapping_output_dir, paste0(output_file_name, ".bam"))
  ### out file name bam sort file
  
  output_file_sorted_bam <-
    file.path(mapping_output_dir,
              paste0(output_file_name, ".sorted.bam"))
  ### Si no esta mapepeado, mapear
  
  if (length(list.files(mapping_output_dir)) == 0) {
    print("#### MAPPING...#####")
    
    comando <-
      paste(
        "bwa mem -HMpP -v 3 -t 8",
        fasta_file,
        fastq_full_path_files[1],
        fastq_full_path_files[2],
        ">",
        output_file_sam
      )
    
    
    system(command = comando, intern = T)
    ### Si ya se ha mapeado pero no esta el bam, crearlo
    
    print("#### SAM TO BAM")
    
    command_sam_to_bam <-
      paste("samtools view -S -b -h -@ 8",
            output_file_sam,
            "-o",
            output_file_bam)
    
    print(command_sam_to_bam)
    system(command_sam_to_bam, intern = T)
    ### Si ya esta el bam, sortearlo
    
    print("### BAM to sorted BAM")
    command_out_bam_sorted <-
      paste(
        "~/tools/gatk-4.3.0.0/gatk SortSam -CREATE_INDEX true -INPUT",
        output_file_bam,
        "-OUTPUT",
        output_file_sorted_bam,
        "-SORT_ORDER coordinate -VALIDATION_STRINGENCY STRICT"
      )
    
    print(command_out_bam_sorted)
    system(command_out_bam_sorted, intern = T)
    ### AHORA SE PROCEDERIA A MERGE LOS BAMS EN CASO DE TRIO
    
  } else if (file.exists(output_file_bam) &&
             !file.exists(output_file_sam) &&
             !file.exists(output_file_sorted_bam)) {
    print("#### SAM TO BAM")
    
    command_sam_to_bam <-
      paste("samtools view -S -b -@ 8 ",
            output_file_sam,
            "-o",
            output_file_bam)
    
    print(command_sam_to_bam)
    system(command_sam_to_bam, intern = T)
    ### Si ya esta el bam, sortearlo
    
    print("### BAM to sorted BAM")
    command_out_bam_sorted <-
      paste(
        "~/tools/gatk-4.3.0.0/gatk SortSam -CREATE_INDEX true -R",
        fasta_file,
        "-INPUT",
        output_file_bam,
        "-OUTPUT",
        output_file_sorted_bam,
        "-SORT_ORDER coordinate -VALIDATION_STRINGENCY STRICT"
      )
    
    print(command_out_bam_sorted)
    system(command_out_bam_sorted, intern = T)
  }
  else if (file.exists(output_file_sam) &&
           file.exists(output_file_bam) &&
           !file.exists(output_file_sorted_bam)) {
    print("### BAM to sorted BAM")
    command_out_bam_sorted <-
      paste(
        "~/tools/gatk-4.3.0.0/gatk SortSam -CREATE_INDEX true -R",
        fasta_file,
        "-INPUT",
        output_file_bam,
        "-OUTPUT",
        output_file_sorted_bam,
        "-SORT_ORDER coordinate -VALIDATION_STRINGENCY STRICT"
      )
    
  } else if (file.exists(output_file_sam) &&
             file.exists(output_file_bam) &&
             file.exists(output_file_sorted_bam)) {
    print("YA SE HA MAPEADO")
    
  }  else{
    stop(message("No se ha mapeado bien", call = T))
  }
  
  
}





markdups <- function(output_dir = output_dir,
                     fastq_dir = fastq_dir) {
  fastq_files <- list.files(fastq_dir, full.names = F)
  
  ## Creamos directorio de mapeo, renombramos al archivo
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  ## quitamos la extension del archivo
  output_file_name <- file_path_sans_ext(output_file_name)
  ## el directorio del mapeo esta creado pero lo tenemos que guardar de nuevo
  mapping_output_dir <- file.path(output_dir, "mapping_output")
  ## nombramos al arhivo bam sorteado
  bam_file <-
    file.path(mapping_output_dir,
              paste0(output_file_name, ".sorted.bam"))
  ## nombramos al archivo bam marcado con duplciados
  mark_file <-
    file.path(mapping_output_dir,
              paste0(output_file_name, ".sorted.mark_dup.bam"))
  ## nombramos al archivo con las meteicas
  metrics_file <-
    file.path(mapping_output_dir,
              paste0(output_file_name, ".sorted.mark_dup.txt"))
  ## Si no existe el archivo crearlo
  if (!file.exists(mark_file)) {
    command <- paste(
      "~/tools/gatk-4.3.0.0/gatk MarkDuplicates -CREATE_INDEX true -INPUT",
      bam_file,
      "-VALIDATION_STRINGENCY STRICT --REMOVE_DUPLICATES true -OUTPUT",
      mark_file,
      "-M",
      metrics_file
    )
    
    print(command)
    system(command = command, intern = T)
    
  } else{
    message("Ya se han marcado los duplicados")
  }
  ## llamamos al comando
  
  
}


create_dict <- function(folder_fasta) {
  ## llamamos al archivo fasta
  fasta_file <- fn_exists_fasta(folder_fasta)
  fasta_dict <- paste0(file_path_sans_ext(fasta_file), ".dict")
  if (!file.exists(fasta_dict)) {
    ## creamos el diccionario del fasta
    command <-
      paste("~/tools/gatk-4.3.0.0/gatk CreateSequenceDictionary -R",
            fasta_file)
    system(command, intern = T)
    
  } else{
    message("Ya se ha creado el diccionario fasta")
  }
  
  
}


creacion_readgroup <-
  function(output_dir = output_dir,
           fastq_dir = fastq_dir) {
    ## el directorio del mapeo esta creado pero lo tenemos que guardar de nuevo
    mapping_output_dir <- file.path(output_dir, "mapping_output")
    fastq_files <- list.files(fastq_dir, full.names = F)
    output_file_name <-
      unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
    ## quitamos la extension del archivo
    output_file_name <- file_path_sans_ext(output_file_name)
    ## obtenemos el archivo marcado con duplicados
    mark_file <-
      file.path(mapping_output_dir,
                paste0(output_file_name, ".sorted.mark_dup.bam"))
    ## quitamos la extension y renombramos al archivo de salida
    out_file <- paste0(file_path_sans_ext(mark_file), "_RG.bam")
    
    ## si ya se ha creado el archivo con grupo
    
    if (!file.exists(out_file)) {
      command <-
        paste(
          "~/tools/gatk-4.3.0.0/gatk AddOrReplaceReadGroups I=",
          mark_file,
          "O=",
          out_file,
          "RGID=1 RGLB=lib2 RGPL=ILLUMINA RGPU=unit1 RGSM=1"
        )
      system(command = command, intern = T)
      
    } else{
      message("Ya estan los grupos")
    }
    
    
    
    
  }




base_recalibrator <-
  function(folder_fasta,
           output_dir,
           folder_data_gatk,
           fastq_dir) {
    ## llamamos al knwon sites file vcf
    known_sites_file <-
      list.files(folder_data_gatk,
                 pattern = ".vcf$",
                 full.names = T)
    ## llamamos al archivo fasta
    fasta_file <- fn_exists_fasta(folder_fasta)
    ## el directorio del mapeo esta creado pero lo tenemos que guardar de nuevo
    mapping_output_dir <- file.path(output_dir, "mapping_output")
    fastq_files <- list.files(fastq_dir, full.names = F)
    output_file_name <-
      unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
    ## quitamos la extension del archivo
    output_file_name <- file_path_sans_ext(output_file_name)
    ## obtenemos el archivo marcado con duplicados
    RG_file <-
      file.path(mapping_output_dir,
                paste0(output_file_name, ".sorted.mark_dup_RG.bam"))
    ## quitamos la extension y renombramos al archivo de salida
    out_file <- file.path(dirname(RG_file), "recal_data.table")
    ## si no existe la tabla de recalibracion, la calculamos
    
    if (!file.exists(out_file)) {
      command <-
        paste(
          "~/tools/gatk-4.3.0.0/gatk BaseRecalibrator -I",
          RG_file,
          " -R",
          fasta_file,
          " --known-sites",
          known_sites_file,
          " -O"  ,
          out_file
        )
      print(command)
      system(command = command, intern = T)
    } else{
      message("ya existe la tabla de base recalculator")
    }
    
    
  }


applybqsr <- function(folder_fasta, output_dir, fastq_dir) {
  ## llamamos al archivo fasta
  fasta_file <- fn_exists_fasta(folder_fasta)
  ## el directorio del mapeo esta creado pero lo tenemos que guardar de nuevo
  mapping_output_dir <- file.path(output_dir, "mapping_output")
  fastq_files <- list.files(fastq_dir, full.names = F)
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  ## quitamos la extension del archivo
  output_file_name <- file_path_sans_ext(output_file_name)
  ## obtenemos el archivo marcado con duplicados
  RG_file <-
    file.path(mapping_output_dir,
              paste0(output_file_name, ".sorted.mark_dup_RG.bam"))
  
  recal_data.table <-
    file.path(dirname(RG_file), "recal_data.table")
  
  out_file <-
    file.path(
      dirname(recal_data.table),
      paste0(output_file_name, ".sorted.mark_dup_RG_bqsr.bam")
    )
  
  if (!file.exists(out_file)) {
    command <-
      paste(
        " ~/tools/gatk-4.3.0.0/gatk ApplyBQSR -I",
        RG_file ,
        "-R",
        fasta_file,
        " --bqsr-recal-file",
        recal_data.table,
        " -O",
        out_file
      )
    print(command)
    system(command = command, intern = T)
  } else{
    message("Ya se ha aplicado el bsqr")
  }
  
}


bam_statistics <- function(folder_fasta, fastq_dir, output_dir) {
  ## el directorio del mapeo esta creado pero lo tenemos que guardar de nuevo
  mapping_output_dir <- file.path(output_dir, "mapping_output")
  fastq_files <- list.files(fastq_dir, full.names = F)
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  output_file_name <- file_path_sans_ext(output_file_name)
  
  ## bam file
  bam_file <-
    file.path(mapping_output_dir,
              paste0(output_file_name, ".sorted.mark_dup_RG_bqsr.bam"))
  ## directorio de salida
  out_dir <-
    paste0(file_path_sans_ext(bam_file), ".CollectMultipleMetrics")
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
  out_file <- file.path(out_dir, "CollectMultipleMetrics")
  ## obtenemos el archivo fasta
  fasta_file <- fn_exists_fasta(folder_fasta)
  ## verificacion de archivos en el directorio
  verificacion <- length(list.files(out_dir))
  ## comando estadistico
  
  if (verificacion < 2) {
    command <-
      paste(
        "~/tools/gatk-4.3.0.0/gatk CollectMultipleMetrics R=",
        fasta_file,
        "I=",
        bam_file,
        "O=",
        out_file
      )
    
    system(command = command, intern = T)
    ## ahora juntamos las estadisticas
    command <- paste("multiqc", out_dir, "-o", out_dir)
    
    system(command = command, intern = T)
    
  } else{
    message("Ya se ha calculado la estadistica bam")
  }
  
  
}



haplotype_caller <- function(output_dir, folder_fasta, fastq_dir) {
  ## recuperamos el archivo fasta
  fasta_file <- fn_exists_fasta(folder_fasta)
  ## recuperamos el ultimo archivo bam
  mapping_output_dir <- file.path(output_dir, "mapping_output")
  fastq_files <- list.files(fastq_dir, full.names = F)
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  output_file_name <- file_path_sans_ext(output_file_name)
  ## bam file
  bam_file <-
    file.path(mapping_output_dir,
              paste0(output_file_name, ".sorted.mark_dup_RG_bqsr.bam"))
  out_dir <- file.path(output_dir, "variantCalling")
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
  out_file <-
    file.path(out_dir, paste0(basename(output_file_name), ".g.vcf.gz"))
  
  if (!file.exists(out_file)) {
    command <-
      paste(
        "~/tools/gatk-4.3.0.0/gatk HaplotypeCaller -I",
        bam_file,
        "-R",
        fasta_file,
        "-ERC GVCF -O",
        out_file
      )
    system(command, intern = T)
  } else{
    message("Ya se han llamado a las variantes")
  }
  
  
}

genotypeGVCF <- function(folder_fasta, output_dir, fastq_dir) {
  ## recuperamos el archivo fasta
  fasta_file <- fn_exists_fasta(folder_fasta)
  ## obtenemos el nombre del archivo
  fastq_files <- list.files(fastq_dir, full.names = F)
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  output_file_name <- file_path_sans_ext(output_file_name)
  ## indicamos el archivo de entrada
  file_in <-
    file.path(output_dir,
              "variantCalling",
              paste0(output_file_name, ".g.vcf.gz"))
  ## indciamos el archivo de salida
  file_out <-
    file.path(
      output_dir,
      "variantCalling",
      paste0(output_file_name, "_apply_genotypeGVCF.g.vcf.gz")
    )
  
  if (!file.exists(file_out)) {
    command <-
      paste(
        "~/tools/gatk-4.3.0.0/gatk GenotypeGVCFs -R",
        fasta_file,
        "-V",
        file_in,
        "-O",
        file_out
      )
    system(command = command, intern = T)
  } else{
    message("Ya se ha calculado la probabilidad posterior de alelo no referente")
  }
  
  
}

variantRecallibrator <-
  function(fastq_dir,
           folder_fasta,
           folder_data_gatk,
           output_dir) {
    ## recuperamos el archivo fasta
    fasta_file <- fn_exists_fasta(folder_fasta)
    ## obtenemos el nombre del archivo
    fastq_files <- list.files(fastq_dir, full.names = F)
    output_file_name <-
      unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
    output_file_name <- file_path_sans_ext(output_file_name)
    
    in_file <-
      file.path(
        output_dir,
        "variantCalling",
        paste0(output_file_name, "_apply_genotypeGVCF.g.vcf.gz")
      )
    
    snps_recal_file <-
      file.path(
        output_dir,
        "variantCalling",
        paste0(output_file_name, "_apply_genotypeGVCF.recal")
      )
    
    tranches_file <-
      file.path(
        output_dir,
        "variantCalling",
        paste0(output_file_name, "_apply_genotypeGVCF.g.tranches")
      )
    
    if (!file.exists(snps_recal_file) |
        !file.exists(tranches_file)) {
      command <-
        paste(
          "~/tools/gatk-4.3.0.0/gatk VariantRecalibrator -V",
          in_file,
          " --resource:hapmap,known=false,training=true,truth=true,prior=15",
          file.path(
            folder_data_gatk,
            "resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz"
          ),
          "--resource:omni,known=false,training=true,truth=true,prior=12",
          file.path(
            folder_data_gatk,
            "resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz"
          ),
          "--resource:1000G,known=false,training=true,truth=false,prior=10",
          file.path(
            folder_data_gatk,
            "resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz"
          ),
          "--resource:dbsnp,known=true,training=false,truth=false,prior=7",
          file.path(
            folder_data_gatk,
            "resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"
          ),
          "-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP",
          "-O",
          snps_recal_file,
          "--tranches-file",
          tranches_file
        )
      system(command = command, intern = T)
      
    } else{
      message("Ya se ha hecho el pre recalibrado de variantes")
    }
    
    
  }

applyVQSR <- function(folder_fasta, fastq_dir, output_dir) {
  ## recuperamos el archivo fasta
  fasta_file <- fn_exists_fasta(folder_fasta)
  ## obtenemos el nombre del archivo
  fastq_files <- list.files(fastq_dir, full.names = F)
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  output_file_name <- file_path_sans_ext(output_file_name)
  
  in_file <-
    file.path(
      output_dir,
      "variantCalling",
      paste0(output_file_name, "_apply_genotypeGVCF.g.vcf.gz")
    )
  
  snps_recal_file <-
    file.path(
      output_dir,
      "variantCalling",
      paste0(output_file_name, "_apply_genotypeGVCF.recal")
    )
  
  tranches_file <-
    file.path(
      output_dir,
      "variantCalling",
      paste0(output_file_name, "_apply_genotypeGVCF.g.tranches")
    )
  
  out_file <-
    file.path(
      output_dir,
      "variantCalling",
      paste0(output_file_name, "_apply_genotypeGVCF.vqsr.vcf")
    )
  
  if (!file.exists(out_file)) {
    command <-
      paste(
        "~/tools/gatk-4.3.0.0/gatk ApplyVQSR -R",
        fasta_file,
        "-V",
        in_file,
        "-O",
        out_file,
        "--truth-sensitivity-filter-level 99.0 --tranches-file",
        tranches_file,
        "--recal-file",
        snps_recal_file,
        "-mode SNP"
      )
    system(command = command, intern = T)
    
  } else{
    message("Ya se ha aplicado VQSR")
  }
  
}
applyVQSR <- function(folder_fasta, fastq_dir, output_dir) {
  ## recuperamos el archivo fasta
  fasta_file <- fn_exists_fasta(folder_fasta)
  ## obtenemos el nombre del archivo
  fastq_files <- list.files(fastq_dir, full.names = F)
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  output_file_name <- file_path_sans_ext(output_file_name)
  
  in_file <-
    file.path(
      output_dir,
      "variantCalling",
      paste0(output_file_name, "_apply_genotypeGVCF.g.vcf.gz")
    )
  
  snps_recal_file <-
    file.path(
      output_dir,
      "variantCalling",
      paste0(output_file_name, "_apply_genotypeGVCF.recal")
    )
  
  tranches_file <-
    file.path(
      output_dir,
      "variantCalling",
      paste0(output_file_name, "_apply_genotypeGVCF.g.tranches")
    )
  
  out_file <-
    file.path(
      output_dir,
      "variantCalling",
      paste0(output_file_name, "_apply_genotypeGVCF.vqsr.vcf")
    )
  
  if (!file.exists(out_file)) {
    command <-
      paste(
        "~/tools/gatk-4.3.0.0/gatk ApplyVQSR -R",
        fasta_file,
        "-V",
        in_file,
        "-O",
        out_file,
        "--truth-sensitivity-filter-level 99.0 --tranches-file",
        tranches_file,
        "--recal-file",
        snps_recal_file,
        "-mode SNP"
      )
    system(command = command, intern = T)
    
  } else{
    message("Ya se ha aplicado VQSR")
  }
  
}


variantFiltration <- function(folder_fasta, output_dir, fastq_dir) {
  ## recuperamos el archivo fasta
  fasta_file <- fn_exists_fasta(folder_fasta)
  ## obtenemos el nombre del archivo
  fastq_files <- list.files(fastq_dir, full.names = F)
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  output_file_name <- file_path_sans_ext(output_file_name)
  
  in_file <-
    file.path(
      output_dir,
      "variantCalling",
      paste0(output_file_name, "_apply_genotypeGVCF.vqsr.vcf")
    )
  
  out_file <-
    file.path(
      output_dir,
      "variantCalling",
      paste0(
        output_file_name,
        "_apply_genotypeGVCF.vqsr.varfilter.vcf"
      )
    )
  
  if (!file.exists(out_file)) {
    command <-
      paste(
        "~/tools/gatk-4.3.0.0/gatk VariantFiltration -R",
        fasta_file,
        "-V",
        in_file,
        "-O",
        out_file,
        "--filter-name LOW_depth1  --filter-expression 'DP< 1'"
      )
    system(command = command, intern = T)
    
    
  } else{
    "Ya se ha filtrado"
  }
  
  
  
}


analysisReady <- function(folder_fasta, output_dir, fastq_dir) {
  ## obtenemos el nombre del archivo
  fastq_files <- list.files(fastq_dir, full.names = F)
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  output_file_name <- file_path_sans_ext(output_file_name)
  
  in_file <-
    file.path(
      output_dir,
      "variantCalling",
      paste0(
        output_file_name,
        "_apply_genotypeGVCF.vqsr.varfilter.vcf"
      )
    )
  out_file <-
    file.path(
      output_dir,
      "variantCalling",
      paste0(
        output_file_name,
        "_apply_genotypeGVCF.vqsr.varfilter.pass.vcf"
      )
    )
  
  if (!file.exists(out_file)) {
    command <-
      paste("bcftools view -f 'PASS,.' -O v -o", out_file, in_file)
    system(command = command, intern = T)
    
  } else{
    "ya se ha hecho el PASS filter"
  }
  
  
}


anotation <-
  function(folder_fasta,
           path_snpeff,
           output_dir,
           fastq_dir) {
    ## obtenemos el nombre del archivo
    fastq_files <- list.files(fastq_dir, full.names = F)
    output_file_name <-
      unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
    output_file_name <- file_path_sans_ext(output_file_name)
    anotacion_dir <- file.path(output_dir, "anotation")
    if (!dir.exists(anotacion_dir)) {
      dir.create(anotacion_dir)
    }
    
    in_file <-
      file.path(
        output_dir,
        "variantCalling",
        paste0(
          output_file_name,
          "_apply_genotypeGVCF.vqsr.varfilter.pass.vcf"
        )
      )
    output_file_anno1 <-
      file.path(anotacion_dir, paste0(output_file_name, "anno_snpeff.vcf"))
    output_file_anno1_1 <-
      file.path(anotacion_dir,
                paste0(output_file_name, "anno_snpeff2.vcf"))
    output_file_anno2 <-
      file.path(anotacion_dir,
                paste0(output_file_name, "anno_snpeff2_clinvar.vcf"))
    
    output_file_anno3 <-
      file.path(anotacion_dir,
                paste0(output_file_name, "anno_snpeff2_clinvar_freqs.vcf"))
    
    output_file_anno4 <-
      file.path(anotacion_dir,
                paste0(output_file_name, "anno_snpeff2_clinvar_freqs_gwas.vcf"))
    output_file_anno5 <-
      file.path(
        anotacion_dir,
        paste0(output_file_name, "anno_snpeff2_clinvar_freqs_gwas_2.vcf")
      )
    
    
    
    comando1 <-
      paste(
        "java -Xmx32g -jar",
        file.path(path_snpeff, "snpEff.jar"),
        "hg38 -v ",
        in_file,
        " > ",
        output_file_anno1
      )
    comando2 <-
      paste(
        "java -Xmx32g -jar",
        file.path(path_snpeff, "SnpSift.jar"),
        " varType -v ",
        output_file_anno1,
        " > ",
        output_file_anno1_1
      )
    comando3 <-
      paste(
        "java -Xmx32g -jar",
        file.path(path_snpeff, "SnpSift.jar"),
        " annotate -v ~/datos_exomas/datos_clinvar/clinvar.vcf.gz",
        output_file_anno1_1,
        " > ",
        output_file_anno2
      )
    comando4 <-
      paste(
        "java -Xmx32g -jar",
        file.path(path_snpeff, "SnpSift.jar"),
        " dbnsfp -v -db ~/datos_exomas/datos_dbsnp/dbNSFP5.1a_grch38.gz",
        output_file_anno2,
        ">",
        output_file_anno3
      )
    comando5 <-
      paste(
        "java -Xmx32g -jar",
        file.path(path_snpeff, "SnpSift.jar"),
        " gwasCat -db ~/datos_exomas/gwas/gwascatalog.txt",
        output_file_anno3,
        ">",
        output_file_anno4
      )
    campos <-
      "aaref,aaalt,rs_dbSNP,HGVSc_snpEff,HGVSp_snpEff,APPRIS,M-CAP_pred,CADD_phred,clinvar_OMIM_id,clinvar_Orphanet_id,clinvar_MedGen_id,AlphaMissense_pred,Reliability_index"
    comando6 <-
      paste(
        "java -Xmx32g -jar",
        file.path(path_snpeff, "SnpSift.jar"),
        " dbnsfp  -v -db ~/datos_exomas/datos_dbsnp/dbNSFP5.1a_grch38.gz -f",
        campos,
        output_file_anno4,
        ">",
        output_file_anno5
      )
    
    if (!file.exists(output_file_anno1)) {
      print(comando1)
      system(command = comando1, intern = T)
      
    }
    if (!file.exists(output_file_anno1_1)) {
      print(comando2)
      system(command = comando2, intern = T)
      
      
    }
    
    
    if (!file.exists(output_file_anno2)) {
      print(comando3)
      system(command = comando3, intern = T)
      
      
    }
    if (!file.exists(output_file_anno3)) {
      print(comando4)
      system(command = comando4, intern = T)
    }
    
    if (!file.exists(output_file_anno4)) {
      print(comando5)
      system(command = comando5, intern = T)
    }
    
    if (!file.exists(output_file_anno5)) {
      print(comando6)
      system(command = comando6, intern = T)
    }
    
  print("YA ESTA ANOTADO")
    
  }


obtener_exoma_overlap <- function(bd_list, exoma_df, cromosomas, muestra) {
  # 1) Unir bd_list en un solo data.frame, convirtiendo todo a caracteres y filtrando cromosomas
  cat("\nCheckpoint 1: Uniendo bd_list. Cols bd_list[[1]]:", paste(names(bd_list[[1]]), collapse = ", "), "\n")
  
  bd_combined <-
    lapply(bd_list, function(x) {
      df <- as.data.frame(lapply(x, as.character), stringsAsFactors = FALSE)
      df[df$Chr %in% cromosomas, ]
    }) %>%
    bind_rows()
  cat("Checkpoint 2: bd_combined columns:", paste(names(bd_combined), collapse = ", "), "\n")
  
  # *** Aquí verifica la existencia de Gene.refGene ***
  if (!"Gene.refGene" %in% names(bd_combined)) {
    stop("La columna 'Gene.refGene' NO existe. Las columnas disponibles son:\n", 
         paste(names(bd_combined), collapse = ", "))
  }
  
  # 2) Agrupar por Chr, Start, End, Gene.refGene; concatenar "codigo" y contar N
    
  bd_grouped <- bd_combined %>%
  group_by(Chr, Start, End, Gene.refGene) %>%
  summarise(
    n = n(),
    gene_name = dplyr::first(Gene.refGene), # Aquí
    paste_m = toString(codigo),
    .groups = "drop"
  ) %>%
  mutate(
    paste_m = sapply(strsplit(paste_m, ",\\s*"), function(v) paste(unique(v), collapse = ",")),
    N = sapply(strsplit(paste_m, ","), length)
  )
# Ya no necesitas rename(gene_name = Gene.refGene)
cat("Checkpoint 3: bd_grouped columns:", paste(names(bd_grouped), collapse = ", "), "\n")

  cat("Checkpoint 3: bd_grouped columns:", paste(names(bd_grouped), collapse = ", "), "\n")
  
  # 3) Filtrar únicamente las filas que contienen la muestra de interés
  bd_filtrado_muestra <- bd_grouped[grep(muestra, bd_grouped$paste_m), ]
  cat("Checkpoint 4: Filtrado muestra. Rows:", nrow(bd_filtrado_muestra), "\n")
  
  # 4) Asegurar que Start y End sean numéricos
  bd_filtrado_muestra <- bd_filtrado_muestra %>%
    mutate(Start = as.numeric(Start), End   = as.numeric(End))
  cat("Checkpoint 5: bd_filtrado_muestra columns:", paste(names(bd_filtrado_muestra), collapse = ", "), "\n")
  
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
subject_ids <- paste0(seqnames(gr_exoma)[subjectHits(hits)], ":", start(gr_exoma)[subjectHits(hits)], "-", end(gr_exoma)[subjectHits(hits)])

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

cat("Duplicados eliminados en hits:", nrow(hits_df) - nrow(hits_df_unique), "\n")

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

cat("Checkpoint 7: resultados columns:", paste(names(resultados), collapse = ", "), "\n")
  resultados_unicos <- resultados[!duplicated(resultados[, c("gene_name", "Start", "End", "CHROM")]), ]
  resultados_unicos$N <- sapply(strsplit(resultados_unicos$samples, ","), length)
  colnames(resultados_unicos)[colnames(resultados_unicos) == "Start"]   <- "POS"
  colnames(resultados_unicos)[colnames(resultados_unicos) == "End"]     <- "END"
  
  primeros <- c("CHROM", "POS", "END", "gene_name", "N", "samples")
  restantes <- setdiff(names(resultados_unicos), primeros)
  bd_exoma_overlap <- resultados_unicos[, c(primeros, restantes)]
  cat("Checkpoint 8: bd_exoma_overlap columns:", paste(names(bd_exoma_overlap), collapse = ", "), "\n")
  
  return(bd_exoma_overlap)
}

process_vcf_to_table <- function(folder_fasta,
                                 output_dir,
                                 fastq_dir,
                                 muestra,
                                 db) {
  print("=== Iniciando process_vcf_to_table ===")
  fastq_files <- list.files(fastq_dir, full.names = F)
  print(paste("FASTQ files encontrados:", paste(fastq_files, collapse = ", ")))
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  output_file_name <- file_path_sans_ext(output_file_name)
  anotacion_dir <- file.path(output_dir, "anotation")
  output_file_anno5 <-
    file.path(anotacion_dir, paste0(output_file_name, "anno_snpeff2_clinvar_freqs_gwas_2.vcf"))
  
  post_process_dir <- file.path(output_dir, "post_process_results")
  if(!dir.exists(post_process_dir)){
    dir.create(post_process_dir,recursive = T)
    print(paste("Directorio creado:", post_process_dir))
  } else {
    print(paste("Directorio ya existe:", post_process_dir))
  }
  
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
  final_vcf_df$ALT <- lapply(final_vcf_df$ALT, function(x) as.character(x[[1]]))
  final_vcf_df <- as.data.frame(lapply(final_vcf_df, as.character))
  print("Unión de dataframes realizada.")
  
  # --- Eliminar columnas específicas
  print("Eliminando columnas no deseadas (1ra ronda)...")
  remove_columns <- c(
    "strand", "FILTER", "BaseQRankSum...14", "ExcessHet...17",
    "FS...18", "MLEAC...20", "MLEAF...21", "MQ...22", "MQRankSum...23",
    "NEGATIVE_TRAIN_SITE...24", "POSITIVE_TRAIN_SITE...25", "RAW_MQandDP...27",
    "ReadPosRankSum...28", "SOR...29", "VQSLOD...30", "culprit...31",
    "ANN", "HOM"
  )
  busca_exoma <- function(terms, cols)
    unique(unlist(sapply(terms, function(x)
      grep(x, cols, ignore.case = TRUE))))
  final_vcf_df <- final_vcf_df[, -busca_exoma(remove_columns, colnames(final_vcf_df)), drop = FALSE]
  print("Columnas eliminadas (1ra ronda).")
  
  # --- Columnas dbNSFP y frecuencias
  print("Extrayendo y limpiando columnas dbNSFP y de frecuencias...")
  cols_dbnsfp <- c(
    "dbNSFP_rs_dbSNP", "dbNSFP_clinvar_OMIM_id", "dbNSFP_clinvar_MedGen_id",
    "dbNSFP_HGVSc_snpEff", "dbNSFP_HGVSp_snpEff", "dbNSFP_clinvar_Orphanet_id",
    "dbNSFP_Reliability_index", "dbNSFP_AlphaMissense_pred", "dbNSFP_SIFT_pred",
    "dbNSFP_Polyphen2_HVAR_pred", "dbNSFP_MutationTaster_pred",
    "dbNSFP_Polyphen2_HDIV_pred", "dbNSFP_CADD_phred", "dbNSFP_Uniprot_acc",
    "dbNSFP_Interpro_domain"
  )
  cols_freqs  <- c(
    "dbNSFP_1000Gp3_SAS_AF", "dbNSFP_1000Gp3_AFR_AF", "dbNSFP_1000Gp3_EUR_AF",
    "dbNSFP_1000Gp3_EAS_AF", "dbNSFP_1000Gp3_AF", "dbNSFP_1000Gp3_AMR_AF"
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
    "END...16", "InbreedingCoeff...19", "SNP...36", "INS...38", "DEL...39",
    "MC...62", "AF_EXAC...68", "AF_ESP...70", "AF_TGP...75", "GWASCAT_TRAIT...107"
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
    "BaseQRankSum", "ExcessHet", "END", "FS", "InbreedingCoeff", "MLEAC", "MLEAF",
    "MQ", "MQRankSum", "NEGATIVE_TRAIN_SITE", "POSITIVE_TRAIN_SITE", "RAW_MQandDP",
    "ReadPosRankSum", "SOR", "VQSLOD", "culprit", "SNP", "INS", "DEL", "MC",
    "AF_EXAC", "AF_ESP", "AF_TGP"
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
  colnames(geno_df_) <- c("GT", "AD", "DP", "GQ", "MIN_DP", "PGT", "PID", "PL", "PS", "RGQ")
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
    "alterno_quitar", "effect", "impact", "gene_name", "gene_name_quitar",
    "effect_quitar", "annotation_id", "gene_biotype", "exon_intron_rank", "nt_change",
    "aa_change", "cDNA_position.cDNA_len", "nt_position", "aa_position", "distance_to_feature",
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
    "CHROM", "START", "END", "width", "paramRangeID", "REF", "ALT", "QUAL", "GT", "AD", "DP",
    "GQ", "MIN_DP", "PGT", "PID", "PL", "PS", "RGQ", "freqs", "AC_1", "AF_1", "AN_1", "DP_1", "QD_1",
    "LOF_1", "NMD_1", "VARTYPE_1", "MNP_1", "MIXED_1", "HET_1", "DBVARID_1", "SCISCV_1", "ALLELEID_1",
    "rs_dbSNP", "HGVSc_snpEff", "HGVSp_snpEff", "gene_name", "annotation_id", "gene_biotype", "exon_intron_rank",
    "nt_change", "aa_change", "cDNA_position.cDNA_len", "nt_position", "aa_position", "distance_to_feature",
    "CLNSIG_1", "CLNVCSO_1", "SCIDNINCL_1", "CLNREVSTAT_1", "ONCREVSTAT_1", "CLNDNINCL_1", "ONC_1",
    "CLNSIGSCV_1", "ORIGIN_1", "ONCINCL_1", "ONCDNINCL_1", "ONCDISDB_1", "SCIREVSTAT_1", "ONCDISDBINCL_1",
    "ONCSCV_1", "CLNDN_1", "ONCCONF_1", "CLNVC_1", "SCIDISDB_1", "CLNVI_1", "ONCDN_1", "CLNSIGINCL_1",
    "CLNDISDB_1", "GENEINFO_1", "CLNDISDBINCL_1", "CLNSIGCONF_1", "SCIDISDBINCL_1", "CLNHGVS_1", "SCIINCL_1",
    "SCIDN_1", "SCI_1", "clinvar_OMIM_id", "clinvar_MedGen_id", "clinvar_Orphanet_id", "Reliability_index",
    "AlphaMissense_pred", "SIFT_pred", "Polyphen2_HVAR_pred", "MutationTaster_pred", "Polyphen2_HDIV_pred",
    "CADD_phred", "Uniprot_acc", "Interpro_domain", "effect", "impact", "errors"
  )
  
  
  #df_limpio <- df_limpio[, orden_columnas[orden_columnas %in% names(df_limpio)], drop = FALSE]
  print("#### columnas df_limpio ########")
  print(colnames(df_limpio))
  print("#### ahora vemos cuales de orden_columnas estan en df_limpio ###")
  print(length(orden_columnas))
  print(orden_columnas)
  print(length(colnames(df_limpio)))
  w <- which(colnames(df_limpio) %in% orden_columnas))
  print(w)
  length(w)
  
  df_limpio <- df_limpio[,w]
  
  print(colnames(df_limpio))
  print(orden_columnas)
  names(df_limpio) <- sub("_1$", "", names(df_limpio))
  df_limpio <- df_limpio[, names(df_limpio) != "DP.1", drop = FALSE]
  print("Columnas finales reordenadas y renombradas.")
  
  # --- Limpieza y formateo final
  print("Limpieza y filtrado final de errores y formatos...")
  df_limpio <- df_limpio %>%
    mutate(
      errors = gsub('["\'\\\\]', '', errors),
      errors = trimws(errors)
    ) %>%
    filter(is.na(errors) | errors == "" | errors == "INFO_REALIGN_3_PRIME")
  df_limpio$AD <- sapply(df_limpio$AD, function(x) paste(x, collapse = "|"))
  df_limpio$PL <- sapply(df_limpio$PL, function(x)
    if (all(is.na(x))) NA_character_ else paste(x, collapse = "|"))
  df_limpio <- bind_cols(SAMPLE = muestra, df_limpio)
  
  print("Primeras filas del dataframe final antes de overlap:")
  print(head(df_limpio))
  
  # --- Overlap con exoma
  print("Aplicando overlap con regiones de exoma...")
  bd_list_ <- readRDS(db)
  cromosomas_ <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
  
  df_limpio <- obtener_exoma_overlap(
    bd_list = bd_list_,
    exoma_df = df_limpio,
    cromosomas = cromosomas_,
    muestra = muestra
  )
  
  print("Overlap con exoma realizado.")
  print("=== Proceso COMPLETO: process_vcf_to_table ===")
  write.csv(df_limpio,"./df_limpio.csv")
  #return(df_limpio)
}









#
# compute_depth <- function(fastq_dir, output_dir) {
#   ## el directorio del mapeo esta creado pero lo tenemos que guardar de nuevo
#   mapping_output_dir <- file.path(output_dir, "mapping_output")
#   fastq_files <- list.files(fastq_dir, full.names = F)
#   output_file_name <-
#     unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
#   ## quitamos la extension del archivo
#   output_file_name <- file_path_sans_ext(output_file_name)
#
#   bam_file <-
#     file.path(mapping_output_dir,
#               paste0(output_file_name, ".sorted.mark_dup_RG_bqsr.bam"))
#
#   dir_coverage <- file.path(output_dir, "coverage_and_stats")
#   if (!dir.exists(dir_coverage)) {
#     dir.create(dir_coverage)
#   }
#   outfile_coverage <- file.path(dir_coverage, "coverage.txt")
#   comando <- paste("./compute_depth.sh", bam_file, ">", outfile_coverage)
#   if (!file.exists(outfile_coverage)) {
#     print("computando cobertura...")
#     print(comando)
#     system(comando, intern = T)
#
#   } else{
#     print("la cobertura ya esta computada")
#   }
#
# }
#
# compute_stats <- function(fastq_dir, output_dir, muestra) {
#   dir_coverage <- file.path(output_dir, "coverage_and_stats")
#   if (!file.exists(file.path(dir_coverage, "stats.csv"))) {
#     dir_coverage <- file.path(output_dir, "coverage_and_stats")
#     cov_file <- file.path(dir_coverage, "coverage.txt")
#     process_dir <- file.path(output_dir, "post_process_results")
#     exoma_file <- file.path(process_dir, "file_ready_analysis.csv")
#
#     cov_data <- read.delim(cov_file, header = F)
#     exoma <- read.csv(exoma_file, na.strings = ".")
#     cromosomas <- c(paste0("chr", 1:22), "chrX", "chrY")
#     exoma <- exoma[exoma$CHROM %in% cromosomas, ]
#     interest.data <- cov_data[cov_data$V1 !=
#                                 "all", "V7"]
#     mean_coverage <- round(mean(interest.data, na.rm = T), 2)
#     mean_coverage20 <- round(mean(interest.data[interest.data > 20], na.rm =
#                                     T), 2)
#     cov.data <- cov_data
#
#     cov.data$V4 <- gsub("_.*", "", cov.data$V4)
#     cov.data <- cov.data[cov.data$V1 != "all", ]
#
#     cov.data <- aggregate(V7 ~ V4, data = cov.data, mean)
#
#     genes_totales <- length(unique(cov.data$V4))
#
#     genes_20 <- length(cov.data[which(cov.data$V7 > 20), "V4"])
#
#     genes_por_20 <- round(100 * genes_20 / genes_totales, 2)
#
#
#     genes_mayor_media <- length(cov.data[which(cov.data$V7 > mean_coverage), "V4"])
#
#     genes_por_media <- round(100 * genes_mayor_media / genes_totales, 2)
#
#
#     variantes_totales.df <- as.data.frame(lapply(exoma[, c("POS", "DP")], as.numeric))
#
#     variantes_dp <- aggregate(DP ~ POS, data = variantes_totales.df, mean)
#
#     variantes_totales <- length(unique(variantes_dp$POS))
#     variantes_20 <- length(unique(variantes_totales.df[which(variantes_totales.df$DP >
#                                                                20), "POS"]))
#     variantes_20_x <- round(100 * variantes_20 / variantes_totales, 2)
#
#     variantes_media <- length(unique(variantes_totales.df[which(variantes_totales.df$DP >
#                                                                   mean_coverage), "POS"]))
#     variantes_media_x <- round(100 * variantes_media / variantes_totales, 2)
#
#     res <- as.data.frame(
#       c(
#         mean_coverage = mean_coverage,
#         mean_coverage20 = mean_coverage20,
#         genes_totales = genes_totales,
#         genes_20 = genes_20,
#         genes_por_20 = genes_por_20,
#         genes_mayor_media = genes_mayor_media,
#         genes_por_media = genes_por_media,
#         variantes_totales = variantes_totales,
#         variantes_20 = variantes_20,
#         variantes_20_x = variantes_20_x,
#         variantes_media = variantes_media,
#         variantes_media_x
#       )
#     )
#     res$Descripcion <- c(
#       "Cobertura media",
#       "Cobertura media > 20X",
#       "Genes totales",
#       "Genes > 20 X",
#       "Genes % > 20X",
#       "Genes > media",
#       "Genes % > media",
#       "Variantes totales",
#       "Variantes > 20X",
#       "Variantes % >20X",
#       "Variantes > media X",
#       "Variantes % > media X"
#     )
#     res <- res[, c(2, 1)]
#     colnames(res) <- c("Descripcion", "Stats")
#
#
#     gcov <- cov_data[cov_data[, 1] == 'all', ]
#     ###
#     longitud <- 300
#     datos.pre <-
#       data.frame(
#         X = gcov[1:longitud, 2],
#         Y = 100 * (1 - cumsum(gcov[1:longitud, 5])),
#         Z = 100 * gcov[1:longitud, 5],
#         relleno = gcov[1:longitud, 1]
#       )
#     datos.pre$relleno <- as.factor(datos.pre$relleno)
#     p1 <-
#       ggplot(data = datos.pre, aes(X, Y, fill = relleno)) + geom_line(color =
#                                                                         "steelblue", linewidth = 2) + xlab("Profundidad de Cobertura") + ylab("Porcentaje de la region >= Profunidad") + theme(legend.position = "none") +
#       theme_classic()
#     p2 <-
#       ggplot(data = datos.pre, aes(X, Z, fill = relleno)) + geom_col() + scale_fill_discrete(type =
#                                                                                                "steelblue") + xlab("Profundidad de Cobertura") + ylab("Porcentaje de la region") + theme(legend.position = "none") + theme_classic()
#
#     p3 <- ggpubr::ggarrange(p1, p2, ncol = 2)
#
#     figura_file_name <- file.path(dir_coverage, "cobertura.jpeg")
#
#     p4 <-
#       ggpubr::annotate_figure(p3,
#                               top = paste(
#                                 muestra,
#                                 "Profundidad media de cobertura:",
#                                 mean_coverage,
#                                 "X"
#                               ))
#     ggsave(filename = figura_file_name, plot = p4)
#     res_file <- file.path(dir_coverage, "stats.csv")
#     write.csv(res, res_file)
#
#     blah <- paste(
#       "<p>Esta muestra se ha estudiado por el metodo de secuenciacion masiva en paralelo del exoma completo. Se analizaron",
#       res[3, 2],
#       "genes. La sensibilidad y la especificidad del metodo son superiores al 98% (SNV< 20 bp INDELS). El porcentaje de genes con una cobertura mayor a 20X es de",
#       round(as.numeric(res[5, 2]), 2),
#       "%. De todas las variantes identificadas, que son un total de",
#       res[8, 2],
#       ",",
#       res[9, 2],
#       "tienen una cobertura mayor a 20X, esto significa un",
#       round(as.numeric(res[10, 2], 2)),
#       "%.</p>"
#     )
#
#     res <- print(xtable(res), type = "html")
#
#     fileconn <- "./aux1.html"
#     writeLines(blah, fileconn)
#     fileconn <- "./aux2.html"
#     writeLines(res, fileconn)
#
#     command <- paste("cat aux1.html aux2.html > ",
#                      file.path(dir_coverage, "doc.html"))
#
#     system(command)
#
#     command <- paste(
#       "pandoc --output",
#       file.path(dir_coverage, "reporte.docx"),
#       file.path(dir_coverage, "doc.html")
#     )
#
#     system(command)
#
#
#
#
#
#   } else{
#     print("ya se computaron las estadisticas")
#   }
# }






start <- Sys.time()
args <- commandArgs(trailingOnly = T)
args <- unlist(strsplit(args, " "))
muestra <- args
# for (muestra in muestras) {
  output_dir <- file.path("~/pipeline/exomas", muestra, "output_dir")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  hpo_file <- "~/datos_exomas/data_pipeline/genes_to_phenotype.txt"
  fastq_dir <- file.path("~/pipeline/exomas", muestra, "fastqfiles")
  folder_fasta <-
    file.path("~/datos_exomas/datos_gatk/hg38")
  folder_data_gatk <- file.path("~/datos_exomas/datos_gatk")
  path_snpeff <- "~/tools/snpEff/"
  bd_data <- "./bd.rds"
  
  control_calidad(fastq_dir, output_dir)
  index_bwa(folder_fasta)
  bwamem(fastq_dir, folder_fasta)
  markdups(output_dir = output_dir, fastq_dir = fastq_dir)
  ## creamos diccionario
  create_dict(folder_fasta)
  ## anadimos reaad group
  creacion_readgroup(output_dir, fastq_dir)
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
  variantRecallibrator(fastq_dir, folder_fasta, folder_data_gatk, output_dir)
  ## apply VQSR
  applyVQSR(folder_fasta, fastq_dir, output_dir)
  ## primer filtraje
  variantFiltration(folder_fasta, output_dir, fastq_dir)
  ## preparamos el archivo listo para elanalisis
  analysisReady(folder_fasta, output_dir, fastq_dir)
  ##
  anotation(folder_fasta, path_snpeff, output_dir, fastq_dir)
  
  valorar_Df <- process_vcf_to_table(folder_fasta = folder_fasta,
                                     output_dir = output_dir,
                                     fastq_dir = fastq_dir,
                                     muestra = muestra,
                                     db = bd_data)
  
  # process(output_dir, fastq_dir, hpo_file, muestra)
  #
  # compute_depth(fastq_dir = fastq_dir, output_dir)
  #
  # unique_variants(output_dir, muestra, bd_data)
  #
  # compute_stats(fastq_dir, output_dir, muestra)
  
  
# }
print(Sys.time() - start)
