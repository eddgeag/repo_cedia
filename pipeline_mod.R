load.libs <- c(
  "data.table",
  "limma",
  "grid",
  "egg",
  "Rbwa",
  "ggpubr",
  "tools",
  "DOSE",
  "GO.db",
  "GSEABase",
  "org.Hs.eg.db",
  "clusterProfiler",
  "dplyr",
  "tidyr",
  "ggplot2",
  "stringr",
  "RColorBrewer",
  "rWikiPathways",
  "RCy3",
  "RISmed",
  "ggplot2",
  "dplyr",
  "xtable"
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

	command <- paste("~/tools/gatk-4.3.0.0/gatk MarkDuplicates -CREATE_INDEX true -INPUT",
  			bam_file,"-VALIDATION_STRINGENCY STRICT --REMOVE_DUPLICATES true -OUTPUT",
			mark_file,"-M",
			metrics_file)

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
          "RGID=1 RGLB=lib2 RGPL=illumina RGPU=unit1 RGSM=1"
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
        " dbnsfp -v -db ~/datos_exomas/datos_dbsnp/dbNSFP4.1a.txt.gz",
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
      "aaref,aaalt,rs_dbSNP151,HGVSc_snpEff,HGVSp_snpEff,APPRIS,M-CAP_pred,CADD_phred,GTEx_V8_gene,GTEx_V8_gene,GTEx_V8_tissue,Geuvadis_eQTL_target_gene,Reliability_index"
    comando6 <-
      paste(
        "java -Xmx32g -jar",
        file.path(path_snpeff, "SnpSift.jar"),
        " dbnsfp  -v -db ~/datos_exomas/datos_dbsnp/dbNSFP4.1a.txt.gz -f",
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
    
    ### convertir a tsv
    
    comando <-
      paste(
        "vk vcf2tsv wide --print-header --ANN",
        output_file_anno5,
        ">",
        file.path(anotacion_dir, "anotacion_completa.tsv")
      )
    if (!file.exists(file.path(anotacion_dir, "anotacion_completa.tsv"))) {
      system(command = comando, intern = T)
      
    } else{
      print("Ya esta anotado todo, falta procesar")
    }
    
    
    
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

funcion_insilico <- function(X.) {
  X. <- gsub("[.]", "", X.)
  X. <- gsub(",", "", X.)
  X. <- (lapply(X., function(X.)
    unique(strsplit(X., "")[[1]])))
  X. <- unlist(lapply(X., function(X.)
    paste0(X., collapse = "")))
  X_conf <- ifelse(nchar(X.) > 1 , "Conflictive", X.)
  X_conf <- ifelse(!X. == "NA", X_conf, NA)
  return(X.)
}
process <- function(output_dir, fastq_dir, hpo_file, muestra) {
  fastq_files <- list.files(fastq_dir, full.names = F)
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  output_file_name <- file_path_sans_ext(output_file_name)
  
  process_dir <- file.path(output_dir, "post_process_results")
  if (!dir.exists(process_dir)) {
    dir.create(process_dir)
  }
  if (!file.exists(file.path(process_dir, "file_raw_analysis.csv"))) {
    ##
    print("procesando vcf a csv...")
    anotacion_dir <- file.path(output_dir, "anotation")
    
    vcf_file <-
      file.path(
        anotacion_dir,
        paste0(
          output_file_name,
          "anno_snpeff2_clinvar_freqs_gwas_2.vcf"
        )
      )
    
    tsv_file <-  file.path(anotacion_dir, "anotacion_completa.tsv")
    
    ## escaneamos archivo tsv
    
    X <- read.delim(tsv_file, na.strings = ".")
    
    X  <- X[, colSums(is.na(X)) < nrow(X)]
    
    ## escaneamos vcf
    
    vcf_scan <- readLines(vcf_file)
    header <- vcf_scan[grep("^#", vcf_scan)]
    data <- vcf_scan[grep("^#", invert = T, vcf_scan)]
    first <- as.data.table(read.delim(vcf_file,
                                      skip = length(header),
                                      header = F))
    
    looklof <- first[grep("LOF", first$V8), c("V2", "V8")]
    looklof <- cbind(looklof$V2, strsplit2(looklof[, V8], "LOF"))
    looklof <- looklof[,-2]
    looklof <- looklof[, 1:2]
    looklof <-
      cbind(looklof[, 1], strsplit2(looklof[, 2], ";"))[, 1:2]
    looklof[, 2] <- gsub("=", "", looklof[, 2])
    colnames(looklof) <- c("POS", "LOF")
    ## errores
    errors <- X$LOF
    X$LOF <- NA
    ## procesamos el archivo tsv
    X <- as.data.frame(X)
    X$POS <- as.numeric(X$POS)
    looklof <- as.data.frame(looklof)
    looklof$POS <- as.numeric(looklof$POS)
    X <- left_join(X, looklof, by = "POS")
    X <- X[,-grep("LOF.x", colnames(X))]
    colnames(X)[which(colnames(X) %in% c("protein_position", "distance_to_feature"))] <-
      c("nt_position", "protein_position")
    colnames(X)[which(colnames(X) == "error")] <-
      "distance_to_feature"
    X <- X[, -grep("_AC$", colnames(X))]
    X <- X[, -which(colnames(X) %in% c("SNP", "DEL", "INS"))]
    X <- X[, -which(colnames(X) == "X1_DP")]
    ## computamos frecuencias
    frecuencias <-
      as.data.frame(lapply(X[, c(grep("^AF_", colnames(X)), grep("_AF$", colnames(X)))], as.numeric))
    freqs <- rowMeans(frecuencias, na.rm = T)
    X <- X[, -which(colnames(X) %in% colnames(frecuencias))]
    X$freq <- freqs
    
    predicciones <- X[, grep("_pred$", colnames(X))]
    
    predicciones <-
      bind_cols(lapply(predicciones, funcion_insilico))
    
    X <- X[, -grep("_pred$", colnames(X))]
    X <- bind_cols(X, predicciones)
    X$END <- ifelse(nchar(X$ALT) > 1, X$POS + nchar(X$ALT), X$POS)
    
    ## buscamos el campo AD dentro del vcf para el computo AB
    
    colnames(first) <-
      c("CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "ANN",
        "colGT",
        "GT")
    
    colGT.unique <- unique(first$colGT)
    
    colGT.unique <- paste0(colGT.unique, collapse = ":")
    
    colGT.unique <- unique(unlist(strsplit(colGT.unique, ":")))
    
    gtfiedls <- strsplit2(first$GT, ":")
    
    indx  <- which(apply(gtfiedls, 1, function(x)
      all(x[6:8] == "")) == T)
    
    gtfiedls[-indx,] <- gtfiedls[-indx, c(1:4, 7, 5, 6, 8)]
    
    colnames(gtfiedls) <- colGT.unique
    
    gtfields <- bind_cols(POS = first$POS, AD = gtfiedls[, 2])
    X$errors <- errors
    X <- X[which(!duplicated(X$POS)),]
    gtfields <- gtfields[which(!duplicated(gtfields$POS)),]
    X <-
      inner_join(X, gtfields, "POS")
    
    ad <- strsplit2(X$AD, ",")
    ad[, 3] <- ifelse(ad[, 3] == "", 0, ad[, 3])
    ad <- apply(ad, 2, as.numeric)
    ab <- ad[, 1] / rowSums(ad) ## computamos alelic balance
    X$AB <- ab
    
    ## obtenemos campos de interes
    
    INFO <-
      c(
        "SAMPLE",
        "errors",
        "freq",
        "AB",
        "CHROM",
        "POS",
        "END",
        "HOM",
        "HET",
        "FILTER",
        "VARTYPE",
        "DP",
        "AF",
        "LOF",
        "NMD",
        "RS",
        "REF",
        "ALT",
        "ORIGIN",
        "gene_name",
        "gene_id",
        "feature_type",
        "feature_id",
        "allele",
        "effect",
        "impact",
        "transcript_biotype",
        "exon_intron_rank",
        "nt_change",
        "distance_to_feature",
        "aa_change",
        "nt_position",
        "cDNA_position.cDNA_len",
        "protein_position",
        "ID",
        "AC",
        "AD",
        "AN",
        "QD",
        "MQ",
        "MC",
        "BaseQRankSum",
        "MQRankSum",
        "SOR",
        "VQSLOD",
        "MLEAF",
        "MLEAC",
        "FS",
        "ExcessHet",
        "ReadPosRankSum",
        "NEGATIVE_TRAIN_SITE",
        "POSITIVE_TRAIN_SITE",
        "GENEINFO"
      )
    
    gwas <- colnames(X)[grep("GWAS", colnames(X))]
    clinvar <-
      c(colnames(X)[44], colnames(X)[grep("CLN", colnames(X))])
    inutiles <- colnames(X)[61:65]
    colnames(X)[grep("LOF", colnames(X))] <- "LOF"
    # interes <- colnames(X)[72:81]
    insilico <- colnames(X)[c(grep("_pred$", colnames(X)),
                              grep("CADD_phred", colnames(X)))]
    eqtl <- colnames(X)[grep("eqtl",ignore.case=T,x=colnames(X))]
    # genotipo <- colnames(X)[c(82,90,41:43,1,2,101,3:39,83:88,89,40,103)]
    
    genotipo <- colnames(X)[grep("X1_", colnames(X))]
    
    X <- X[, c(INFO, genotipo, gwas, clinvar, inutiles, insilico,eqtl)]
    
    colnames(X)[grep("X1_", colnames(X))] <-
      gsub("X1_", "", colnames(X[grep("X1_", colnames(X))]))
    
    hpo <- read.delim(hpo_file,
                      skip = 1,
                      header = F)[, c(2, 4)]
    
    hpo_ <- aggregate(V4 ~ V2, hpo, FUN = paste, collapse = ";")
    colnames(hpo_) <- c("gene_name", "hpo")
    X <- left_join(X, hpo_, by = "gene_name")
    
    
    X <- buscar_herencia(X)
    
    X$SAMPLE <- muestra
    X <- X[!duplicated(X[,c("gene_name","POS","END")]),]
    Y <- X[which(X$errors==""|X$errors=="INFO_REALIGN_3_PRIME"),]

    write.csv(Y, file = file.path(process_dir, "file_ready_analysis.csv"))
    write.csv(X, file=file.path(process_dir,"file_raw_analysis.csv"))
  } else{
    "Ya estamos listo para el analisis, me saliste cubano wey"
  }
  
  
}

compute_depth <- function(fastq_dir, output_dir) {
  ## el directorio del mapeo esta creado pero lo tenemos que guardar de nuevo
  mapping_output_dir <- file.path(output_dir, "mapping_output")
  fastq_files <- list.files(fastq_dir, full.names = F)
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  ## quitamos la extension del archivo
  output_file_name <- file_path_sans_ext(output_file_name)
  
  bam_file <-
    file.path(mapping_output_dir,
              paste0(output_file_name, ".sorted.mark_dup_RG_bqsr.bam"))
  
  dir_coverage <- file.path(output_dir,"coverage_and_stats")
  if(!dir.exists(dir_coverage)){
    dir.create(dir_coverage)
  }
  outfile_coverage <- file.path(dir_coverage,"coverage.txt")
  comando <- paste("./compute_depth.sh", bam_file, ">",outfile_coverage)
  if(!file.exists(outfile_coverage)){
    print("computando cobertura...")
    print(comando)
    system(comando,intern = T)
    
  }else{
    print("la cobertura ya esta computada")
  }
  
}

compute_stats <- function(fastq_dir, output_dir,muestra){
  
dir_coverage <- file.path(output_dir,"coverage_and_stats")
if(!file.exists(file.path(dir_coverage,"stats.csv"))){
  
  dir_coverage <- file.path(output_dir,"coverage_and_stats")
  cov_file <- file.path(dir_coverage,"coverage.txt")
  process_dir <- file.path(output_dir, "post_process_results")
  exoma_file <- file.path(process_dir, "file_ready_analysis.csv")
  
  cov_data <- read.delim(cov_file,header = F)
  exoma <- read.csv(exoma_file,na.strings = ".")
  interest.data <- cov_data[cov_data$V1 !=
                              "all", "V7"]
  mean_coverage <- round(mean(interest.data,na.rm=T), 2)
  mean_coverage20 <- round(mean(interest.data[interest.data>20],na.rm=T), 2)
  cov.data <- cov_data
  
  cov.data$V4 <- gsub("_.*","",cov.data$V4)
  cov.data <- cov.data[cov.data$V1 !="all",]
  
  cov.data <- aggregate(V7 ~ V4, data=cov.data, mean)
  
  genes_totales <- length(unique(cov.data$V4))
  
  genes_20 <- length(cov.data[which(cov.data$V7>20),"V4"])
  
  genes_por_20 <- round(100*genes_20/genes_totales,2)
  
  
  genes_mayor_media <- length(cov.data[which(cov.data$V7>mean_coverage),"V4"])
  
  genes_por_media <- round(100*genes_mayor_media/genes_totales,2)
  
  
  variantes_totales.df <- as.data.frame(lapply(exoma[,c("POS","DP")],as.numeric))
  
  variantes_dp <- aggregate(DP ~ POS, data=variantes_totales.df, mean)
  
  variantes_totales <- length(unique(variantes_dp$POS))
  variantes_20 <- length(unique(variantes_totales.df[which(variantes_totales.df$DP>20),"POS"]))
  variantes_20_x <- round(100*variantes_20/variantes_totales,2)
  
  variantes_media <- length(unique(variantes_totales.df[which(variantes_totales.df$DP>mean_coverage),"POS"]))
  variantes_media_x <- round(100*variantes_media/variantes_totales,2)
  
  res <- as.data.frame(c(mean_coverage=mean_coverage,
                         mean_coverage20=mean_coverage20,
                         genes_totales=genes_totales,
                         genes_20=genes_20,
                         genes_por_20=genes_por_20,
                         genes_mayor_media=genes_mayor_media,
                         genes_por_media=genes_por_media,
                         variantes_totales=variantes_totales,
                         variantes_20=variantes_20,
                         variantes_20_x=variantes_20_x,
                         variantes_media=variantes_media,
                         variantes_media_x))
  res$Descripcion <- c("Cobertura media",
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
                       "Variantes % > media X")
  res <- res[,c(2,1)]
  colnames(res) <- c("Descripcion","Stats")
  
  
  gcov <- cov_data[cov_data[, 1] == 'all',]
  ###
  longitud <- 300
  datos.pre <-
    data.frame(
      X = gcov[1:longitud, 2],
      Y = 100 * (1 - cumsum(gcov[1:longitud, 5])),
      Z = 100 * gcov[1:longitud, 5],
      relleno = gcov[1:longitud, 1]
    )
  datos.pre$relleno <- as.factor(datos.pre$relleno)
  p1 <-
    ggplot(data = datos.pre, aes(X, Y, fill = relleno)) + geom_line(color =
                                                                      "steelblue", linewidth = 2) + xlab("Profundidad de Cobertura") + ylab("Porcentaje de la region >= Profunidad") + theme(legend.position = "none") +
    theme_classic()
  p2 <-
    ggplot(data = datos.pre, aes(X, Z, fill = relleno)) + geom_col() + scale_fill_discrete(type =
                                                                                             "steelblue") + xlab("Profundidad de Cobertura") + ylab("Porcentaje de la region") + theme(legend.position = "none") + theme_classic()
  
  p3 <-ggpubr::ggarrange(p1, p2, ncol = 2)
  
  figura_file_name <- file.path(dir_coverage,"cobertura.jpeg")
  
  p4 <-
    ggpubr::annotate_figure(p3,
                    top = paste(
                      muestra,
                      "Profundidad media de cobertura:",
                      mean_coverage,
                      "X"
                    ))
  ggsave(filename = figura_file_name,plot=p4)
  res_file <- file.path(dir_coverage,"stats.csv")
  write.csv(res,res_file)
  
  blah <- paste("<p>Esta muestra se ha estudiado por el metodo de secuenciacion masiva en paralelo del exoma completo. Se analizaron",res[3,2],"genes. La sensibilidad y la especificidad del metodo son superiores al 98% (SNV< 20 bp INDELS). El porcentaje de genes con una cobertura mayor a 20X es de",round(as.numeric(res[5,2]),2),"%. De todas las variantes identificadas, que son un total de",res[8,2],",",res[9,2],"tienen una cobertura mayor a 20X, esto significa un",round(as.numeric(res[10,2],2)),"%.</p>")

  res <- print(xtable(res),type="html")
  
  fileconn <- "./aux1.html"
  writeLines(blah,fileconn)
  fileconn <- "./aux2.html"
  writeLines(res,fileconn)
  
  command <- paste("cat aux1.html aux2.html > ",file.path(dir_coverage,"doc.html"))

  system(command)

  command <- paste("pandoc --output",file.path(dir_coverage,"reporte.docx"),file.path(dir_coverage,"doc.html"))

  system(command)


  

  
}else{
  print("ya se computaron las estadisticas")
}
}

unique_variants <- function(output_dir,muestra,bd_data){
if(!file.exists(file.path(output_dir,"post_process_results","file_ready_analysis_optimized.csv"))){

  
  process_dir <- file.path(output_dir, "post_process_results")
  exoma <- file.path(process_dir, "file_raw_analysis.csv")
  bd <- readRDS(bd_data)
  exoma <- read.csv(exoma,na.strings = ".")
  exoma_pos <- exoma[,c("CHROM","POS","END","gene_name","SAMPLE")]
  colnames(exoma_pos) <- c("Chr","Start","End","Gene.refGene","codigo")
  cromosomas <- c(paste0("chr",1:22),"chrX","chrY")
  
  if (!any(names(bd) == muestra)) {
    bd[[muestra]] <- exoma_pos
    saveRDS(bd, bd_data)
  }
  
  ## convertimos a df
  
  bd <-   lapply(bd, function(x) as.data.frame(lapply(x,as.character)))
  bd <- lapply(bd, function(x) x[x$Chr %in% cromosomas,])
  bd <- bd %>% bind_rows %>% group_by(Chr,Start,End,Gene.refGene) %>% summarise(n=n(),
                                                                                         paste_m=toString(codigo))
  bd$paste_m <- unlist(lapply(bd$paste_m,function(x) paste0(unique(unlist(strsplit(x,", "))),collapse = ",")))
  
  
  ## buscamos que este en el exoma

  posiciones_start <- bd[grep(muestra,bd$paste_m),c("Start")]
  posiciones_start <- as.numeric(unlist(posiciones_start))
  bd_exoma <- bd[which(as.numeric(bd$Start) %in% posiciones_start),]
  colnames(bd_exoma) <- c("CHROM","POS","END","gene_name","N","samples")
  bd_exoma$POS <- as.numeric(bd_exoma$POS)
bd_exoma$END <-	as.numeric(bd_exoma$END)
exoma$POS <- as.numeric(exoma$POS)
exoma$END <- as.numeric(exoma$END)
exoma_optimized <- dplyr::right_join(bd_exoma, exoma, by = c("CHROM", "gene_name", "POS","END"))
exoma_optimized <- exoma_optimized[!duplicated(exoma_optimized[, c("gene_name", "POS", "CHROM","END")]), ]
for (n in 1:nrow(exoma_optimized)) {
  exoma_optimized$N[n] <- length(unlist(strsplit(exoma_optimized$samples[n], ",")))

}
   
  write.csv(exoma_optimized,file.path(output_dir,"post_process_results","file_raw_analysis_optimized.csv"))
  

  exoma_optimized<- exoma_optimized[!duplicated(exoma_optimized[,c("gene_name","POS","END")]),]
  exoma_optimized_ready <- exoma_optimized[which(exoma_optimized$errors==""|exoma_optimized$errors=="INFO_REALIGN_3_PRIME"),]
  write.csv(exoma_optimized_ready,file.path(output_dir,"post_process_results","file_ready_analysis_optimized.csv"))

  unique_variants <- exoma_optimized_ready[exoma_optimized_ready$N==1,]
  acmg <- read.delim("./acmg.txt",header=T)[,1]
  secondary <- exoma_optimized_ready[exoma_optimized_ready$gene_name %in% acmg,]
  secondary_unique <- secondary[secondary$N==1,]
  write.csv(unique_variants,file.path(output_dir,"post_process_results","file_unique_variants_analysis_optimized.csv"))
  write.csv(secondary,file.path(output_dir,"post_process_results","file_secondary_findings_analysis_optimized.csv"))
  write.csv(secondary_unique,file.path(output_dir,"post_process_results","file_unique_variants_secondary_analysis_optimized.csv"))
 ## low coverage
  acmg_low <- secondary[secondary$DP<20,]
  acmg_low_unique <- secondary_unique[secondary_unique$DP<20,]
  unique_low <- unique_variants[unique_variants$DP<20,]
  raw_low <- exoma_optimized_ready[exoma_optimized_ready$DP<20,]
 ## high coverage >=20
  acmg_high <- secondary[secondary$DP<20,]
  acmg_high_unique <- secondary_unique[secondary_unique$DP>=20,]
  unique_high <- unique_variants[unique_variants$DP>=20,]
  raw_high <- exoma_optimized_ready[exoma_optimized_ready$DP>=20,]
 ## DIRECTORIOS
  dir_low <- file.path(output_dir,"post_process_results","low_coverage")
  if(!dir.exists(dir_low)){
	dir.create(dir_low)
	}
  dir_high <- file.path(output_dir,"post_process_results","high_coverage")
 if(!dir.exists(dir_high)){
        dir.create(dir_high)
	}

  ### write low coverage
  write.csv(acmg_low,file.path(dir_low,"acmg_low.csv"))
  write.csv(acmg_low_unique,file.path(dir_low,"acmg_low_unique.csv"))
  write.csv(unique_low,file.path(dir_low,"unique_low.csv"))
  write.csv(raw_low,file.path(dir_low,"raw_low.csv"))
  ### write high coverage
  write.csv(acmg_high,file.path(dir_high,"acmg_high.csv"))
  write.csv(acmg_high_unique,file.path(dir_high,"acmg_high_unique.csv"))
  write.csv(unique_high,file.path(dir_high,"unique_high.csv"))
  write.csv(raw_high,file.path(dir_high,"raw_high.csv"))


}else{
print("ya estan las unicas variantes")
}
  ## buscamos variantes unicas
  
  
#  posiciones_start <- bd[grep(paste0(paste0("^",muestra,"$")),bd$paste_m),c("Start")]
#  posiciones_start <- as.numeric(unlist(posiciones_start))
  
#  exoma_filtrado <- exoma[exoma$POS %in% posiciones_start,]
#
 # exoma_filtrado <- exoma_filtrado[!duplicated(exoma_filtrado[,c("gene_name","POS","END")]),]






# exoma_ready <- exoma_filtrado[which(exoma_filtrado$errors==""|exoma_filtrado$errors=="INFO_REALIGN_3_PRIME"),]
  
#  write.csv(exoma_filtrado,file = file.path(process_dir,paste0("variantes_unicas",muestra,".csv")))
#  write.csv(exoma_ready,file = file.path(process_dir,paste0("variantes_unicas_ready",muestra,".csv")))
#
#  X.raw <- read.csv(file.path(process_dir,paste0("file_raw_analysis.csv")))
#  X <- read.csv(file.path(process_dir,paste0("file_ready_analysis.csv")))
#  
#  X.raw <- right_join(exoma_filtrado,X.raw,by=c("CHROM","POS","gene_name"))
#  X <- right_join(exoma_ready,X,by=c("CHROM","POS","gene_name"))
#	
 # X.raw <- X.raw[!duplicated(X.raw[,c("CHROM","POS","gene_name")]),]
#  X <- X[!duplicated(X[,c("CHROM","POS","gene_name")]),]

#  write.csv(X,file=file.path(process_dir,paste0("file_ready_analysis_optimized.csv")))
#  write.csv(X.raw,file=file.path(process_dir,paste0("file_ray_analysis_optimized.csv")))

}





start <- Sys.time()
args <- commandArgs(trailingOnly = T)
args <- unlist(strsplit(args, " "))
muestras <- args
for (muestra in muestras) {
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
  
  process(output_dir, fastq_dir, hpo_file, muestra)
  
  compute_depth(fastq_dir = fastq_dir,output_dir)
  
  unique_variants(output_dir,muestra,bd_data)

  compute_stats(fastq_dir, output_dir,muestra)
  
  
}
print(Sys.time() -start)
