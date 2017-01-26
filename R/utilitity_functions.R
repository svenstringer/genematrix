#' Download source url if it does not exist in destination path
#'
#' Checks if data source already exists in destination path.
#' If not it downloads the data source and saves to destionation path.
#'
#' @param source_url url to source file to download
#' @param dest_path path to save source file in
#'
#' @return None
download_source <- function(source_url, dest_path) {

  if (!file.exists(dest_path)) {
    message("Downloading ", source_url, " to ", dest_path)
    download_ok <- !download.file(source_url, dest_path,mode="wb")
    stopifnot(download_ok)
    message("Download OK")
  } else {
    message("File ", dest_path, " already exist; no need to download")
  }
}


#' Install stand-alone magma if not available yet
install_magma <- function(settings=gm_settings){

  stopifnot( c("magma_url","magma_dir","magma_executable") %in%  names(settings))

  magma_url <- settings$magma_url
  magma_ref_url <- settings$magma_ref_url
  magma_dir <- settings$magma_dir
  magma_ref_prefix <- settings$magma_ref_prefix

  magmazip_path <- file.path(settings$magma_dir,basename(settings$magma_url))
  magmarefzip_path <- file.path(settings$magma_dir,basename(settings$magma_ref_url))

  create_dir(magma_dir)
  #Download magma executable
  if (!file.exists(settings$magma_executable)){
    message("Downloading platform-specific magma binary...")
    download_ok <- !download.file(magma_url, magmazip_path,mode="wb")
    stopifnot(download_ok)
    message("Download OK")
    unzip(magmazip_path,exdir=magma_dir)
    stopifnot(file.exists(settings$magma_executable))
  }

  #Download magma 1000 genomes reference
  if (!file.exists(paste0(magma_ref_prefix,".bed"))){
    message("Downloading magma reference file...")
    download_ok <- !download.file(magma_ref_url, magmarefzip_path,mode="wb")
    stopifnot(download_ok)
    message("Download OK")
    unzip(magmarefzip_path,exdir=magma_dir)
    stopifnot(file.exists(paste0(magma_ref_prefix,".bed")))
    stopifnot(file.exists(paste0(magma_ref_prefix,".bim")))
    stopifnot(file.exists(paste0(magma_ref_prefix,".fam")))
  }

  #Make magma executable on linux/max
  if (get_magma_osversion() %in% c("static","mac")){
    system(paste0("chmod +x ",gm_settings$magma_executable))
  }
}

#' Create gene annotation file for MAGMA (and install MAGMA if not available)
#' @param settings list with magma installation settings
#'
#' @export
annotate_magma <- function(gene_matrix,settings=gm_settings){

  stopifnot( c("magma_ref_prefix","magma_dir","magma_executable") %in%  names(settings))

  magma_ref_prefix <- settings$magma_ref_prefix
  magma_executable <- settings$magma_executable
  magma_geneloc_file <- settings$magma_geneloc_file
  magma_snpmap_file <- settings$magma_snpmap_file
  magma_annot_prefix <- settings$magma_annot_prefix

  if(!file.exists(magma_executable)){install_magma(settings)}

  # Create snploc  map
  #if(!file.exists(magma_snpmap_file)){
  #message("Create magma SNP location map... (can take a couple of minutes)")
  #snpmap <- fread(paste0(magma_ref_prefix, ".bim"))
  #setnames(snpmap, c("CHR", "SNP", "CM", "POS", "A1", "A2"))
  #snpmap[, `:=`(SNP, paste0(CHR, ":", POS, ":", ifelse(A1 < A2, A1, A2), ":", ifelse(A1 >= A2, A1, A2)))]
  #snpmap <- snpmap[, .(SNP, CHR, POS)]
  #write.table(snpmap, file = magma_snpmap_file, sep = " ", quote = F, row.names = F, col.names = T)
  #}

  # Create a gene loc
  if(!file.exists(magma_geneloc_file)){
  message("Create a magma gene location file from core matrix...")
  magma_geneloc <- gene_matrix[, c("entrez_id", "chr", "start", "end", "strand"), with = F]
  write.table(magma_geneloc, file = magma_geneloc_file, sep = " ", quote = F, row.names = F, col.names = F)
  }

  #Annotate genes
  if(!file.exists(paste0(magma_annot_prefix,".genes.annot"))){
  message("Map snps to genes")
  cmd <- paste0(magma_executable, " --annotate --snp-loc ", magma_ref_prefix, ".bim --gene-loc ", magma_geneloc_file,
                " --out ", magma_annot_prefix)
  system(cmd)
  }
}

#' Create directory if it does not exist
#'
#' @param dir_name directory name to check
#'
#' @return None
create_dir <- function(dir_name) {
  dir.create(dir_name, showWarnings = FALSE)
  if (!file.exists(dir_name)) {
    stop("ERROR: directory ", dir_name, " could not be created")
  }
}

#' Checks if path exists
#'
#' Checks if path exists. If not, stop execution.
#'
#' @param path_name to check
#'
#' @return None
check_path <- function(path_name) {
  if (!file.exists(path_name))  stop("ERROR: file ", path_name, " does not exist")
}
