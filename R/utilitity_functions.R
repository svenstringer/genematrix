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
#' @export
install_magma <- function(settings=gm_settings){

  stopifnot( c("magma_url","magma_dir","magma_executable") %in%  names(settings))

  magma_url <- settings$magma_url
  magmaref_url <- settings$magma_ref_url
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
    stopifnot(!file.exists(settings$magma_executable))
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
}

#' Install stand-alone magma if not available yet
annotate_magma <- function(settings=gm_settings){

  stopifnot( c("magma_ref_prefix","magma_dir","magma_executable") %in%  names(settings))

  magma_ref_prefix <- settings$magma_ref_prefix
  magma_executable <- settings$magma_executable

  if(!file.exists(magma_executable)){install_magma(settings)}

  # Create snploc to rsid map
  snpmap <- fread(paste0(magma_ref_prefix, ".bim"))
  setnames(snpmap, c("CHR", "SNP", "CM", "POS", "A1", "A2"))
  snpmap[, `:=`(snpid, paste0(CHR, ":", POS, ":", ifelse(A1 < A2, A1, A2), ":", ifelse(A1 >= A2, A1, A2)))]

  # Create a gene loc
  magma_geneloc <- core[, c("entrez_id", "chr", "start", "end", "strand"), with = F]
  write.table(magma_geneloc, file = magma_geneloc_file, sep = " ", quote = F, row.names = F, col.names = F)

  # Code is format summary data (summaryfile SPECIFIC!!)
  df <- fread(summary_file, select = c("SNP", "CHR", "BP", "P", "A1", "A2"))
  df[, `:=`(snpid, paste0(CHR, ":", BP, ":", ifelse(A1 < A2, A1, A2), ":", ifelse(A1 >= A2, A1, A2)))]

  df2 <- merge(df, snpmap, by.x = c("snpid", "CHR", "BP"), by.y = c("snpid", "CHR", "POS"), suffixes = c("", ".snpmap"))
  setnames(df2, "SNP", "SNP2")
  setnames(df2, "SNP.snpmap", "SNP")

  write.table(df2[, .(SNP, CHR, BP, P, A1, A2)], file = magma_summary_file, sep = " ", quote = F, row.names = F, col.names = T)

  #Annotation magma
  cmd <- paste0(magma_executable, " --annotate --snp-loc ", magma_summary_file, " --gene-loc ", magma_geneloc_file,
                " --out ", magma_annot_prefix)
  system(cmd)
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
