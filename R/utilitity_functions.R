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

#' Get compiled magma platform version to download: windows (win), mac (mac osx), static (linux fully static )
get_magma_osversion <- function(){
  os <- Sys.info()[["sysname"]]
  return(switch(os,
         Windows="win",
         Darwin="mac",
         Linux="static"))
}

#' Install stand-alone magma if not available yet
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
