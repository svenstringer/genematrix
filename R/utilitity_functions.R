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
    download_ok <- !download.file(source_url, dest_path)
    stopifnot(download_ok)
    message("Download OK")
  } else {
    message("File ", dest_path, " already exist; no need to download")
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
