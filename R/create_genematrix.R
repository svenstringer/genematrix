#Assumes R3.3.1+

#' Save annotated Gene Matrix
#'
#' Main function to generate a default Gene Matrix and store it in a data.table.
#' The first time this may take a long time (up to an hour). Will be much quicker due to caching later.
#'
#' @param settings a named list with user specified settings. Defaults to global defaults in gm_setting
#'     stored in config.R
#' @export
generate_genematrix <- function(settings=gm_settings){

stopifnot(c("gencode_version","cache_dir") %in% names(settings) )

#If no gene_matrix_path is found in settings list a default filename is save in the current working directory.
if(! "gene_matrix_prefix" %in% names(settings) ){
  gene_matrix_prefix <- paste0("genematrix_gencodev", settings$gencode_version, "_v", settings$gmversion)
  gene_matrix_path <- file.path(getwd(), gene_matrix_prefix)
}

create_dir(settings$cache_dir)
check_path(settings$cache_dir)

message("Gene matrix file will be saved as: ", gene_matrix_path)

# Create core gene matrix by processing and merging gencode, hgnc, and entrez data
core <- get_core_matrix(settings)

# Create mapping from alias to official gene symbols
gene_translation_table <- get_symbol_table(core,  settings)

# Add custom annotation to create final gene matrix
message("Add annotations...")
gene_matrix <- add_annotations(core, gene_translation_table, settings)

#Save specified columns to final output
message("Save gene matrix to file...")
publish_genematrix(gene_matrix, gene_matrix_prefix)

message("Gene matrix saved as ", gene_matrix_prefix)
}


#' Add annotations to core gene matrix
#'
#' Adds several types of gene annotation to the core matrix, such as constraint scores,
#' gene-based pvalues, GWAS catalog associations, and OMIM associations.
#'
#' @param core core gene matrix to add annotation to
#' @param gene_translation_table list mapping aliases to official symbol used in core gene matrix
#'
#' @return returns a data.table with annotated genes
#'
#' @export
add_annotations <- function(core, gene_translation_table,settings) {

  if(is.null(settings$output_cols)){
    gene_matrix <- core
  } else {
    gene_matrix <- core[,settings$output_cols,with=F]
  }


  #restrict columns in core matrix


  # Add pli scores from exac
  gene_matrix <- merge_exacpli("fullexac", gene_matrix, gene_translation_table, settings)
  gene_matrix <- merge_exacpli("nonpsychexac", gene_matrix, gene_translation_table, settings)

  #Add omim
  if(file.exists(file.path(settings$cache_dir,settings$omim_morbidmap_file))){
    gene_matrix <- merge_omim(gene_matrix, gene_translation_table, settings)
  }else{
    message("Omim file ",settings$omim_morbidmap_file," could not be found. This file requires a license and needs to be manually placed in the cachedir to be used for annotation.")
  }
  #Add gwas catalog
  gene_matrix <- merge_gwascatalog(gene_matrix, gene_translation_table, settings)
  #Add brain expression
  if(file.exists(file.path(settings$cache_dir,settings$commonmind_DLPFC_diffexpr_file))){
    gene_matrix <- merge_commonmindDLPFC(gene_matrix, gene_translation_table, settings)
  }else{
    message("Common mind DLPFC diff expression file ",commonmind_DLPFC_diffexpr_file," could not be found. This file requires a license and needs to be manually placed in the cachedir to be used for annotation.")
  }

  return(gene_matrix)
}


#' Save customized gene matrix based on settings
#'
#' Publish gene matrix in .csv and .Rdata format
#'
#' @param gene_matrix a data.table with gene annotation to save
#' @param settings a list including which columns to save and destination path
#'
#' @return None
#'
#' @export
publish_genematrix <- function(gene_matrix, gene_matrix_prefix=NULL,output_cols=gm_settings$output_colnames) {

  if(is.null(gene_matrix_prefix)){
    #If no gene_matrix_path is found in settings list a default filename is save in the current working directory.
    gene_matrix_prefix <- paste0("genematrix_gencodev", settings$gencode_version, "_v", settings$gmversion)
    gene_matrix_path <- file.path(getwd(), gene_matrix_prefix)
  }

  setkey(gene_matrix,chr,start)

  write.table(gene_matrix,file=paste0(gene_matrix_prefix,".tsv"),quote=F,sep='\t',row.names=F)
  write.table(gene_matrix,file=paste0(gene_matrix_prefix,".csv"),quote=T,sep=',',row.names=F)
  write.foreign(df=gene_matrix, datafile=paste0(gene_matrix_prefix,"_sas.csv"), codefile=paste0(gene_matrix_prefix,".sas"), package="SAS")
  save(gene_matrix,file=paste0(gene_matrix_prefix,".Rdata"))

  message("Gene matrix saved")
}
