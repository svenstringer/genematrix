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

stopifnot(c("gencode_version","cache_dir") %in% settings )
attach(settings)
on.exit(detach(settings), add = T)

#If no gene_matrix_path is found in settings list a default filename is save in the current working directory.
if(! gene_matrix_path %in% names(settings) ){
  gene_matrix_file <- paste0("genematrix_core_gencode", gencode_version, "_", Sys.Date(), ".csv")
  gene_matrix_path <- file.path(getwd(), gene_matrix_file)
}

create_dir(cache_dir)
check_path(cache_dir)

message("Gene matrix file will be saved as: ", gene_matrix_path)

# Create core gene matrix by processing and merging gencode, hgnc, and entrez data
core <- get_core_matrix(settings)

# Create mapping from alias to official gene symbols
gene_translation_table <- get_symbol_table(core,  settings)


# Add custom annotation to create final gene matrix
gene_matrix <- add_annotations(core, gene_translation_table, settings)

#Save specified columns to final output
publish_genematrix(gene_matrix, settings)

message("Gene matrix saved as ", gene_matrix_path)

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
  gene_matrix <- core

  # Add pli scores from exac
  gene_matrix <- merge_exacpli(settings$fullexacpli_url, gene_matrix, gene_translation_table, settings$gene_bp_dif)
  gene_matrix <- merge_exacpli(settings$nonpsychexacpli_url, gene_matrix, gene_translation_table, settings$gene_bp_dif)

  return(gene_matrix)
}

#' Save customized gene matrix based on settings
#'
#' Publish gene matrix in .csv format including user-set columns
#'
#' @param gene_matrix a data.table with gene annotation to save
#' @param settings a list including which columns to save and destination path
#'
#' @return None
#'
#' @export
publish_genematrix <- function(gene_matrix, settings) {
  output_cols <- c('chr','chr_name','chr_plink','start','end','symbol','aliases','strand','merge_trial','tag','source','feature','remap_status','hgnc_id')
  output_df <- gene_matrix[,output_cols,with=F]
  write.table(output_df,file=gene_matrix_path,quote=T,sep='\t',row.names=F)

  message("Publish gene matrix")
}






