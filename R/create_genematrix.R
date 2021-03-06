#Assumes R3.3.1+

#' Save annotated Gene Matrix
#'
#' Main function to generate a default Gene Matrix and store it in a data.table.
#' The first time this may take a long time (up to an hour). Will be much quicker due to caching later.
#'
#' @param settings a named list with user specified settings. Defaults to global defaults in gm_setting
#'     stored in config.R
#' @return data.table with gene matrix
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
gene_matrix
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


  #restrict columns in core matrix

  # Add genebased pvalues
  gene_matrix <- merge_genebased_pvalues(gene_matrix,settings)


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
  gene_matrix <- merge_gwascatalog(gene_matrix, settings)

  #Add brain expression
  if(file.exists(file.path(settings$cache_dir,settings$commonmind_DLPFC_diffexpr_file))){
    gene_matrix <- merge_commonmindDLPFC(gene_matrix, settings)
  }else{
    message("Common mind DLPFC diff expression file ",commonmind_DLPFC_diffexpr_file," could not be found. This file requires a license and needs to be manually placed in the cachedir to be used for annotation.")
  }


  #Final renames
  setnames(gene_matrix,
           old=c("symbol","chr_prefix","chr_name","hgnc_full_name",
                 "gene_id","hgnc_ucsc_id","hgnc_refseq_accession","hgnc_uniprot_ids","hgnc_ccds_id"),
           new=c("gene_name","chr_full","chr_char","hgnc_gene_product",
                 "ensembl_id_long","ucsc_id","refseq_accession","uniprot_ids","ccds_id"))

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
publish_genematrix <- function(gene_matrix, gene_matrix_prefix=NULL,settings=gm_settings,utility_version=T) {

  output_col_file=settings$output_col_file

  if(is.null(gene_matrix_prefix)){
    suffix <- ""
    if(!utility_version) suffix <-  "_extended"
    #If no gene_matrix_path is found in settings list a default filename is save in the current working directory.
    gene_matrix_prefix <- paste0("genematrix_gencodev", settings$gencode_version, "_v", settings$gmversion,suffix)
    gene_matrix_path <- file.path(getwd(), gene_matrix_prefix)
  }

  #Set column order and rename column names
  stopifnot(file.exists(output_col_file))
  output_cols <- fread(output_col_file)

  genep_lab_prefix <- "gene_based_P_"
  genep_row_idx <- which(startsWith(output_cols$label,genep_lab_prefix))

  genep_columns <- names(gene_matrix)[startsWith(names(gene_matrix),genep_lab_prefix)]
  out_df <- output_cols[1:(genep_row_idx-1),]
  if(length(genep_columns)>0){
    out_df <- rbind(out_df,data.table(label=genep_columns,utility_version=output_cols[genep_row_idx,utility_version],description="gene-based pvalue",source="magma"))
  }
  out_df <- rbind(out_df,output_cols[(genep_row_idx+1):nrow(output_cols),])

  setcolorder(gene_matrix, out_df$label)
  gm <- gene_matrix
  if(utility_version){
    out_df <- subset(out_df, utility_version=="Y")
    gm <- gene_matrix[,(out_df$label),with=F]
  }

  setkey(gm,chr_plink,start)

  # Tab separated for scripting
  write.table(gm,file=paste0(gene_matrix_prefix,".tsv"),quote=F,sep='\t',row.names=F)
  #SAS conversion
  write.foreign(df=gm, datafile=paste0(gene_matrix_prefix,"_sas.csv"), codefile=paste0(gene_matrix_prefix,".sas"), package="SAS")
  # R binary format
  save(gm,file=paste0(gene_matrix_prefix,".Rdata"))

  gene_matrix_excel <- gm

  for(col in names(gene_matrix_excel)){
    if(endsWith(col,"_link")){
      excel_link <- paste0('=HYPERLINK("',gm[[col]],'","',col,'")')
      gene_matrix_excel[,(col):=excel_link]
    }
  }

  #.csv for easy import in Excel (add hyperlink)
  write.table(gene_matrix_excel,file=paste0(gene_matrix_prefix,".csv"),quote=T,sep=',',row.names=F,qmethod="double")

  labels <- c("publication_date", names(settings))
  values <- c(date(),unlist(settings))
  log_df <- data.table(label=labels,value=values)
  write.table(log_df,file=paste0(gene_matrix_prefix,".settings"),sep="\t",row.names = F,col.names=T)

  message("Gene matrix saved")
}
