#' Process omim morbidmap.txt file (source data.omim.org/downloads, registration/license required)
process_commonmindDLPFC <- function(settings=gm_settings) {

  message("Process commonmind DLPF differential expression SCZ vs control data...")

  expr_file <- file.path(settings$cache_dir,settings$commonmind_DLPFC_diffexpr_file)

  cm <- fread(expr_file)


  return(cm)
}

#' Merge omim phenotypes to gene matrix
merge_commonmindDLPFC <- function(gene_matrix, settings=gm_settings) {

  setkey(gene_matrix, symbol)

  annot_label <- "cm_DLPFC_"

  cm <- process_commonmindDLPFC(settings)
  setnames(cm,names(cm),paste0(annot_label,names(cm)))
  merged_df <- merge(gene_matrix,cm,by.x=c("ensembl_id","symbol"),by.y=paste0(annot_label,c("genes","MAPPED_genes")),all.x=T)
  merged_df
}
