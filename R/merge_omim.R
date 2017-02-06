#' Process omim morbidmap.txt file (source data.omim.org/downloads, registration/license required)
#' @export
process_omim <- function(settings=gm_settings) {

  message("Process OMIM file...")

  omim_file <- file.path(settings$cache_dir,settings$omim_morbidmap_file)
  stopifnot(file.exists(omim_file))
  omim <- readLines(omim_file)
  omim <- omim[!startsWith(omim,"#")]
  pheno <- sapply(strsplit(omim,split="\t"),function(x){x[[1]]})
  gene_symbols <- sapply(strsplit(omim,split="\t"),function(x){x[[2]]})
  omim <- data.table(pheno=pheno,gene_symbols=gene_symbols)

  omim <- omim[endsWith(pheno,"(3)") & !startsWith(pheno,"?") ]
  omim[,gene:=sapply(strsplit(omim$gene_symbols,split=","),function(x)x[[1]])]
  omim[,pheno:=sapply(omim$pheno,function(x)substr(x[[1]],0,nchar(x[[1]])-4))]
  omim[,pheno_id:=sapply(strsplit(omim$pheno,split=", "),function(x)ifelse(length(x)>=2,x[[length(x)]],NA))]

  omim <- omim[,.SD,by=.(gene,pheno_id)]

  #Take the first row of each gene-pheno id combi
  setkey(omim,gene,pheno_id)
  omim<- omim[J(unique(gene,pheno_id)),mult="first"]
  omim <- unique(omim[,.(pheno,gene)])

  setkey(omim,gene)

  omim <- omim[,paste0(pheno,collapse=settings$value_sep),by=gene]

  setnames(omim,c("gene","omim_phenos"))
  stopifnot(length(omim$gene)==length(unique(omim$gene)))

  return(omim)
}

#' Merge omim phenotypes to gene matrix
#' @export
merge_omim <- function(gene_matrix, gene_translation_table, settings=gm_settings) {

  annot_label <- "omim"
  map_suffix <- settings$annot_match_suffix

  setkey(gene_matrix, symbol)

  annot_df <- process_omim(settings)

  ##map gene column to official symbol
  symbols <- lookup_symbol(annot_df$gene, gene_translation_table)
  stopifnot(length(symbols) == nrow(annot_df))

  #Number of official gene symbols a gene name could be mapped
  #if gene name is already official symbol choose that one out of several alternative mappings)
  n_genes <- sapply(symbols,length)

  mapped_gene <-as.character(sapply(names(symbols),
                                    function(al)ifelse(al %in% symbols[al],al,
                                                       ifelse(length(symbols[[al]])>=1,symbols[[al]][1],NA))))
  annot_df[,gene:=mapped_gene]
  annot_df[,`:=`(paste0(annot_label,map_suffix),(n_genes==1))]

  annot_df <- annot_df[!is.na(gene),]

  annot_df <- annot_df[,.(paste0(omim_phenos,collapse=settings$value_sep),all(omim_single_gene_match)),by=gene,]
  setnames(annot_df,c("gene",paste0(annot_label,"_phenos"),paste0(annot_label,map_suffix)))


  merge_annotation_by_genename(annot_label,annot_df,gene_matrix,gene_translation_table,settings)

}
