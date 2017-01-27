# wget http://omim.org/static/omim/data/mim2gene.txt -O omim.mim2gene.txt wget
# http://data.omim.org/downloads/MbiFy5hvQ92TJUMFEpE_7w/mimTitles.txt -O omim.mimTitles.txt wget
# http://data.omim.org/downloads/MbiFy5hvQ92TJUMFEpE_7w/genemap.txt -O omim.genemap.txt wget
# http://data.omim.org/downloads/MbiFy5hvQ92TJUMFEpE_7w/morbidmap.txt -O omim.morbidmap.txt wget
# http://data.omim.org/downloads/MbiFy5hvQ92TJUMFEpE_7w/genemap2.txt -O omim.genemap2.txt

# === NHGRI/EBI gwas catalog === use new 'alternative version wget http://www.genome.gov/admin/gwascatalog.txt
# wget https://www.ebi.ac.uk/gwas/api/search/downloads/full -O gwascatalog.txt wget
# https://www.ebi.ac.uk/gwas/api/search/downloads/alternative -O gwascatalog.txt


#https://data.omim.org/downloads/MbiFy5hvQ92TJUMFEpE_7w/morbidmap.txt

#' Process omim morbidmap.txt file (source data.omim.org/downloads, registration/license required)
process_omim <- function(settings=gm_settings) {

  message("Process OMIM file (warning can be ignored)...")

  omim_file <- file.path(settings$cache_dir,settings$omim_morbidmap_file)
  stopifnot(file.exists(omim_file))
  omim <- fread(omim_file)
  setnames(omim,c("pheno","gene_symbols","MIMnr","location"))

  omim <- omim[endsWith(pheno,"(3)") & !startsWith(pheno,"?") ]
  omim[,gene:=sapply(strsplit(omim$gene_symbols,split=","),function(x)x[[1]])]
  omim <- unique(omim[,.(pheno,gene)])


  omim <- omim[,paste(pheno,collapse=settings$value_sep),by=gene]
  setnames(omim,c("gene","omim_phenos"))

  return(omim)
}

#' Merge omim phenotypes to gene matrix
merge_omim <- function(gene_matrix, gene_translation_table, settings=gm_settings) {

  setkey(gene_matrix, symbol)

  annot_label <- "omim"

  merge_annotation_type(annot_label,process_omim,gene_matrix,gene_translation_table,settings)

}
