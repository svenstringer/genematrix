#' @import data.table
#' @import R.utils
#' @importFrom stats end start
#' @importFrom utils download.file write.table
NULL

#' Get compiled magma platform version to download: windows (win), mac (mac osx), static (linux fully static )
#' @export
get_magma_osversion <- function(){
  os <- Sys.info()[["sysname"]]
  return(switch(os,
                Windows="win",
                Darwin="mac",
                Linux="static"))
}


gencode_version <- 24
magma_version <- "1.06"
cache_dir <- file.path(getwd(),"datacache")
magma_dir <- file.path(getwd(),"magma")
sumstats_dir <- file.path(getwd(),"sumstats")
magma_ref_url <- "http://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip"

magma_osversion <- get_magma_osversion()

#' Global settings available after loading GeneMatrix package
#' @export
gm_settings <- list(gencode_version = gencode_version,  #v24 used in UCSC browser, although v25 available (August 2016)
                 gencode_url = paste0("ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_", gencode_version,
                      "/GRCh37_mapping/gencode.v", gencode_version,
                      "lift37.annotation.gtf.gz"), # Source: http://www.gencodegenes.org/releases/grch37_mapped_releases.html
                  hgnc_url = "ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt", # HGNC source
                  entrez_url = "ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz", # Entrez gene info
                  fullexacpli_url = "ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/functional_gene_constraint/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt", # Additional annotation files
                  nonpsychexacpli_url = "ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/functional_gene_constraint/fordist_cleaned_nonpsych_z_pli_rec_null_data.txt",
                  gwascatalog_url = "https://www.ebi.ac.uk/gwas/api/search/downloads/alternative",
                  omim_url = "http://omim.org/static/omim/data/mim2gene.txt",
                  magma_url = paste0("http://ctg.cncr.nl/software/MAGMA/prog/magma_v",magma_version,"_",magma_osversion,".zip"),
                  magma_ref_url = magma_ref_url,
                  sumstats_dir = sumstats_dir,
                  cache_dir = cache_dir,
                  gtt_path = file.path(cache_dir,"gene_symbol_table.Rdata"),
                  core_path = file.path(cache_dir,"core_matrix.Rdata"),
                  value_sep = "//",
                  add_excel_collisions = T,
                  gmversion = "1.0.0",
                  gene_bp_dif = 5e+05,
                  snp_gene_bp_dif = 25000,
                  magma_dir = magma_dir,
                  magma_executable = file.path(magma_dir,ifelse(magma_osversion=="win","magma.exe","magma")),
                  magma_ref_prefix = file.path(magma_dir,substr(basename(magma_ref_url),0,nchar(basename(magma_ref_url))-4)),
                  magma_geneloc_file = file.path(magma_dir,"magma_gencode_geneloc.txt"),
                  magma_annot_prefix = file.path(magma_dir,"magma_gencode")
                 )
