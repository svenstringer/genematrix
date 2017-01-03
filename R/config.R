#' @import data.table
#' @import R.utils
#' @importFrom stats end start
#' @importFrom utils download.file write.table
NULL

gencode_version <- 24
cache_dir <- "./datacache"

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
                  sumstats_dir = "./sumstats",
                  cache_dir = cache_dir,
                  gtt_path = file.path(cache_dir,"gene_symbol_table.Rdata"),
                  core_path = file.path(cache_dir,"core_matrix.Rdata"),
                  value_sep = "//",
                  add_excel_collisions = T,
                  annotations = c("exacpli", "omim", "gwascatalog", "genebasedp","expression"),
                  gmversion = 1,
                  gene_bp_dif = 5e+05)
