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
                Windows="win_static",
                Darwin="mac",
                Linux="static"))
}


gencode_version <- 24
magma_version <- "1.06"
cache_dir <- file.path(".","datacache")
magma_dir <- file.path(".","magma")
sumstats_dir <- file.path(".","sumstats")
magma_ref_url <- "http://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip"

magma_osversion <- get_magma_osversion()

output_cols <-  NULL #c('symbol','chr','chr_name','chr_plink','start','end',
                     #'symbol','aliases','strand','merge_trial',
                     #'tag','source','feature','remap_status','hgnc_id','entrez_id')

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
                  omim_morbidmap_file = "morbidmap.txt",
                  commonmind_DLPFC_diffexpr_file = "CMC_MSSM-Penn-Pitt_DLPFC_mRNA_IlluminaHiSeq2500_gene-adjustedSVA-differentialExpression-includeAncestry-DxSCZ-DE.tsv",
                  magma_url = paste0("http://ctg.cncr.nl/software/MAGMA/prog/magma_v",magma_version,"_",magma_osversion,".zip"),
                  magma_ref_url = magma_ref_url,
                  sumstats_dir = sumstats_dir,
                  cache_dir = cache_dir,
                  gtt_path = file.path(cache_dir,"gene_symbol_table.Rdata"),
                  core_path = file.path(cache_dir,"core_matrix.Rdata"),
                  value_sep = "//",
                  add_excel_collisions = T,
                  gmversion = "1.0",
                  gene_bp_dif = 5e+05,
                  snp_gene_bp_dif = 25000,
                  magma_dir = magma_dir,
                  magma_executable = file.path(magma_dir,ifelse(magma_osversion=="win_static","magma.exe","magma")),
                  magma_ref_prefix = file.path(magma_dir,substr(basename(magma_ref_url),0,nchar(basename(magma_ref_url))-4)),
                  magma_geneloc_file = file.path(magma_dir,paste0("magma_gencode",gencode_version,"_geneloc.txt")),
                  magma_annot_prefix = file.path(magma_dir,paste0("magma_gencode",gencode_version)),
                  magma_snpmap_file = file.path(magma_dir,paste0("magma_gencode",gencode_version,"_snpmap.txt")),
                  magma_model = "snp-wise=mean", #can be "snp-wise=mean, snp-wise=top, or multi=snp-wise"
                  example_sumstat_url = "https://www.med.unc.edu/pgc/files/resultfiles/daner_PGC_SCZ49_1000G-frq.sh2_mds10.gz",
                  example_sumstat_path = file.path(sumstats_dir,"daner_PGC-SCZ49_1000G-frq_mds10.txt.gz"),
                  annot_match_suffix="_single_gene_match",
                  output_cols = output_cols
                 )
