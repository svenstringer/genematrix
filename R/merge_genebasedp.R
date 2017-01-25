#' Compute gene-based pvalue based on daner summary stats files available
#'
#' @param settings a list containing all magma-related settings
#'
#' @export
compute_genebased_pvalues <- function(settings=gm_settings){

sum_files <- Sys.glob(file.path(settings$sumstats_dir,"daner_*"))

#Prepare example daner file
if(length(sum_files)==0){
  stopifnot(!is.null(gm_settings$example_sumstat_path))
  summary_file <- gm_settings$example_sumstat_path
  message("Downloading example daner sumstat file...")
  download_ok <- !download.file(gm_settings$example_sumstat_url, summary_file,mode="wb")
  stopifnot(download_ok)
  message("Download OK")
}

for(f in sum_files){
  compute_genebased_p(f)
}
}

#' Compute gene-based pvalue based on daner summary stats files available
#'
#' @param settings a list containing all magma-related settings
#'
#' @export
merge_genebased_pvalues <- function(gene_matrix,settings=gm_settings){

  genep_files <- Sys.glob(file.path(settings$cache_dir,"*.genes.out"))

  for(f in genep_files){
    res<- merge_genebased_p(gene_matrix,f)
  }
  res
}




#' Compute gene-based pvalue based on daner summary stats file
#'
#' @param settings a list containing all
#'
#' @return Data table with one p-value per gene
#'
#' @export
compute_genebasedp <- function(summary_file=NULL,output_prefix=NULL,annot_prefix=NULL,ref_file=NULL,snpmap_file=NULL,magma_executable=NULL){


  stopifnot(startsWith(basename(summary_file),"daner_"))

  if(is.null(output_prefix)) output_prefix <- file.path(gm_settings$cache_dir,strsplit(basename(summary_file),split="_")[[1]][2])
  if(is.null(ref_file)){
    stopifnot(!is.null(gm_settings$magma_ref_prefix))
    ref_file <- gm_settings$magma_ref_prefix
  }
  if(is.null(annot_prefix)){
    stopifnot(!is.null(gm_settings$magma_annot_prefix))
    annot_prefix <- gm_settings$magma_annot_prefix
  }
  if(is.null(snpmap_file)){
    stopifnot(!is.null(gm_settings$magma_snpmap_file))
    snpmap_file <- gm_settings$magma_snpmap_file
  }
  if(is.null(magma_executable)){
    stopifnot(!is.null(gm_settings$magma_executable))
    magma_executable <- gm_settings$magma_executable
  }

  snpmap <- load(snpmap_file)
  magma_sum_file <- paste0(f,".magma")
  if(endsWith(f,".gz")){
    unzip_name <- substr(f,0,nchar(f)-3)
    if(get_magma_osversion() == "win_static"){
      if(!file.exists(unzip_name)) gunzip(f,remove=F)
      df <- fread(unzip_name,select = c("SNP", "CHR", "BP", "P", "A1", "A2"))
    }else{
      df <- fread(paste0("gunzip -c ",f),select = c("SNP", "CHR", "BP", "P", "A1", "A2"))
    }
  }else{
    df <- fread(f,select = c("SNP", "CHR", "BP", "P", "A1", "A2"))
  }
  setnames(df,old=colnames(df),new=c("SNP", "CHR", "POS", "P", "A1", "A2"))

  df[, `:=`(SNP, paste0(CHR, ":", POS, ":", ifelse(A1 < A2, A1, A2), ":", ifelse(A1 >= A2, A1, A2)))]
  df <- merge(snpmap,df,by=c("SNP","CHR","POS"))

  write.table(df[,.(SNP,CHR,POS,P)],file=magma_sum_file, sep = " ", quote = F, row.names = F, col.names = T)

  stopifnot(file.exists(paste0(ref_file,".bim")))
  stopifnot(file.exists(paste0(ref_file,".bed")))
  stopifnot(file.exists(paste0(ref_file,".fam")))
  stopifnot(file.exists(magma_sum_file))
  stopifnot(file.exists(paste0(annot_prefix,".genes.annot")))
  stopifnot(file.exists(magma_executable))

  message("Compute genebased p-value with MAGMA")
  #N is irrelevant for outcome in this magma analysis, set arbitrary
  N<-10000
  cmd <- paste0(magma_executable, " --bfile ", ref_file, " --pval ", magma_sum_file, " N=",N," --gene-model multi --gene-annot ",
                annot_prefix, ".genes.annot --out ", output_prefix)
  system(cmd)

  genes <- fread(paste0(output_prefix, ".genes.out"))

  setorder(genes, P, ZSTAT)
  genes
}

merge_genebasedp <- function(gene_matrix,f){
  message("Merge genebased p file", f)
  gene_matrix
}

# MAGMA analysis example on ADNI
# magma_executable <- "./magma_v1.05b_static/magma"
# summary_file <- "./IGAP_base_stage1_stats.txt"  #'/home/sven/lisa/ctgukbio/datasets/adni/qc/final/IGAP_base_stage1_stats.txt'
# magma_geneloc_file <- "./magma_gencode_geneloc.txt"
# magma_ref_file <- "./magma_v1.05b_static/g1000_eur/g1000_eur"
# magma_summary_file <- "./magma_igap_gwas.txt"
# magma_gene_prefix <- "./magma_igap"
# magma_annot_prefix <- "./magma_gencode"
#
#
# output_prefix <- "gencode"
#
# # Create snploc to rsid map
# snpmap <- fread(paste0(magma_ref_file, ".bim"))
# setnames(snpmap, c("CHR", "SNP", "CM", "POS", "A1", "A2"))
# snpmap[, `:=`(snpid, paste0(CHR, ":", POS, ":", ifelse(A1 < A2, A1, A2), ":", ifelse(A1 >= A2, A1, A2)))]
#
# # Create a gene loc
# magma_geneloc <- core[, c("entrez_id", "chr", "start", "end", "strand"), with = F]
# write.table(magma_geneloc, file = magma_geneloc_file, sep = " ", quote = F, row.names = F, col.names = F)
#
#
# # Code is format summary data (summaryfile SPECIFIC!!)
# df <- fread(summary_file, select = c("SNP", "CHR", "BP", "P", "A1", "A2"))
# df[, `:=`(snpid, paste0(CHR, ":", BP, ":", ifelse(A1 < A2, A1, A2), ":", ifelse(A1 >= A2, A1, A2)))]
#
# df2 <- merge(df, snpmap, by.x = c("snpid", "CHR", "BP"), by.y = c("snpid", "CHR", "POS"), suffixes = c("", ".snpmap"))
# setnames(df2, "SNP", "SNP2")
# setnames(df2, "SNP.snpmap", "SNP")
#
# write.table(df2[, .(SNP, CHR, BP, P, A1, A2)], file = magma_summary_file, sep = " ", quote = F, row.names = F, col.names = T)
#
# # Annotation magma
# cmd <- paste0(magma_executable, " --annotate --snp-loc ", magma_summary_file, " --gene-loc ", magma_geneloc_file,
#               " --out ", magma_annot_prefix)
# system(cmd)
#
# # Gene-based test (N=74046) Based on mean z-value in gene
# cmd <- paste0(magma_executable, " --bfile ", magma_ref_file, " --pval ", magma_summary_file, " N=74046 --gene-annot ",
#               magma_annot_prefix, ".genes.annot --out ", magma_gene_prefix, "_mean")
# system(cmd)
#
# # based on min z-value in gene
