#' Compute gene-based pvalue based on daner summary stats files available
#'
#' @param settings a list containing all magma-related settings
#'
#' @export
compute_genebased_pvalues <- function(settings=gm_settings){

create_dir(settings$sumstats_dir)
sum_files <- Sys.glob(file.path(settings$sumstats_dir,"daner_*"))

#Download example daner file if sumstat dir is empty
if(length(sum_files)==0){
  stopifnot(!is.null(gm_settings$example_sumstat_path))
  if(get_magma_osversion()=="win_static"){
    message("Cannot process example daner file on windows. Prepare your own summary files and put in sumstats directory to be processed")
    return()
  }else{
  message("Downloading and preparing example daner sumstat file...")

  summary_file <- gm_settings$example_sumstat_path

  download_ok <- !download.file(gm_settings$example_sumstat_url, summary_file,mode="wb")
  stopifnot(download_ok)
  message("Download OK")
  gunzip(summary_file)
  summary_file <- substr(summary_file,0,nchar(summary_file)-3)
  stopifnot(file.exists(summary_file))
  system(paste0("awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' ",summary_file," > ",summary_file,".temp"))
  system(paste0("mv ",summary_file,".temp ",summary_file))
  sum_files <- summary_file
}}

for(f in sum_files){

   output_prefix <- file.path(gm_settings$cache_dir,strsplit(basename(f),split="_")[[1]][2])
  if(!file.exists(paste0(output_prefix, ".genes.out"))){
    if(endsWith(f,".gz")){
      stop("Error: file ",f," is .gz file and will not be processed. Make sure you have all summary files prepared unzipped daner files")
    }else{
    compute_genebased_p(f,output_prefix=output_prefix,magma_model=settings$magma_model)
    }
  }
}
}

#' Compute gene-based pvalue based on daner summary stats files available
#'
#' @param settings a list containing all magma-related settings
#'
#' @export
merge_genebased_pvalues <- function(gene_matrix,settings=gm_settings){

  annotate_magma(gene_matrix,settings)

  res <- merge_snpnrs(gene_matrix,settings)

  compute_genebased_pvalues(settings)

  genep_files <- Sys.glob(file.path(settings$cache_dir,"*.genes.out"))

  for(f in genep_files){
    res<- merge_genebased_p(res,f)
  }
  res
}


#' Merge nr of SNPs per gene used in gene-based test
#'
#' @param gene_matrix a data.table with gene annotation
#' @param settings a list containing all magma-related settings
#'
#' @return a data.table with nr of SNPs added
#'
#' @export
merge_snpnrs <- function(gene_matrix,settings=gm_settings){
  message("Merge nr of snps per gene...")
  annot_file <- paste0(settings$magma_annot_prefix,".genes.annot")
  stopifnot(file.exists(annot_file))
  #automatically split on : (dirty hack)
  snp_annot <- fread(annot_file)
  entrez_ids <- as.numeric(sapply(strsplit(snp_annot$V1,split="\t"),function(x){x[[1]]}))
  snpn <- sapply(strsplit(snp_annot$V3,split="\t"),function(x){length(x)-1})

  res <- merge(gene_matrix,data.table(entrez_id=entrez_ids,N_snps_1000g=snpn),all.x=T,by.x="entrez_id",by.y="entrez_id")
  res
}


#' Compute gene-based pvalue based on daner summary stats file
#'
#' Computes gene-based pvalues with magma and saves them to file.
#'
#' @param settings a list containing all
#'
#' @return filename with gene-based pvalues
#'
#' @export
compute_genebased_p <- function(summary_file=NULL,output_prefix=NULL,annot_prefix=NULL,ref_file=NULL,snpmap_file=NULL,magma_executable=NULL,magma_model="snp-wise=mean"){
  stopifnot(magma_model %in% c("multi","multi=snp-wise","snp-wise=top","snp-wise=top,1","snp-wise=mean"))
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
  if(is.null(magma_executable)){
    stopifnot(!is.null(gm_settings$magma_executable))
    magma_executable <- gm_settings$magma_executable
  }

  message("Process summary file ",summary_file," ...")


  stopifnot(file.exists(paste0(ref_file,".bim")))
  stopifnot(file.exists(paste0(ref_file,".bed")))
  stopifnot(file.exists(paste0(ref_file,".fam")))
  stopifnot(file.exists(paste0(annot_prefix,".genes.annot")))
  stopifnot(file.exists(magma_executable))

  message("Compute genebased p-value with MAGMA")
  #N is irrelevant for outcome in this magma analysis, set arbitrary
  N<-10000
  cmd <- paste0(magma_executable, " --bfile ", ref_file, " --pval ", summary_file, " N=",N," --gene-model ",magma_model," --gene-annot ",
                annot_prefix, ".genes.annot --out ", output_prefix)
  system(cmd)
  return(paste0(output_prefix, ".genes.out"))
}

merge_genebased_p <- function(gene_matrix,f){
  message("Merge genebased p file", f)
  genes <- fread(f,colClasses = c(CHR="character"))
  pval_colname <- substr(basename(f),start=1,stop=nchar(basename(f))-10)
  if("P" %in% names(genes)){
    genes <- genes[,.(GENE,CHR,P)]
    genes[,`:=`(paste0("gene_based_P_",pval_colname),P)]
    genes[,P:=NULL]
  }else{
    genes <- genes[,.(GENE,CHR,P_JOINT)]
    genes[,`:=`(paste0("gene_based_P_",pval_colname),P_JOINT)]
    genes[,P_JOINT:=NULL]
  }
  gene_matrix <- merge(gene_matrix,genes,by.x=c("entrez_id","chr"),by.y=c("GENE","CHR"),all.x=T)
}
