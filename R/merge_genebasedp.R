# The program cant start because libgcc_s_seh-1.dll is missing from your computer. Tyr reinstalling the program to fix thsi problem.
merge_genebasedpvalues <- function(gene_matrix,gene_translation_table,settings=gm_settings){

  annotate_magma(gene_matrix,settings=settings)

}


#' Compute gene-based pvalue based on daner summary stats file
#'
#' @param settings a list containing all
#'
#' @return Data table with one p-value per gene
#'
#' @export
compute_genebasedp <- function(summary_file=NULL,output_prefix=NULL,annot_prefix=NULL,ref_file=NULL,magma_executable=NULL){

  #Check if magma executable exists, if not install
  if(is.null(magma_executable)){
    stopifnot(!is.null(gm_settings$magma_executable))
    magma_executable <- gm_settings$magma_executable
    if(!file.exists(magma_executable))
    {

    }
  }


  if(is.null(summary_file)){
    stopifnot(!is.null(gm_settings$example_sumstat_path))
    summary_file <- gm_settings$example_sumstat_path
  }
  stopifnot(startsWith(basename(summary_file),"daner_"))

  if(is.null(output_prefix)) output_prefix <- file.path(gm_settings$cache_dir,strsplit(basename(summary_file),split="_")[[1]][2])
  if(is.null(ref_file)){
    stopifnot(!is.null(gm_settings$magma_ref_prefix))
    ref_file <- gm_settings$magma_ref_prefix
  }
  if(is.null(annot_prefix)){
    stopifnot(!is.null(gm_settings$magma_annot_prefix))
    ref_file <- gm_settings$magma_annot_prefix
  }



  stopifnot(file.exists(ref_file))
  stopifnot(file.exists(summary_file))
  stopifnot(file.exists(annot_prefix))

  message("Compute genebased p-value with MAGMA")
  #N is irrelevant for outcome in this magma analysis, set arbitrary
  N<-10000
  cmd <- paste0(magma_executable, " --bfile ", ref_file, " --pval ", summary_file, " N=",N," --gene-model multi --gene-annot ",
                annot_prefix, ".genes.annot --out ", output_prefix)
  system(cmd)

  genes <- fread(paste0(output_prefix, ".genes.out"))

  setorder(genes, P, ZSTAT)
  genes
}

merge_genebasedp <- function(settings=gm_settings){
  message("Merge genebased p")
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
#DAner file HEADER (snp IS RS ID)
#CHR     SNP     BP      A1      A2      FRQ_A_33640     FRQ_U_43456     INFO    OR      SE      P       ngt     Direction       HetISqt HetChiSq        HetDf   HetPVa




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





# Code is format summary data SPECIFIC!!
#message("Merge summary file info with magma snpmap")
#df <- fread(summary_file, select = c("SNP", "CHR", "BP", "P", "A1", "A2"))
#df[, `:=`(snpid, paste0(CHR, ":", BP, ":", ifelse(A1 < A2, A1, A2), ":", ifelse(A1 >= A2, A1, A2)))]

##df2 <- merge(df, snpmap, by.x = c("snpid", "CHR", "BP"), by.y = c("snpid", "CHR", "POS"), suffixes = c("", ".snpmap"))
#setnames(df2, "SNP", "SNP2")
#setnames(df2, "SNP.snpmap", "SNP")

#write.table(df2[, .(SNP, CHR, BP, P, A1, A2)], file = magma_summary_file, sep = " ", quote = F, row.names = F, col.names = T)
