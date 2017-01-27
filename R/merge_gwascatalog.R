# === NHGRI/EBI gwas catalog === use new 'alternative version wget http://www.genome.gov/admin/gwascatalog.txt
# wget https://www.ebi.ac.uk/gwas/api/search/downloads/full -O gwascatalog.txt wget
# https://www.ebi.ac.uk/gwas/api/search/downloads/alternative -O gwascatalog.txt

#' Process exac pli constraint file
process_gwascatalog <- function(settings=gm_settings) {
  stopifnot("gwascatalog_url" %in% names(settings))

  gc_url <- settings$gwascatalog_url
  gc_path <- file.path(settings$cache_dir, "gwascatalog_alternative.txt")
  stopifnot(!is.null(gc_path))
  download_source(gc_url, gc_path)

  message("Process gwas catalog file")

  df <- fread(gc_path,colClasses=list(character=c("P-VALUE","CHR_ID","SNP_GENE_IDS","SNP_ID_CURRENT")))

  df[,P := as.double(df$`P-VALUE`)]
  df[,`P-VALUE` := NULL]

  df <- subset(df,CNV=="N")
  df <- subset(df,P >= 0 & P<=5E-8)
  df <- subset(df,startsWith(`STRONGEST SNP-RISK ALLELE`,"rs"))
  df[,SNPID := sapply(strsplit(`STRONGEST SNP-RISK ALLELE`,split="-"),function(x)x[1])]

  col_selection <- c("SNPID","PUBMEDID", "P","MAPPED_TRAIT")
  df <- df[, col_selection, with = F]

  df <-  df[, .SD[which.min(P)], by=.(MAPPED_TRAIT,SNPID)]

  message("Read snp reference file...")
  bim <- fread(paste0(settings$magma_ref_prefix,".bim"))
  setnames(bim,c("CHR","SNPID","CM","POS","A1","A2"))

  message("Merge SNP position to gwas catalog data")
  df <- merge(df,bim[,.(SNPID,CHR,POS)],by="SNPID")

  setnames(df,c("gwascatalog_rsid","gwascatalog_mappedtrait","gwascatalog_pubmedid","gwascatalog_p","chr","pos"))

  return(df)
}

is_snp_in_gene <- function(chr_gwc, pos_gwc, chr_gene, start_gene, end_gene){
  as.character(chr_gene) == as.character(chr_gwc) & pos_gwc >= start_gene & pos_gwc <= end_gene
}


#' Merge exac pli constraint data into gene matrix
merge_gwascatalog <- function(gene_matrix, gene_translation_table, settings=gm_settings) {

  gwc_file <- file.path(settings$cache_dir,"gwascatalog_mapped.Rdata")

  if(!file.exists(gwc_file)){
  setkey(gene_matrix, symbol)

  gwc <- process_gwascatalog(settings)
  setkey(gwc, chr,pos)

  gwc[,symbol:=as.character(NA)]

  snp_to_gene_map <- sapply(1:nrow(gwc),function(i){
    snp=gwc[i]
    gene_matrix[is_snp_in_gene(snp$chr,snp$pos,gene_matrix$chr_plink,gene_matrix$start,gene_matrix$end),symbol]
  })

  stopifnot(length(snp_to_gene_map)==nrow(gwc))

  for(i in 1:length(snp_to_gene_map)){
    snp <- snp_to_gene_map[[i]]
    if(length(snp)>0){
      gwc[i,symbol:= snp[1]]
      if(length(snp)>1){
        for(s in snp[-1]){
          nrow <- gwc[i,]
          nrow[1,symbol:= s]
          rbind(gwc,nrow)
        }
      }}}

  gwc <- gwc[,.(paste0(gwascatalog_rsid,collapse=settings$value_sep),
                paste0(gwascatalog_mappedtrait,collapse=settings$value_sep),
                paste0(gwascatalog_p,collapse=settings$value_sep),
                paste0(gwascatalog_pubmedid,collapse=settings$value_sep)),
                by=c("symbol")]
  setnames(gwc,old=c("V1","V2","V3","V4"),new=c("gwascatalog_rsid","gwascatalog_mappedtrait","gwascatalog_p","gwascatalog_pubmedid"))
  gwc <- subset(gwc,!is.na(symbol))

  save(gwc,file=gwc_file)
  }
  load(gwc_file)

  merged_df <- merge(gene_matrix, gwc[,.(symbol,gwascatalog_rsid,gwascatalog_mappedtrait,gwascatalog_p,gwascatalog_pubmedid)], by= "symbol", all.x = T,all.y=F)
  stopifnot(nrow(merged_df)==nrow(gene_matrix))
  n_matches <- sum(!is.na(merged_df$gwascatalog_rsid),na.rm=T)
  message(n_matches, " genes out of ", nrow(merged_df), " are uniquely matched to mapped genes in gwascatalog data")
  merged_df
}



