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
# cmd <- paste0(magma_executable, " --bfile ", magma_ref_file, " --pval ", magma_summary_file, " N=74046 --gene-model multi --gene-annot ",
#               magma_annot_prefix, ".genes.annot --out ", magma_gene_prefix)
# system(cmd)
#
# magma_gene_prefix <- "magma_igap"
#
# genes <- fread(paste0(magma_gene_prefix, ".genes.out"))
#
# setorder(genes, P, ZSTAT)
# genes
