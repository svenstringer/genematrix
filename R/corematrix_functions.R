#' Preprocess gencode file
#'
#' Filter gencode file to include only genes on chromosome 1-22,X,Y,M and reformat before returning
#' as \code{data.table}
#'
#' Gencode files are downloaded from
#' \link{http://www.gencodegenes.org/releases/grch37_mapped_releases.html}
#'
#' Assumed input format .gtf.gz:
#' \itemize{
#'  \item{1) }{chromosome name chr\{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,M\}
#'            or GRC accession}
#'  \item{2) }{annotation source \{ENSEMBL,HAVANA\}}
#'  \item{3) }{feature type \{gene,transcript,exon,CDS,UTR,start_codon,stop_codon,Selenocysteine\}}
#'  \item{4) }{genomic start location integer-value (1-based)}
#'  \item{5) }{genomic end location integer-value}
#'  \item{6) }{score (not used)}
#'  \item{7) }{genomic strand \{+,-\}}
#'  \item{8) }{genomic phase (for CDS features) \{0,1,2,.\}}
#'  \item{9) }{additional information as key-value pairs}
#'}
#' First, the entries are filtered on \code{feature_type == 'gene'} and \code{status == 'KNOWN'}.
#' This mostly excludes transcripts. The 'chr' prefix is removed from chromosome values and any
#' chromosomes other than 1-2,X,Y,M are removed. The resulting chromosome values are cast into
#' an ordered factor (ordering: 1-22,X,Y,M). Then additional columns are extracted from the
#' key,value pairs in the info column. Any genes with gene_types in \code{c('misc_RNA','snoRNA','snRNA')}
#' are removed. Finally the redundant columns score, phase, and info are removed and a
#' new column ensembl_gene_id is created from gene_id that does not contain subnumbering
#' (i.e. id is x instead of x.y). The resulting file still contains duplicate gene names,
#' but these will be removed after the merge with the canonical hgnc data.

#' @param gencode_path file path of downloaded .gtf.gz file
#' @return processed gencode table as \code{data.table}
#' @export
#' @examples
#' \dontrun{
#' gencode_data <- process_gencodefile(gencode_path)
#' }
process_gencodefile <- function(gencode_path) {

    message("Filter gencode file on feature type == 'gene' and status KNOWN...")

    #
    #gencode <- fread(paste("gunzip -c", gencode_path, "| awk '$3 == \"gene\" {print $0}' - | grep KNOWN"))
    if(.Platform$OS.type=="windows"){
      unzip_name <- substr(gencode_path,0,nchar(gencode_path)-3)
      if(!file.exists(unzip_name)) gunzip(gencode_path,remove=F)
      gencode <- fread(unzip_name)
    }else{
      gencode <- fread(paste("gunzip -c", gencode_path, "| awk '$3 == \"gene\" {print $0}' - | grep KNOWN"))
    }

    setnames(gencode, c("chr", "source", "feature", "start", "end", "score", "strand", "phase", "info"))

    if(.Platform$OS.type=="windows"){
      gencode <- gencode[feature=="gene",]
      gencode <- gencode[sapply(1:nrow(gencode),function(i)length(grep("KNOWN",gencode[i,"info",with=T],fixed=T))>0),]
    }

    message("Gencode table contains ", nrow(gencode), " gene entries")


    ## Process chromosome column
    message("Add chr_prefix column with 'chr' prefix removed")
    # strip chr from chromosome id
    gencode[, chr_prefix := chr]
    gencode[, chr := gsub("chr", "", chr, fixed = T)]

    # Remove chr ids other than 1-22,X,Y,M (1-25) and order on chr,end
    message("Remove chromosomes other than 1-22,X,Y,M...")

    chr_ids <- c(as.character(1:22), "X", "Y", "M")
    gencode <- subset(gencode, chr %in% chr_ids)
    gencode[, chr_id := factor(chr, levels = chr_ids, ordered = T)]

    # Map genes on Y-PAR to X
    message("Annotate PAR genes")
    gencode[, is_par := as.numeric(is_par(chr, start, end))]

    # Add plink chromosome convention column
    message("Add plink chromosome name convention chr_plink column (1-22,23=X,24,=Y,25=XY/PAR,26=M)")
    gencode[, chr_plink := as.numeric(chr_id)]
    gencode[chr_plink == 25, ]$chr_plink <- 26 #Map M to 26
    gencode[is_par == 1, ]$chr_plink <- 25 #Map PAR region (XY) to 25


    # Add chromosome name column
    message("Add chromosome name convention chr_name column (1-22,X,Y,XY,M)")
    gencode[, chr_name := factor(chr_plink,levels=c(1:22,"X","Y","XY","M"))]
    gencode[chr_plink == 23, ]$chr_name <- "X"
    gencode[chr_plink == 24, ]$chr_name <- "Y"
    gencode[chr_plink == 25, ]$chr_name <- "XY"
    gencode[chr_plink == 26, ]$chr_name <- "M"

    # identify additional implicit info columns
    message("Extract additional columns from info columns...")
    s <- strsplit(gencode$info, split = "; ?", fixed = F)
    info_colnames <- unique(sapply(strsplit(unlist(s), split = " "), function(x) x[1]))

    # Add additional columns with info from info column
    info_cols <- replicate(length(info_colnames), rep(NA, nrow(gencode)), simplify = F)
    names(info_cols) <- info_colnames
    for (colname in info_colnames) {
        message("Adding column ", colname)
        info_cols[[colname]] <- sapply(s, function(x) ifelse(any(startsWith(x, colname)),
                                                      gsub("\"", "", strsplit(x[startsWith(x,
                                                      colname)], split = " ")[[1]][2]), NA))
    }
    gencode[, `:=`((info_colnames), info_cols)]

    #process 'level' variable and rename into verification_level_gencode
    level_labs <- c("Verified","Manually annotated","Automatically annotated")
    gencode[, gencode_verification_level:=level_labs[as.numeric(gencode$level)]]
    gencode[, level:=NULL]

    # Remove genes with multiple mappings
    message("Deleting genes with >1 mapping or remap_status not 'full_contig' or 'automatic_gene'")
    del_tags <- c("fragmented_locus", "reference_genome_error", "semi-processed")
    message("Deleting genes with remap tag: ", paste(del_tags, collapse = ", "))
    gencode <- gencode[(remap_num_mappings == 1 | is.na(remap_num_mappings)) &
                         remap_status %in% c("full_contig",
                          "automatic_gene") & !(tag %in% del_tags), ]

    # Filter on gene_type
    del_gene_types <- c("misc_RNA", "snoRNA", "snRNA")
    message("Deleting genes with gene_type: ", paste(del_gene_types, collapse = ", "))
    gencode <- gencode[!gene_type %in% del_gene_types, ]

    # Deleting uninformative columns
    message("Deleting score column")
    gencode[, `:=`(score, NULL)]
    message("Deleting phase column")
    gencode[, `:=`(phase, NULL)]
    message("Deleting info column")
    gencode[, `:=`(info, NULL)]
    message("Deleting feature column")
    gencode[, `:=`(feature, NULL)]
    message("Deleting gene_status column")
    gencode[, `:=`(gene_status, NULL)]
    message("Deleting remap_substituted_missing_target column")
    gencode[, `:=`(remap_substituted_missing_target, NULL)]
    message("Deleting remap_num_mappings column")
    gencode[, `:=`(remap_num_mappings, NULL)]

    message("Add ensembl_gene_id column")
    gencode[, `:=`(ensembl_id, sapply(strsplit(gene_id, split = ".", fixed = T), function(x) x[1]))]

    message("Add vega_id column")
    gencode[, `:=`(vega_id, sapply(strsplit(havana_gene, split = ".", fixed = T), function(x) x[1]))]

    message("Order on chr,start position")
    setkey(gencode, chr, start)
    message("Gencode table contains ", nrow(gencode), " genes")

    setnames(gencode,old=c("havana_gene","gene_name","source","tag","gene_type"),
             c("havana_id","gencode_gene_name","gencode_source","gencode_tag","gencode_gene_type"))

    return(gencode)
}

#' Preprocess hgnc file as downloaded from
#' http://www.genenames.org/cgi-bin/statistics > Complete hgnc_dataset (.txt)
#'
#' Assumed input format .txt:
#'
#' Processing steps:
#'
#' 1) Restrict to relevant columns (as specified in the function)
#'
#' 2) Remove genes that are withdrawn
#'
#' 3) Add alias column which merges gene symbols from alias_symbol and prev_symbol in a single column
#'
#' 4) Change '|' value separater into ';'
#'
#' @param hgnc_path file path of downloaded .txt file
#' @param value_sep character string to separate multi-value fields
#' @return processed hgnc table as data.table
#' @export
process_hgncfile <- function(hgnc_path, value_sep) {
    columns <- c("hgnc_id", "symbol", "name", "locus_group", "locus_type", "gene_family",
                 "gene_family_id", "alias_symbol", "alias_name", "prev_symbol", "location",
                 "date_symbol_changed", "date_modified", "entrez_id", "ensembl_gene_id",
                 "vega_id", "ucsc_id", "refseq_accession", "ccds_id", "uniprot_ids",
                 "pubmed_id", "omim_id")



    hgnc <- fread(hgnc_path, colClasses = "character")

    hgnc$entrez_id <- as.integer(hgnc$entrez_id)


    message("Keep the following columns: ", paste(columns, collapse = ", "))
    message("Remove hgnc withdrawn genes")
    hgnc <- hgnc[locus_group != "withdrawn", columns, with = F]

    message("Add alias column wih all aliases for a gene")

    s1 <- strsplit(hgnc$alias_symbol, split = "|", fixed = T)
    alias_list1 <- sapply(s1, function(x) ifelse(length(x) == 0, "", paste(x, collapse = value_sep)))
    s2 <- strsplit(hgnc$prev_symbol, split = "|", fixed = T)
    alias_list2 <- sapply(s2, function(x) ifelse(length(x) == 0, "", paste(x, collapse = value_sep)))

    hgnc[, aliases := gsub(paste0("^", value_sep, "|", value_sep, "$"),
                              "", paste(alias_list1, alias_list2, sep = value_sep))]

    message("Change field separator string into custom separator ", value_sep)

    for (col_name in names(hgnc)) {
        if (class(hgnc[, col_name, with = F][[1]]) == "character")
           hgnc[, (col_name):= gsub("|", value_sep, hgnc[, col_name, with = F][[1]], fixed = T),
              with = T]
    }

    message("Add chr column based on location and delete genes with locations not in 1-22,X,Y,M")
    chr_ids <- c(as.character(1:22), "X", "Y", "M")
    chr <- toupper(sapply(strsplit(hgnc$location, split = "[ a-ln-wA-LN-WzZ_]", fixed = F),
                          function(x) x[1]))
    hgnc[, `:=`(chr, factor(chr, levels = chr_ids, ordered = T))]
    hgnc <- subset(hgnc, !is.na(chr))

    #Remove HGNC columns
    hgnc[, alias_symbol := NULL]
    hgnc[, date_symbol_changed := NULL]
    hgnc[, date_modified := NULL]
    #hgnc[, date_symbol_changed := NULL]
    #hgnc[, date_modified := NULL]
    hgnc[, prev_symbol := NULL]
    hgnc[, omim_id := NULL]


    message("HGNC table contains ", nrow(hgnc), " gene entries")

    return(hgnc)
}

#' Preprocess entrez file as downloaded from
#' ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz > gene info file from ncbi
#'
#' Assumed input format .gz:
#'
#' Processing steps:
#'
#' 1) Change | separator into specified value separator
#'
#'
#' @param entrez gene info file path of downloaded .gz file
#' @param value_sep value separator used for multi-value columns
#' @return processed gene info table as data.table
#' @export
#' @examples
#' \dontrun{
#' entrez_file <- process_entrez(entrez_path)
#' }
process_entrezfile <- function(entrez_path, value_sep) {
    message("Process entrez file...")
    if(.Platform$OS.type=="windows"){
      unzip_name <- substr(entrez_path,0,nchar(entrez_path)-3)
      if(!file.exists(unzip_name)) gunzip(entrez_path,remove=F)
      entrez <- fread(unzip_name)
    }else{
      entrez <- fread(paste("gunzip -c", entrez_path))
    }

    old_columns <- c("GeneID", "Symbol", "chromosome",
                     "description", "type_of_gene",
                     "Other_designations")

    new_names <- c("entrez_id", "entrez_gene_name", "chr",
                    "entrez_description", "entrez_type_of_gene",
                   "entrez_other_designations")

    entrez <- entrez[,old_columns,with=F]

    setnames(entrez,old=old_columns, new=new_names)

    entrez[chr=="X|Y",chr:="XY"]
    entrez[chr=="MT",chr:="M"]


    message("Add chr column based on location and delete genes with locations not in 1-22,X,Y,XY,M")
    chr_ids <- c(as.character(1:22), "X", "Y", "XY","M")
    entrez[, chr_factor := factor(chr, levels = chr_ids, ordered = T)]
    entrez <- subset(entrez, (chr %in% chr_ids))


    message("Change field separator string into custom separator ", value_sep)

    for (col_name in names(entrez)) {
        if (class(entrez[, col_name,with=F][[1]])[1] == "character")
            entrez[, (col_name) :=  gsub("|", value_sep, entrez[, col_name, with = F][[1]],
                                         fixed = T)]
    }


    message("Entrez table contains ", nrow(entrez), " gene entries")

    return(entrez)
}


#' Merge gencode and hgnc file
#' @export
merge_gencode_hgnc <- function(gencode, hgnc) {

    gencode_hgnc_suffixes <- c("_gencode", "_hgnc")

    # First merge on gene_name,ensembl_id(30594 unique symbols/ensemble ids)
    df <- merge(gencode, hgnc, by.x = c("ensembl_id", "gencode_gene_name", "chr_id"),
                by.y = c("ensembl_gene_id", "symbol",
        "chr"), suffixes = gencode_hgnc_suffixes)
    df[, merge_trial :=1 ]
    df[, symbol := gencode_gene_name]

    message("First merge on gene symbol, ensembl id, and chromosome: ", nrow(df), " genes")

    # Second merge on gene_name,chr
    gencode_remain <- gencode[!ensembl_id %in% df$ensembl_id | !gencode_gene_name %in% df$symbol, ]
    hgnc_remain <- hgnc[!symbol %in% df$symbol | !ensembl_gene_id %in% df$ensembl_gene_id, ]

    df2 <- merge(gencode_remain, hgnc_remain, by.x = c("gencode_gene_name", "chr_id"), by.y = c("symbol", "chr"),
                 suffixes = gencode_hgnc_suffixes)
    df2[, merge_trial := 2]
    df2[, symbol := gencode_gene_name]
    df2[, ensembl_gene_id := NULL]

    #gencode_temp_col <- paste0("ensembl_gene_id", gencode_hgnc_suffixes[1])
    #hgnc_temp_col <- paste0("ensembl_gene_id", gencode_hgnc_suffixes[2])
    #df[, ensembl_id := paste0("ensembl_id", gencode_hgnc_suffixes[1])]

    df <- rbind(df, df2)
    message("Second merge on gene symbol, and chromosome: ", nrow(df2), " additional genes")

    # Third merge on ensembl_id,chr
    gencode_remain <- gencode[!ensembl_id %in% df$ensembl_id | !gencode_gene_name %in% df$symbol, ]
    hgnc_remain <- hgnc[!symbol %in% df$symbol | !ensembl_gene_id %in% df$ensembl_gene_id, ]

    df2 <- merge(gencode_remain, hgnc_remain, by.x = c("ensembl_id", "chr_id"),
                  by.y = c("ensembl_gene_id", "chr"),
                  suffixes = gencode_hgnc_suffixes)
    df2[, merge_trial := 3]

    df <- rbind(df, df2)

    message("Third merge on ensemble gene id, and chromosome: ", nrow(df2), " additional genes")

    # Check and fix vega_id matching between ENCODE and hgnc
    vega_id_nomatch <- as.numeric(df[, paste0("vega_id",
                                  gencode_hgnc_suffixes[1]), with = F][[1]] != df[, paste0("vega_id",
                                  gencode_hgnc_suffixes[2]), with = F][[1]])
    n_nonmatching_vegaid <- sum(vega_id_nomatch, na.rm = T)
    if (n_nonmatching_vegaid > 0) {
        message("Warning: vega_id is different between GENCODE and HGNC file for ",
                n_nonmatching_vegaid, " genes ")
        message("Warning: vega_id of GENCODE file is leading downstream")
    }


    gencode_temp_col <- paste0("vega_id", gencode_hgnc_suffixes[1])
    hgnc_temp_col <- paste0("vega_id", gencode_hgnc_suffixes[2])

    df[, vega_id := gencode_temp_col ]
    df[, (gencode_temp_col) := NULL]
    df[, (hgnc_temp_col) := NULL]
    df[, no_vegaid_match := vega_id_nomatch]
    df[, location := NULL]


    message("Rename entrez headers...")
    setnames(df,old=c("name","locus_group","locus_type",
                      "gene_family","gene_family_id","alias_name",
                      "ucsc_id", "refseq_accession", "ccds_id",
                      "uniprot_ids","pubmed_id","merge_trial"),
             new=c("hgnc_full_name","hgnc_locus_group","hgnc_locus_type",
                   "hgnc_gene_family","hgnc_gene_family_id","hgnc_alias_names",
                    "hgnc_ucsc_id", "hgnc_refseq_accession", "hgnc_ccds_id",
                    "hgnc_uniprot_ids","hgnc_pubmed_id","hgnc_merge_trial"))

    message("Total number of genes in gene matrix after merge gencode and hgnc: ", nrow(df))


    return(df)
}

#' Merge gencode and hgnc file
#' @export
merge_core_entrez <- function(core, entrez) {
    df <- merge(core, entrez, by.x = c("entrez_id", "symbol", "chr_name"),
                by.y = c("entrez_id", "entrez_gene_name", "chr_factor"),suffixes=c("","_entrez"))
    df[,chr_entrez:=NULL]
    uniq_entrez_ids <- as.numeric(names(table(df$entrez_id)[table(df$entrez_id) == 1]))

    message("Removing ", nrow(df) - length(uniq_entrez_ids), " gene entries due to duplicate entrez id")
    df <- df[entrez_id %in% uniq_entrez_ids, ]

    message("Total number of genes in core matrix after merge with entrez info: ", nrow(df))

    return(df)
}

#' Check whether a genomic position is in a PAR regions
#' Assumes chr as factor 1-22,X,Y,M
is_par <- function(chr, start_pos, end_pos, build = "b37") {
    if (build == "b37") {
        return((chr == "Y" & (end_pos <= 2649520 | start_pos >= 59034050)) |
               (chr == "X" & (end_pos <= 2699520 | start_pos >= 154931044)))
    } else {
        stop("Error: Genome Build ", build, " is not implemented")
    }
}

#' Core matrix is a merge between gencode, entrez and hgnc data.
#'
#' Gene definitions are based on gencode. Entrez ids and official hgnc symbols are merged in
#' Gene symbols in the resulting file are unique and only contains gencode genes that could be
#' reliably matched with entrez and hgnc data. All three data sources are preprocessed and filtered
#' before matching (see repsective preprocessing and merge functions for details).
#' @param gencode_url url to download gencode file from.
#' @param hgnc_url url to download hgnc file from.
#' @param entrez_url url to download entrez file from.
#' @param value_sep character string separating values in multi-value fields.
#' @param download_dir directory to save gencode, hgnc, and entrez source files.
#' @return a data.table with core gene matrix
#' @export
create_core_matrix <- function(gencode_url,hgnc_url,entrez_url,value_sep, download_dir) {

    gencode_path <- file.path(download_dir, basename(gencode_url))
    hgnc_path <- file.path(download_dir, basename(hgnc_url))
    entrez_path <- file.path(download_dir, basename(entrez_url))

    # Download gencode hgnc and entrez files if they do not already exist
    download_source(gencode_url, gencode_path)
    download_source(hgnc_url, hgnc_path)
    download_source(entrez_url, entrez_path)

    gencode <- process_gencodefile(gencode_path)
    hgnc <- process_hgncfile(hgnc_path, value_sep = value_sep)
    entrez <- process_entrezfile(entrez_path, value_sep = value_sep)
    core <- merge_gencode_hgnc(gencode, hgnc)
    core <- merge_core_entrez(core, entrez)
    core <- merge_url_links(core)
    return(core)
}


#' Merge url links to common resources
#'
#' @param core data.table with gene definition
#' @return data.table with core gene matrix including url links to relevant resources
#' @export
merge_url_links <- function(core){

  core[,ucsc_link := paste0("http://genome.ucsc.edu/cgi-bin/hgTracks?&org=Human&db=hg19&position=",chr_prefix,"%3A",format(start,scientific=F,trim=T),"-",format(end,scientific=F,trim=T))]
  core[,genecards_link := paste0("http://www.genecards.org/cgi-bin/carddisp.pl?gene=",symbol)]
  core[,wikipedia_link := paste0("https://en.wikipedia.org/wiki/",symbol)]
  core[,pubmed_link := paste0("https://www.ncbi.nlm.nih.gov/pubmed?term=",symbol)]
  core[,wikigenes_link := paste0("http://www.wikigenes.org/?search=",symbol)]

  return(core)
}


#' Load core gene matrix and create it first if it does not exist
#'
#' A cache wrapper around \code{create_core_matrix()}
#'
#' @param settings list with settings (defaults to global genematrix settings gm_settings)
#' @return data.table with core gene matrix
#' @export
get_core_matrix <- function(settings=gm_settings){
  for (setting in c("gencode_url","hgnc_url","entrez_url","value_sep","cache_dir","core_path"))
  {
    if(!setting %in% names(settings)) stop("ERROR: ",setting, " not in settings")
  }


  if(file.exists(settings$core_path)){
    load(settings$core_path)
  }else{
    core_matrix <- create_core_matrix(settings$gencode_url, settings$hgnc_url,
                               settings$entrez_url, settings$value_sep, settings$cache_dir)
    save(core_matrix,file=settings$core_path)
  }
  return(core_matrix)

}

