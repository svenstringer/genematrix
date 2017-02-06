#' Process exac pli constraint file
#' @export
process_exacpli <- function(pli_suffix,settings=gm_settings) {
    stopifnot(pli_suffix %in% c("fullexac","nonpsychexac"))

    if (pli_suffix == "fullexac") {
      pli_url <- settings$fullexacpli_url
    }
    if (pli_suffix == "nonpsychexac") {
      pli_url <- settings$nonpsychexacpli_url
    }

    message("Process exac pLi file", pli_suffix)

    col_selection <- c("transcript", "gene", "chr", "cds_start", "pLI", "pRec", "pNull")
    pli_path <- file.path(settings$cache_dir, basename(pli_url))

    stopifnot(!is.null(pli_path))

    download_source(pli_url, pli_path)

    pli <- fread(pli_path, colClasses = list(character = "chr"))

    if (pli_url == settings$nonpsychexacpli_url) {
        old_names <- c("tx_start", "pRecessive")
        new_names <- c("cds_start", "pRec")
        setnames(pli, old_names, new_names)
    }
    pli <- pli[, col_selection, with = F]
    return(pli)
}

is_pli_match <- function(chr_pli, chr, start, cds_start, gene_bp_dif) {
    !is.na(chr_pli) & as.character(chr) == chr_pli & abs(start - cds_start) < gene_bp_dif
}


#' Merge exac pli constraint data into gene matrix
#' @export
merge_exacpli <- function(pli_suffix, gene_matrix, gene_translation_table, settings=gm_settings) {

    gene_bp_dif <- settings$gene_bp_dif

    stopifnot(pli_suffix %in% c("fullexac","nonpsychexac"))

    if (pli_suffix == "fullexac") {
      pli_url <- settings$fullexacpli_url
    }
    if (pli_suffix == "nonpsychexac") {
        pli_url <- settings$nonpsychexacpli_url
    }

    setkey(gene_matrix, symbol)

    pli <- process_exacpli(pli_suffix,settings)
    setkey(pli, gene)

    symbols <- lookup_symbol(pli$gene, gene_translation_table)
    stopifnot(length(symbols) == nrow(pli))

    max_n_genes <- max(as.vector(sapply(symbols, function(x) length(x))))

    unmatched_df <- gene_matrix


    # Loop over possible symbol aliases and retry merge gene_matrix and pli score matrix
    for (i in 1:max_n_genes) {
        pli[, lookup_symbol := as.vector(sapply(symbols, function(x) x[i]))]
        merged_df <- merge(unmatched_df[, names(gene_matrix), with = F], pli, by.x = "symbol", by.y = "lookup_symbol",
            suffixes = c("", "_pli"), all.x = T)
        unmatched_df <- merged_df[!is_pli_match(chr_pli, chr, start, cds_start, gene_bp_dif), ]
        new_matches <- merged_df[is_pli_match(chr_pli, chr, start, cds_start, gene_bp_dif), ]
        stopifnot(nrow(unmatched_df) + nrow(new_matches) == nrow(merged_df))
        if (i == 1) {
            final_df <- new_matches
        } else {
            final_df <- rbind(final_df, new_matches)
        }
    }
    pli_colnames <- setdiff(names(unmatched_df), names(gene_matrix))

    unmatched_df[, (pli_colnames) := NA]

    final_df <- rbind(final_df, unmatched_df)

    # Remove pli matches that map to multiple genes
    duplicate_genes <- names(table(final_df$symbol)[table(final_df$symbol) > 1])
    final_df[symbol %in% duplicate_genes, (pli_colnames) := NA]

    final_df <- unique(final_df)

    stopifnot(nrow(final_df) == nrow(gene_matrix))

    n_matches <- nrow(final_df[is_pli_match(chr_pli, chr, start, cds_start, gene_bp_dif), ])
    message(n_matches, " genes out of ", nrow(final_df), " are uniquely matched to genes in pli_", pli_suffix, " data")

    final_df[, c("gene", "chr_pli", "cds_start") := NULL]

    # Rename
    old_names <- c("transcript", "pLI", "pRec", "pNull")
    setnames(final_df, old_names, paste0(old_names, "_", pli_suffix))

    final_df
}
