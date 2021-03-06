#' Standardize gene names to upper case and remove punctuations (other than .)
#' @param genenames vector of genenames
#' @return a vector of sanitized gene names (capital letters, most punctuations removed)
uniform_genenames <- function(genenames) {
    unif_symbols <- toupper(genenames)
    unif_symbols <- gsub("[^A-Z0-9\\.]", "", unif_symbols)
    return(unif_symbols)
}

#' Load gene name translation table and create if it does not exist
#' @param core_matrix data.frame with core matrix
#' @param settings list with global settings
#' @return a data.frame with gene_translation_table
#' @export
get_symbol_table <- function(core_matrix, settings) {
  if(file.exists(settings$gtt_path)){
    load(settings$gtt_path)
  }else{
    gene_translation_table <- create_symbol_table(core_matrix,settings$value_sep,settings$add_excel_collisions)
    save(gene_translation_table,file=settings$gtt_path)
  }
  return(gene_translation_table)
}

#' Create gene symbol translation table
#'
create_symbol_table <- function(core, value_sep, excel_col = T) {

    message("Creating gene translation table... (may take a couple of minutes)")

    stopifnot(all(c("aliases", "symbol") %in% names(core)))

    symbol_table <- list()

    symbols <- core$symbol
    aliases <- core$aliases

    for (i in 1:length(symbols)) {
        if (i%%1000 == 0)
            message(i)
        # add official symbol, convert to uppercase, remove nondigit/nonletter/nondot chars, remove duplicates
        current_aliases <- uniform_genenames(c(symbols[i], strsplit(aliases[i], split = value_sep)[[1]]))

        for (al in current_aliases) {
            if (is.null(symbol_table[[al]])) {
                symbol_table[[al]] <- symbols[i]
            } else {
                symbol_table[[al]] <- unique(c(symbol_table[[al]], symbols[i]))
            }
        }
    }

    if (excel_col)
        symbol_table <- add_excel_collisions(symbol_table)

    if(any(is.na(symbol_table[[al]]))){stop(paste(current_aliases,al,i))}

    return(symbol_table)
}

#' Symbol
add_to_symbol_table <- function(symbol_alias_df, symbol_table, ignore_new_symbols = F) {
    df <- symbol_alias_df[ignore_new_symbols == F | alias_in_table(df$alias, symbol_table), ]
    stopifnot(sort(names(df)) == c("alias", "symbol"))

    for (i in 1:nrow(df)) {
        pair <- df[i, ]
        if (ignore_new_symbols == F | symbol_in_table(pair$symbol, symbol_table)) {
            symbol_table <- add_alias_to_symbol_table(pair$alias, pair$symbol, symbol_table)
        }
    }
    return(symbol_table)
}

symbol_in_table <- function(symbol, symbol_table) {
    return(symbol %in% unlist(symbol_table))
}

alias_in_table <- function(alias, symbol_table) {
    return(alias %in% names(symbol_table))
}


add_alias_to_symbol_table <- function(alias, symbol, symbol_table) {
    if (is.null(symbol_table[[alias]])) {
        symbol_table[[alias]] <- symbol
    } else {
        symbol_table[[alias]] <- unique(c(symbol_table[[alias]], symbol))
    }
    return(symbol_table)
}

add_excel_collisions <- function(symbol_table) {
    aliases <- names(symbol_table)

    months_full <- c("JANUARY", "FEBRUARY", "MARCH", "APRIL", "JUNE", "JULY", "AUGUST", "SEPTEMBER", "OCTOBER", "NOVEMBER",
        "DECEMBER")
    # gene aliases starting with letters ending with digits (potential collisions)
    df <- data.frame(alias = aliases[grepl("^[A-Z]+[0-9]+$", aliases)], stringsAsFactors = F)
    df$alias_prefixes <- as.character(sapply(df$alias, function(x) substr(x, 1, attr(regexpr("^[A-Z]+", x), "match.length"))))
    df$alias_digits <- as.integer(sapply(df$alias, function(x) substr(x, 1 + nchar(x) - attr(regexpr("[0-9]+$", x),
        "match.length"), nchar(x))))

    df$matches <- pmatch(df$alias_prefixes, months_full, duplicates.ok = T)
    df <- subset(df, sapply(1:length(df$alias_prefixes), function(i) !is.na(df$matches[i]) & nchar(df$alias_prefixes[i]) >
        2))

    df <- subset(df, nchar(df$alias_prefixes) > 2 & as.numeric(df$alias_digits) <= 31)
    df$matches <- NULL

    df$conversion <- paste0(df$alias_digits, substr(df$alias_prefixes, 1, 3))

    new_table <- symbol_table
    message("Add excel conversion genes to gene translation table")

    for (i in 1:nrow(df)) {
        alias <- df[i, "conversion"]
        symbols <- unique(as.character(lookup_symbol(df[i, "alias"], symbol_table)[[1]]))
        length(symbols)
        for (j in 1:length(symbols)) {
            message("Add ", alias, " -> ", symbols[j], " to gene translation table")
            new_table <- add_alias_to_symbol_table(alias, symbols[j], new_table)
        }
    }

    return(new_table)
}

#' Lookup gene symbol in gene-translation-table (gtt)
#'
#' @param symbols character of gene symbols to look up in gtt
#' @param symbol_table list containing gene_translation_table (mappging from aliases to official gene name)
#' @return character vector with official gene names
#' @export
lookup_symbol <- function(symbols, symbol_table) {
    names_std <- uniform_genenames(symbols)
    res <- sapply(names_std, function(nm) ifelse(is.null(symbol_table[[nm]]), list(NA), list(symbol_table[[nm]])))
    return(res)
}
