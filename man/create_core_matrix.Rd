% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corematrix_functions.R
\name{create_core_matrix}
\alias{create_core_matrix}
\title{Core matrix is a merge between gencode, entrez and hgnc data.}
\usage{
create_core_matrix(gencode_url, hgnc_url, entrez_url, value_sep, download_dir)
}
\arguments{
\item{gencode_url}{url to download gencode file from.}

\item{hgnc_url}{url to download hgnc file from.}

\item{entrez_url}{url to download entrez file from.}

\item{value_sep}{character string separating values in multi-value fields.}

\item{download_dir}{directory to save gencode, hgnc, and entrez source files.}
}
\value{
a data.table with core gene matrix
}
\description{
Gene definitions are based on gencode. Entrez ids and official hgnc symbols are merged in
Gene symbols in the resulting file are unique and only contains gencode genes that could be
reliably matched with entrez and hgnc data. All three data sources are preprocessed and filtered
before matching (see repsective preprocessing and merge functions for details).
}

