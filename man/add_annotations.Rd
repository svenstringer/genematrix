% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_genematrix.R
\name{add_annotations}
\alias{add_annotations}
\title{Add annotations to core gene matrix}
\usage{
add_annotations(core, gene_translation_table, settings)
}
\arguments{
\item{core}{core gene matrix to add annotation to}

\item{gene_translation_table}{list mapping aliases to official symbol used in core gene matrix}
}
\value{
returns a data.table with annotated genes
}
\description{
Adds several types of gene annotation to the core matrix, such as constraint scores,
gene-based pvalues, GWAS catalog associations, and OMIM associations.
}

