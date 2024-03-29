% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seqwalker.R
\name{go_term_analysis}
\alias{go_term_analysis}
\title{gprofiler go term analysis}
\usage{
go_term_analysis(gene_list, organism)
}
\arguments{
\item{gene_list}{list of mutiple samples with genes of vectors}

\item{organism}{eg. "hg38"}
}
\value{
list of outputs related to gprofiler and gost
}
\description{
gprofiler go term analysis
}
\examples{
go <- go_term_analysis(list_rnaseqkd_de_shC1_c1509_shS_c1509, "hsapiens")
}
