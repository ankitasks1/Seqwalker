% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seqwalker.R
\name{diffbind_norm}
\alias{diffbind_norm}
\title{diffbind normalize counts and perform DE analysis}
\usage{
diffbind_norm(obj)
}
\arguments{
\item{obj}{diffbind object}
}
\value{
Normalize and prepare contrast
}
\description{
diffbind normalize counts and perform DE analysis
}
\examples{
norm <- diffbind_norm(atacseqkd_diffbind_list$counts$dba_obj)
}
