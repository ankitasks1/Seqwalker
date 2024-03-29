% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seqwalker.R
\name{maplot}
\alias{maplot}
\title{ggmaplot maplot}
\usage{
maplot(df, columnbyorder, fdr, fc)
}
\arguments{
\item{df}{dataframe of output from deseq2 or other ways}

\item{columnbyorder}{since columns we renamed again the exact order need to be specified}

\item{fdr}{p.adj-value threshold}

\item{fc}{fold change threshold}
}
\value{
MA plot
}
\description{
ggmaplot maplot
}
\examples{
ma <- maplot(atacseqkd_diffbind_limma_list$dar_analysis$dars$coef_shSUZ12_shControl$all, c(2,1,5,7), 0.05, 2)
}
