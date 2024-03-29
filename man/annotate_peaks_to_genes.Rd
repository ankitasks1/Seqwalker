% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seqwalker.R
\name{annotate_peaks_to_genes}
\alias{annotate_peaks_to_genes}
\title{bedtools annotation}
\usage{
annotate_peaks_to_genes(
  peakfileslist,
  genefile,
  path,
  assaytype,
  software,
  outformat
)
}
\arguments{
\item{peakfileslist}{list of peak file}

\item{genefile}{datframe of gene file}

\item{path}{location of peaks files}

\item{assaytype}{eg. "atacseq"}

\item{software}{eg. "deseq2"}

\item{outformat}{eg. ".txt"}
}
\value{
distance calculated by peakfile and gene using bedtools closest  with -d option
}
\description{
bedtools annotation
}
\examples{
anno <- annotate_peaks_to_genes(public_list$peakfilelist, public_list$gene, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration","chipseq", "encode", "txt")
}
