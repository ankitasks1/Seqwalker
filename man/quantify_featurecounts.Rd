% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seqwalker.R
\name{quantify_featurecounts}
\alias{quantify_featurecounts}
\title{featurecounts quantify per feature in a given interval}
\usage{
quantify_featurecounts(
  peaks_path,
  assaytype,
  bamfiles,
  sites_files,
  sites_type,
  sites_column_to_rearrange,
  pairedend = FALSE,
  refgenome,
  delim,
  merge_sites_files = FALSE
)
}
\arguments{
\item{peaks_path}{location of peaks files}

\item{assaytype}{eg. "atacseq"}

\item{bamfiles}{eg. "\\.bam$"}

\item{sites_files}{eg. "\\_peaks_id.bed$"}

\item{sites_type}{eg. "histone_marks"}

\item{sites_column_to_rearrange}{This is suitable for 12-bed format as suggested in Footnotes
eg. c(12,4,1:3,7). If you have specific format of peaks files please rearrange appropriately
but makes sure you give correct index}

\item{pairedend}{boolean {TRUE, FALSE}}

\item{refgenome}{eg. "hg38}

\item{delim}{eg. "_", the separator you want to use for assigning new id}

\item{merge_sites_files}{boolean {TRUE, FALSE} keep it FALSE}
}
\value{
It will out put the quantification of bam files withing a given peak(s)/region(s) of interest
}
\description{
featurecounts quantify per feature in a given interval
}
\examples{
featurecounts_out <- quantify_featurecounts("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/", "atacseqkd","\\\\.bam$", "\\\\_peaks_id.bed$","histone_marks", c(12,4,1:3,7), pairedend=TRUE, "hg38", "_",merge_sites_files=FALSE)
}
