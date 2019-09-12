\name{dmrcatedata}
\alias{dmrcatedata}
\docType{data}
\title{Accompanying data package for DMRcate}
\description{
This package contains a series of objects a) used to filter Illumina array data based on SNP, cross-hybridising and sex chromosome information, b) used to annotate bisulfite-sequencing generated DMRs with their overlapping genes, and c) used to plot DMRs in \code{DMR.plot()}.
}
\usage{data(dmrcatedata)}
\format{character vector (\code{crosshyb}), data.frame (\code{snpsall}), GRanges (\code{hg19.generanges}, \code{hg38.generanges} and \code{mm10.generanges}) and GeneRegionTrack (\code{hg19.grt}, \code{hg38.grt} and \code{mm10.grt}).}
\source{
\code{snpsall} sourced from \url{http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1066-1}, accessed September 2019, \url{http://www.sickkids.ca/MS-Office-Files/Research/Weksberg Lab/48640-polymorphic-CpGs-Illumina450k.xlsx}, accessed September 2019
\code{crosshyb} sourced from \url{http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1066-1}, accessed September 2019, \url{http://www.sickkids.ca/MS-Office-Files/Research/Weksberg Lab/48639-non-specific-probes-Illumina450k.xlsx}, accessed September 2019;
\code{hg19.generanges}, \code{hg19.grt}, \code{hg38.generanges}, \code{hg38.grt}, \code{mm10.generanges} and \code{mm10.grt} generated from inst/scripts/make-Data.R.}
\keyword{datasets}