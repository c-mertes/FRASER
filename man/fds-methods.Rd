% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics-definitions.R, R/AllGenerics.R
\name{samples}
\alias{samples}
\alias{samples<-}
\alias{condition}
\alias{condition<-}
\alias{bamFile}
\alias{bamFile<-}
\alias{name}
\alias{name<-}
\alias{strandSpecific}
\alias{strandSpecific<-}
\alias{pairedEnd}
\alias{pairedEnd<-}
\alias{workingDir}
\alias{workingDir<-}
\alias{scanBamParam}
\alias{scanBamParam<-}
\alias{nonSplicedReads}
\alias{nonSplicedReads<-}
\alias{fds-methods}
\alias{samples,FraserDataSet-method}
\alias{samples<-,FraserDataSet-method}
\alias{condition,FraserDataSet-method}
\alias{condition<-,FraserDataSet-method}
\alias{bamFile,FraserDataSet-method}
\alias{bamFile<-,FraserDataSet-method}
\alias{name,FraserDataSet-method}
\alias{name<-,FraserDataSet-method}
\alias{workingDir,FraserDataSet-method}
\alias{workingDir<-,FraserDataSet-method}
\alias{strandSpecific,FraserDataSet-method}
\alias{strandSpecific<-,FraserDataSet-method}
\alias{pairedEnd,FraserDataSet-method}
\alias{pairedEnd<-,FraserDataSet-method}
\alias{scanBamParam,FraserDataSet-method}
\alias{scanBamParam<-,FraserDataSet-method}
\alias{nonSplicedReads,FraserDataSet-method}
\alias{nonSplicedReads<-,FraserDataSet-method}
\alias{FRASER.mcols.get}
\alias{FRASER.rowRanges.get}
\alias{mapSeqlevels}
\title{Getter/Setter methods for the FraserDataSet}
\usage{
samples(object)

samples(object) <- value

condition(object)

condition(object) <- value

bamFile(object)

bamFile(object) <- value

name(object)

name(object) <- value

strandSpecific(object)

strandSpecific(object) <- value

pairedEnd(object)

pairedEnd(object) <- value

workingDir(object)

workingDir(object) <- value

scanBamParam(object)

scanBamParam(object) <- value

nonSplicedReads(object)

nonSplicedReads(object) <- value

\S4method{samples}{FraserDataSet}(object)

\S4method{samples}{FraserDataSet}(object) <- value

\S4method{condition}{FraserDataSet}(object)

\S4method{condition}{FraserDataSet}(object) <- value

\S4method{bamFile}{FraserDataSet}(object)

\S4method{bamFile}{FraserDataSet}(object) <- value

\S4method{name}{FraserDataSet}(object)

\S4method{name}{FraserDataSet}(object) <- value

\S4method{workingDir}{FraserDataSet}(object)

\S4method{workingDir}{FraserDataSet}(object) <- value

\S4method{strandSpecific}{FraserDataSet}(object)

\S4method{strandSpecific}{FraserDataSet}(object) <- value

\S4method{pairedEnd}{FraserDataSet}(object)

\S4method{pairedEnd}{FraserDataSet}(object) <- value

\S4method{scanBamParam}{FraserDataSet}(object)

\S4method{scanBamParam}{FraserDataSet}(object) <- value

\S4method{nonSplicedReads}{FraserDataSet}(object)

\S4method{nonSplicedReads}{FraserDataSet}(object) <- value

FRASER.mcols.get(x, type = NULL, ...)

FRASER.rowRanges.get(x, type = NULL, ...)

mapSeqlevels(fds, style = "UCSC", ...)
}
\arguments{
\item{object}{A FraserDataSet object.}

\item{value}{The new value that should replace the current one.}

\item{x}{A FraserDataSet object.}

\item{type}{The psi type (psi3, psi5 or theta)}

\item{...}{Further parameters. For mapSeqLevels: further parameters 
passed to GenomeInfoDb::mapSeqlevels().}

\item{fds}{FraserDataSet}

\item{style}{The style of the chromosome names.}
}
\value{
Getter method return the respective current value.
}
\description{
The following methods are getter and setter methods to extract 
or set certain values of a FraserDataSet object. 

\code{samples} sets or gets the sample IDs; \code{condition} ;
\code{}
\code{nonSplicedReads} return a RangedSummarizedExperiment object 
containing the counts for the non spliced reads overlapping splice 
sites in the fds.
\code{}
}
\examples{
fds <- createTestFraserDataSet()
samples(fds)
samples(fds) <- 1:dim(fds)[2]
condition(fds)
condition(fds) <- 1:dim(fds)[2]
bamFile(fds) # file.paths or objects of class BamFile
bamFile(fds) <- file.path("bamfiles", samples(fds), "rna-seq.bam")
name(fds)
name(fds) <- "My Analysis"
workingDir(fds)
workingDir(fds) <- tempdir()
strandSpecific(fds)
strandSpecific(fds) <- TRUE
strandSpecific(fds) <- "reverse"
strandSpecific(fds)
scanBamParam(fds)
scanBamParam(fds) <- ScanBamParam(mapqFilter=30)
nonSplicedReads(fds)
rowRanges(fds)
rowRanges(fds, type="theta")
mcols(fds, type="psi5")
mcols(fds, type="theta")
seqlevels(fds)
seqlevels(mapSeqlevels(fds, style="UCSC"))
seqlevels(mapSeqlevels(fds, style="Ensembl"))
seqlevels(mapSeqlevels(fds, style="dbSNP"))

}
\author{
Christian Mertes \email{mertes@in.tum.de}

Ines Scheller \email{scheller@in.tum.de}
}
