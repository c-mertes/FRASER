% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotationOfRanges.R
\name{annotateRanges}
\alias{annotateRanges}
\alias{annotateRangesWithTxDb}
\title{Annotates the given FraserDataSet with the HGNC symbol with biomaRt}
\usage{
annotateRanges(
  fds,
  feature = "hgnc_symbol",
  featureName = feature,
  biotype = list("protein_coding"),
  ensembl = NULL,
  GRCh = 37
)

annotateRangesWithTxDb(
  fds,
  feature = "SYMBOL",
  featureName = "hgnc_symbol",
  keytype = "ENTREZID",
  txdb = NULL,
  orgDb = NULL,
  filter = list()
)
}
\arguments{
\item{fds}{FraserDataSet}

\item{feature}{Defines which feature (default is HGNC symbol) should be 
annotated. Has to be the \code{biomaRt} feature name or a 
column name in \code{orgDb}.}

\item{featureName}{The column name of the feature in the FraserDataSet mcols.}

\item{biotype}{The biotype for \code{biomaRt}.}

\item{ensembl}{The ensembl that should be used. If NULL, the default one is 
used (hsapiens_gene_ensembl, GRCh37).}

\item{GRCh}{GRCh version to connect to. If this is NULL, then the current 
GRCh38 is used. Otherwise, this can only be 37 (default) 
at the moment (see \code{\link[biomaRt]{useEnsembl}}).}

\item{keytype}{The keytype or column name of gene IDs in the \code{TxDb}
object (see 
\code{\link[AnnotationDbi:AnnotationDb-class]{keytypes}}
for a list of available ID types).}

\item{txdb}{A \code{TxDb} object. If this is NULL, then the default 
one is used, currently this is 
\code{TxDb.Hsapiens.UCSC.hg19.knownGene}.}

\item{orgDb}{An \code{orgDb} object or a data table to map the feature names.
If this is NULL, then \code{org.Hs.eg.db} is used as the default.}

\item{filter}{A named list specifying the filters which should be applied to 
subset to e.g. only protein-coding genes for annotation. 
\code{names(filter)} needs to be column names in the given 
orgDb object (default: no filtering).}
}
\value{
FraserDataSet
}
\description{
Annotates the given FraserDataSet with the HGNC symbol with biomaRt
}
\examples{

fds <- createTestFraserDataSet()

### Two ways to annotage ranges with gene names: 
# either using biomart with GRCh38
try({
  fds <- annotateRanges(fds, GRCh=38)
  rowRanges(fds, type="j")[,c("hgnc_symbol")]
})

# either using biomart with GRCh37
try({
  fds <- annotateRanges(fds, featureName="hgnc_symbol_37", GRCh=37)
  rowRanges(fds, type="j")[,c("hgnc_symbol_37")]
})
 
# or with a provided TxDb object
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
require(org.Hs.eg.db)
orgDb <- org.Hs.eg.db
fds <- annotateRangesWithTxDb(fds, txdb=txdb, orgDb=orgDb)
rowRanges(fds, type="j")[,"hgnc_symbol"]

}
