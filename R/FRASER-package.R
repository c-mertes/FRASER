#' 
#' FRASER: A package providing a workflow to detect aberrant splicing events
#'     in RNA sequencing data in a rare disease cohort
#'
### PACKAGE NAME
#'
#' @name FRASER
#' @noRd
#' 
#' @import data.table
#'
#' @importFrom Biobase rowMax
#'
### Import of generic functions
#' 
#' @importMethodsFrom OUTRIDER results
#' 
### Parallel computing
#'
#' @import BiocParallel
#' @importFrom pcaMethods pca loadings
#'
### GRange/Experiment/bamFile packages
#' @importFrom BiocGenerics updateObject counts counts<- strand strand<- which
#' @importFrom GenomicFeatures makeTxDbFromGFF intronsByTranscript genes exons
#'          fiveUTRsByTranscript threeUTRsByTranscript
#' @importFrom GenomicAlignments junctions readGAlignments summarizeJunctions
#'          readGAlignmentPairs
#' @importFrom SummarizedExperiment assay assay<- assays assays<- assayNames
#'          colData colData<- rowData rowRanges rowRanges<- SummarizedExperiment
#'          rbind Assays
#' @importFrom GenomicRanges findOverlaps granges GRanges GRangesList
#'          makeGRangesFromDataFrame invertStrand start end start<- end<-
#' @importFrom IRanges subsetByOverlaps from to IRanges ranges nearest distance
#' @importFrom Rsamtools ScanBamParam scanBamHeader bamMapqFilter
#'          bamWhich bamWhich<- BamFile idxstatsBam
#' @importFrom Rsubread featureCounts
#' @importFrom BSgenome getBSgenome
#'
### Annotation
#'
#' @importFrom biomaRt useEnsembl getBM
#' @importFrom AnnotationDbi select
#'
### Plotting
#'
#' @importFrom plotly plot_ly subplot layout add_trace ggplotly
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom cowplot theme_cowplot
#' @importFrom ggrepel geom_text_repel
#'
#'
### Data handling
#'
#' @importFrom HDF5Array writeHDF5Array path HDF5Array 
#'          saveHDF5SummarizedExperiment loadHDF5SummarizedExperiment
#' @importFrom rhdf5 h5ls H5Fopen H5Fclose H5Pclose H5Dget_create_plist
#'          H5Pget_layout H5Pget_chunk h5delete
#'
### P-Value calculation
#'
#' @importFrom stats sd rbinom na.omit p.adjust ppoints qbeta rnorm
#'          predict cor cutree dbinom dist hclust lm optim optimize pbinom
#'          rlnorm rnbinom pnorm
#' @importFrom VGAM rbetabinom vglm Coef pbetabinom pbetabinom.ab betabinomial
#'          dbetabinom.ab dbetabinom
#'
### Miscelenious functions
#'
#' @importFrom BBmisc isScalarCharacter isScalarLogical chunk %nin% seq_row
#'          seq_col isScalarInteger isFALSE is.error isScalarValue 
#'          isScalarNumeric
#' @importFrom R.utils renameFile withTimeout
#' @importFrom tools file_path_as_absolute
#' @importFrom methods as callNextMethod is new show slot slot<- validObject
#' @importFrom utils capture.output packageVersion
#'
#'
#'
### To be added into the functions above
#'
#' @importFrom S4Vectors DataFrame metadata Rle SimpleList mcols mcols<-
#'          start end metadata metadata<- subjectHits queryHits elementMetadata
#' @importFrom grDevices colorRampPalette
#' @importFrom GenomeInfoDb keepStandardChromosomes seqlevels<- seqlevels
#'          seqlengths seqlengths<- seqlevelsStyle<- seqlevelsStyle seqnames 
#'          seqinfo standardChromosomes dropSeqlevels keepSeqlevels
#' @importFrom DelayedArray rowMaxs rowMeans path<- cbind plogis qlogis 
#'          DelayedArray
#' @importFrom DelayedMatrixStats colSds rowMedians rowSds colMeans2 rowMeans2
#'          rowQuantiles
#' @importFrom matrixStats colMaxs colMedians colMins colAnys
#' @importFrom stats runif median quantile
#' @importFrom extraDistr rdirmnom dbbinom pbbinom
#' @importFrom PRROC pr.curve
#' @importFrom ggplot2 ggtitle xlab ylab ggplot geom_point geom_line 
#'          geom_smooth aes geom_line geom_hline geom_vline geom_abline 
#'          geom_segment geom_ribbon scale_color_manual scale_x_log10 
#'          scale_y_log10 scale_color_gradientn labs theme_bw theme
#'          scale_color_brewer scale_color_discrete scale_linetype_manual 
#'          annotate geom_histogram scale_fill_manual xlim scale_colour_manual
#'          element_blank annotation_logticks
#'
#' @importFrom tibble as_tibble
#'
#' @useDynLib FRASER
#'
NULL

#' TODO: Should be removed in the end!
#' @noRd
globalVariables(c(".", "J", ".N", ".asDataFrame", "End", "first_feature", 
        "FN", "HTML", "other_features", "Start", 
        "TP", "deltaPsi", "curgr", "gene", "lty", "hgnc_symbol", "id",
        "ldat", "mapped", "p.adj", "pval", "pValue", "sampleID", "sampleGroup", 
        "chr", "isOptimalQ",
        "symbol", "type", "pseudocount", "known_intron", "known_start", 
        "known_end", "..columns2Write", "maxVal", "idxGroup", "idxInCount", 
        "idxInGroup", "lower", "nOk", "nPerGroup", "nsubset", "obs", "p", 
        "plotPoint", "rn", "upper", 
        "..sum_cols", "groupID", "groupSize", "Aggregation",
        "aroc", "dpsi", "exMask", "fact",
        "fdrByFeature", "fds_inputK", "fds_inputN", "fds_new", "fds_rhoIn",
        "featureID", "features", "flog.debug", "flog.error", "flog.fatal",
        "flog.info", "flog.namespace", "flog.trace", "flog.warn", "geneID",
        "idx", "iteration", "Iteration", "IterSteps", "junctionID", "k", "Loss",
        "model", "mu", "n", ",nsubset", "o3", "o5", "obsPsi", "os", "pa",
        "padj", "passed", "pByFeature", "pointNr", "predPsi", "psi3", "psi5",
        "psiType", "psiValue", "seqlength", "seqlevel", "Step", "traceNr",
        "uniqueID", "V1", "value", "zscore", "maxDTheta", 
        "blacklist", "potentialImpact", "causesFrameshift", "annotatedJunction", 
        "distNearestGene", "UTR_overlap", "meanCount", "medianCount", 
        "potentialImpact2"),
        package="FRASER")
