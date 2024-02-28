#######################
## FraserDataSet
## ====================


#' FraserDataSet
#'
#' This class is designed to store the whole FRASER data set
#' needed for an analysis of a disease cohort
#'
#' @author Christian Mertes \email{mertes@@in.tum.de}
setClass("FraserDataSet",
    contains="RangedSummarizedExperiment",
    slots = list(
        name            = "character",
        bamParam        = "ScanBamParam",
        workingDir      = "character",
        nonSplicedReads = "RangedSummarizedExperiment"
    ),
    prototype = list(
        name            = "Data Analysis",
        bamParam        = ScanBamParam(mapqFilter=0),
        workingDir      = "FRASER_output",
        nonSplicedReads = SummarizedExperiment(rowRanges=GRanges())
    )
)

## Validity
## ========

#
# check sample annotation within the colData slot of the SE object
#
validateSampleAnnotation <- function(object) {
    sampleData <- as.data.table(colData(object))
    if(nrow(sampleData) == 0){
        return(NULL)
    }
    if(!"sampleID" %in% colnames(sampleData)){
        return("Please provide a 'sampleID' column with a ID for each sample.")
    }
    if(any(duplicated(sampleData$sampleID))){
        return("The 'sampleID' column needs to be unique.")
    }
    if(any(samples(object) != rownames(colData(object)))){
        return("Please set the rownames of your colData to the sampleIDs")
    }
    NULL
}

validateName <- function(object){
    if(!isScalarCharacter(object@name)){
        return("The name of the experiment should be of type 'character'.")
    }
    if(length(object@name) == 0){
        return("The experiment name can not be empty.")
    }
    if(!grep("^[a-zA-Z0-9 ._-]+$", object@name, perl=TRUE)){
        return(paste("For readabilty the name of the experiment should only ",
                "contain the following characters: 'a-zA-Z0-9 ._-'"
        ))
    }
    NULL
}

validateBamParam <- function(object) {
    if(!is(scanBamParam(object), "ScanBamParam")) {
        return("The 'bamParam' option must be a ScanBamParam object.")
    }
    NULL
}

validateStrandSpecific <- function(object) {
    sampleData <- as.data.table(colData(object))
    if("strand" %in% colnames(sampleData) && 
            (!is.integer(sampleData[,strand]) 
                || any(!sampleData[,strand] %in% 0:2))){
        return(paste("The 'strand' column in the sample annotation in",
                "'colData(fds)' must be empty or contain only integers:",
                "0L == 'no', 1L == 'yes', 2L == 'reverse'."))
    }
    NULL
}

validatePairedEnd <- function(object) {
    sampleData <- as.data.table(colData(object))
    if("pairedEnd" %in% colnames(sampleData) && 
            any(!is.logical(sampleData[,pairedEnd]))){
        return(paste("The 'pairedEnd' column in the sample annotation in",
                "'colData(fds)' must only contain logical values ",
                "(TRUE or FALSE)."))
    }
    NULL
}

validateWorkingDir <- function(object) {
    if(!isScalarCharacter(object@workingDir)){
        return(paste("The path to the working directory needs",
                "to be a set as a character."))
    }
    if(object@workingDir == ""){
        return("The working directory can not be empty.")
    }
    if(dir.exists(object@workingDir)){
        if(file.access(object@workingDir, mode = 4) != 0){
            return(paste(
                    "Make sure we can write to the given working directory '",
                    object@workingDir, "'."))
        }
    }
    NULL
}

validateNonSplicedReadsType <- function(object) {
    if(!is(object@nonSplicedReads, "RangedSummarizedExperiment")) {
        return("'nonSplicedReads' must be a RangedSummarizedExperiment object")
    }
    if(length(object) != 0 && dim(object@nonSplicedReads)[2] != dim(object)[2]){
        return(paste("The nonSplicedReads dimensions are not correct.",
                "This is a internal error!"))
    }
    ans <- validObject(object@nonSplicedReads)
    if(!isScalarLogical(ans) || ans == FALSE){
        return(ans)
    }
    NULL
}

validateAssays <- function(object){
    if(length(assayNames(object)) > 1){
        if(any(duplicated(assayNames(object)))){
            return(
                "FRASER enforces unique assay names! Please provie such names."
            )
        }
    }
    NULL
}

# for non-empty fds objects check if non-spliced reads are overlapping 
# with at least 1 donor/acceptor site
validateNonSplicedReadsSanity <- function(object){
    # fds object must have samples and junctions
    if(all(dim(object) > c(0,0))){

        # fds object must be annotated with start/end/spliceSite indexes
        if(all(c("startID", "endID") %in% names(mcols(object))) &&
                        "spliceSiteID" %in% names(mcols(object, type="theta"))){

            # check that every spliceSiteID matches either a start or end index
            if(nrow(nonSplicedReads(object)) != length(intersect(
                    mcols(object, type="theta")$spliceSiteID, 
                    unlist(mcols(object)[, c("startID", "endID")])))){
                return(paste("The nonSplicedReads do not have corresponding",
                        "splitReads. This is probably the result of merging"))
            }
        }
    }
    NULL
}


## general validate function
validateFraserDataSet <- function(object) {
    c(
        validateSampleAnnotation(object),
        validatePairedEnd(object),
        validateName(object),
        validateBamParam(object),
        validateStrandSpecific(object),
        validateWorkingDir(object),
        validateNonSplicedReadsType(object),
        validateNonSplicedReadsSanity(object),
        validateAssays(object)
    )
}
setValidity("FraserDataSet", validateFraserDataSet)


## Cosmetics (the show function)
## =============================

## show method for FraserDataSet
showFraserDataSet <- function(object) {
    if(dim(object)[2] < 1){
        cat("This Fraser object does not contain any sample! Please add one.")
        return()
    }

    # taken from SummarizedExperiment show function
    scat <- function(fmt, vals=character(), exdent=2, ...){
        vals <- ifelse(nzchar(vals), vals, "''")
        lbls <- paste(S4Vectors:::selectSome(vals), collapse=" ")
        txt <- sprintf(fmt, length(vals), lbls)
        cat(strwrap(txt, exdent=exdent, ...), sep="\n")
    }
    cat("-------------------- Sample data table -----------------\n")
    sampleData <- as.data.table(colData(object))
    if("bamFile" %in% sampleData & all(
                vapply(sampleData$bamFile, isScalarCharacter, logical(1)))){
        sampleData$bamFile <- gsub("... [^/]+/", ".../",
            vapply(sampleData$bamFile, function(str){
                if(nchar(str) <= 29) return(str)
                paste("...", substr(str, nchar(str) - 25, nchar(str)))
            }, FUN.VALUE="")
        )
    }
    show(as_tibble(sampleData))
    cat("\n")
    
    if(length(object) > 0){
        cat(paste0("Number of samples:      ", dim(object)[2]), "\n")
        cat(paste0("Number of junctions:    ", length(object)), "\n")
        cat(paste0("Number of splice sites: ", 
                length(nonSplicedReads(object))), "\n")
        scat("assays(%d):    %s\n", assayNames(object))
        cat("\n")
    }

    cat("----------------------- Settings -----------------------\n")
    cat(paste0("Analysis name:               ", name(object)), "\n")
    cat(paste0("Analysis is strand specific: ", getStrandString(object)), "\n")
    cat(paste0("Working directory:           '", workingDir(object), "'"), "\n")
    cat("\n")
    
    cat("-------------------- BAM parameters --------------------\n")
    # if(identical(scanBamParam(FraserDataSet()), scanBamParam(object))){
    #     cat(paste0("Default used with: ",
    #             "bamMapqFilter=", bamMapqFilter(scanBamParam(object))
    #     ))
    # } else {
        show(scanBamParam(object))
    # }
    cat("\n\n")
}

setMethod("show", "FraserDataSet", function(object) {
    showFraserDataSet(object)
})

## Constructor
## ==========

#'
#' The FRASER dataset object
#'
#' Constructs an FRASER object based on the given input. It can take only the
#' annotation (colData) or count tables (junctions/spliceSites).
#'
#' @param colData A DataFrame containing the annotation of the samples
#' @param junctions,spliceSites A data.frame like object containing the 
#'                raw counts for each junction or splice site.
#'                It requires the columns \code{startID} and \code{endID} for 
#'                the junctions and \code{spliceSiteID} and \code{type} for the
#'                splice sites. Those columns identifies the corresponding
#'                splice site for the given junction and map to the splice site.
#'                For each sample the counts are saved in a corresponding 
#'                column with the same name. It can also be a GRange object.
#' @param ... Any parameters corresponding to the slots and their possible
#'                values. See \linkS4class{FraserDataSet}
#' @return A FraserDataSet object.
#' @author Christian Mertes \email{mertes@@in.tum.de}
#' @export
#' @examples
#'   fraser <- FraserDataSet()
#'   
#'   # example sample annoation
#'   sampleTable <- fread(system.file("extdata",
#'           "sampleTable_countTable.tsv", package="FRASER", mustWork=TRUE))
#'  
#'   # get raw counts 
#'   junctionCts   <- fread(system.file("extdata", 
#'           "raw_junction_counts.tsv.gz", package="FRASER", mustWork=TRUE))
#'   spliceSiteCts <- fread(system.file("extdata", 
#'           "raw_site_counts.tsv.gz", package="FRASER", mustWork=TRUE))
#'   
#'   # create FRASER object
#'   fds <- FraserDataSet(colData=sampleTable, junctions=junctionCts,
#'           spliceSites=spliceSiteCts, name="Example Dataset")
#'   
FraserDataSet <- function(colData=NULL, junctions=NULL, spliceSites=NULL, ...) {
    if(!is.null(colData)){
        if(is.data.table(colData)){
            colData <- DataFrame(colData)
        }
        if(is.null(rownames(colData))){
            rownames(colData) <- colData[['sampleID']]
        }
        if(is.null(junctions) & is.null(spliceSites))
            return(new("FraserDataSet", colData=colData, ...))
        if(is.null(junctions)){
            stop("Please provide junctions counts if you provide ",
                    "spliceSite counts.")
        }
        if(is.null(spliceSites)){
            stop("Please provdie splice site counts if you provide ",
                    "junction counts.")
        }
        if(!"sampleID" %in% colnames(colData)){
            stop("Please provide a column in colData specifying the sampleID.")
        }
        if(is.character(junctions)){
            if(!file.exists(junctions))
                stop("Junction file '", junctions, "' does not exists")
            junctions <- fread(junctions)
        }
        if(is.character(spliceSites)){
            if(!file.exists(spliceSites))
                stop("SpliceSite file '", spliceSites, "' does not exists")
            spliceSites <- fread(spliceSites)
        }
        if(is.data.frame(junctions)){
            junctions <- makeGRangesFromDataFrame(junctions,
                    keep.extra.columns=TRUE)
        }
        if(is.data.frame(spliceSites)){
            spliceSites <- makeGRangesFromDataFrame(spliceSites,
                    keep.extra.columns=TRUE)
        }
        if(!is(junctions, "GRanges")){
            stop("Provided junction object is not a GRanges object.")
        }
        if(!is(spliceSites, "GRanges")){
            stop("Provided spliceSite object is not a GRanges object.")
        }
        nsr <- SummarizedExperiment(
                rowRanges=spliceSites[,c("spliceSiteID", "type")],
                assays=SimpleList(rawCountsSS=as.data.frame(
                        mcols(spliceSites)[colData[,"sampleID"]]), a=NULL)[1])
        se <- SummarizedExperiment(
                rowRanges=junctions[,c("startID", "endID")],
                colData=colData,
                assays=SimpleList(rawCountsJ=as.data.frame(
                        mcols(junctions)[colData[,"sampleID"]]), a=NULL)[1])
        obj <- new("FraserDataSet", se, nonSplicedReads=nsr, ...)
    } else {
        obj <- new("FraserDataSet", ...)
    }
    
    metadata(obj)[["version"]] <- packageVersion("FRASER")

    validObject(obj)
   
    return(obj)
}



