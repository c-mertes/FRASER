asSE <- function(x){
    if(is(x, "RangedSummarizedExperiment")){
        return(as(x, "RangedSummarizedExperiment"))
    }
    as(x, "SummarizedExperiment")
}

asFDS <- function(x){
    return(as(x, "FraserDataSet"))
}

#'
#' @title Getter/Setter methods for the FraserDataSet
#'
#' @description The following methods are getter and setter methods to extract 
#' or set certain values of a FraserDataSet object. 
#' 
#' \code{samples} sets or gets the sample IDs; \code{condition} ;
#' \code{}
#' \code{nonSplicedReads} return a RangedSummarizedExperiment object 
#' containing the counts for the non spliced reads overlapping splice 
#' sites in the fds.
#' \code{}
#'
#' @param object A FraserDataSet object.
#' @param value The new value that should replace the current one.
#' @param x A FraserDataSet object.
#' @param type The psi type (psi3, psi5 or theta)
#' @return Getter method return the respective current value.
#' @examples
#' fds <- createTestFraserDataSet()
#' samples(fds)
#' samples(fds) <- 1:dim(fds)[2]
#' condition(fds)
#' condition(fds) <- 1:dim(fds)[2]
#' bamFile(fds) # file.paths or objects of class BamFile
#' bamFile(fds) <- file.path("bamfiles", samples(fds), "rna-seq.bam")
#' name(fds)
#' name(fds) <- "My Analysis"
#' workingDir(fds)
#' workingDir(fds) <- tempdir()
#' strandSpecific(fds)
#' strandSpecific(fds) <- TRUE
#' strandSpecific(fds) <- "reverse"
#' strandSpecific(fds)
#' scanBamParam(fds)
#' scanBamParam(fds) <- ScanBamParam(mapqFilter=30)
#' nonSplicedReads(fds)
#' rowRanges(fds)
#' rowRanges(fds, type="theta")
#' mcols(fds, type="psi5")
#' mcols(fds, type="theta")
#' seqlevels(fds)
#' seqlevels(mapSeqlevels(fds, style="UCSC"))
#' seqlevels(mapSeqlevels(fds, style="Ensembl"))
#' seqlevels(mapSeqlevels(fds, style="dbSNP"))
#' 
#' @author Christian Mertes \email{mertes@@in.tum.de}
#' @author Ines Scheller \email{scheller@@in.tum.de}
#'
#' @rdname fds-methods
#' @name fds-methods
NULL


#' @rdname fds-methods
#' @export
setMethod("samples", "FraserDataSet", function(object) {
    return(as.character(colData(object)[,"sampleID"]))
})

#' @rdname fds-methods
#' @export
setReplaceMethod("samples", "FraserDataSet", function(object, value) {
    colData(object)[,"sampleID"] <- as.character(value)
    rownames(colData(object)) <- colData(object)[,"sampleID"]
    validObject(object)
    return(object)
})


#' @export
#' @rdname fds-methods
setMethod("condition", "FraserDataSet", function(object) {
    if("condition" %in% colnames(colData(object))){
        return(colData(object)[,"condition"])
    }
    return(seq_len(dim(object)[2]))
})

#' @export
#' @rdname fds-methods
setReplaceMethod("condition", "FraserDataSet", function(object, value) {
    colData(object)[,"condition"] <- value
    validObject(object)
    return(object)
})


#' @export
#' @rdname fds-methods
setMethod("bamFile", "FraserDataSet", function(object) {
    bamFile <- colData(object)[,"bamFile"]
    if(all(vapply(bamFile, is, class2="BamFile", FUN.VALUE=logical(1)))){
        bamFile <- vapply(bamFile, path, "")
    }
    return(bamFile)
})

#' @export
#' @rdname fds-methods
setReplaceMethod("bamFile", "FraserDataSet", function(object, value) {
    colData(object)[,"bamFile"] <- value
    validObject(object)
    return(object)
})


#' @export
#' @rdname fds-methods
setMethod("name", "FraserDataSet", function(object) {
    return(slot(object, "name"))
})

#' @export
#' @rdname fds-methods
setReplaceMethod("name", "FraserDataSet", function(object, value) {
    slot(object, "name") <- value
    validObject(object)
    return(object)
})

#' @export
#' @rdname fds-methods
setMethod("workingDir", "FraserDataSet", function(object) {
    return(slot(object, "workingDir"))
})

#' @export
#' @rdname fds-methods
setReplaceMethod("workingDir", "FraserDataSet", function(object, value) {
    slot(object, "workingDir") <- value
    validObject(object)
    return(object)
})

#' @export
#' @rdname fds-methods
setMethod("strandSpecific", "FraserDataSet", function(object) {
    return(slot(object, "strandSpecific"))
})

#' @export
#' @rdname fds-methods
setReplaceMethod("strandSpecific", "FraserDataSet", function(object, value) {
    if(is.logical(value)){
        value <- as.integer(value)
    }
    if(is.character(value)){
        value <- switch(tolower(value),
                        'no' = 0L,
                        'unstranded' = 0L,
                        'yes' = 1L,
                        'stranded' = 1L,
                        'reverse' = 2L,
                        -1L)
    }
    slot(object, "strandSpecific") <- value
    validObject(object)
    return(object)
})


#' @export
#' @rdname fds-methods
setMethod("pairedEnd", "FraserDataSet", function(object) {
    if(!("pairedEnd" %in% colnames(colData(object)))){
        return(rep(FALSE, length(samples(object))))
    }
    pairedEnd <- colData(object)[,"pairedEnd"]
    return(pairedEnd)
})

#' @export
#' @rdname fds-methods
setReplaceMethod("pairedEnd", "FraserDataSet", function(object, value) {
    if(any(!is.logical(value))){
        warning("Value need to be logical. Converting to logical.")
    }
    colData(object)[,"pairedEnd"] <- as.logical(value)
    validObject(object)
    return(object)
})


#' @export
#' @rdname fds-methods
setMethod("scanBamParam", "FraserDataSet", function(object) {
    return(slot(object, "bamParam"))
})

#' @export
#' @rdname fds-methods
setReplaceMethod("scanBamParam", "FraserDataSet", function(object, value) {
    slot(object, "bamParam") <- value
    validObject(object)
    return(object)
})


#' @export
#' @rdname fds-methods
setMethod("nonSplicedReads", "FraserDataSet", function(object){
    return(slot(object, "nonSplicedReads"))
})

#' @export
#' @rdname fds-methods
setReplaceMethod("nonSplicedReads", "FraserDataSet", function(object, value){
    slot(object, "nonSplicedReads") <- value
    validObject(object)
    return(object)
})

#'
#' Subsetting by indices for junctions
#'
#' Providing subsetting by indices through the single-bracket operator
#'
#' @param x A \code{FraserDataSet} object
#' @param i A integer vector to subset the rows/ranges
#' @param j A integer vector to subset the columns/samples
#' @param by a character (j or ss) definig if we subset by
#'             junctions or splice sites
#' @return A subsetted \code{FraserDataSet} object
#' @examples
#'     fds <- createTestFraserDataSet()
#'     fds[1:10,2:3]
#'     fds[,samples(fds) %in% c("sample1", "sample2")]
#'     fds[1:10,by="ss"]
#'
#' @rdname subset
subset.FRASER <- function(x, i, j, by=c("j", "ss")){
    if(length(by) == 1){
        by <- whichReadType(x, by)
    }
    by <- match.arg(by)

    if(missing(i) && missing(j)){
        return(x)
    }
    if(missing(i)){
        i <- TRUE
    }
    if(missing(j)){
        j <- TRUE
    }

    if(is(i, "RangedSummarizedExperiment") | is(i, "GRanges")){
        if(by=="ss"){
            i <- unique(to(findOverlaps(i, nonSplicedReads(x))))
        } else {
            i <- unique(to(findOverlaps(i, x)))
        }
    }

    nsrObj <- nonSplicedReads(x)
    if(length(nsrObj) == 0 & by=="ss"){
        stop("Cannot subset by splice sites, if you not counted them!")
    }
    ssIdx <- NULL
    if(by == "ss"){
        ssIdx <- unlist(rowData(nonSplicedReads(x)[i])["spliceSiteID"])
        i <- as.vector(unlist(rowData(x, typ="psi3")["startID"]) %in% ssIdx |
            unlist(rowData(x, typ="psi3")["endID"]) %in% ssIdx)
    }

    # first subset non spliced reads if needed
    if(length(nsrObj) > 0){
        if(is.null(ssIdx)){
            # get selected splice site IDs
            ssIdx <- unique(unlist(
                rowData(x, type="j")[i, c("startID", "endID")]
            ))
        }

        # get the selection vector
        idxNSR <- rowData(x, type="ss")[['spliceSiteID']] %in% ssIdx

        # subset it
        nsrObj <- nsrObj[idxNSR,j]
    }

    # subset the inheritate SE object
    if(length(x) == 0){
        i <- NULL
    }
    subX <- as(as(x, "RangedSummarizedExperiment")[i,j], "FraserDataSet")

    # create new FraserDataSet object
    newx <- new("FraserDataSet",
            subX,
            name            = name(x),
            bamParam        = scanBamParam(x),
            strandSpecific  = strandSpecific(x),
            workingDir      = workingDir(x),
            nonSplicedReads = nsrObj
    )
    validObject(newx)
    return(newx)
}
#' @rdname subset
#' @export
setMethod("[", c("FraserDataSet", "ANY", "ANY"), subset.FRASER)


#'
#' Returns the assayNames of FRASER
#' 
#' @param x FraserDataSet
#' 
#' @return Character vector
#' @export
setMethod("assayNames", "FraserDataSet", function(x) {
    return(c(
        assayNames(asSE(x)),
        assayNames(nonSplicedReads(x))
    ))
})


FraserDataSet.assays.replace_pre <-
            function(x, withDimnames=TRUE, HDF5=TRUE, type=NULL, ..., value){
    if(any(names(value) == "")) stop("Name of an assay can not be empty!")
    if(any(duplicated(names(value)))) stop("Assay names need to be unique!")
    if(is.null(type)){
        type <- names(value)
    }
    # make sure all slots are HDF5
    if(isTRUE(HDF5)){
        for(i in seq_along(value)){
            if(!any(class(value[[i]]) %in% c("HDF5Matrix", "DelayedMatrix")) ||
                tryCatch(!is.character(path(value[[i]])), 
                            error=function(e){TRUE})){
                
                aname <- names(value)[i]
                h5obj <- saveAsHDF5(x, aname, object=value[[i]])
                value[[i]] <- h5obj
            }
        }
    }

    # first replace all existing slots
    nj <- names(value) %in% assayNames(asSE(x))
    ns <- names(value) %in% assayNames(nonSplicedReads(x))
    jslots <- value[nj]
    sslots <- value[ns]

    # add new slots if there are some
    if(sum(!(nj | ns)) > 0){
        type <- vapply(type, checkReadType, fds=x, FUN.VALUE="")
        jslots <- c(jslots, value[!(nj | ns) & type=="j"])
        sslots <- c(sslots, value[!(nj | ns) & type=="ss"])
    }

    # assign new assays only for non split
    # and do the rest in the specific functions
    assays(nonSplicedReads(x), withDimnames=withDimnames, ...) <- sslots
    
    return(list(x=x, value=jslots))
}

FraserDataSet.assays.replace_r40 <- function(x, withDimnames=TRUE, HDF5=TRUE, 
                type=NULL, ..., value){
    ans <- FraserDataSet.assays.replace_pre(x, withDimnames=withDimnames, 
            HDF5=HDF5, type=type, ..., value=value)
    
    # retrieve adapted objects and set final assays on main SE object
    value <- ans[['value']]
    x <- ans[['x']]
    x <- callNextMethod()
    
    # validate and return
    validObject(x)
    return(x)
}

FraserDataSet.assays.replace_r36 <- function(x, ..., HDF5=TRUE, type=NULL, 
                withDimnames=TRUE, value){
    ans <- FraserDataSet.assays.replace_pre(x, withDimnames=withDimnames, 
            HDF5=HDF5, type=type, ..., value=value)
    
    # retrieve adapted objects and set final assays on main SE object
    value <- ans[['value']]
    x <- ans[['x']]
    x <- callNextMethod()
    
    # validate and return
    validObject(x)
    return(x)
}

FraserDataSet.assays.set_r40 <- function(x, withDimnames=TRUE, ...){
    return(c(
        assays(asSE(x), withDimnames=withDimnames, ...),
        assays(nonSplicedReads(x), withDimnames=withDimnames, ...)
    ))
}

FraserDataSet.assays.set_r36 <- function(x, ..., withDimnames=TRUE){
    FraserDataSet.assays.set_r40(x=x, ..., withDimnames=withDimnames)
}

if(compareVersion(package.version("SummarizedExperiment"), "1.17-2") < 0){
    FraserDataSet.assays.set <- FraserDataSet.assays.set_r36
    FraserDataSet.assays.replace <- FraserDataSet.assays.replace_r36
} else {
    FraserDataSet.assays.set <- FraserDataSet.assays.set_r40
    FraserDataSet.assays.replace <- FraserDataSet.assays.replace_r40
}

#'
#' Returns the assay for the given name/index of the FraserDataSet
#'
#' @param x FraserDataSet
#' @param ... Parameters passed on to SummarizedExperiment::assays() 
#' @param withDimnames Passed on to SummarizedExperiment::assays() 
#' @param HDF5 Logical value indicating whether the assay should be stored as 
#' a HDF5 file. 
#' @param type The psi type. 
#' @param value The new value to which the assay should be set. 
#'
#' @return (Delayed) matrix.
#' @export
setMethod("assays", "FraserDataSet", FraserDataSet.assays.set)

#' @rdname assays-FraserDataSet-method
#' @export
setReplaceMethod("assays", c("FraserDataSet", "SimpleList"),
        FraserDataSet.assays.replace)

#' @rdname assays-FraserDataSet-method
#' @export
setReplaceMethod("assays", c("FraserDataSet", "list"),
        FraserDataSet.assays.replace)

#' @rdname assays-FraserDataSet-method
#' @export
setReplaceMethod("assays", c("FraserDataSet", "DelayedMatrix"),
        FraserDataSet.assays.replace)


#'
#' retrieve the length of the object (aka number of junctions)
#' 
#' @param x FraserDataSet
#'
#' @return Length of the object.
#' @export
setMethod("length", "FraserDataSet", function(x) callNextMethod())

#' @rdname fds-methods 
#' @export
FRASER.mcols.get <- function(x, type=NULL, ...){
    type <- checkReadType(x, type)
    if(type=="j"){
        return(mcols(asSE(x), ...))
    }
    mcols(nonSplicedReads(x), ...)
}
FRASER.mcols.replace <- function(x, type=NULL, ..., value){
    type <- checkReadType(x, type)
    if(type=="j") {
        return(callNextMethod())
        # se <- asSE(x)
        # mcols(se, ...) <- value
        # return(asFDS(x))
    }
    mcols(nonSplicedReads(x), ...) <- value
    return(x)
}
setMethod("mcols", "FraserDataSet", FRASER.mcols.get)
setReplaceMethod("mcols", "FraserDataSet", FRASER.mcols.replace)

#' @rdname fds-methods
#' @export
FRASER.rowRanges.get <- function(x, type=NULL, ...){
    type <- checkReadType(x, type)
    if(type=="j")  return(callNextMethod())
    if(type=="ss") return(rowRanges(nonSplicedReads(x), ...))
}
FRASER.rowRanges.replace <- function(x, type=NULL, ..., value){
    type <- checkReadType(x, type)
    if(type=="j") return(callNextMethod())
    rowRanges(nonSplicedReads(x), ...) <- value
    return(x)
}

setMethod("rowRanges", "FraserDataSet", FRASER.rowRanges.get)
setReplaceMethod("rowRanges", "FraserDataSet", FRASER.rowRanges.replace)

#'
#' Getter/setter for count data
#'
#' @param fds,object FraserDataSet
#' @param type The psi type.
#' @param side "ofInterest" for junction counts, "other" for sum of counts of 
#' all other junctions at the same donor site (psi5) or acceptor site (psi3), 
#' respectively. 
#' @param value An integer matrix containing the counts.
#' @param ... Further parameters that are passed to assays(object,...)
#' 
#' @return FraserDataSet
#' @rdname counts
#' @examples 
#'  fds <- createTestFraserDataSet()
#'  
#'  counts(fds, type="psi5", side="ofInterest")
#'  counts(fds, type="psi5", side="other")
#'  head(K(fds, type="psi3"))
#'  head(N(fds, type="theta"))
#'  
setMethod("counts", "FraserDataSet", function(object, type=NULL,
            side=c("ofInterest", "otherSide")){
    side <- match.arg(side)
    if(side=="ofInterest"){
        type <- checkReadType(object, type)
        aname <- paste0("rawCounts", toupper(type))
        if(!aname %in% assayNames(object)){
            stop(paste0("Missing rawCounts for '", type, "'. ",
                    "Please count your data first. And then try again."))
        }
        return(assays(object)[[aname]])
    }

    # extract psi value from type
    type <- whichPSIType(type)
    if(length(type) == 0 | length(type) > 1){
        stop(paste0("Please provide a correct psi type: psi5, psi3, or ",
                    "theta. Not the given one: '", type, "'."))
    }
    aname <- paste0("rawOtherCounts_", type)
    if(!aname %in% assayNames(object)){
        stop(paste0("Missing rawOtherCounts for type '", type, "'.",
                " Please calculate PSIValues first. And then try again."))
    }
    return(assays(object)[[aname]])
})

#'
#' setter for count data
#' 
#' @rdname counts
setReplaceMethod("counts", "FraserDataSet", function(object, type=NULL,
                    side=c("ofInterest", "otherSide"), ..., value){
    side <- match.arg(side)

    if(side=="ofInterest"){
        type <- checkReadType(object, type)
        aname <- paste0("rawCounts", toupper(type))
    } else {
        type <- whichPSIType(type)
        aname <- paste0("rawOtherCounts_", type)
    }
    assays(object, ...)[[aname]] <- value
    validObject(value)
    return(object)
})


setAs("DelayedMatrix", "data.table", function(from){ 
    as.data.table(from) })
setAs("DataFrame",     "data.table", function(from){ 
    as.data.table(from) })
setAs("DelayedMatrix", "matrix", function(from){ 
    as.matrix(as(from, "data.table")) })
setAs("DataFrame", "matrix", function(from){
    as.matrix(as(from, "data.table")) })

#'
#' retrieve a single sample result object
#' @noRd
resultsSingleSample <- function(sampleID, gr, pvals, padjs, zscores, psivals,
                rawCts, rawTotalCts, deltaPsiVals, muPsi, psiType, fdrCut,
                zscoreCut, dPsiCut, rowMeansK, rowMeansN, minCount,
                additionalColumns){

    zscore  <- zscores[,sampleID]
    dpsi    <- deltaPsiVals[,sampleID]
    pval    <- pvals[,sampleID]
    padj    <- padjs[,sampleID]

    goodCut <- !logical(length(zscore))
    if(!is.na(zscoreCut)){
        goodCut <- goodCut & na2default(abs(zscore) >= zscoreCut, TRUE)
    }
    if(!is.na(dPsiCut)){
        goodCut <- goodCut & na2default(abs(dpsi) >= dPsiCut, TRUE)
    }
    if(!is.na(fdrCut)){
        goodCut <- goodCut & na2false(padj <= fdrCut)
    }
    if(!is.na(minCount)){
        goodCut <- goodCut & rawTotalCts[,sampleID] >= minCount
    }

    ans <- granges(gr[goodCut])

    if(!any(goodCut)){
        return(ans)
    }

    if(!"hgnc_symbol" %in% colnames(mcols(gr))){
        mcols(gr)$hgnc_symbol <- NA_character_
    }

    # extract data
    mcols(ans)$sampleID        <- Rle(sampleID)
    if("hgnc_symbol" %in% colnames(mcols(gr))){
        mcols(ans)$hgncSymbol <- Rle(mcols(gr[goodCut])$hgnc_symbol)
    }
    if("other_hgnc_symbol" %in% colnames(mcols(gr))){
        mcols(ans)$addHgncSymbols <- Rle(mcols(gr[goodCut])$other_hgnc_symbol)
    }
    mcols(ans)$type            <- Rle(psiType)
    mcols(ans)$pValue          <- signif(pval[goodCut], 5)
    mcols(ans)$padjust         <- signif(padj[goodCut], 5)
    mcols(ans)$zScore          <- Rle(round(zscore[goodCut], 2))
    mcols(ans)$psiValue        <- Rle(round(psivals[goodCut,sampleID], 2))
    mcols(ans)$deltaPsi        <- Rle(round(dpsi[goodCut], 2))
    mcols(ans)$meanCounts      <- Rle(round(rowMeansK[goodCut], 2))
    mcols(ans)$meanTotalCounts <- Rle(round(rowMeansN[goodCut], 2))
    mcols(ans)$counts          <- Rle(rawCts[goodCut, sampleID])
    mcols(ans)$totalCounts     <- Rle(rawTotalCts[goodCut, sampleID])
    
    if(!is.null(additionalColumns)){
        for(column in additionalColumns){
            mcols(ans)[,column] <- Rle(mcols(gr[goodCut])[,column])
        }
    }

    return(ans[order(mcols(ans)$pValue)])
}

FRASER.results <- function(object, sampleIDs, fdrCutoff, zscoreCutoff, 
                    dPsiCutoff, psiType, BPPARAM=bpparam(), maxCols=20, 
                    minCount, additionalColumns=NULL){

    # check input
    checkNaAndRange(fdrCutoff,    min=0, max=1,   scalar=TRUE, na.ok=TRUE)
    checkNaAndRange(dPsiCutoff,   min=0, max=1,   scalar=TRUE, na.ok=TRUE)
    checkNaAndRange(zscoreCutoff, min=0, max=100, scalar=TRUE, na.ok=TRUE)
    checkNaAndRange(minCount,     min=0, max=Inf, scalar=TRUE, na.ok=TRUE)

    stopifnot(is(object, "FraserDataSet"))
    stopifnot(all(sampleIDs %in% samples(object)))
    
    if("annotatedJunction" %in% colnames(mcols(object, type="j")) && 
            !("annotatedJunction" %in% additionalColumns)){
        additionalColumns <- c(additionalColumns, "annotatedJunction")
    }

    resultsls <- bplapply(psiType, BPPARAM=BPPARAM, function(type){
        message(date(), ": Collecting results for: ", type)
        currentType(object) <- type
        gr <- rowRanges(object, type=type)

        # first get row means
        rowMeansK <- rowMeans(K(object, type=type))
        rowMeansN <- rowMeans(N(object, type=type))

        # then iterate by chunk
        chunkCols <- getMaxChunks2Read(fds=object, assayName=type, max=maxCols)
        sampleChunks <- getSamplesByChunk(fds=object, sampleIDs=sampleIDs,
                chunkSize=chunkCols)

        ans <- lapply(seq_along(sampleChunks), function(idx){
            message(date(), ": Process chunk: ", idx, " for: ", type)
            sc <- sampleChunks[[idx]]
            tmp_x <- object[,sc]

            # extract values
            rawCts       <- as.matrix(K(tmp_x))
            rawTotalCts  <- as.matrix(N(tmp_x))
            pvals        <- as.matrix(pVals(tmp_x))
            padjs        <- as.matrix(padjVals(tmp_x))
            zscores      <- as.matrix(zScores(tmp_x))
            psivals      <- as.matrix(assay(tmp_x, type))
            muPsi        <- as.matrix(predictedMeans(tmp_x))
            psivals_pc   <- (rawCts + pseudocount()) /
                                (rawTotalCts + 2*pseudocount())
            deltaPsiVals <- psivals_pc - muPsi

            if(length(sc) == 1){
                colnames(pvals) <- sc
                colnames(padjs) <- sc
                colnames(zscores) <- sc
                colnames(deltaPsiVals) <- sc
            }
            
            # create result table
            sampleRes <- lapply(sc,
                    resultsSingleSample, gr=gr, pvals=pvals, padjs=padjs,
                    zscores=zscores, psiType=type, psivals=psivals,
                    deltaPsiVals=deltaPsiVals, muPsi=muPsi, rawCts=rawCts,
                    rawTotalCts=rawTotalCts, fdrCut=fdrCutoff,
                    zscoreCut=zscoreCutoff, dPsiCut=dPsiCutoff,
                    rowMeansK=rowMeansK, rowMeansN=rowMeansN, 
                    minCount=minCount, additionalColumns=additionalColumns)

            # return combined result
            return(unlist(GRangesList(sampleRes)))
        })

        unlist(GRangesList(ans))
    })

    # merge results
    ans <- unlist(GRangesList(resultsls))

    # sort it if existing
    if(length(ans) > 0){
        ans <- ans[order(ans$pValue)]
    }

    # return only the results
    return(ans)
}

#'
#' Extracting results and aberrant splicing events
#'
#' The result function extracts the results from the given analysis object
#' based on the given options and cutoffs. The aberrant function extracts 
#' aberrant splicing events based on the given cutoffs.
#'
#' @param object A \code{\link{FraserDataSet}} object
#' @param sampleIDs A vector of sample IDs for which results should be 
#' retrieved
#' @param padjCutoff The FDR cutoff to be applied or NA if not requested.
#' @param zScoreCutoff The z-score cutoff to be applied or NA if not requested.
#' @param deltaPsiCutoff The cutoff on delta psi or NA if not requested.
#' @param minCount The minimum count value of the total coverage of an intron 
#' to be considered as significant.
#' result
#' @param psiType The psi types for which the results should be retrieved.
#' @param additionalColumns Character vector containing the names of additional 
#' columns from mcols(fds) that should appear in the result table 
#' (e.g. ensembl_gene_id). Default is \code{NULL}, so no additional columns 
#' are included. 
#' @param BPPARAM The BiocParallel parameter.
#' @param res Result as created with \code{results()}
#' @param geneColumn The name of the column in \code{mcols(res)} that contains 
#'     the gene symbols.   
#' @param method The p.adjust method that is being used to adjust p values per
#'     sample.
#' @param type Splicing type (psi5, psi3 or theta)
#' @param by By default \code{none} which means no grouping. But if 
#'              \code{sample} or \code{feature} is specified the sum by 
#'              sample or feature is returned
#' @param aggregate If TRUE the returned object is based on the grouped 
#'              features
#' @param ... Further arguments can be passed to the method. If "zscores", 
#'              "padjVals" or "dPsi" is given, the values of those arguments
#'              are used to define the aberrant events.
#'
#' @return For \code{results}: GRanges object containing significant results.
#'     For \code{aberrant}: Either a of logical values of size 
#'     introns/genes x samples if "by" is NA or a vector with the 
#'     number of aberrant events per sample or feature depending on 
#'     the vaule of "by"
#' 
#' @rdname results
#' @examples
#' # get data, fit and compute p-values and z-scores
#' fds <- createTestFraserDataSet()
#' 
#' # extract results: for this example dataset, z score cutoff of 2 is used to
#' # get at least one result and show the output
#' res <- results(fds, padjCutoff=NA, zScoreCutoff=3, deltaPsiCutoff=0.05)
#' res
#' 
#' # aggregate the results by genes (gene symbols need to be annotated first 
#' # using annotateRanges() function)
#' resultsByGenes(res)
#'
#' # get aberrant events per sample: on the example data, nothing is aberrant
#' # based on the adjusted p-value
#' aberrant(fds, type="psi5", by="sample")
#' 
#' # get aberrant events per gene (first annotate gene symbols)
#' fds <- annotateRangesWithTxDb(fds)
#' aberrant(fds, type="psi5", by="feature", zScoreCutoff=2, padjCutoff=NA,
#'         aggregate=TRUE)
#'         
#' # find aberrant junctions/splice sites
#' aberrant(fds, type="psi5")
#' @export
setMethod("results", "FraserDataSet", function(object, 
                    sampleIDs=samples(object), padjCutoff=0.05,
                    zScoreCutoff=NA, deltaPsiCutoff=0.3,
                    minCount=5, psiType=c("psi3", "psi5", "theta"),
                    additionalColumns=NULL, BPPARAM=bpparam(), ...){
    FRASER.results(object=object, sampleIDs=sampleIDs, fdrCutoff=padjCutoff,
            zscoreCutoff=zScoreCutoff, dPsiCutoff=deltaPsiCutoff,
            minCount=minCount, psiType=match.arg(psiType, several.ok=TRUE),
            additionalColumns=additionalColumns, BPPARAM=BPPARAM)
})

#' @rdname results
#' @export
resultsByGenes <- function(res, geneColumn="hgncSymbol", method="BY"){
    # sort by pvalue
    res <- res[order(res$pValue)]

    # extract subset
    if(is(res, "GRanges")){
        ans <- as.data.table(mcols(res)[,c(geneColumn, "pValue", "sampleID")])
        colnames(ans) <- c("features", "pval", "sampleID")
    } else {
        ans <- featureNames <- res[,.(
                features=get(geneColumn), pval=pValue, sampleID=sampleID)]
    }

    # remove NAs
    naIdx <- ans[,is.na(features)]
    ansNoNA <- ans[!is.na(features)]

    # compute pvalues by gene
    ansNoNA[,pByFeature:=min(p.adjust(pval, method="holm")),
            by="sampleID,features"]

    # subset to lowest pvalue by gene
    dupIdx <- duplicated(ansNoNA[,.(features,sampleID)])
    ansGenes <- ansNoNA[!dupIdx]

    # compute FDR
    ansGenes[,fdrByFeature:=p.adjust(pByFeature, method=method),
            by="sampleID"]

    # get final result table
    finalAns <- res[!naIdx][!dupIdx]
    finalAns$pValueGene  <- ansGenes$pByFeature
    finalAns$padjustGene <- ansGenes$fdrByFeature
    finalAns
}

#'
#' Mapping of chromosome names
#'
#' @param fds FraserDataSet
#' @param style The style of the chromosome names.
#' @param ... Further parameters. For mapSeqLevels: further parameters 
#'     passed to GenomeInfoDb::mapSeqlevels().
#' 
#' @rdname fds-methods
#' @export
mapSeqlevels <- function(fds, style="UCSC", ...){

    mappings <- na.omit(GenomeInfoDb::mapSeqlevels(seqlevels(fds), style, ...))
    # fix missing names() when fds has only a single chromosome
    if(is.null(names(mappings))){ 
        names(mappings) <- seqlevels(fds)
    }

    if(length(mappings) != length(seqlevels(fds))){
        message(date(), ": Drop non standard chromosomes for compatibility.")
        fds <- keepStandardChromosomes(fds)
        nonSplicedReads(fds) <- keepStandardChromosomes(nonSplicedReads(fds))
        validObject(fds)
    }
    fds <- fds[as.vector(seqnames(fds)) %in% names(mappings)]

    seqlevels(fds) <- as.vector(mappings)
    seqlevels(nonSplicedReads(fds)) <- as.vector(mappings)

    return(fds)
}


aberrant.FRASER <- function(object, type=currentType(object), padjCutoff=0.05,
                    deltaPsiCutoff=0.3, zScoreCutoff=NA, minCount=5,
                    by=c("none", "sample", "feature"), aggregate=FALSE, ...){

    checkNaAndRange(zScoreCutoff,   min=0, max=Inf, na.ok=TRUE)
    checkNaAndRange(padjCutoff,     min=0, max=1,   na.ok=TRUE)
    checkNaAndRange(deltaPsiCutoff, min=0, max=1,   na.ok=TRUE)
    by <- match.arg(by)

    dots <- list(...)
    if("n" %in% names(dots)){
        n <- dots[['n']]
    } else {
        n <- N(object, type=type)
    }
    if("zscores" %in% names(dots)){
        zscores <- dots[['zscores']]
    } else {
        zscores <- zScores(object, type=type)
    }
    if("padjVals" %in% names(dots)){
        padj <- dots[['padjVals']]
    } else {
        padj <- padjVals(object, type=type)
    }
    if("dPsi" %in% names(dots)){
        dpsi <- dots[['dPsi']]
    } else {
        dpsi <- deltaPsiValue(object, type=type)
    } 
    
    
    # create cutoff matrix
    goodCutoff <- matrix(TRUE, nrow=nrow(zscores), ncol=ncol(zscores),
            dimnames=dimnames(zscores))
    if("hgnc_symbol" %in% colnames(mcols(object, type=type)) &
                nrow(mcols(object, type=type)) == nrow(goodCutoff)){
        rownames(goodCutoff) <- mcols(object, type=type)[,"hgnc_symbol"]
    } else if(isTRUE(aggregate)){
        stop("Please provide hgnc symbols to compute gene p values!")
    }
    
    # check each cutoff if in use (not NA)
    if(!is.na(minCount)){
        goodCutoff <- goodCutoff & as.matrix(n >= minCount)
    }
    if(!is.na(zScoreCutoff)){
        goodCutoff <- goodCutoff & as.matrix(abs(zscores) >= zScoreCutoff)
    }
    if(!is.na(deltaPsiCutoff)){
        goodCutoff <- goodCutoff & as.matrix(abs(dpsi) >= deltaPsiCutoff)
    }
    if(!is.na(padjCutoff)){
        goodCutoff <- goodCutoff & as.matrix(padj <= padjCutoff)
    }
    goodCutoff[is.na(goodCutoff)] <- FALSE
    
    # check if we should go for aggregation
    # TODO to speed it up we only use any hit within a feature
    # but should do a holm's + BY correction per gene and genome wide
    if(isTRUE(aggregate)){
        goodCutoff <- as.matrix(data.table(goodCutoff, keep.rownames=TRUE)[,
                as.data.table(t(colAnys(as.matrix(.SD)))), by=rn][,-1])
        rownames(goodCutoff) <- unique(mcols(object, type=type)[,"hgnc_symbol"])
        colnames(goodCutoff) <- colnames(zscores)
    }
    
    # return results
    if(by == "feature"){
        return(rowSums(goodCutoff))
    }
    if(by == "sample"){
        return(colSums(goodCutoff))
    }
    return(goodCutoff)
}

#' @rdname results
#' @export
setMethod("aberrant", "FraserDataSet", aberrant.FRASER)
