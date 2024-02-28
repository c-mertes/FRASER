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
    if(!is.null(colnames(object))){
        return(colnames(object))
    }
    return(as.character(colData(object)[,"sampleID"]))
})

#' @rdname fds-methods
#' @export
setReplaceMethod("samples", "FraserDataSet", function(object, value) {
    colData(object)[,"sampleID"] <- as.character(value)
    rownames(colData(object)) <- colData(object)[,"sampleID"]
    colnames(object) <- as.character(value)
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
    if(!"strand" %in% colnames(colData(object))){
        warning("Strand is not specified. Please set the used RNA-seq",
                " protocol by using 'strandSpecific(object) <- c(...)'.",
                "\n\nWe assume as default a non stranded protocol.")
        return(0)
    }
    return(colData(object)$strand)
})

#' @export
#' @rdname fds-methods
setReplaceMethod("strandSpecific", "FraserDataSet", function(object, value) {
    if(is.logical(value)){
        value <- as.integer(value)
    }
    if(is.character(value)){
        val_chr <- tolower(value)
        value <- sapply(val_chr, switch, 'no' = 0L, 'unstranded' = 0L, 
                'yes' = 1L, 'stranded' = 1L, 'reverse' = 2L, -1L)
        if(any(value < 0)){
          stop("Please specify correct strandness of the samples.", 
                  " It needs to be one of: 'no', 'unstranded', 'yes',",
                  " 'stranded' or 'reverse'.", "\n\nIt was specified: ", 
                  paste(collapse = ", ", val_chr[value < 0]))
        }
    }
    colData(object)$strand <- as.integer(value)
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
#' @param by a character (j or ss) defining if we subset by
#'             junctions or splice sites
#' @param ... Parameters currently not used or passed on
#' @param drop No dimension reduction is done. And the \code{drop}
#'             parameter is currently not used at all.
#' @return A subsetted \code{FraserDataSet} object
#' @examples
#'     fds <- createTestFraserDataSet()
#'     fds[1:10,2:3]
#'     fds[,samples(fds) %in% c("sample1", "sample2")]
#'     fds[1:10,by="ss"]
#'
#' @rdname subset
subset.FRASER <- function(x, i, j, by=c("j", "ss"), ..., drop=FALSE){
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
        nsrObj <- nsrObj[idxNSR,j,drop=FALSE]
    }

    # subset the inheritate SE object
    if(length(x) == 0){
        i <- NULL
    }
    subX <- as(as(x, "RangedSummarizedExperiment")[i,j,drop=FALSE], "FraserDataSet")

    # create new FraserDataSet object
    newx <- new("FraserDataSet",
            subX,
            name            = name(x),
            bamParam        = scanBamParam(x),
            workingDir      = workingDir(x),
            nonSplicedReads = nsrObj
    )
    validObject(newx)
    return(newx)
}
#' @rdname subset
#' @export
setMethod("[", c("FraserDataSet", "ANY", "ANY", drop="ANY"), subset.FRASER)


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
#'  counts(fds, side="ofInterest")
#'  counts(fds, type="jaccard", side="other")
#'  head(K(fds))
#'  head(K(fds, type="psi5"))
#'  head(K(fds, type="psi3"))
#'  head(N(fds, type="theta"))
#'  
setMethod("counts", "FraserDataSet", function(object, type=currentType(object),
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
        stop(paste0("Please provide a correct psi type: psi5, psi3, ",
                    "theta or jaccard. Not the given one: '", 
                    type, "'."))
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
setReplaceMethod("counts", "FraserDataSet", function(object, 
                    type=currentType(object),
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

#'
#' retrieve a single sample result object
#' @noRd
resultsSingleSample <- function(sampleID, gr, pvals, padjs, 
                                psivals, rawCts, rawTotalCts, rawNonsplitCts, 
                                rawNsProportion, nsProportion_99quantile,
                                deltaPsiVals, psiType, rowMeansK, rowMeansN, 
                                aberrant, aggregate, rho, 
                                pvalsGene=NULL, padjsGene=NULL, 
                                aberrantGene, additionalColumns,
                                geneColumn="hgnc_symbol"){
    mcols(gr)$idx <- seq_along(gr)
    # if gene level results, find the most aberrant junction per gene first
    if(isTRUE(aggregate)){
        goodGenes <- rownames(aberrantGene)[aberrantGene[,sampleID] & 
                                        !is.na(aberrantGene[,sampleID])]
        geneJunctions <- findJunctionsForAberrantGenes(gr=gr, 
                            aberrantGenes=goodGenes, 
                            pvals=pvals[,sampleID], 
                            dpsi=deltaPsiVals[,sampleID],
                            geneColumn=geneColumn,
                            aberrantJunctions=aberrant[,sampleID])    
        goodCut <- rep(FALSE, nrow(pvals))
        goodCut[geneJunctions] <- TRUE
    } else{
        goodCut <- aberrant[,sampleID]
    }
    
    ans <- granges(gr[goodCut])
    
    if(!any(goodCut)){
        return(ans)
    }
    mcols(ans)$idx <- mcols(gr)$idx[goodCut]
    
    if(!geneColumn %in% colnames(mcols(gr))){
        mcols(gr)[,geneColumn] <- NA_character_
    }
    
    # extract data
    mcols(ans)$sampleID        <- Rle(sampleID)
    if("hgnc_symbol" %in% colnames(mcols(gr))){
        mcols(ans)$hgncSymbol <- Rle(mcols(gr[goodCut])[,geneColumn])
    }
    
    mcols(ans)$type            <- Rle(psiType)
    mcols(ans)$pValue          <- signif(pvals[goodCut,sampleID], 5)
    mcols(ans)$padjust         <- signif(padjs[goodCut,sampleID], 5)
    mcols(ans)$psiValue        <- Rle(round(psivals[goodCut,sampleID], 2))
    mcols(ans)$deltaPsi        <- round(deltaPsiVals[goodCut,sampleID], 2)
    mcols(ans)$counts          <- Rle(rawCts[goodCut, sampleID])
    mcols(ans)$totalCounts     <- Rle(rawTotalCts[goodCut, sampleID])
    mcols(ans)$meanCounts      <- Rle(round(rowMeansK[goodCut], 2))
    mcols(ans)$meanTotalCounts <- Rle(round(rowMeansN[goodCut], 2))
    
    if(psiType == "jaccard"){
        mcols(ans)$nonsplitCounts <- 
            Rle(round(rawNonsplitCts[goodCut, sampleID], 2))
        mcols(ans)$nonsplitProportion <- 
            Rle(round(rawNsProportion[goodCut, sampleID], 2))
        mcols(ans)$nonsplitProportion_99quantile  <- 
            Rle(round(nsProportion_99quantile[goodCut], 2))
    }
    
    if(!is.null(additionalColumns)){
        for(column in additionalColumns){
            mcols(ans)[,column] <- Rle(mcols(gr[goodCut])[,column])
        }
    }
    
    if(isTRUE(aggregate)){
        # report junction more than once if it is significant for several genes
        nrGenesPerJunction <- table(geneJunctions)
        ans <- rep(ans, nrGenesPerJunction[as.character(mcols(ans)$idx)])
        mcols(ans)$hgncSymbol <- 
            as.data.table(ans)[, names(geneJunctions)[geneJunctions == idx], 
                                by = eval(colnames(mcols(ans)))][,V1]
        
        # add gene level pvalue    
        mcols(ans)$pValueGene  <- 
            signif(pvalsGene[mcols(ans)$hgncSymbol,sampleID], 5)
        mcols(ans)$padjustGene <- 
            signif(padjsGene[mcols(ans)$hgncSymbol,sampleID], 5)
        mcols(ans)$hgncSymbol <- Rle(mcols(ans)$hgncSymbol)
    }
    
    # remove helper column
    mcols(ans)$idx <- NULL
    
    
    return(ans[order(mcols(ans)$pValue, -abs(mcols(ans)$deltaPsi))])
}

FRASER.results <- function(object, sampleIDs, fdrCutoff, 
                            dPsiCutoff, minCount, rhoCutoff, psiType, 
                            maxCols=20, aggregate=FALSE, collapse=FALSE,
                            geneColumn="hgnc_symbol", BPPARAM=bpparam(), 
                            subsetName=NULL, all=all, additionalColumns=NULL){
    
    stopifnot(is(object, "FraserDataSet"))
    stopifnot(all(sampleIDs %in% samples(object)))
    
    if("annotatedJunction" %in% colnames(mcols(object, type="j")) && 
            !("annotatedJunction" %in% additionalColumns)){
        additionalColumns <- c(additionalColumns, "annotatedJunction")
    }
    
    # only extract results for requested psiTypes if pvals exist for them
    stopifnot(all(psiType %in% psiTypes))
    if(is.na(rhoCutoff)){
        rhoCutoff <- 1
    }
    pvalsAvailable <- checkPadjAvailableForFilters(object, type=psiType,
                                                   filters=list(rho=rhoCutoff), 
                                                   aggregate=aggregate,
                                                   subsetName=subsetName)
    psiType <- psiType[pvalsAvailable]
    if(all(isFALSE(pvalsAvailable))){
        stop("For the splice metric(s), pvalues are not yet computed. \n", 
             "Please compute them first by running the ",
             "calculatePadjValues function.")
    }
    
    resultsls <- bplapply(psiType, BPPARAM=BPPARAM, function(type){
        message(date(), ": Collecting results for: ", type,
                ifelse(is.null(subsetName), " (transcriptome-wide)",
                       paste0(" (", subsetName, ")")))
        currentType(object) <- type
        gr <- rowRanges(object, type=type)
        
        # first get row means
        rowMeansK <- rowMeans(K(object, type=type))
        rowMeansN <- rowMeans(N(object, type=type))
        
        # get proportion of nonsplitCounts among all counts (N) for each intron
        if(type == "jaccard"){
            rawNonsplitCts  <- as.matrix(assay(object, "rawCountsJnonsplit"))
            rawNsProportion <- rawNonsplitCts / as.matrix(N(object))
            nsProportion_99quantile <- 
                rowQuantiles(rawNsProportion, probs=0.99)
        } else{
            rawNonsplitCts  <- NULL
            rawNsProportion <- NULL
            nsProportion_99quantile <- NULL
        }
        
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
            pvals        <- as.matrix(pVals(tmp_x, 
                                            filters=list(rho=rhoCutoff)))
            padjs        <- as.matrix(padjVals(tmp_x, 
                                            subsetName=subsetName,
                                            filters=list(rho=rhoCutoff)))
            psivals      <- as.matrix(assay(tmp_x, type))
            muPsi        <- as.matrix(predictedMeans(tmp_x))
            psivals_pc   <- (rawCts + pseudocount()) /
                (rawTotalCts + 2*pseudocount())
            deltaPsiVals <- deltaPsiValue(tmp_x, type)
            rho          <- rho(tmp_x, type)
            aberrant     <- aberrant.FRASER(tmp_x, type=type, 
                                            padjCutoff=ifelse(isTRUE(aggregate), 
                                                        NA, fdrCutoff), 
                                            deltaPsiCutoff=dPsiCutoff, 
                                            minCount=minCount, 
                                            rhoCutoff=rhoCutoff,
                                            aggregate=FALSE,
                                            all=all,
                                            geneColumn=geneColumn,
                                            subsetName=subsetName)
            if(isTRUE(aggregate)){
                pvalsGene    <- as.matrix(pVals(tmp_x, level="gene", 
                                                filters=list(rho=rhoCutoff)))
                padjsGene    <- as.matrix(padjVals(tmp_x, level="gene",
                                                   subsetName=subsetName, 
                                                filters=list(rho=rhoCutoff)))
                aberrantGene <- aberrant.FRASER(tmp_x, type=type, 
                                                padjCutoff=fdrCutoff, 
                                                deltaPsiCutoff=dPsiCutoff, 
                                                minCount=minCount, 
                                                rhoCutoff=rhoCutoff,
                                                aggregate=TRUE,
                                                all=all,
                                                geneColumn=geneColumn,
                                                subsetName=subsetName)
            } else{
                pvalsGene    <- NULL
                padjsGene    <- NULL
                aberrantGene <- NULL
            }
            
            if(length(sc) == 1){
                colnames(pvals) <- sc
                colnames(padjs) <- sc
                colnames(deltaPsiVals) <- sc
            }
            
            # create result table
            sampleRes <- lapply(sc,
                        resultsSingleSample, gr=gr, pvals=pvals, 
                        padjs=padjs, psiType=type, 
                        psivals=psivals, deltaPsiVals=deltaPsiVals, 
                        rawCts=rawCts, rawTotalCts=rawTotalCts, 
                        rawNonsplitCts=rawNonsplitCts[,sc,drop=FALSE], 
                        rawNsProportion=rawNsProportion[,sc,drop=FALSE],
                        nsProportion_99quantile=nsProportion_99quantile,
                        rowMeansK=rowMeansK, rowMeansN=rowMeansN, 
                        aberrant=aberrant, aggregate=aggregate,
                        rho=rho, geneColumn=geneColumn,
                        pvalsGene=pvalsGene, padjsGene=padjsGene,
                        aberrantGene=aberrantGene, 
                        additionalColumns=additionalColumns)
            
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
        if(is.null(subsetName)){
            mcols(ans)[["FDR_set"]] <- "transcriptome-wide"
        } else{
            mcols(ans)[["FDR_set"]] <- subsetName
        }
    }
    
    # collapse into one row per gene if requested
    if(isTRUE(aggregate) && isTRUE(collapse)){
        ans <- collapseResTablePerGene(ans)
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
#' @param deltaPsiCutoff The cutoff on delta psi or NA if not requested.
#' @param minCount The minimum count value of the total coverage of an intron 
#' to be considered as significant.
#' result
#' @param rhoCutoff The cutoff value on the fitted rho value 
#' (overdispersion parameter of the betabinomial) above which 
#' junctions are filtered 
#' @param psiType The psi types for which the results should be retrieved.
#' @param additionalColumns Character vector containing the names of additional 
#' columns from mcols(fds) that should appear in the result table 
#' (e.g. ensembl_gene_id). Default is \code{NULL}, so no additional columns 
#' are included. 
#' @param BPPARAM The BiocParallel parameter.
#' @param type Splicing type (psi5, psi3 or theta)
#' @param by By default \code{none} which means no grouping. But if 
#'              \code{sample} or \code{feature} is specified the sum by 
#'              sample or feature is returned
#' @param aggregate If TRUE the returned object is aggregated to the feature 
#'              level (i.e. gene level).
#' @param collapse Only takes effect if \code{aggregate=TRUE}. 
#'              If TRUE, collapses results across the different psi 
#'              types to return only one row per feature (gene) and sample.
#' @param geneColumn The column name of the column that has the gene annotation 
#'              that will be used for gene-level pvalue computation.
#' @param all By default FALSE, only significant introns (or genes) are listed 
#'             in the results. If TRUE, results are assembled for all 
#'             samples and introns/genes regardless of significance. 
#' @param returnTranscriptomewideResults If FDR corrected pvalues for subsets 
#'              of genes of interest have been calculated, this parameter 
#'              indicates whether additionally the transcriptome-wide results 
#'              should be returned as well (default), or whether only results 
#'              for those subsets should be retrieved.
#' @param subsetName The name of a subset of genes of interest for which FDR 
#'             corrected pvalues were previously computed. Those FDR values 
#'             on the subset will then be used to determine aberrant status. 
#'             Default is NULL (using transcriptome-wide FDR corrected pvalues).
#' @param ... Further arguments can be passed to the method. If "n",  
#'              "padjVals", "dPsi" or "rhoVals" are given, the values of those 
#'              arguments are used to define the aberrant events.
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
#' # extract results: for this example dataset, no cutoffs are used to
#' # show the output of the results function
#' res <- results(fds, all=TRUE)
#' res
#' 
#' # aggregate the results by genes (gene symbols need to be annotated first 
#' # using annotateRanges() function)
#' results(fds, padjCutoff=NA, deltaPsiCutoff=0.1, aggregate=TRUE)
#'
#' # aggregate the results by genes and collapse over all psi types to obtain 
#' # only one row per gene in the results table 
#' results(fds, padjCutoff=NA, deltaPsiCutoff=0.1, aggregate=TRUE, 
#'         collapse=TRUE)
#'
#' # get aberrant events per sample: on the example data, nothing is aberrant
#' # based on the adjusted p-value
#' aberrant(fds, type="jaccard", by="sample")
#' 
#' # get aberrant events per gene (first annotate gene symbols)
#' fds <- annotateRangesWithTxDb(fds)
#' aberrant(fds, type="jaccard", by="feature", padjCutoff=NA, aggregate=TRUE)
#'         
#' # find aberrant junctions/splice sites
#' aberrant(fds, type="jaccard")
#' 
#' # retrieve results limiting FDR correction to only a subset of genes
#' # first, we need to create a list of genes per sample that will be tested
#' geneList <- list('sample1'=c("TIMMDC1"), 'sample2'=c("MCOLN1"))
#' fds <- calculatePadjValues(fds, type="jaccard", 
#'                  subsets=list("exampleSubset"=geneList))
#' results(fds, all=TRUE, returnTranscriptomewideResults=FALSE)
#' 
#' @export
setMethod("results", "FraserDataSet", function(object, 
                    sampleIDs=samples(object), padjCutoff=0.1,
                    deltaPsiCutoff=0.1,
                    rhoCutoff=NA, aggregate=FALSE, collapse=FALSE,
                    minCount=5, psiType=psiTypes,
                    geneColumn="hgnc_symbol", all=FALSE,
                    returnTranscriptomewideResults=TRUE,
                    additionalColumns=NULL, BPPARAM=bpparam()){
    psiType <- match.arg(psiType, several.ok=TRUE)
    FDRsets <- availableFDRsubsets(object)
    
    if(isFALSE(returnTranscriptomewideResults) && is.null(FDRsets)){
        warning("Retrieving transcriptome-wide results as no other ",
                "FDR subsets are available in the fds object.")
        returnTranscriptomewideResults <- TRUE
    } 
    if(isTRUE(returnTranscriptomewideResults)){
        res <- FRASER.results(object=object, sampleIDs=sampleIDs, 
                fdrCutoff=padjCutoff, dPsiCutoff=deltaPsiCutoff, 
                rhoCutoff=rhoCutoff, minCount=minCount, 
                psiType=psiType, all=all,
                aggregate=aggregate, collapse=collapse, geneColumn=geneColumn,
                subsetName=NULL, additionalColumns=additionalColumns, 
                BPPARAM=BPPARAM)
    }
    
    # add results for FDR_subsets if requested
    if(!is.null(FDRsets)){
        resls_subsets <- lapply(FDRsets, function(setName){
            res_sub <- FRASER.results(object=object, sampleIDs=sampleIDs, 
                fdrCutoff=padjCutoff, dPsiCutoff=deltaPsiCutoff, 
                rhoCutoff=rhoCutoff, minCount=minCount, 
                psiType=psiType, all=all,
                aggregate=aggregate, collapse=collapse, geneColumn=geneColumn,
                subsetName=setName, additionalColumns=additionalColumns, 
                BPPARAM=BPPARAM)
        })
        
        if(isTRUE(returnTranscriptomewideResults)){
            res <- unlist(GRangesList(unlist(list(res, resls_subsets))))
        } else{
            res <- unlist(GRangesList(unlist(resls_subsets)))
        }
        
        # sort it if existing
        if(length(res) > 0){
            res <- res[order(res$pValue)]
            if(isTRUE(aggregate)){
                res <- res[!is.na(res$pValueGene)]
            }
        }
    }
    return(res)
})

aberrant.FRASER <- function(object, type=fitMetrics(object), 
                                padjCutoff=0.1, deltaPsiCutoff=0.1, 
                                minCount=5, rhoCutoff=NA,
                                by=c("none", "sample", "feature"), 
                                aggregate=FALSE, geneColumn="hgnc_symbol", 
                                subsetName=NULL, all=FALSE, ...){
    
    checkNaAndRange(padjCutoff,     min=0, max=1,   scalar=TRUE,   na.ok=TRUE)
    checkNaAndRange(deltaPsiCutoff, min=0, max=1,   scalar=TRUE,   na.ok=TRUE)
    checkNaAndRange(rhoCutoff,      min=0, max=1,   scalar=TRUE,   na.ok=TRUE)
    checkNaAndRange(minCount,       min=0, max=Inf, scalar=TRUE,   na.ok=TRUE)
    by <- match.arg(by)
    type <- match.arg(type)
    
    if(is.na(rhoCutoff)){
        rhoCutoff <- 1
    }
    
    dots <- list(...)
    if("n" %in% names(dots)){
        n <- dots[['n']]
    } else {
        n <- N(object, type=type)
    }
    if("padjVals" %in% names(dots)){
        padj <- dots[['padjVals']]
    } else {
        # check if padj values are available for the given filters
        pvalsAvailable <- checkPadjAvailableForFilters(object, type=type,
                                                filters=list(rho=rhoCutoff), 
                                                aggregate=aggregate,
                                                subsetName=subsetName)
        if(isFALSE(pvalsAvailable)){
            stop("For the given filters, pvalues are not yet computed. \n", 
                    "Please compute them first by running the ",
                    "calculatePadjValues function with the requested filters.")
        }
        padj <- padjVals(object, type=type, level="site", subsetName=subsetName,
                        filters=list(rho=rhoCutoff))
    }
    if("dPsi" %in% names(dots)){
        dpsi <- dots[['dPsi']]
    } else {
        dpsi <- deltaPsiValue(object, type=type)
    } 
    if("rhoVals" %in% names(dots)){
        rho <- dots[['rhoVals']]
    } else {
        rho <- matrix(rho(object, type=type), 
                        nrow=nrow(dpsi), ncol=ncol(dpsi))
    } 
    if(isTRUE(aggregate)){
        if("padjGeneVals" %in% names(dots)){
            padj_gene <- dots[['padjGeneVals']]
        } else{
            padj_gene <- padjVals(object, type=type, level="gene", 
                                subsetName=subsetName,
                                filters=list(rho=rhoCutoff))    
        }
        
    }
    
    if(is.na(padjCutoff)){
        padjCutoff <- 1
    }
    if(isTRUE(aggregate)){
        padjCutoffGene <- padjCutoff
        padjCutoff <- 1
    }
    
    if(isTRUE(all)){
        aberrantEvents <- matrix(TRUE, nrow=nrow(object), ncol=ncol(object))
        colnames(aberrantEvents) <- colnames(object)
    } else{
        aberrantEvents <- as.matrix(padj) <= padjCutoff
        
        # check each cutoff if in use (not NA)
        if(!is.na(minCount)){
            aberrantEvents <- aberrantEvents & as.matrix(n >= minCount)
        }
        if(!is.na(deltaPsiCutoff)){
            aberrantEvents <- aberrantEvents & 
                as.matrix(abs(dpsi) >= deltaPsiCutoff)
        }
        if(!is.na(rhoCutoff)){
            aberrantEvents <- aberrantEvents & as.matrix(rho <= rhoCutoff)
        }
        aberrantEvents[is.na(aberrantEvents)] <- FALSE
    }
        
    if(isTRUE(aggregate)){
        if(is.null(rownames(padj_gene))){
            stop("Missing rownames for gene-level padj values.")
        }
        # reduce aberrant matrix to one row per gene 
        # (TRUE if any junction is aberrant for each sample)
        ab_dt <- data.table(geneID=getGeneIDs(object, type=type, unique=FALSE,
                                                geneColumn=geneColumn), 
                            aberrantEvents)
        ab_dt[, dt_idx:=seq_len(.N)]
        dt_tmp <- ab_dt[!is.na(geneID), splitGenes(geneID), by="dt_idx"]
        ab_dt <- ab_dt[dt_tmp$dt_idx]
        ab_dt[,`:=`(geneID=dt_tmp$V1, dt_idx=NULL)]
        ab_dt <- ab_dt[,lapply(.SD, any), by="geneID"]
        aberrantEvents <- as.matrix(ab_dt[,-1])
        rownames(aberrantEvents) <- ab_dt[,geneID]
        
        if(isFALSE(all)){
            aberrantEvents <- aberrantEvents & as.matrix(
                    padj_gene[rownames(aberrantEvents),colnames(aberrantEvents)]
                ) <= padjCutoffGene
        }
    }
    
    return(switch(match.arg(by),
                    none = aberrantEvents,
                    sample = colSums(aberrantEvents, na.rm=TRUE),
                    feature = rowSums(aberrantEvents, na.rm=TRUE)
    ))
}

#' @rdname results
#' @export
setMethod("aberrant", "FraserDataSet", aberrant.FRASER)

