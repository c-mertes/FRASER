#'
#' Input check functions
#'
#' Checks all user input and returns corresponding messages
#' checkFraserDataSet
#' 
#' @return logical(1)
#' @rdname checkInputFunctions
#' @noRd
checkFraserDataSet <- function(fds){
    if(!is(fds, "FraserDataSet")){
        stop("Please provide a FraserDataSet object.")
    }
    return(invisible(TRUE))
}

#' @rdname checkInputFunctions
#' @noRd
checkCountData <- function(fds, stop=TRUE){
    checkFraserDataSet(fds)
    if(!all(c("rawCountsJ", "rawCountsSS") %in% assayNames(fds))){
        if(isFALSE(stop)) return(invisible(FALSE))
        stop("No counts detected! Please provide counts first.")
    }
    if(!all(paste0("rawOtherCounts_", psiTypes) %in% assayNames(fds))){
        if(isFALSE(stop)) return(invisible(FALSE))
        stop("Please compute first the total expression at each junction.")
    }
    return(invisible(TRUE))
}

#'
#' clear the files in the cache to start fresh
#'
#' @return nothing
#' @examples
#'     fds <- createTestFraserSettings()
#'     cleanCache(fds)
#' @noRd
cleanCache <- function(fds, all=FALSE, cache=TRUE, assays=FALSE, results=FALSE){
    stopifnot(is(fds, "FraserDataSet"))

    dirs2delete <- c()
    fdsDirName <- nameNoSpace(fds)
    if(cache == TRUE || all == TRUE){
        dirs2delete <- "cache"
    }
    if(assays == TRUE || all == TRUE){
        dirs2delete <- c(dirs2delete, file.path("savedObjects", fdsDirName))
    }
    if(results == TRUE || all == TRUE){
        dirs2delete <- c(dirs2delete, file.path("results", fdsDirName))
    }

    # clean cache
    for(d in dirs2delete){
        full_dir <- file.path(workingDir(fds), d)
        if(dir.exists(full_dir)){
            message(date(), ": Remove directory: '", d, "'.")
            unlink(full_dir, recursive=TRUE)
        }
    }
}

#'
#' checks if the given type is part of the correct category
#' to map correctly to a read type
#'
#' @noRd
checkReadType <- function(fds, type){

    # check if type is null or missing
    if(missing(type) | is.null(type)){
        if(verbose(fds) > 3){
            warning("Read type was not specified!",
                    "We will assume the default: 'j'")
        }
        return("j")
    }
    type <- unique(type)
    stopifnot(isScalarCharacter(type))
    correctTypes <- c(psi3="j", psi5="j", theta="ss", jaccard="j")

    # check if it is already the correct type
    if(type %in% correctTypes) return(type)

    # check if psitype is given
    if(type %in% names(correctTypes)) return(correctTypes[type])

    # check assay names
    atype <- whichReadType(fds, type)
    if(!is.na(atype)) return(atype)

    # regex on the psi type
    atype <- correctTypes[vapply(names(correctTypes), FUN=grepl, type, 
                                FUN.VALUE=logical(1))]
    if(length(atype) == 1){
        return(atype)
    }

    stop("Given read type: '", type, "' not recognized. ",
            "It needs to be 'j' (junction) or 'ss' (splice sites)",
            "\nor an existing assay name within the given object."
    )
}

#'
#' returns the corresonding PSI type to the given name
#'
#' @noRd
whichPSIType <- function(type){
    unlist(regmatches(type, gregexpr("psi(3|5)|theta|jaccard", type, perl=TRUE)))
}

#'
#' returns the read type based on the given assay name
#'
#' @noRd
whichReadType <- function(fds, name){
    stopifnot(isScalarCharacter(name))

    # check writing
    if(name == "ss" | endsWith(name, "theta"))
        return("ss")
    if(name == "j"  | endsWith(name, "psi5") | endsWith(name, "psi3") | 
            endsWith(name, "jaccard"))
        return("j")

    # check assay names
    fdsNames <- assayNames(fds)
    if(name %in% fdsNames){
        nsrNamesL <- length(assayNames(nonSplicedReads(fds)))
        fdsNamesL <- length(fdsNames)

        return(ifelse(
            which(fdsNames == name) <= fdsNamesL - nsrNamesL,
            "j",
            "ss"
        ))
    }

    stop("Could not find read type: ", name)
}

#'
#' Removes the white spaces to have a cleaner file path
#'@noRd
nameNoSpace <- function(name){
    if(is(name, "FraserDataSet")) name <- name(name)
    stopifnot(isScalarCharacter(name))
    gsub("\\s+", "_", name, perl=TRUE)
}

#'
#' Convert to default values
#'
#' Convert all NA's of a input vector or of a
#' single dimension matrix/data.table to FALSE
#'
#' Convert NULL to NA or to another default value
#'
#' @return vector
#' @examples
#'   a <- c(TRUE, FALSE, NA, TRUE, NA)
#'   na2false(a)
#'
#'   dt <- data.table(a)
#'   na2false(dt)
#'
#'   null2na(NULL)
#'   null2na(1:10)
#'
#' @rdname na2default
#' @aliases na2false na2zero null2na null2default
#' @noRd
na2default <- function(x, default=FALSE){
    if(any(class(x) %in% c("DataFrame", "matrix", "data.frame"))){
        stopifnot(dim(x)[2] == 1)
        x <- as.vector(as.matrix(x)[,1])
    }
    x[is.na(x)] <- default
    return(x)
}

#' @rdname na2default
#' @noRd
null2default <- function(x, default=NA){
    if(is.null(x)){
        return(default)
    }
    return(x)
}

#' @rdname na2default
#' @noRd
na2false <- function(x){
    na2default(x, FALSE)
}

#' @rdname na2default
#' @noRd
na2zero <- function(x){
    na2default(x, 0)
}

#' @rdname na2default
#' @noRd
null2na <- function(x){
    null2default(x, NA)
}

#'
#' the qq plot function with confidence band of 5%
#' @noRd
fraserQQplotPlotly <- function(pvalues, ci=TRUE, reducePoints=FALSE,
                    sampleWise=TRUE, main="FRASER QQ-Plot"){
    if(isTRUE(reducePoints)){
        reducePoints <- c(50, 10)
    }

    # convert it to matrix if its a vector
    if(any(is(pvalues, "numeric"))){
        pvalues <- matrix(pvalues)
        colnames(pvalues) <- "observed pvalues"
    } else if(!sampleWise){
        pvalues <- t(pvalues)
    }
    if(!is.matrix(pvalues)){
        pvalues <- as.matrix(pvalues)
    }

    # convert NA to 1
    pvalues[is.na(pvalues)] <- 1

    # length of pvalues
    len_pval <- dim(pvalues)[1]

    # check colnames
    if(is.null(colnames(pvalues))){
        colnames(pvalues) <- seq_len(dim(pvalues)[2])
    }

    # my observerd and expected values
    zeroOffset <- 10e-100
    observ <- -log10(pvalues + zeroOffset)
    expect <- -log10(ppoints(len_pval))

    # create main plot object
    p <- plot_ly(type="scatter", mode="lines")

    # add theoretical trace
    p <- add_trace(p, x=expect, y=expect, mode="lines",
            line=list(color="#FF3030"), name="theoretical-line")
    p <- layout(p, title=main,
            xaxis=list(title="Expected -log<sub>10</sub>(<i>P</i>-value)"),
            yaxis=list(title="Observed -log<sub>10</sub>(<i>P</i>-value)"))

    # add confidence interval
    if(ci){
        if(FALSE){
            # confidence qnorm based from car::qqPlot
            o  <- abs(rnorm(10000)%%1*10^-seq(0, 3, length.out = 10000))
            o  <- sort(o)
            P  <- ppoints(o)
            n  <- length(P)
            zz <- qnorm(1-(1-0.95)/2)
            SE <- (1/dnorm(P))*sqrt(P*(1 - P)/n)
            lower <- zz*SE + P
            upper <- P*P/lower
            lower[lower==1] <- 0.999
            upper[upper==1] <- 0.999
            if(FALSE){
                plot(-log10(P), -log10(o),type="n")
                grid()
                points(-log10(P), -log10(o))
                abline(0,1, col="red")
                lines(-log10(P), -log10(upper), col="blue")
                lines(-log10(P), -log10(lower), col="green")
            }
        } else {
            # confidence qbeta based from GWASTools::qqPlot
            a <- seq_along(expect)
            upper <- -log10(qbeta(0.025, rev(a), a))
            lower <- -log10(qbeta(0.975, rev(a), a))
        }

        path  <- paste("L", c(rev(expect), expect), c(upper, rev(lower)))
        p <- layout(p, shapes=list(list(
            type="path", fillcolor="grey", opacity = 0.3,
            path=paste("M 0 0", paste(path, collapse = " "), "Z")
        )))
    }

    for(idx in seq_len(dim(pvalues)[2])){
        dat <- data.table(
            expect=expect,
            observ=sort(observ[,idx], decreasing=TRUE, na.last=TRUE)
        )
        if(length(reducePoints) > 0 & is.numeric(reducePoints)){
            nEdge <- 50
            nBy   <- 10
            if(is.numeric(reducePoints) & reducePoints[1] > 0 &
                        reducePoints[1] <= len_pval){
                nEdge <- reducePoints[1]
                if(length(reducePoints) == 2 && reducePoints[2] > 0
                            && reducePoints[2] <= 100){
                    nBy <- reducePoints[2]
                }
            }
            dat <- dat[sort(unique(c(
                seq_len(nEdge), -(nEdge-1):0+ldat, seq(1, ldat, nBy)
            )))]
        }
        p <- add_trace(p, data=dat, mode="markers",
                x=~expect, y=~observ, name=colnames(pvalues)[idx], opacity=0.3
        )
    }

    # return object
    return(p)
}

#'
#' logger function for internal use only
#' @noRd
logger <- function(type="INFO", name=flog.namespace(), ...){
    stopifnot(isScalarCharacter(type))
    type <- toupper(type)
    stopifnot(type %in% c("TRACE", "DEBUG", "INFO", "WARN", "ERROR", "FATAL"))

    fun <- switch(type,
            TRACE=flog.trace,
            DEBUG=flog.debug,
            INFO=flog.info,
            WARN=flog.warn,
            ERROR=flog.error,
            FATAL=flog.fatal
    )
    fun(name=name, ...)
}

#'
#' check if the given assay already exists within the object
#' @noRd
assayExists <- function(fds, assayName){
    stopifnot(isScalarCharacter(assayName))
    stopifnot(is(fds, "FraserDataSet"))

    aexists <- assayName %in% assayNames(fds)
    if(aexists){
        message(date(), ": The ", assayName, " are already computed and will ",
                "be used now. If you want to recompute them, please remove ",
                "the following assay: ", assayName, " by issuing following ",
                "command: assays(fds)[['", assayName, "']] <- NULL")
    }
    return(aexists)
}

getAssayAsVector <- function(fds, prefix, psiType=currentType(fds), sampleID){
    as.vector(assay(fds, paste0(prefix, psiType))[,sampleID])
}


variableJunctions <- function(fds, type=currentType(fds), minDeltaPsi=0.1){
    psi <- K(fds, type=type)/N(fds, type=type)
    j2keep <- rowMaxs(abs(psi - rowMeans(psi, na.rm=TRUE)), na.rm=TRUE)
    j2keep >= minDeltaPsi
}

subsetKMostVariableJunctions <- function(fds, type=currentType(fds), n){
    curX <- x(fds, type=type, all=TRUE, center=FALSE, noiseAlpha=NULL)
    xsd <- colSds(curX)
    nMostVarJuncs <- which(xsd >= sort(xsd, TRUE)[min(length(xsd), n*2)])
    ans <- logical(length(xsd))
    ans[sample(nMostVarJuncs, min(length(xsd), n))] <- TRUE
    ans
}

getSubsetVector <- function(fds, type=currentType(fds), minDeltaPsi=0.1, 
                            nSubset=15000){
    # get any variable intron
    ans <- variableJunctions(fds, type, minDeltaPsi=minDeltaPsi)

    # subset most variable intron
    fds_sub <- fds[ans,,by=type]
    ans_sub <- subsetKMostVariableJunctions(fds_sub, type, nSubset)

    # set correct exclusion mask for x computation
    ans[ans] <- ans_sub
    featureExclusionMask(fds) <- exMask
}

pasteTable <- function(x, ...){
    tab <- table(x, ...)
    paste(names(tab), tab, collapse="\t", sep=": ")
}

#'
#' Map between individual seq level style and dataset common one
#' for counting and aggregating the reads
#' @noRd
checkSeqLevelStyle <- function(gr, fds, sampleID, sampleSpecific=FALSE,
                    coldata=colData(fds)){
    if(!"SeqLevelStyle" %in% colnames(coldata)){
        return(gr)
    }
    style <- coldata[sampleID,"SeqLevelStyle"]
    if(isFALSE(sampleSpecific)){
        style <- names(sort(table(coldata[,"SeqLevelStyle"]), TRUE)[1])
        if(length(unique(coldata[,"SeqLevelStyle"])) > 1){
            gr <- keepStandardChromosomes(gr, pruning.mode="coarse")
        }
    }

    gr <- keepSeqlevels(gr, unique(seqnames(gr)))
    seqlevelsStyle(gr) <- style
    gr
}

uniformSeqInfo <- function(grls){
    tmpSeqlevels <- unique(data.table(
        seqlevel  = unlist(lapply(grls, seqlevels)),
        seqlength = unlist(lapply(grls, seqlengths))
    )[order(seqlevel)])

    chromosomeNames <- tmpSeqlevels[,seqlevel]
    if(any(duplicated(chromosomeNames))){
        nonUniqueChromosomes <- chromosomeNames[duplicated(chromosomeNames)]
        warning("There are non unique chromosomes in this dataset (", 
                nonUniqueChromosomes, ")!\n This means that these chromosomes ",
                "are associated with > 1 chromosome length!\n These ",
                "chromosomes will be removed from the dataset as mapping ",
                "between intron\n positions between samples is not possible.")
        tmpSeqlevels <- tmpSeqlevels[seqlevel != nonUniqueChromosomes]
    }

    ans <- lapply(grls, function(x){
        seqlevels(x, pruning.mode="coarse")  <- tmpSeqlevels[,seqlevel]
        seqlengths(x) <- tmpSeqlevels[,seqlength]
        x
    })
    ans
}

getHDF5ChunkSize <- function(fds, assayName){
    # taken from here : https://github.com/grimbough/rhdf5/commit/52af7840c7
    # To be backwards compatible to R version 3.5.0
    # 
    h5file <- getFraserHDF5File(fds, assayName)
    h5obj <- H5Fopen(h5file, flags="H5F_ACC_RDONLY")
    h5dataset <- h5obj&assayName
    
    pid <- H5Dget_create_plist(h5dataset)
    
    on.exit(H5Pclose(pid), add=TRUE)
    on.exit(H5Fclose(h5obj), add=TRUE)
    
    if (H5Pget_layout(pid) != "H5D_CHUNKED")
        return(NULL)
    else 
        return(rev(H5Pget_chunk(pid)))
}

getMaxChunks2Read <- function(fds, assayName, max=15, axis=c("col", "row")){
    axis <- match.arg(axis)
    if(!any(c("DelayedArray", "DelayedMatrix") %in%
            class(assay(fds, assayName)))){
        if(axis == "col"){
            return(ceiling(ncol(assay(fds, assayName))/bpnworkers(bpparam())))
        }
        return(ceiling(nrow(assay(fds, assayName))/bpnworkers(bpparam())))
    }

    dims <- getHDF5ChunkSize(fds, assayName)
    if(axis == "col"){
        ans <- dims[2]
    } else {
        ans <- dims[1]
    }
    max(1, ans/ceiling(ans/max))
}

getSamplesByChunk <- function(fds, sampleIDs, chunkSize){
    chunks <- trunc(0:(ncol(fds)-1)/chunkSize)
    ans <- lapply(0:max(chunks), function(x){
        intersect(sampleIDs, samples(fds)[chunks == x])
    })
    ans[vapply(ans, length, integer(1)) >0]
}

checkNaAndRange <- function(x, min=-Inf, max=Inf, scalar=TRUE, na.ok=FALSE){
    xname <- deparse(substitute(x))
    if(isTRUE(scalar) & !isScalarValue(x)){
        stop(xname, " should be a scalar value!")
    }
    if(any(is.na(x)) && isFALSE(na.ok)){
        stop(xname, " contains NA values, which is not allowed.")
    }
    if(sum(!is.na(x)) == 0){
        return(invisible(TRUE))
    }
    if(!is.numeric(x[!is.na(x)])){
        stop(xname, " should be numeric!")
    }
    if(x[!is.na(x)] < min){
        stop(xname, " should be bigger than ", min)
    }
    if(x[!is.na(x)] > max){
        stop(xname, " should be smaller than ", max)
    }
    invisible(TRUE)
}

putCounts2Memory <- function(fds, type=currentType(fds)){
    counts(fds, type=type, side="other", HDF5=FALSE) <-
            as.matrix(counts(fds, type=type, side="other"))
    counts(fds, type=type, side="ofInterest", HDF5=FALSE) <- 
            as.matrix(counts(fds, type=type, side="ofInterest"))
    fds
}


plotBasePlot <- function(ggplot, basePlot=FALSE){
    if(isFALSE(basePlot)){
        ggplot$labels <- lapply(ggplot$labels, function(x){
                if(typeof(x) == "expression"){
                    warning("Found expression for plotly. Please adapt it!")
                    return(as.character(x))
                }
                x})
        return(plotly::ggplotly(ggplot, tooltip="text"))
    }
    ggplot
}

getBPParam <- function(worker, tasks=0, ...){
    if(worker < 2){
        return(SerialParam(...))
    }
    worker <- min(worker, multicoreWorkers())
    if(.Platform$OS.type != "unix") {
        ans <- SnowParam(workers=worker, tasks=tasks, ...)
    } else {
        ans <- MulticoreParam(workers=worker, tasks, ...)
    }
    ans
}

getStrandString <- function(fds){
    strand <- sapply(strandSpecific(fds)+1L, switch, "no", "yes", "reverse")
    strand <- paste(collapse=", ", strand)
    return(strand)
}


#'
#' Check if adjusted pvalues have been computed for a given set of filters.
#' @noRd
checkPadjAvailableForFilters <- function(fds, type=currentType(fds), 
                    filters=list(), dist="BetaBinomial", aggregate=FALSE,
                    subsetName=NULL){
    dist <- match.arg(dist, choices=c("BetaBinomial", "Binomial", "Normal"))
    aname <- paste0("padj", dist)
    aname <- ifelse(isTRUE(aggregate), paste0(aname, "_gene"), aname)
    aname <- ifelse(!is.null(subsetName), paste0(aname, "_", subsetName), aname)

    # add information on used filters
    for(n in sort(names(filters))){
        aname_new <- paste0(aname, "_", n, filters[[n]])
        if(n == "rho" && filters[[n]] == 1){
            if(any(grepl(aname_new, assayNames(fds))) ||
                    any(grepl(aname_new, names(metadata(fds))))){
                aname <- aname_new
            }
        }else{
            aname <- aname_new
        }
    }
    aname <- paste(aname, type, sep="_")
    if(isTRUE(aggregate)){
        pvalsAvailable <- aname %in% names(metadata(fds))
    } else{
        pvalsAvailable <- aname %in% assayNames(fds)
    }
    return(pvalsAvailable)
}

#'
#' Find most aberrant junction for each aberrant gene
#' 
#' @param gr GRanges object with information about junctions.
#' @param aberrantGenes Significant genes for which the corresponding junction 
#'     should be extracted. 
#' @param pvals Vector of pvalues (for one sample).
#' @param dpsi Vector of delta psi values (for one sample).
#' @param aberrantJunctions Vector indicating which junctions are considered 
#'     aberrant. 
#' @param geneColumn Name of the column in mcols(fds) that has gene annotation.
#' @noRd
findJunctionsForAberrantGenes <- function(gr, aberrantGenes, pvals, dpsi, 
                                aberrantJunctions, geneColumn="hgnc_symbol"){
    dt <- data.table(idx=mcols(gr)$idx, 
                        geneID=mcols(gr)[,geneColumn],
                        pval=pvals, 
                        dpsi=abs(dpsi),
                        aberrant=aberrantJunctions)
    dt[, dt_idx:=seq_len(.N)]
    dt_tmp <- dt[, splitGenes(geneID), by="dt_idx"]
    dt <- dt[dt_tmp$dt_idx,]
    dt[,`:=`(geneID=dt_tmp$V1, dt_idx=NULL)]
    dt <- dt[geneID %in% aberrantGenes,]
    dt <- dt[!is.na(aberrant) & aberrant == TRUE,]
    
    # sort per gene by lowest pvalue / highest deltaPsi and return index
    dt <- dt[order(geneID, -aberrant, pval, -dpsi)]
    dt <- dt[!duplicated(dt, by="geneID"),]
    
    # remove gene-level significant result if no junction in that gene passed
    # the filters
    dt <- dt[!is.na(pval),]
    
    junctionsToReport <- dt[,idx]
    names(junctionsToReport) <- dt[,geneID]
    junctionsToReport <- sort(junctionsToReport)
    return(junctionsToReport)
}

collapseResTablePerGene <- function(res, geneColumn="hgncSymbol"){
    if(length(res) == 0){
        return(res)
    }
    if(!is.data.table(res)){
        res <- as.data.table(res)
    }
    
    if(any(!c("pValue", "pValueGene", geneColumn) %in% colnames(res))){
        stop("For collapsing per gene, the results table needs to contain ",
            "the columns pValue, pValueGene and ", geneColumn, ".")
    }
    
    res <- res[order(res$pValueGene, res$pValue)]
    naIdx <- is.na(res[, get(geneColumn)])
    ansNoNA <- res[!is.na(res[, get(geneColumn)]),]
    
    # get final result table
    dupIdx <- duplicated(data.table(as.vector(ansNoNA[, get(geneColumn)]), 
                                    as.vector(ansNoNA$sampleID)))
    ans <- res[!naIdx,][!dupIdx,]
    return(ans)
}

#' ignores NA in unique if other values than NA are present
#' @noRd
uniqueIgnoreNA <- function(x){
    uniq <- unique(x)
    if(length(uniq) > 1) uniq <- uniq[!is.na(uniq)]
    return(uniq)
}

#' split string of gene names into vector
#' @noRd
splitGenes <- function(x, sep=";"){
    return(unlist(strsplit(as.character(x), sep, fixed=TRUE)))
}

#' cap string of gene names to show max 3 gene names
#' @noRd
limitGeneNamesList <- function(gene_names, maxLength=3){
    gene_names <- as.character(gene_names)
    numFeatures <- unlist(lapply(gene_names, function(x) length(splitGenes(x))))
    gene_names[numFeatures > maxLength] <- 
        unlist(lapply(gene_names[numFeatures > maxLength], function(x){
            paste(c(splitGenes(x)[seq_len(maxLength)], "..."), 
                    collapse=";") 
        } ))
    return(gene_names)
}

checkForAndCreateDir <- function(object, dir){
    verbose <- 0
    if(is(object, "FraserDataSet")){
        verbose <- verbose(object)
        if(missing(dir)){
            dir <- workingDir(object)
        }
    }
    if(!dir.exists(dir)){
        if(verbose > 1){
            message(date(), ": The given working directory '", 
                    dir, "' does not exists. We will create it.")
        }
        dir.create(dir, recursive=TRUE)
    }
    if(!dir.exists(dir)){
        stop("Can not create workding directory: ", dir)
    }
    return(TRUE)
}
