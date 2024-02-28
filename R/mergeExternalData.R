#' 
#' Merge external data
#' 
#' To boost its own sequencing data, one can download existing and precounted 
#' data. This function merges the existing \code{FraserDataSet} with
#' external count data.
#' 
#' For more details on existing datasets have a look at:
#' <https://github.com/gagneurlab/drop#datasets>
#' 
#' Since FRASER can not hand NA values, the merge will return only the 
#' intersecting regions and will drop any non overlapping features. This has to
#' be kept in mind when analysing rare disease samples.
#' 
#' @param fds A \code{FraserDataSet}
#' @param countFiles A character vector of file names pointing to the external 
#'           count data. The vector has to be names or the files have to start 
#'           with \code{k_j}, \code{k_theta}, \code{n_psi3}, \code{n_psi5}, 
#'           \code{n_theta}.
#' @param sampleIDs The samples to be merged from the external data. 
#' @param annotation A sample annotation of the external data (optional).
#' 
#' @return Merged \code{FraserDataSet} object.
#' 
#' @examples
#' anno <- fread(system.file("extdata", "externalCounts", 
#'         "annotation.tsv.gz", package="FRASER"))
#' ctsFiles <- list.files(full.names = TRUE, pattern="counts",
#'         system.file("extdata", "externalCounts", package="FRASER"))
#' 
#' fds <- createTestFraserDataSet()
#' fds_merged <- mergeExternalData(fds, ctsFiles, anno[,sampleID], anno)
#' 
#' K(fds, "psi5")
#' K(fds_merged, "psi5")
#' 
#' @export
mergeExternalData <- function(fds, countFiles, sampleIDs, annotation=NULL){

    # check input
    checkFraserDataSet(fds)
    checkCountData(fds)

    if(length(countFiles) != 5){
        stop("You have to provide 5 files, but only ", length(countFiles),
                " were provided.")
    }
    if(is.null(names(countFiles))){
        names(countFiles) <- gsub("(_counts)?.tsv.gz", "", basename(countFiles))
    }
    reqNames <- c("k_j", "k_theta", "n_psi3", "n_psi5", "n_theta")
    if(any(!reqNames %in% unique(names(countFiles)))){
        stop("Please name the input or the files correctly. We are missing: ",
                paste(collapse=", ",
                        reqNames[!reqNames %in% names(countFiles)]))
    }
    if(any(!file.exists(countFiles))){
        stop("Provided files are missing! We are missing: ",
                paste(collapse=", ", countFiles[!file.exists(countFiles)]))
    }
    if(any(unique(sampleIDs) != sampleIDs)){
        stop("Provided sampleIDs have to be unique!")
    }
    sampleIDs <- as.character(sampleIDs)

    # load external counts
    message("Loading provided counts")
    names(reqNames) <- reqNames
    extCts <- lapply(reqNames, function(id){
        gr <- makeGRangesFromDataFrame(fread(countFiles[id]),
                keep.extra.columns=TRUE)
        seqlevelsStyle(gr) <- seqlevelsStyle(fds) #force fds style onto external counts
        if(any(!sampleIDs %in% colnames(mcols(gr)))){
            stop("Can not find provided sampleID in count data. Missing IDs: ",
                    paste(collapse=", ",
                            sampleIDs[!sampleIDs %in% colnames(mcols(gr))]))
        }
        gr[,sampleIDs]
    })
    stopifnot(all(granges(extCts[['k_j']]) == granges(extCts[['n_psi3']])))
    stopifnot(all(granges(extCts[['k_j']]) == granges(extCts[['n_psi5']])))
    stopifnot(all(granges(extCts[['k_theta']]) == granges(extCts[['n_theta']])))

    #
    # merging colData
    #
    message(date(), ": Merging data ...")
    if(!is.null(annotation)){
        annotation <- DataFrame(annotation)
    } else {
        annotation <- DataFrame(sampleID=as.character(sampleIDs))
    }
    rownames(annotation) <- annotation[,"sampleID"]
    newColData <- DataFrame(rbind(fill=TRUE,
            as.data.table(colData(fds)),
            as.data.table(annotation[sampleIDs,])))
    rownames(newColData) <- newColData[,"sampleID"]

    #
    # merge psi5/psi3 data
    #
    extractExtData <- function(fds, countFun, type, ov, extData, extName){
        ctsOri <- as.matrix(countFun(fds, type=type)[from(ov),,drop=FALSE])
        ctsExt <- as.matrix(mcols(extData[[extName]])[to(ov),,drop=FALSE])
        ans <- cbind(ctsOri, ctsExt)
        mode(ans) <- "integer"
        ans
    }

    # find overlap
    if(all(strand(rowRanges(fds, type="j")) == "*")){
        for(id in reqNames){
            strand(extCts[[id]]) <- "*"
        }
    }
    ov <- findOverlaps(rowRanges(fds, type="j"), extCts[['k_j']], type="equal")

    newCtsK_J    <- extractExtData(fds, K, "j",    ov, extCts, "k_j")
    newCtsN_psi5 <- extractExtData(fds, N, "psi5", ov, extCts, "n_psi5")
    newCtsN_psi3 <- extractExtData(fds, N, "psi3", ov, extCts, "n_psi3")

    # get ranges after merging
    SR_ranges <- rowRanges(fds)[from(ov),c("startID", "endID")]

    #
    # merge theta data
    #
    # find overlap
    ovss <- findOverlaps(rowRanges(fds, type="theta"),
            extCts[['k_theta']], type="equal")

    newCtsK_theta <- extractExtData(fds, K, "theta", ovss, extCts, "k_theta")
    newCtsN_theta <- extractExtData(fds, N, "theta", ovss, extCts, "n_theta")
    NSR_ranges <- rowRanges(fds, type="theta")[from(ovss),c("spliceSiteID", "type")]

    # Find the overlaps of the NSR with SR after merging/filtering
    NSR_index <- which(NSR_ranges$spliceSiteID %in% c(SR_ranges$startID, SR_ranges$endID))

    # Only take NSR that have at least 1 split read over the same junction.
    NSR_ranges <- NSR_ranges[NSR_index]
    newCtsK_theta <- newCtsK_theta[NSR_index,]
    newCtsN_theta <- newCtsN_theta[NSR_index,]

    #
    # finalize merged FraserDataObject
    #
    nsr <- SummarizedExperiment(
            colData = newColData,
            assays = SimpleList(
                    rawCountsSS = newCtsK_theta,
                    rawOtherCounts_theta = (newCtsN_theta - newCtsK_theta)),
            rowRanges= NSR_ranges
    )

    ans <- new("FraserDataSet",
            name = name(fds),
            bamParam = scanBamParam(fds),
            workingDir = workingDir(fds),
            colData = newColData,
            assays = Assays(SimpleList(
                    rawCountsJ = newCtsK_J,
                    rawOtherCounts_psi5 = newCtsN_psi5 - newCtsK_J,
                    rawOtherCounts_psi3 = newCtsN_psi3 - newCtsK_J)),
            nonSplicedReads = nsr,
            rowRanges = SR_ranges,
            elementMetadata = DataFrame(newCtsK_J[,integer(0)]),
            metadata = metadata(fds)
    )

    #
    # compute new psi values
    #
    ans <- calculatePSIValues(ans)

    ans
}
