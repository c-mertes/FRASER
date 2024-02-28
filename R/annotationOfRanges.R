#'
#' Annotates the given FraserDataSet with the HGNC symbol with biomaRt
#' 
#' @param fds FraserDataSet
#' @param feature Defines which feature (default is HGNC symbol) should be 
#'             annotated. Has to be the \code{biomaRt} feature name or a 
#'             column name in \code{orgDb}.
#' @param featureName The column name of the feature in the FraserDataSet mcols.
#' @param biotype The biotype for \code{biomaRt}.
#' @param ensembl The ensembl that should be used. If NULL, the default one is 
#'             used (hsapiens_gene_ensembl, GRCh37).
#' @param GRCh GRCh version to connect to. If this is NULL, then the current 
#'             GRCh38 is used. Otherwise, this can only be 37 (default) 
#'             at the moment (see \code{\link[biomaRt]{useEnsembl}}).
#' @param txdb A \code{TxDb} object. If this is NULL, then the default 
#'             one is used, currently this is 
#'             \code{TxDb.Hsapiens.UCSC.hg19.knownGene}.
#' @param orgDb An \code{orgDb} object or a data table to map the feature names.
#'             If this is NULL, then \code{org.Hs.eg.db} is used as the default.
#' @param filter A named list specifying the filters which should be applied to 
#'             subset to e.g. only protein-coding genes for annotation. 
#'             \code{names(filter)} needs to be column names in the given 
#'             orgDb object (default: no filtering).
#' @param keytype The keytype or column name of gene IDs in the \code{TxDb}
#'             object (see 
#'             \code{\link[AnnotationDbi:AnnotationDb-class]{keytypes}}
#'             for a list of available ID types).
#' 
#' @return FraserDataSet
#' 
#' @examples
#'
#' fds <- createTestFraserDataSet()
#' 
#' ### Two ways to annotage ranges with gene names: 
#' # either using biomart with GRCh38
#' try({
#'   fds <- annotateRanges(fds, GRCh=38)
#'   rowRanges(fds, type="j")[,c("hgnc_symbol")]
#' })
#' 
#' # either using biomart with GRCh37
#' try({
#'   fds <- annotateRanges(fds, featureName="hgnc_symbol_37", GRCh=37)
#'   rowRanges(fds, type="j")[,c("hgnc_symbol_37")]
#' })
#'  
#' # or with a provided TxDb object
#' require(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' require(org.Hs.eg.db)
#' orgDb <- org.Hs.eg.db
#' fds <- annotateRangesWithTxDb(fds, txdb=txdb, orgDb=orgDb)
#' rowRanges(fds, type="j")[,"hgnc_symbol"]
#' 
#' @rdname annotateRanges
#' @export
annotateRanges <- function(fds, feature="hgnc_symbol", featureName=feature,
            biotype=list("protein_coding"), ensembl=NULL, GRCh=37){

    # check input
    stopifnot(is(fds, "FraserDataSet"))
    if(length(fds) == 0) return(fds)
    
    # useEnsembl only understands GRCh=37 or GRCh=NULL (uses 38 then)
    if(GRCh == 38){
        GRCh <- NULL
    }

    if(is.null(ensembl)){
        tryCatch({
            ensemblOutput <- capture.output(ensembl <- useEnsembl(
                    biomart="ensembl", dataset="hsapiens_gene_ensembl", 
                    GRCh=GRCh))
        },
        error=function(e){
                message("\nCheck if we have a internet connection!",
                        " Could not connect to ENSEMBL. With error: `", e, "`.")
        })
        if(is.null(ensembl)){
            message("Nothing was annotated!")
            return(fds)
        }
    }

    # check if USCS naming scheme is used
    useUSCS <- all(startsWith(seqlevels(fds), "chr"))

    # get gene annotation
    annotation <- getFeatureAsGRange(ensembl, feature, featureName,
            biotype, useUSCS)

    # annotate splice sites first
    gr <- rowRanges(fds, type="theta")
    if(any(strand(gr) == "*")){
        strand(annotation) <- "*"
    }
    annos <- getAnnotationFeature(data=gr, featureName, annotation)
    mcols(fds, type="theta")[[featureName]] <- annos
    
    # annotate junctions with genes at donor and acceptor sites
    fds <- annotateFeatureFromSpliceSite(fds, featureName)

    return(fds)
}

#' @rdname annotateRanges
#' @export
annotateRangesWithTxDb <- function(fds, feature="SYMBOL", 
                    featureName="hgnc_symbol", keytype="ENTREZID",
                    txdb=NULL, orgDb=NULL, filter=list()){
    gene_id <- NULL
    
    # check input
    stopifnot(is(fds, "FraserDataSet"))
    if(length(fds) == 0) return(fds)
    
    if(is.null(txdb)){
        if(requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene")){
            txdb <- 
        TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene 
        } else{
            stop("Please provide a TxDb object as input.")
        }
        
    }
    if(is.null(orgDb)){
        if(requireNamespace("org.Hs.eg.db")){
            orgDb <- org.Hs.eg.db::org.Hs.eg.db
        } else{
            stop("Please provide an OrgDb object to extract gene symbols")
        }
    }
    
    # get GRanges object with the splice sites which should be annotated
    gr <- rowRanges(fds, type="theta")
    
    # get the annotation to compare to
    anno <- genes(txdb)
    if(is.data.table(orgDb)){
        tmp <- merge(x=as.data.table(anno)[,.(gene_id)], y=orgDb, 
                by.y=keytype, by.x="gene_id", all.x=TRUE, sort=FALSE)[,
                            c("gene_id", feature, names(filter)), with=FALSE]
    } else {
        tmp  <- as.data.table(select(orgDb, keys=mcols(anno)[,"gene_id"], 
                columns=c(feature, names(filter)), keytype=keytype))
    }
        
    # filter genes as specified by user (e.g. only protein_coding)
    tmp[, include:=TRUE]
    if(!is.null(filter) & length(filter) > 0 & !is.null(names(filter))){
        for(n in names(filter)){
            stopifnot(n %in% colnames(tmp))
            tmp[!(get(n) %in% filter[[n]]), include:=FALSE]
        }
    }
    
    # add the new feature to the annotation
    tmp[, uniqueID := .GRP, by=keytype]
    tmp <- tmp[include == TRUE,]
    anno <- anno[tmp[,uniqueID]]
    mcols(anno)[[featureName]] <- tmp[,get(feature)]
        
    # clean up of NA and "" ids
    anno <- anno[!is.na(mcols(anno)[,featureName]),]
    anno <- anno[mcols(anno)[,featureName] != "",]
    if(any(strand(gr) == "*")){
        strand(anno) <- "*"
    }
    
    # retrieve the feature of interest for the splice sites
    annos <- getAnnotationFeature(data=gr, featureName, anno)
    mcols(fds, type="theta")[[featureName]] <- annos
    
    # transfer annoated features for splice sites to junctions
    fds <- annotateFeatureFromSpliceSite(fds, featureName)
    
    return(fds)   
}

#'
#' use biomart to extract the current feature annotation
#'
#' @noRd
getFeatureAsGRange <- function(ensembl, feature, featureName,
                    biotype, useUSCS=TRUE){

    # contact biomaRt to retrive hgnc symbols
    ensemblResults <- getBM(attributes=c(feature, "chromosome_name",
                    "start_position", "end_position", "strand"),
            filters=c("biotype"), values=c(biotype), mart=ensembl)
    setnames(ensemblResults, feature, featureName)

    # remove emtpy symbols or non standard chromosomes
    ensemblResults <- ensemblResults[!is.na(ensemblResults[[featureName]]),]
    ensemblResults <- ensemblResults[ensemblResults[[featureName]] != "",]
    ensemblResults <- ensemblResults[
            !grepl("_|\\.", ensemblResults$chromosome_name),]
    ensemblResults[,"strand"]  <- ifelse(
            ensemblResults[,"strand"] < 0, "-", "+")
    
    # create a grange object
    results <- makeGRangesFromDataFrame(ensemblResults,
            start.field="start_position", end.field="end_position", 
            strand.field="strand", keep.extra.columns=TRUE)

    if(useUSCS){
        seqlevels(results) <- paste0("chr", seqlevels(results))
    }

    return(results)
}


#'
#' merge the retrieved annotations from a biomart with the ranges of our data
#'
#' @noRd
getAnnotationFeature <- function(data, feature, annotation){

    # find overlap with the query
    # suppress seqnames not in anothers object!
    # TODO: is there a nicer way to do this?
    suppressWarnings(hits <- findOverlaps(data, annotation))

    # find overlap and extract it
    featureDT <- data.table(
            from=from(hits),
            feature=mcols(annotation[to(hits)])[[feature]]
    )
    missingValues <- setdiff(seq_along(data), unique(from(hits)))
    if(length(missingValues) > 0){
        featureDT <- rbind(
                featureDT, data.table(from=missingValues, feature=NA)
        )
    }
    
    # extract only the feature and group them with a ";"
    featureDT <- featureDT[,feature:=paste(unique(feature), collapse = ";"),
                            by="from"]
    featureDT <- featureDT[!duplicated(featureDT),]
    featureDT[feature == "NA", feature:=NA]

    return(featureDT[order(from),feature])
}


#'
#' Find annotated junctions
#' 
#' Annotate the object with a given annotation. 
#' 
#' @return FraserDataSet
#' @examples 
#'     TODO <- 1
#' @noRd
findAnnotatedJunction <- function(fds, annotation, annotateNames=TRUE,
                    seqLevelStyle=seqlevelsStyle(fds)[1],
                    stranded=strandSpecific(fds), ...){
    
    # import annotation if given as file
    if(isScalarCharacter(annotation)){
        txdb <- makeTxDbFromGFF(annotation, ...)
    }
    
    # extract introns
    intronsBT <- intronsByTranscript(txdb, use.names=TRUE)
    introns <- unlist(intronsBT)
    
    # check if strandspecific data is used
    gr <- rowRanges(fds, type="psi5")
    
    # 
    if(0 %in% stranded){
        if(length(unique(stranded)) > 1){
            warning("You use a mixture of stranded and unstranded data.", 
                    " This is highly discouraged and ",
                    " can lead to unwanted results. We treat all as stranded.")
        } else {
            strand(gr) <- "*"
        }
    }
    
    # set correct seq level names
    seqlevelsStyle(introns) <- seqLevelStyle
    seqlevelsStyle(gr) <- seqLevelStyle
    
    # find overlaps with introns
    ov <- findOverlaps(gr, introns, type="equal")
    newAnno <- data.table(idx=seq_row(fds), known_intron=FALSE, 
            known_start=FALSE, known_end=FALSE, TranscriptNames=NA_character_)
    newAnno[idx %in% from(ov), known_intron:=TRUE]
    
    if(isTRUE(annotateNames)){
        ovanno <- as.data.table(ov)
        ovanno[,name:=names(introns)[subjectHits]]
        ovanno <- ovanno[,.(
                TranscriptNames=paste(unique(name), collapse=",")),
                by=queryHits]
        newAnno <- merge(newAnno, ovanno, by.x="idx", 
                by.y="queryHits", all.x=TRUE)
    }
    
    # overlap start/stop
    ov <- findOverlaps(gr, introns, type="start")
    ovdt <- data.table(from=from(ov), to=to(ov), 
            strand=as.vector(strand(introns)[to(ov)]))[,
                    .(strand=strand[1]),by=from]
    newAnno[ovdt$from, known_start:=known_start | ovdt$strand %in% c("+", "*")]
    newAnno[ovdt$from, known_end  :=known_end   | ovdt$strand %in% c("-")]
    
    ov <- findOverlaps(gr, introns, type="end")
    ovdt <- data.table(from=from(ov), to=to(ov), 
            strand=as.vector(strand(introns)[to(ov)]))[,
                    .(strand=strand[1]),by=from]
    newAnno[ovdt$from, known_start:=known_start | ovdt$strand %in% c("-")]
    newAnno[ovdt$from, known_end  :=known_end   | ovdt$strand %in% c("+", "*")]
    
    # update fds object
    columns2Write <- c("known_intron", "known_start", "known_end", 
            "TranscriptNames")
    col2Take <- !colnames(mcols(fds, type="psi5")) %in% columns2Write
    mcols(fds, type="psi5") <- cbind(
            mcols(fds, type="psi5")[,col2Take], newAnno[,..columns2Write])
    
    # return annotated object
    fds
}

#' annotate junctions with genes at donor and acceptor sites
#' @noRd
annotateFeatureFromSpliceSite <- function(fds, featureName){
    ssdt <- data.table(spliceSiteID=mcols(fds, type="theta")$spliceSiteID,
                        genes=mcols(fds, type="theta")[[featureName]] 
    )
    junction_dt <- data.table(startID=mcols(fds, type="psi3")$startID,
                                endID=mcols(fds, type="psi3")$endID
    )
    junction_dt <- merge(junction_dt, ssdt, all.x=TRUE, 
                        by.x="startID", by.y="spliceSiteID", sort=FALSE)
    setnames(junction_dt, "genes", "genes_donor")
    junction_dt <- merge(junction_dt, ssdt, all.x=TRUE, 
                        by.x="endID", by.y="spliceSiteID", sort=FALSE)
    setnames(junction_dt, "genes", "genes_acceptor")
    
    junction_dt[,genes:=paste(uniqueIgnoreNA(
                    c(splitGenes(genes_donor), splitGenes(genes_acceptor))), 
                    collapse=";"),
                by="startID,endID"]
    junction_dt[genes == "NA", genes:=NA]
    mcols(fds, type="j")[[featureName]] <- junction_dt[,genes]
    return(fds)
}
