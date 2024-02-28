#'
#' Create an simulated example data set for FRASER
#'
#' Simulates a data set based on random counts following a
#' beta binomial (or Dirichlet-Multinomial) distribution.
#'
#' @param j Number of simulated junctions
#' @param m Number of simulated samples
#' @param q number of simulated latent variables.
#' @param distribution Either "BB" for a beta-binomial simulation or "DM" for a 
#' dirichlet-multinomial simulation.
#' @param ... Further arguments used to construct the FraserDataSet.
#'
#' @return An FraserDataSet containing an example dataset based on
#'         simulated data
#'
#' @examples
#' # A generic dataset
#' fds1 <- makeSimulatedFraserDataSet()
#' fds1
#'
#' # A generic dataset with specificed sample size and injection method
#' fds2 <- makeSimulatedFraserDataSet(m=10, j=100, q=3)
#' fds2
#'
#' @rdname makeSimulatedFraserDataSet
#' @aliases makeSimulatedFraserDataSet
#' @export
makeSimulatedFraserDataSet <- function(m=100, j=500, q=10,
                distribution=c("BB", "DM"), ...){
    distribution <- match.arg(distribution)
    fds <- switch (distribution,
        BB = makeSimulatedFraserDataSet_BetaBinomial(m=m, j=j, q=q, ...),
        DM = makeSimulatedFraserDataSet_Multinomial(m=m, j=j, q=q, ...))
    
    dontWriteHDF5(fds) <- TRUE
    fds
}

makeSimulatedFraserDataSet_BetaBinomial <- function(m=200, j=10000, q=10,
                nbMeanLog=8, nbSize=0.55, rhoMeanlog=-5, rhoSdlog=1,
                Hmean=0, Hsd=100, Dmean=0, Dsd=0.007, ...){
    
    # Simulate total coverage N = nonSplit + n at each junction using lognormal
    junction_n <- rlnorm(j, meanlog=nbMeanLog, sdlog=nbSize)
    N <- t(vapply(junction_n, function(x){
        # draw counts for every sample from negative binomial
        n <- rnbinom(m, mu=x, size=nbSize)},
        double(m)))
    # N <- matrix(rnbinom(j*m, mu=nbMean, size=nbSize), nrow=j, ncol=m)
    mode(N) <- 'integer'
    
    #
    # Simulate covariates for SE
    #
    H_true <- matrix(rnorm(m*q, mean=Hmean, sd=Hsd), nrow=m, ncol=q)
    D_true <- matrix(rnorm(j*q, mean=Dmean, sd=Dsd), nrow=j, ncol=q)
    b_true <- sample(c(rnorm(round(j*(1/6)), mean=3.5, sd=1.5),
            rnorm(j - round(j*(1/6)), mean=-2.5, sd=2.5)))
    y_se   <- t(predictYCpp(H_true, D_true, b_true))
    mu_se     <- predictMuCpp(y_se)
    # betaBin dispersion
    rho_se <- rlnorm(j, meanlog=rhoMeanlog, sdlog=rhoSdlog)
    
    # Draw nonSplit reads from BB
    nonSplit <-matrix(rbetabinom(j*m, size=N, prob=mu_se, rho=rho_se),
            nrow=j, ncol=m)
    mode(nonSplit) <- 'integer'
    
    # Set n = N - nonSplit
    n <- N - nonSplit
    
    #
    # Simulate covariates for PSI3=PSI5
    #
    H_true <- matrix(rnorm(m*q, mean=Hmean, sd=Hsd), nrow=m, ncol=q)
    D_true <- matrix(rnorm(j*q, mean=Dmean, sd=Dsd), nrow=j, ncol=q)
    b_true <- sample(c(rnorm(round(j*(1/6)), mean=-2.5, sd=2.5),
            rnorm(j-round(j*(1/6)), mean=3.5, sd=1.5)))
    y_psi  <- t(predictYCpp(H_true, D_true, b_true))
    mu_psi <- predictMuCpp(y_psi)
    # betaBin dispersion
    rho_psi <- rlnorm(j, meanlog=rhoMeanlog, sdlog=rhoSdlog)
    
    # Draw nonSplit reads from BB
    k <- matrix(rbetabinom(j*m, size=n, prob=mu_psi, rho=rho_psi),
            nrow=j, ncol=m)
    mode(k) <- 'integer'
    
    
    #
    # Create FRASER data set
    #
    sampleIDs <- paste0("sample", seq_len(m))
    anno <- data.table(sampleID = sampleIDs, bamFile=rep(NA, m))
    fds <- FraserDataSet(colData=anno, ...)
    
    # put in n as rawcountsJ first so it doesn't complain later
    # when assinging k to it
    junctionData <- SummarizedExperiment(
        colData=colData(fds),
        assays=list(rawCountsJ=k),
        rowRanges=GRanges(seqnames=rep("chr1", j),
                ranges=IRanges(start=seq_len(j), width=1 )))
    
    nonSplitData <- SummarizedExperiment(
        colData=colData(fds),
        assays=list(rawCountsSS=nonSplit),
        rowRanges=GRanges(seqnames=rep("chr1", j),
                ranges=IRanges(start=seq_len(j), width=1 )))
    
    fds <- new("FraserDataSet",
            junctionData,
            name            = name(fds),
            bamParam        = scanBamParam(fds),
            workingDir      = workingDir(fds),
            nonSplicedReads = nonSplitData)
    
    dontWriteHDF5(fds) <- TRUE
    metadata(fds)[['optimalEncDim']]  <- q
    metadata(fds)[['encDimTable']]    <- data.table(
            encodingDimension=q, evaluationLoss=1, evalMethod='simulation')
    
    # Store "other" counts
    counts(fds, type="psi3", side="other", withDimnames=FALSE) <-    n - k
    counts(fds, type="psi5", side="other", withDimnames=FALSE) <-    n - k
    counts(fds, type="theta", side="other", withDimnames=FALSE) <- n
    
    # store information about the simulation in the fds
    # (same values for psi3 and psi5)
    for(type in c("psi3", "psi5")){
        setAssayMatrix(fds=fds, name="truePSI", type=type, 
                withDimnames=FALSE) <- mu_psi
        setAssayMatrix(fds=fds, name="trueLogitPSI", type=type,
                withDimnames=FALSE) <- y_psi
        mcols(fds, type=type)[,"trueRho"] <- rho_psi
        
        # needed so that subsetting the fds works later
        mcols(fds, type=type)[["startID"]] <- 
            seq_len(nrow(mcols(fds, type=type)))
        mcols(fds, type=type)[["endID"]]   <- 
            seq_len(nrow(mcols(fds, type=type)))
    }
    
    # store info for SE
    setAssayMatrix(fds=fds, name="truePSI", type="theta", 
            withDimnames=FALSE) <- mu_se
    setAssayMatrix(fds=fds, name="trueLogitPSI", type="theta", 
            withDimnames=FALSE) <- y_se
    mcols(fds, type="theta")[,"trueRho"] <- rho_se
    
    # needed so that subsetting the fds works later
    siteIDs <- seq_row(mcols(fds, type="theta"))
    mcols(fds, type="theta")[["startID"]]      <- siteIDs
    mcols(fds, type="theta")[["endID"]]        <- siteIDs
    mcols(fds, type="theta")[["spliceSiteID"]] <- siteIDs
    
    return(fds)
    
}

#
# Generate a simulated fds using a Dirichlet-Mulitnomial distribution
#
makeSimulatedFraserDataSet_Multinomial <- function(m=200, j=1000, q=10,
                groups=round(j*0.65), nbMeanLog=8, nbSize=0.55,
                rhoMeanlog=-5, rhoSdlog=1, Hmean=0, Hsd=100, Dmean=0,
                Dsd=0.007, ...){
    
    #
    # Simulate groups of junctions having the same donor/acceptor site
    #
    d <- groups                  # nr of donors/acceptors
    donorGroups <- seq_len(d)    # ids of the groups (1,2,...,d)
    # assign junction to donor/acceptor groups
    junctionGroups <- sort(sample(x = donorGroups, size=j, replace = TRUE))   
    
    # Set d to the actual number of groups (some might not have been assigned 
    # any junction at all)
    donorGroups <- donorGroups[donorGroups %in% unique(junctionGroups)]
    d <- length(donorGroups)
    
    #
    # SE simulated as additional junction for every group
    #
    nrJunctions <- j                      # store for later
    psiJunctionGroups <- junctionGroups   # store for later
    j <- j + d                            # add one junction for every group
    junctionGroups <- sort(c(junctionGroups, donorGroups))
    
    
    #
    # Simulate n for each junction using lognormal and negative binomial
    #
    junction_n <- rlnorm(j, meanlog=nbMeanLog, sdlog=nbSize)
    n <- t(vapply(junction_n, function(x){
        # draw counts for every sample from negative binomial
        n <- rnbinom(m, mu=x, size=nbSize)},
        double(m)))
    # n <- matrix(rnbinom(j*m, mu=nMean, size=theta), nrow=j, ncol=m)
    
    # Sum up the values of n within one donor/acceptor group to get the n value 
    # of the group
    dt_n <- as.data.table(n)
    sum_n <- cbind(dt_n, junctionGroups)
    res <- lapply(colnames(dt_n),function(x){
        sum_n[,c(paste0(x, "_sum")):=sum(get(x)),by=junctionGroups]
    })
    
    # Extract final n matrix (all junctions within a group have the same n 
    # value)
    sum_cols <- paste0(colnames(dt_n), "_sum")
    n <- as.matrix(sum_n[,..sum_cols])
    colnames(n) <- NULL
    mode(n) <- 'integer'
    
    
    
    #
    # Simulate betaBin dispersion for every donor/acceptor group, same 
    # dispersion for all junctions within a group
    #
    group_rho <- rlnorm(d, meanlog=rhoMeanlog, sdlog=rhoSdlog)# group dispersion
    # group_rho <- abs(rnorm(d, mean=0.0001, sd=0.05))     # group dispersion
    dt1 <- data.table(group=donorGroups, rho=group_rho)
    dt2 <- data.table(group=junctionGroups)
    rhoFull <- merge(dt1, dt2, by="group")
    rho <- rhoFull[,rho]                    # set dispersion for every junction
    
    
    
    #
    # Simulate covariates.
    #
    H_true <- matrix(rnorm(m*q, mean=Hmean, sd=Hsd), nrow=m, ncol=q)
    D_true <- matrix(rnorm(j*q, mean=Dmean, sd=Dsd), nrow=j, ncol=q)
    b_true <- double(j)
    y_true <- t(predictYCpp(H_true, D_true, b_true))
    ae_mu     <- predictMuCpp(y_true)
    
    # Use softmax on mu to get mu within one group to sum up to 1
    softmax <- function(x){
        sum <- sum(exp(x))
        return( exp(x) / sum)
    }
    softmax_mu <- lapply(donorGroups, function(groupID){
        group_mu <- ae_mu[which(junctionGroups == groupID),]
        return( apply(group_mu, 2, softmax) )
    })
    mu <- do.call(rbind, softmax_mu)
    
    # Calcluate alpha values
    alpha <- mu * (1-rho)/rho
    
    
    
    #
    # Simulate psi3=psi5 and SE jointly with Dirchilet-Multinomial
    #
    res <- vapply(donorGroups, function(groupID){
        
        # Get the indices of the junctions in this group
        pos <- which(junctionGroups == groupID)
        
        # get alpha and n values for this group
        group_alpha <- alpha[pos,]
        group_n   <- n[pos,]
        
        # draw counts from a dirichlet multinomial distribution 
        # (jointly for PSi and SE)
        counts <- t(rdirmnom(m, group_n[1,], t(group_alpha)))
        
        # First row of resulting counts represents non split reads for SE, 
        # the other rows are the counts for the "real" junctions
        nonSplit <- counts[1,]
        k <- counts[-1,]
        
        # Substract nonSplit reads from total N to get n for PSI (=other for SE)
        group_n   <- group_n - matrix(nonSplit, nrow=nrow(group_n), 
                ncol=ncol(group_n), byrow = TRUE)
        
        # Also get mu and rho for this group (to split it into PSI and SE part 
        # for storing later)
        group_mu  <- mu[pos,]
        group_rho <- rho[pos]
        
        # First row/value is the one for SE
        se_n     <- group_n[1,]
        se_mu    <- group_mu[1,]
        se_alpha <- group_alpha[1,]
        se_rho   <- group_rho[1]
        
        # The other rows/values are the ones relevant for PSI
        psi_n     <- group_n[-1,]
        psi_alpha <- group_alpha[-1,]
        psi_rho   <- group_rho[-1]
        psi_mu <- group_mu[-1,]
        if(is.null(nrow(psi_mu))){  # only one junction in the group, convert 
            # to matrix with one row so that colSums 
            # works
            psi_mu <- matrix(psi_mu, nrow=1)
        }
        # scale mu's so that they sum up to 1 again (after mu for SE is removed)
        psi_mu    <- psi_mu /  matrix(colSums(psi_mu), nrow=nrow(psi_mu), 
                ncol=ncol(psi_mu), byrow = TRUE)
        
        # return everyting relevant for storing counts, n, alpha, mu, 
        # rho split into the parts for PSI and SE
        return(list(k=k, nonSplit=nonSplit, n=psi_n, group_n=se_n, 
                psi_mu=psi_mu, se_mu=se_mu, psi_alpha=psi_alpha, 
                se_alpha=se_alpha, psi_rho=psi_rho, se_rho=se_rho))
    }, FUN.VALUE=list(k=numeric(m), nonSplit=numeric(m), n=numeric(m), 
            group_n=numeric(m), psi_mu=numeric(m), 
            se_mu=numeric(m), psi_alpha=numeric(m), 
            se_alpha=numeric(m), psi_rho=numeric(1), 
            se_rho=numeric(1)))
    
    # Extract k and nonSplit reads
    k <- do.call(rbind, res[1,])
    nonSplit <- do.call(rbind, res[2,])
    mode(k) <- 'integer'
    mode(nonSplit) <- 'integer'
    
    # Extract n value for each junction and group n value for SE
    psi_n <- do.call(rbind, res[3,])
    se_n <- do.call(rbind, res[4,])
    
    # Extract mu values for PSI and SE
    psi_mu <- do.call(rbind, res[5,])
    se_mu <- do.call(rbind, res[6,])
    
    # Extract alpha values for PSI and SE
    psi_alpha <- do.call(rbind, res[7,])
    se_alpha <- do.call(rbind, res[8,])
    
    # Extract rho values for PSI and SE
    psi_rho <- unlist(res[9,])
    se_rho <- unlist(res[10,])
    
    
    
    #
    # Create FRASER data set
    #
    sampleIDs <- paste0("sample", seq_len(m))
    anno <- data.table(sampleID = sampleIDs, bamFile=rep(NA, m))
    fds <- FraserDataSet(colData=anno, ...)
    
    # put in k as rawcountsJ (= nr reads spanning each junction)
    junctionData <- SummarizedExperiment(
        colData=colData(fds),
        assays=list(rawCountsJ=k),
        rowRanges=GRanges(seqnames=rep("chr1", nrJunctions), 
                ranges=IRanges(start=psiJunctionGroups, width=1)))
    
    nonSplitData <- SummarizedExperiment(
        colData=colData(fds),
        assays=list(rawCountsSS=nonSplit),
        rowRanges=GRanges(seqnames=rep("chr1", d), 
                ranges=IRanges(start=donorGroups, width=1)))
    
    fds <- new("FraserDataSet",
            junctionData,
            name            = name(fds),
            bamParam        = scanBamParam(fds),
            workingDir      = workingDir(fds),
            nonSplicedReads = nonSplitData)
    
    dontWriteHDF5(fds) <- TRUE
    metadata(fds)[['optimalEncDim']]  <- q
    
    
    # store other counts (n-k) for psi3=psi5 and SE (other=n)
    counts(fds, type="psi3", side="other", withDimnames=FALSE)    <- psi_n - k
    counts(fds, type="psi5", side="other", withDimnames=FALSE)    <- psi_n - k
    counts(fds, type="theta", side="other", withDimnames=FALSE) <- se_n
    
    for(type in c("psi3", "psi5")){
        # store information about the simulation in the fds 
        # (same values for psi3 and psi5)
        setAssayMatrix(fds=fds, name="truePSI", type=type,
                withDimnames=FALSE) <- psi_mu
        setAssayMatrix(fds=fds, name="trueLogitPSI", type=type,
                withDimnames=FALSE) <- qlogis(psi_mu)
        setAssayMatrix(fds=fds, name="trueAlpha", type=type,
                withDimnames=FALSE) <- psi_alpha
        mcols(fds, type=type)[,"trueRho"] <- psi_rho
        
        # needed so that subsetting the fds works later
        mcols(fds, type=type)[["startID"]] <- psiJunctionGroups
        mcols(fds, type=type)[["endID"]]   <- psiJunctionGroups
    }
    
    # store info for SE
    setAssayMatrix(fds=fds, name="truePSI", type="theta",
            withDimnames=FALSE) <- se_mu
    setAssayMatrix(fds=fds, name="trueLogitPSI", type="theta",
            withDimnames=FALSE) <- qlogis(se_mu)
    setAssayMatrix(fds=fds, name="trueAlpha", type="theta",
            withDimnames=FALSE) <- se_alpha
    mcols(fds, type="theta")[,"trueRho"] <- se_rho
    
    # needed so that subsetting the fds works later
    mcols(fds, type="theta")[["startID"]]      <- donorGroups
    mcols(fds, type="theta")[["endID"]]        <- donorGroups
    mcols(fds, type="theta")[["spliceSiteID"]] <- donorGroups
    
    return(fds)
    
}

#'
#' Inject artificial outliers in an existing fds
#' 
#' @param fds FraserDataSet
#' @param type The psi type
#' @param freq The injection frequency.
#' @param minDpsi The minimal delta psi with which outliers will be injected.
#' @param minCoverage The minimal total coverage (i.e. N) required for a 
#' junction to be considered for injection of an outlier.
#' @param deltaDistr The distribution from which the delta psi value of 
#' the injections is drawn (default: uniform distribution).
#' @param verbose Should additional information be printed during computation?
#' @param method Defines by which method the new psi of injections is computed, 
#' i.e. to which value the delta psi of the injection is added: "meanPSI" for 
#' adding to the mean psi of the junction over all samples or "samplePSI" to 
#' add to the psi value of the junction in the specific sample. "simulatedPSI" 
#' is only possible if a simulated dataset is used.
#' @param BPPARAM A BiocParallel object to run the computation in parallel
#' 
#' @return FraserDataSet
#' 
#' @examples 
#' # A generic dataset
#' fds <- makeSimulatedFraserDataSet()
#' fds <- calculatePSIValues(fds)
#' fds <- injectOutliers(fds, minDpsi=0.2, freq=1E-3)
#' @export
injectOutliers <- function(fds, type=psiTypes,
                    freq=1E-3, minDpsi=0.2, minCoverage=2,
                    deltaDistr="uniformDistr", verbose=FALSE,
                    method=c('samplePSI', 'meanPSI', 'simulatedPSI'),
                    BPPARAM=bpparam()){
    type <- match.arg(type, several.ok=TRUE)
    method <- match.arg(method)
    
    if(length(type) > 1){
        for(t in type){
            fds <- injectOutliers(fds, type=t, freq=freq, minDpsi=minDpsi,
                    minCoverage=minCoverage, deltaDistr=deltaDistr,
                    verbose=verbose, method=method, BPPARAM=BPPARAM)
        }
        return(fds)
    }
    # copy original k and o
    if(type == "theta"){
        setAssayMatrix(fds, type="theta", "originalCounts",
                withDimnames=FALSE) <- 
                        counts(fds, type="theta", side="ofInterest")
        setAssayMatrix(fds, type="theta", "originalOtherCounts",
                withDimnames=FALSE) <- 
                        counts(fds, type="theta", side="other")
    }
    else{
        setAssayMatrix(fds, type=type, "originalCounts",
                withDimnames=FALSE) <- 
                        counts(fds, type=type, side="ofInterest")
        setAssayMatrix(fds, type="psi5", "originalOtherCounts",
                withDimnames=FALSE) <- 
                        counts(fds, type="psi5", side="other")
        setAssayMatrix(fds, type="psi3", "originalOtherCounts",
                withDimnames=FALSE) <- 
                        counts(fds, type="psi3", side="other")
        setAssayMatrix(fds, type="jaccard", "originalOtherCounts",
                       withDimnames=FALSE) <- 
                        counts(fds, type="jaccard", side="other")
    }

    # get infos from the fds
    m <- ncol(fds)
    j <- nrow(mcols(fds, type=type))

    k <- as.matrix(K(fds, type=type))
    n <- as.matrix(N(fds, type=type))
    o <- as.matrix(counts(fds, type=type, side="other"))

    psi <- switch(match.arg(method),
            samplePSI    = (k + pseudocount())/(n + 2*pseudocount()),
            meanPSI      = matrix(nrow=j, ncol=m,
                    rowMeans( (k + pseudocount()) / (n + 2*pseudocount()))),
            simulatedPSI = getAssayMatrix(fds, "truePSI", type=type) )

    # nedded to have comparably many outliers when type is SE
    if(type == "theta"){
        freq <- freq/10
    }

    # data table with information for calculating groups
    dt <- data.table(
            junctionID = seq_len(j), 
            groupID = getSiteIndex(fds, type))
    dt[,groupSize:=.N, by=groupID]

    # Get groups where outlier can be injected
    available_groups <- dt[groupSize > ifelse(type == "theta" | 
                                                type == "jaccard", 0, 1), 
                            unique(groupID)]
    
    # e.g. for psi3/5: no donor/acceptor 
    # groups with at least 2 junctions (e.g in simulationBB)
    if(length(available_groups) == 0){ 
        available_groups <- seq_len(j)
        freq <- freq/10
    }
    
    list_index <- data.frame(row=integer(0), col=integer(0))
    count <- 0
    while(nrow(list_index) < 1){
        if(count == 20){
            stop("Could not inject at least 2 outliers.", 
                    " Please make sure you have enough junctions and samples",
                    " in you dataset.")
        }
        count <- count + 1
        
        indexOut_groups <- matrix(sample(c(0,1,-1), length(available_groups)*m, 
                replace=TRUE, prob=c(1-freq, freq/2, freq/2)), ncol=m)
    
        # positions where outliers will be injected
        list_index <- which(indexOut_groups != 0, arr.ind = TRUE)
    }
    
    # apply injection function to each outlier
    message(date(), ": Injecting ", nrow(list_index), " outliers ...")
    result <- bplapply(seq_len(nrow(list_index)), list_index=list_index,
            indexOut_groups=indexOut_groups, type=type, psi=psi, n=n, dt=dt,
            minDpsi=minDpsi, verbose=verbose, BPPARAM=SerialParam(),
            FUN=function(j, list_index, indexOut_groups, type,
                            psi, n, dt=dt, minDpsi, verbose){
                # extract group, sample and injecetion 
                # direction (i.e +1/up or -1/down)
                row       <- list_index[j,'row']
                group     <- available_groups[row]
                sample    <- list_index[j,'col']
                injDirection    <- indexOut_groups[row, sample]

                # sample one junction from all junction within this group
                group_junctions <- dt[groupID == group, junctionID]
                junction <- if(length(group_junctions)==1){ 
                            group_junctions } # group }
                        else {
                            sample(group_junctions, 1) }

                # get current psi of this junction and calculate
                # maximal possible value of delta psi for the injection
                junction_psi    <- psi[junction, sample]
                maxDpsi <- if(injDirection > 0){ 
                            1 - junction_psi }
                        else{
                            junction_psi }

                # if not possible to inject here with at least minDpsi,
                # change injection direction
                if(maxDpsi < minDpsi){
                    injDirection <- -injDirection
                    indexOut_groups[row, sample] <- injDirection
                    maxDpsi <- if(injDirection > 0){ 
                                1 - junction_psi }
                            else{
                                junction_psi }
                }

                # sample delta psi for injection from uniform
                # distribution between min and max dpsi
                minDpsi <- ifelse(minDpsi < maxDpsi, minDpsi, maxDpsi)
                injDpsi <- injDirection * switch(as.character(deltaDistr),
                        uniformDistr = runif(1, minDpsi, maxDpsi),
                        ifelse(as.double(deltaDistr) > maxDpsi, 
                                maxDpsi, as.double(deltaDistr)) )

                # get N of this junction
                n_ji <- n[junction,sample]

                # new counts after injection
                newKs     <- integer(length(group_junctions))
                indexDpsi <- double(length(group_junctions))
                indOut  <- integer(length(group_junctions))

                # for all other junctions in this group
                for(i in seq_len(length(group_junctions))){
                    junction_k <- group_junctions[i]
                    # get new_psi and change k based on n and new_psi 
                    # (and take pseudocounts into account): 
                    #     (k+1)/(n+2)=psi -> k = psi*(n+2) - 1
                    if(junction_k == junction){
                        new_psi <- junction_psi + injDpsi
                        new_k <- round( new_psi * 
                                (n_ji + 2*pseudocount()) - pseudocount() )

                        # store position of outlier
                        indOut[i]      <- injDirection
                        indexDpsi[i]   <- injDpsi
                    }
                    else{
                        deltaPSI_k <- - (psi[junction_k,sample] / 
                                (1-junction_psi)) * injDpsi
                        new_psi <- psi[junction_k,sample] + deltaPSI_k
                        new_k <- round( new_psi*(n_ji + 2*pseudocount()) - 
                                pseudocount() )

                        # also store position of other junction used in swap
                        indOut[i]      <- -injDirection * 2
                        indexDpsi[i]   <- deltaPSI_k
                        
                        if(deltaPSI_k > 1 | deltaPSI_k < -1){
                            warning("Calculated a injected |deltaPSI| > 1!")
                        }
                    }
                    # for SE: ensure new_k <= n_ij 
                    # (so that o=n-k is always >= 0) 
                    # (not needed for psi3/5 because o
                    # will recalculated from k's there)
                    if(length(group_junctions)==1){ # if(type == "theta"){
                        new_k <- min(new_k, n_ji)
                    }
                    # ensure new_k >= 0 and assign k_ij <- new_k
                    newKs[i] <- max(0, new_k)
                }

                if(verbose){
                    print(paste("injected outlier", j, "with delta PSI of", 
                            injDpsi, "at junction", junction, "and sample", 
                            sample))
                }

                return(list(newKs = newKs, newOs=sum(newKs)-newKs, 
                        junctions = group_junctions, injDeltaPSI = indexDpsi,
                        injDirections = indOut, 
                        sample = rep(sample, length(group_junctions))))
            })

    # get all junctions, samples, ... where outlier injection changed the counts
    junctions     <- unlist(lapply(result, "[[", 'junctions'))
    samples       <- unlist(lapply(result, "[[", 'sample'))
    newKs         <- unlist(lapply(result, "[[", 'newKs'))
    newOs         <- unlist(lapply(result, "[[", 'newOs'))
    injDirection  <- unlist(lapply(result, "[[", 'injDirections'))
    injDeltaPSI   <- unlist(lapply(result, "[[", 'injDeltaPSI'))

    # matrices to store indices of outliers and their delta psi
    indexOut <- matrix(0, nrow = j, ncol = m)
    indexDeltaPSI <- matrix(0, nrow = j, ncol = m)
    
    # set counts to changed counts after the injection
    replaceIndices                <- matrix(c(junctions,samples), ncol=2)
    k[replaceIndices]             <- newKs
    if(length(available_groups) == j){
        o                           <- n - k
    } else{
        o[replaceIndices]           <- newOs
    }
    indexOut[replaceIndices]      <- injDirection
    indexDeltaPSI[replaceIndices] <- injDeltaPSI

    # store positions of true outliers and their true delta PSI value
    setAssayMatrix(fds=fds, name="trueOutliers", type=type,
            withDimnames=FALSE) <- indexOut
    setAssayMatrix(fds=fds, name="trueDeltaPSI", type=type,
            withDimnames=FALSE) <- indexDeltaPSI

    # store new k counts which include the outlier counts
    counts(fds, type=type, side="ofInterest", withDimnames=FALSE) <- k

    # store modified other counts
    counts(fds, type=type, side="other", withDimnames=FALSE) <- o

    return(fds)
}

removeInjectedOutliers <- function(fds, type){
    
    # copy injected k and o counts
    if(type == "theta"){
        setAssayMatrix(fds, type="theta", "outlierCounts",
            withDimnames=FALSE) <- counts(fds, type="theta", side="ofInteres")
        setAssayMatrix(fds, type="theta", "outlierOtherCounts",
            withDimnames=FALSE) <- counts(fds, type="theta", side="other")
    }
    else{
        setAssayMatrix(fds, type="psi5", "outlierCounts",
            withDimnames=FALSE) <- counts(fds, type="psi5", side="ofInterest")
        setAssayMatrix(fds, type="psi5", "outlierOtherCounts",
            withDimnames=FALSE) <- counts(fds, type="psi5", side="other")
        setAssayMatrix(fds, type="psi3", "outlierOtherCounts", 
            withDimnames=FALSE) <- counts(fds, type="psi3", side="other")
    }
    
    # assign original k and o to rawCountsJ and rawOtherCounts
    if(type == "theta"){
        counts(fds, type="theta", side="ofInterest", withDimnames=FALSE) <- 
                getAssayMatrix(fds, type="theta", "originalCounts")
        counts(fds, type="theta", side="other", withDimnames=FALSE) <- 
                getAssayMatrix(fds, type="theta", "originalOtherCounts")
        
        assays(fds)[['originalCounts_theta']] <- NULL
        assays(fds)[['originalOtherCounts_theta']] <- NULL
    }
    else{
        counts(fds, type="psi5", side="ofInterest", withDimnames=FALSE) <- 
                getAssayMatrix(fds, type="psi5", "originalCounts")
        counts(fds, type="psi5", side="other", withDimnames=FALSE) <-
                getAssayMatrix(fds, type="psi5", "originalOtherCounts")
        counts(fds, type="psi3", side="other", withDimnames=FALSE) <-
                getAssayMatrix(fds, type="psi3", "originalOtherCounts")
        
        assays(fds)[['originalCounts_psi5']] <- NULL
        assays(fds)[['originalOtherCounts_psi5']] <- NULL
        assays(fds)[['originalOtherCounts_psi3']] <- NULL
    }
    
    return(fds)
    
}

restoreInjectedOutliers <- function(fds, type){
    
    # copy original k and o counts
    if(type == "theta"){
        setAssayMatrix(fds, type="theta", "originalCounts") <- 
            counts(fds, type="theta", side="ofInterest")
        setAssayMatrix(fds, type="theta", "originalOtherCounts") <- 
            counts(fds, type="theta", side="other")
    }
    else{
        setAssayMatrix(fds, type="psi5", "originalCounts") <- 
            counts(fds, type="psi5", side="ofInterest")
        setAssayMatrix(fds, type="psi5", "originalOtherCounts") <- 
            counts(fds, type="psi5", side="other")
        setAssayMatrix(fds, type="psi3", "originalOtherCounts") <- 
            counts(fds, type="psi3", side="other")
    }
    
    # assign injected k and o to rawCountsJ and rawOtherCounts
    if(type == "theta"){
        counts(fds, type="theta", side="ofInterest") <- 
            getAssayMatrix(fds, type="theta", "outlierCounts")
        counts(fds, type="theta", side="other") <- 
            getAssayMatrix(fds, type="theta", "outlierOtherCounts")
        
        assays(fds)[['outlierCounts_theta']] <- NULL
        assays(fds)[['outlierOtherCounts_theta']] <- NULL
    }
    else{
        counts(fds, type="psi5", side="ofInterest") <- 
            getAssayMatrix(fds, type="psi5", "outlierCounts")
        counts(fds, type="psi5", side="other") <- 
            getAssayMatrix(fds, type="psi5", "outlierOtherCounts")
        counts(fds, type="psi3", side="other") <-
            getAssayMatrix(fds, type="psi3", "outlierOtherCounts")
        
        assays(fds)[['outlierCounts_psi5']] <- NULL
        assays(fds)[['outlierOtherCounts_psi5']] <- NULL
        assays(fds)[['outlierOtherCounts_psi3']] <- NULL
    }
    
    return(fds)
    
}


# log2 fold change of measured psi vs AE predicted psi
log2fc <- function(realPsi, predictedPsi){
    return( log2(realPsi) - log2(predictedPsi) )
}
