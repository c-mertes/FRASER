#'
#' Main autoencoder fit function
#'
#' @noRd
fitAutoencoder <- function(fds, q, type=currentType(fds), noiseAlpha=1, 
                    minDeltaPsi=0.1,
                    rhoRange=c(-30, 30), lambda=0, convergence=1e-5,
                    iterations=15, initialize=TRUE, control=list(),
                    BPPARAM=bpparam(), verbose=FALSE, nrDecoderBatches=5,
                    weighted=FALSE, nSubset=15000, multiRho=FALSE,
                    latentSpace=c('AE', 'PCA')){

    if(!bpisup(BPPARAM)){
        bpstart(BPPARAM)
    }
    dims <- c(row=nrow(mcols(fds, type=type)), col=ncol(fds))

    # set alpha for noise injection for denoising AE
    # init it corresponding to the size of 'x' not 'fds'
    currentNoiseAlpha(fds) <- noiseAlpha
    if(!is.null(noiseAlpha)){
        noise(fds, type=type) <- matrix(rnorm(prod(dims)), nrow=dims['col'])
    }

    # make sure its only in-memory data for k and n
    currentType(fds) <- type
    counts(fds, type=type, side="other", HDF5=FALSE) <- 
        as.matrix(N(fds) - K(fds))
    counts(fds, type=type, side="ofInterest", HDF5=FALSE) <- as.matrix(K(fds))

    # copy fds object to save original input object
    # and create second object with only the subset to fit the encoder
    copy_fds <- fds
    exMask <- variableJunctions(fds, type, minDeltaPsi=minDeltaPsi)
    fds <- fds[exMask,,by=type]
    exMask2 <- subsetKMostVariableJunctions(fds, type, nSubset)
    fds <- fds[exMask2,,by=type]
    featureExclusionMask(fds) <- TRUE

    # set correct exclusion mask for x computation
    exMask[exMask == TRUE] <- exMask2
    featureExclusionMask(copy_fds) <- exMask

    # initialize E and D using PCA and bias as zeros.
    if(isTRUE(initialize) | is.null(E(fds)) | is.null(D(fds))){
        fds <- initAutoencoder(fds, q, rhoRange, type=type, BPPARAM=BPPARAM)
    }

    # initial loss
    lossList <- lossED(fds, lambda, byRows=TRUE)
    colnames(lossList) <- "init_pca"
    message('Initial PCA loss: ', mean(lossList[,1]))

    if(match.arg(latentSpace) == 'AE'){

        #inital rho values
        if(isTRUE(verbose)){
            print(summary(rho(fds)))
        }

        # Dont use batche runs in E fitting state
        batches4EFit <- nrDecoderBatches
        if(nrow(mcols(fds, type=type)) < nrow(mcols(copy_fds, type=type))){
            batches4EFit <- 1
        }

        # initialize D
        fds <- updateD(fds, type=type, lambda=lambda, control=control,
                        BPPARAM=BPPARAM, verbose=verbose, 
                        nrDecoderBatches=batches4EFit,
                        multiRho=multiRho, weighted=FALSE)
        lossList <- updateLossList(fds, lossList, 'init', 'D', lambda, 
                                    verbose=verbose)

        # initialize rho step
        fds <- updateRho(fds, type=type, rhoRange, BPPARAM=BPPARAM, 
                            verbose=verbose)
        lossList <- updateLossList(fds, lossList, 'init', 'Rho', lambda, 
                                    verbose=verbose)

        # optimize log likelihood
        t1 <- Sys.time()
        currentLoss <- lossED(fds, lambda, byRows=TRUE)
        for(i in seq_len(iterations)){
            t2 <- Sys.time()

            # update E step
            fds <- updateE(fds, control=control, BPPARAM=BPPARAM, 
                            verbose=verbose)
            lossList <- updateLossList(fds, lossList, i, 'E', lambda, 
                                        verbose=verbose)

            # update D step
            fds <- updateD(fds, type=type, lambda=lambda, control=control,
                            BPPARAM=BPPARAM, verbose=verbose, 
                            nrDecoderBatches=batches4EFit,
                            multiRho=multiRho, weighted=FALSE)
            lossList <- updateLossList(fds, lossList, i, 'D', lambda, 
                                        verbose=verbose)

            # update rho step
            fds <- updateRho(fds, type=type, rhoRange, BPPARAM=BPPARAM, 
                            verbose=verbose)
            lossList <- updateLossList(fds, lossList, i, 'Rho', lambda, 
                                        verbose=verbose)

            if(isTRUE(verbose)){
                print(paste('Time for one autoencoder loop:', Sys.time() - t2))
            } else {
                print(paste0(date(), ': Iteration: ', i, ' loss: ',
                            mean(lossList[,ncol(lossList)])))
            }

            # check
            curLossDiff <- rowMax(abs(
                matrix(currentLoss, ncol=3, nrow=length(currentLoss))
                - lossList[,ncol(lossList) - 2:0]))
            if(all(max(curLossDiff) < convergence)){
                message(date(), ': the AE correction converged with: ',
                        mean(lossList[,ncol(lossList)]))
                break
            } else {
                if(isTRUE(verbose)){
                    message(date(), ": Current max diff is: ", max(curLossDiff))
                    message(date(), ": Summary: ", 
                            paste(collapse=", ", sep=": ",
                                    names(summary(curLossDiff)),
                                    signif(summary(curLossDiff), 2)))
                }
            }
            currentLoss <- lossList[,ncol(lossList)]
        }

        print(Sys.time() - t1)
    }

    # TODO when using all features (!=nrow(fds)) in fitting for SE: also stop 
    # here and set copy_fds to fitted fds but with all junctions
    # (at the moment: fds contains only a subset of all junctions, 
    # so copy_fds <- fds doesn't work for SE)
    if(nrow(fds) == nrow(copy_fds) & !isTRUE(weighted)){    
        copy_fds <- fds
    } else {
        # update the D matrix and theta
        print(paste0("Finished with fitting the E matrix. Starting now with ",
                    "the D fit. ..."))

        # set noiseAlpha to 0 or NULL to NOT use noise now 
        # (latent space already fitted)
        currentNoiseAlpha(fds) <- NULL

        copy_fds <- initAutoencoder(copy_fds, q, rhoRange, type=type, 
                                    BPPARAM=BPPARAM)
        E(copy_fds) <- E(fds)
        currentLoss <- lossED(copy_fds, lambda, byRows=TRUE)

        # adapt loss list to full matrix
        newLossList <- matrix(NA_real_, ncol=ncol(lossList), nrow=dims["row"])
        newLossList[featureExclusionMask(copy_fds),] <- lossList
        colnames(newLossList) <- colnames(lossList)
        lossList <- newLossList

        for(i in seq_len(iterations)){
            t2 <- Sys.time()

            # update D step
            copy_fds <- updateD(copy_fds, type=type, lambda=lambda, 
                                control=control, BPPARAM=BPPARAM, 
                                verbose=verbose, 
                                nrDecoderBatches=nrDecoderBatches, 
                                multiRho=multiRho, weighted=weighted)
            lossList <- updateLossList(copy_fds, lossList, paste0("final_", i), 
                                        'D', lambda, verbose=verbose)

            # update rho step
            copy_fds <- updateRho(copy_fds, type=type, rhoRange, 
                                    BPPARAM=BPPARAM, verbose=verbose)
            lossList <- updateLossList(copy_fds, lossList, paste0("final_", i), 
                                        'Rho', lambda, verbose=verbose)

            if(isTRUE(verbose)){
                print(paste('Time for one D & Rho loop:', Sys.time() - t2))
            } else {
                print(paste0(date(), ': Iteration: final_', i, ' loss: ',
                            mean(lossList[,ncol(lossList)]), " (mean); ",
                            max(lossList[,ncol(lossList)]), " (max)"))
            }

            # check
            curLossDiff <- rowMax(abs(
                    matrix(currentLoss, ncol=2, nrow=length(currentLoss))
                    - lossList[,ncol(lossList) - c(1,0)]))
            if(all(curLossDiff < convergence)){
                message(date(), ': the final AE correction converged with:',
                        mean(lossList[,ncol(lossList)]))
                break
            } else {
                if(isTRUE(verbose)){
                    message(date(), ": Current max diff is: ", max(curLossDiff))
                }
            }
            currentLoss <- lossList[, ncol(lossList)]
        }
    }

    print(paste0(i, ' Final betabin-AE loss: ', 
                mean(lossList[, ncol(lossList)])))
    bpstop(BPPARAM)

    # add additional values for the user to the object
    metadata(copy_fds)[[paste0('dim_', type)]] <- dim(copy_fds)
    metadata(copy_fds)[[paste0('loss_', type)]] <- lossList
    metadata(copy_fds)[[paste0('convList_', type)]] <- lossList


    # add correction factors
    predictedMeans <- t(predictMu(copy_fds))
    stopifnot(identical(dim(K(copy_fds)), dim(predictedMeans)))
    predictedMeans(copy_fds, type, withDimnames=FALSE) <- predictedMeans

    # store weights if weighted version
    if(isTRUE(weighted)){
        weights <- weights(copy_fds, type)
        weights(copy_fds, type, withDimnames=FALSE) <- weights
    }

    # validate object
    validObject(copy_fds)
    return(copy_fds)
}

initAutoencoder <- function(fds, q, rhoRange, type, BPPARAM){

    x <- x(fds, all=TRUE, center=FALSE)
    pca <- pca(as.matrix(x(fds, all=TRUE)), nPcs=q, center=FALSE)
    pc  <- pcaMethods::loadings(pca)

    # Set initial values from PCA
    D(fds) <- pc
    E(fds) <- pc[featureExclusionMask(fds),]
    b(fds) <- colMeans2(x)

    # initialize rho
    # rho(fds) <- methodOfMomentsRho(K(fds), N(fds))
    fds <- updateRho(fds, type=type, rhoRange, BPPARAM=BPPARAM, verbose=FALSE)

    # reset counters
    mcols(fds, type=type)[paste0('NumConvergedD_', type)] <- 0

    return(fds)
}

updateLossList <- function(fds, lossList, i, stepText, lambda, verbose){
    currLoss <- lossED(fds, lambda, byRows=TRUE)
    lossList <- cbind(lossList, currLoss)
    colnames(lossList)[ncol(lossList)] <- paste0(i, '_', stepText)
    if(isTRUE(verbose)){
        print(paste0(date(), ': Iteration: ', i, ' ',
                stepText, ' loss: ',
                round(mean(currLoss), 7), " (mean); ",
                round(max(currLoss), 7),  " (max)"))
    }
    return(lossList)
}

lossED <- function(fds, lambda=0, byRows=FALSE,
                    noiseAlpha=currentNoiseAlpha(fds)){
    K <- K(fds)
    N <- N(fds)
    rho <- matrix(rho(fds), ncol=ncol(K), nrow = nrow(K))
    D <- D(fds)

    y <- predictY(fds, noiseAlpha=noiseAlpha)

    return(fullNLL(y, rho, as.matrix(K), as.matrix(N), D, lambda, 
                    byRows=byRows))
}

# lossED <- function(fds){
#   K <- as.matrix(K(fds))
#   N <- as.matrix(N(fds))
#   mu <- predictMu(fds)
#   rho <- matrix(rho(fds), ncol=ncol(fds), nrow = nrow(fds))
#
#   r  <- (1-rho)/rho
#   eps <- 0.5
#   alpha  <- lgamma(mu*r)
#   alphaK <- lgamma(mu*r + K + eps)
#   beta   <- lgamma((mu-1)*(-r))
#   betaNK <- lgamma((mu-1)*(-r) + (N - K + eps))
#
#   #mean negative log likelihood with pseudocounts
#   mean(alpha + beta - alphaK - betaNK - lgamma(N+1+2*eps) + lgamma(K+1+eps) + 
#          lgamma(N-K+1+eps) + lgamma(r + N + 2*eps) - lgamma(r))
# }
