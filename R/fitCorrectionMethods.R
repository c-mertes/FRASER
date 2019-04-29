
fit <- function(fds, correction=c("FraseR", "fullFraseR", "PCA", "PEER", "PEERdecoder", "BB"), q, type="psi3", rhoRange=c(1e-5, 1-1e-5), 
                  noiseAlpha=1, lambda=0, convergence=1e-5, iterations=15, initialize=TRUE,
                  control=list(), BPPARAM=bpparam(), verbose=FALSE, 
                  recommendedPEERq=TRUE){
  
  method <- match.arg(correction)
  
  fds <- switch(method,
         FraseR = fitFraserAE(fds=fds, q=q, type=type, noiseAlpha=noiseAlpha, rhoRange=rhoRange, lambda=lambda,
                               convergence=convergence, iterations=iterations, initialize=initialize,
                               control=control, BPPARAM=BPPARAM, verbose=verbose, subset=TRUE),
         fullFraseR = fitFraserAE(fds=fds, q=q, type=type, noiseAlpha=noiseAlpha, rhoRange=rhoRange, lambda=lambda,
                                convergence=convergence, iterations=iterations, initialize=initialize,
                                control=control, BPPARAM=BPPARAM, verbose=verbose, subset=FALSE),
         PCA    = fitPCA(fds=fds, q=q, psiType=type, rhoRange=rhoRange, BPPARAM=BPPARAM),
         PEER   = fitPEER(fds=fds, q=q, psiType=type, recomendedQ=recommendedPEERq, rhoRange=rhoRange, BPPARAM=BPPARAM),
         PEERdecoder = fitPEERDecoder(fds=fds, q=q, psiType=type, recomendedQ=recommendedPEERq, rhoRange=rhoRange, BPPARAM=BPPARAM),
         BB     = fitBB(fds=fds, psiType=type) )
  
  return(fds)
  
}

fitPCA <- function(fds, q, psiType, rhoRange=c(1e-5, 1-1e-5), BPPARAM=parallel(fds)){
  
  #+ subset fitting
  currentType(fds) <- psiType
  curDims <- dim(K(fds, psiType))
  probE <- max(0.001, min(1,30000/curDims[1]))
  featureExclusionMask(fds) <- sample(c(TRUE, FALSE), curDims[1],
                                          replace=TRUE, prob=c(probE, 1-probE))
  print(table(featureExclusionMask(fds)))
  
  # PCA on subset -> E matrix
  currentNoiseAlpha(fds) <- NULL
  message(date(), " Computing PCA ...")
  pca <- pca(as.matrix(x(fds)), nPcs=q)
  pc <- pcaMethods::loadings(pca)
  E(fds) <- pc
  # linear regression to fit D matrix
  lmFit <- lm(as.matrix(x(fds, all=TRUE, center=FALSE)) ~ as.matrix(H(fds)))
  D(fds) <- t(lmFit$coefficients[-1,])
  b(fds) <- lmFit$coefficients[1,]
  
  # fit rho
  message(date(), " Fitting rho ...")
  counts(fds, type=psiType, side="other", HDF5=FALSE) <- as.matrix(N(fds) - K(fds))
  counts(fds, type=psiType, side="ofInterest", HDF5=FALSE) <- as.matrix(K(fds))
  fds <- updateRho(fds, type=psiType, rhoRange=rhoRange, BPPARAM=BPPARAM, verbose=TRUE) 
  # store corrected logit psi 
  predictedMeans(fds, psiType) <- t(predictMu(fds))
  
  return(fds)
}
  
fitPEER <-function(fds, q, psiType, recomendedQ=TRUE, rhoRange=c(1e-5, 1-1e-5), BPPARAM=parallel(fds)){
  
  # set featureExclusionMask of all junctions to TRUE for peer
  currentType(fds) <- psiType
  featureExclusionMask(fds, type=psiType) <- rep(TRUE, nrow(mcols(fds, type=psiType)))
  
  #+ PEER
  require(peer)
  #+ prepare PEER model
  if(isTRUE(recomendedQ)){
    # recommendation by PEER: min(0.25*n, 100)
    q <- min(as.integer(0.25* ncol(fds)), 100) 
  }
  maxFactors <- q                        # number of known hidden factors
  model <- PEER()
  PEER_setPhenoMean(model, as.matrix(x(fds, all=TRUE, center=FALSE)))
  PEER_setNk(model, maxFactors)          # nr of hidden confounders
  PEER_setNmax_iterations(model, 1000)   # 1000 iterations is default
  # PEER_setAdd_mean(model, TRUE)        # should mean expression be added as additional factor? currently FALSE
  
  #+ run full Peer pipeline
  message(date(), "Fitting PEER model ...")
  PEER_update(model)
  
  #+ extract PEER data
  peerResiduals <- PEER_getResiduals(model)
  peerLogitMu <- t(as.matrix(x(fds, all=TRUE, center=FALSE)) - peerResiduals)
  
  #+ save peer model in fds object
  setAssayMatrix(fds, "peerLogitMu", type=psiType) <- peerLogitMu
  predictedMeans(fds, psiType) <- predictMuCpp(peerLogitMu)
  metadata(fds)[[paste0("PEERmodel_", psiType)]] <- list(
    alpha     = PEER_getAlpha(model),
    residuals = PEER_getResiduals(model),
    W         = PEER_getW(model),
    hiddenSpace = PEER_getX(model))
  
  #+ fit rho
  message(date(), "Fitting rho ...")
  k <- as.matrix(K(fds, psiType))
  n <- as.matrix(N(fds, psiType))
  y <- peerLogitMu
  fitparameters <- bplapply(seq_len(nrow(k)), estRho,
                            k=k, n=n, y=y, rhoRange=rhoRange,
                            BPPARAM=BPPARAM, nll=truncNLL_rho)
  rho(fds) <- vapply(fitparameters, "[[", double(1), "minimum")
  print(summary(rho(fds)))
  
  return(fds)
  
}

fitPEERDecoder <-function(fds, q, psiType, recomendedQ=TRUE, rhoRange=c(1e-5, 1-1e-5), BPPARAM=parallel(fds)){
  
  # set featureExclusionMask of all junctions to TRUE for peer
  currentType(fds) <- psiType
  featureExclusionMask(fds, type=psiType) <- rep(TRUE, nrow(mcols(fds, type=psiType)))
  
  #+ PEER
  require(peer)
  #+ prepare PEER model
  maxFactors <- q                        # number of known hidden factors
  if(isTRUE(recomendedQ)){
    # recommendation by PEER: min(0.25*n, 100)
    maxFactors <- min(as.integer(0.25* ncol(fds)), 100) 
  }
  model <- PEER()
  x <- x(fds, all=TRUE, center=FALSE)
  PEER_setPhenoMean(model, as.matrix(x))
  PEER_setNk(model, maxFactors)          # nr of hidden confounders
  PEER_setNmax_iterations(model, 1000)   # 1000 iterations is default
  # PEER_setAdd_mean(model, TRUE)        # should mean expression be added as additional factor? currently FALSE
  
  #+ run full Peer pipeline
  message(date(), "Fitting PEER model ...")
  PEER_update(model)
  
  #+ extract PEER data
  peerResiduals <- PEER_getResiduals(model)
  peerLogitMu <- t(as.matrix(x(fds, all=TRUE, center=FALSE)) - peerResiduals)
  
  #+ save peer model in fds object
  setAssayMatrix(fds, "peerLogitMu", type=psiType) <- peerLogitMu
  metadata(fds)[[paste0("PEERmodel_", psiType)]] <- list(
    alpha     = PEER_getAlpha(model),
    residuals = PEER_getResiduals(model),
    W         = PEER_getW(model),
    hiddenSpace = PEER_getX(model))
  
  #+ fit rho
  message(date(), "Fitting rho ...")
  k <- as.matrix(K(fds, psiType))
  n <- as.matrix(N(fds, psiType))
  y <- peerLogitMu
  fitparameters <- bplapply(seq_len(nrow(k)), estRho,
                            k=k, n=n, y=y, rhoRange=rhoRange,
                            BPPARAM=BPPARAM, nll=truncNLL_rho)
  rho(fds) <- vapply(fitparameters, "[[", double(1), "minimum")
  print(summary(rho(fds)))
  
  #+ fit decoder on peer hidden space
  D <- matrix(rnorm(prod(dim(fds))), nrow=nrow(fds), ncol=q)
  b <- colMeans2(x)
  rho <- rho(fds)
  alphas <- PEER_getAlpha(model)[,1]
  H <- PEER_getX(model)[,order(alphas, decreasing = TRUE)[1:q]]
  
  # run fits
  fitls <- bplapply(seq_len(nrow(k)), singleDFit, D=D, b=b, k=k, n=n, H=H,
                    rho=rho, lambda=0, control=list(), BPPARAM=parallel(fds))
  
  # extract infos
  parMat <- vapply(fitls, '[[', double(ncol(D) + 1), 'par')
  
  # update b and D
  b(fds) <- parMat[1,]
  D(fds) <- t(parMat)[,-1]
  
  # predict and save (logit) psi
  y <- predictYCpp(as.matrix(H), D(fds), b(fds))
  predictedMeans(fds, psiType) <- predictMuCpp(y)
  
  return(fds)
  
}
  
fitBB <- function(fds, psiType){
  currentType(fds) <- psiType
  fds <- pvalueByBetaBinomialPerType(fds=fds, aname=paste0("pvalue_", psiType), psiType=psiType, pvalFun=betabinVglmTest) 
  
  return(fds)
}

fitFraserAE <- function(fds, q, type, noiseAlpha, rhoRange, lambda, convergence, iterations, initialize,
                                  control, BPPARAM, verbose, subset=TRUE){
  #+ subset fitting
  curDims <- dim(K(fds, type))
  if(isTRUE(subset)){
    probE <- max(0.001, min(1,30000/curDims[1]))
    featureExclusionMask(fds) <- sample(c(TRUE, FALSE), curDims[1],
                                        replace=TRUE, prob=c(probE, 1-probE))
  } else{
    featureExclusionMask(fds) <- rep(TRUE,curDims[1])
  }
  print(table(featureExclusionMask(fds)))
  
  fds <- fitAutoencoder(fds=fds, q=q, type=type, noiseAlpha=noiseAlpha, rhoRange=rhoRange, lambda=lambda,
                 convergence=convergence, iterations=iterations, initialize=initialize,
                 control=control, BPPARAM=BPPARAM, verbose=verbose)
  return(fds)
  
}