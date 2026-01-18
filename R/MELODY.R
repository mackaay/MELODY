#' Main functions
#'
#' The Main function of MELODY
#' @param sig_matrix  sig_matrix file path to gene expression from isolated cells, or a matrix of expression profile of cells.
#'
#' @param mixture_file mixture_file file path to heterogenous mixed expression file, or a matrix of heterogenous mixed expression
#'
#' @param perm Number of permutations
#' @param maxSize maximum size for the computation, to be passed to the future.global.maxSize. Default to 500 MB
#' @import utils
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom stats sd
#' @export
#' @examples
#' \dontrun{
#'   ## example 1
#'   sig_matrix <- system.file("extdata", "blood_signature.txt", package = "MELODY")
#'   mixture_file <- system.file("extdata", "example.txt", package = "MELODY")
#'   results <- MELODY(sig_matrix, mixture_file)
#'   ## example 2
#'   data(blood_signature.txt)
#'   data(example.txt)
#'   results <- MELODY(sig_matrix = blood_signature, mixture_file = example)
#' }

library(e1071)
library(furrr)

svm_core <- function(X, y, maxSize=500){
  #try different values of nu
  svm_itera <- 3
  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
    model
  }
  enableParallel(maxSize=maxSize)
  if (Sys.info()['sysname'] == 'Windows') {
    out <- future_map(1:svm_itera, res)
  } else {
    if (svm_itera <= availableCores() - 2) {
      enableParallel(nThreads = svm_itera, maxSize=maxSize)
    } else {
      enableParallel(maxSize=maxSize)
    }
    out <- future_map(1:svm_itera, res)
  }
  nusvm <- rep(0,svm_itera)
  corrv <- rep(0,svm_itera)
  #do MELODY
  t <- 1
  while(t <= svm_itera) {
    weights = t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights<0)]<-0
    w <- weights/sum(weights)
    u <- sweep(X,MARGIN=2,w,'*')
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }
  #pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]
  #get and normalize coefficients
  q <- t(model$coefs) %*% model$SV
  q[which(q<0)]<-0
  w <- (q/sum(q))
  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]
  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
}


MELODY <- function(sig_matrix, mixture_file, maxSize=500, perm = 0){
  #read in data
  if (is.character(sig_matrix)) {
    signature.mat <- read.delim(sig_matrix, header=T, sep="\t", row.names=1, check.names = F)
    signature.mat <- data.matrix(signature.mat)
  } else {
    signature.mat <- sig_matrix
  }

  if (is.character(mixture_file)) {
    mixture <- read.delim(mixture_file, header=T, sep="\t", row.names=1, check.names = F)
    mixture <- data.matrix(mixture)
  } else {
    mixture <- mixture_file
  }

  #order
  signature.mat <- signature.mat[order(rownames(signature.mat)),]
  mixture <- mixture[order(rownames(mixture)),]

  rownames(mixture) <- gsub("-", "_", rownames(mixture))

  #check log transformed,  if max < 50 in mixture file
  if(max(mixture) < 50) {mixture <- 2^mixture}

  #intersect genes
  keep <- intersect(row.names(signature.mat), row.names(mixture))
  signature.mat <- signature.mat[keep,]
  mixture <- mixture[keep,]

  #Standardize sig matrix
  signature.mat <- (signature.mat - mean(signature.mat)) / sd(as.vector(signature.mat))

  #empirical null distribution of correlation coefficients
  if (perm > 0) {
    nulldist <- sort(doPerm(perm, signature.mat, mixture)$dist)
  }

  header <- c('Mixture',colnames(signature.mat),"P-value","Correlation","RMSE")

  output <- matrix()
  itera <- 1
  mixtures <- dim(mixture)[2]
  pval <- 0.9999

  #iterate through mixtures
  while (itera <= mixtures) {
    query <- mixture[,itera]
    #standardize mixture
    query <- (query - mean(query)) / sd(query)

    #run SVR core algorithm
    result <- svm_core(signature.mat, query, maxSize)

    #get results
    mix_w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse
    #calculate p-value
    if (perm > 0) {
      pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))
    }

    #print output
    out <- c(colnames(mixture)[itera],mix_w,pval,mix_r,mix_rmse)
    if(itera == 1) {
      output <- out
    } else {
      output <- rbind(output, out)
    }

    itera <- itera + 1
  }

  #return matrix object containing all results
  obj <- rbind(header,output)
  obj <- obj[,-1]
  obj <- obj[-1,]
  obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
  rownames(obj) <- colnames(mixture)
  colnames(obj) <- c(colnames(signature.mat),"P-value","Correlation","RMSE")
  obj
}


