## numerically stable backsolving via cholesky decomp
bsi <- function(CH,v) backsolve(CH,backsolve(CH, v, transpose = TRUE))
##' An Ensemble Kalman Filter
##'
##' Content to be written
##' @title Ensemble Kalman-Filter
##' @param X 
##' @param H 
##' @param d 
##' @param eps 
##' @return X2 a matrix of conditioned data
##' @author Pete Dodd
enkf <- function(X,H,d,eps=1){
    if(!is.matrix(X) | !is.matrix(H))stop('X and H must be matrices!')
    nens <- ncol(X)
    nobs <- length(d)
    ndim <- nrow(X)
    if(!all(dim(H)==c(nobs,ndim)))stop('dim mismatch between X,H,d!')
    if(!(length(eps)==1 | length(eps)==nobs)){cat('nobs =',nobs,'; length(eps) =',length(eps),'\n');stop('length(eps) should be 1 or len(d)!');}
    if(length(eps)==1) eps <- rep(eps,nobs)
    D <- matrix(d,ncol=nens,nrow=nobs)
    E <- matrix(rnorm(nobs*nens)*eps,ncol=nens,nrow=nobs) 
    D <- D + E
    EX <- rowMeans(X)
    A <- X - matrix(EX,ncol=nens,nrow=ndim)
    C <- A %*% t(A) / (nens-1)
    HA <- H %*% A
    R <- diag(eps^2)
    U <- chol(HA %*% t(HA)/(nens-1) + R)
    X2 <- X + C %*% t(H) %*% bsi(U, D - H %*% X)
    return(X2)
}
