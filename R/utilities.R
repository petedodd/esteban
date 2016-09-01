x2list <- function(x,nmz=NULL){
    if(is.null(nmz) & is.null(names(x))){
        names(x) <- rownames(rngs)
        if(length(x)!=nrow(rngs)) stop('If no names supplied, x must be lenth nrow(rngs)!')
    }
    if(is.null(nmz)) nmz <- names(x)  #can be shorter
    if(is.null(names(x))) names(x) <- nmz;
    ans <- as.list(rowMeans(rngs))
    for(nm in nmz) ans[[nm]] <- x[nm]
    ans[rownames(rngs)]
}


toparm <- function(X,iset=1:nrow(rngs)){
  X <- ilogit(X)
  X <- X * (rngs[iset,2]-rngs[iset,1])
  X <- X + rngs[iset,1]
  X
}


fromparm <- function(x,iset=1:nrow(rngs)){
  x <- x - rngs[iset,1]
  x <- x / (rngs[iset,2]-rngs[iset,1])
  logit(x)
}

toparmM <- function(X,iset=1:nrow(rngs)){
  X <- exp(X)
  X <- X/(1+X)
  X <- t(X)                           #col=dim
  X <- X %*% diag(rngs[iset,2]-rngs[iset,1])
  X <- X + matrix(rngs[iset,1],ncol=length(iset),nrow=nrow(X),byrow=TRUE)
  X
}

logit <- function(x) log(x/(1-x))
ilogit <- function(x) exp(x)/(1+exp(x))

logprior <- function(x,var=NULL){                   #take a list
    if(!is.null(var)){
        ans <- 0                        #default value
        if(var=='bet')ans <- dlnorm(2*x,mean=1.678,sdlog=0.371,log=TRUE)
        if(var=='eps')ans <- dlnorm(x,2.303,0.329,log=TRUE)
        if(var=='pr')ans <- dbeta(x,13.730,123.568,log=TRUE)
        if(var=='v')ans <- dbeta(x,20.70,77.88,log=TRUE)
        if(var=='foi')ans <- -x*100
        if(var=='dcdr')ans <- -x*100
        if(var=='theta')ans <- dlnorm(x,-6.502,0.605,log=TRUE)
        if(var=='rho')ans <- dlnorm(x,-1.022,0.167,log=TRUE)
        if(var=='HRa')ans <- dlnorm(x,-1.050,0.115,log=TRUE)
        ## CFRs
        if(var=='HCFR00')ans <- dbeta(x,23.676,6.678,log=TRUE)
        if(var=='HCFR01')ans <- dbeta(x,7.776,78.621,log=TRUE)
        if(var=='HCFRi0')ans <- dbeta(x,9.541,5.848,log=TRUE)
        if(var=='HCFRi1')ans <- dbeta(x,1.347,35.369,log=TRUE)
        if(var=='HCFRe0')ans <- dbeta(x,11.881,12.366,log=TRUE)
        if(var=='HCFRe1')ans <- dbeta(x,3.924,71.172,log=TRUE)
        if(var=='CFRd')ans <- dbeta(x,4.369,109.921,log=TRUE)
        if(var=='CFRu')ans <- dbeta(x,25.482,33.779,log=TRUE)
    } else {
        ans <- sum(mapply(FUN=logprior,unlist(x),names(x)))
    }
    ans
}

LLfun <- function(x,LLP,iset=1:nrow(rngs)){
    if(any(is.nan(x))) stop('NaN into log-likelihood!')
    ans <- (x - 2*log(1+exp(x)))  - x^2/18
    ans <- sum(ans)
    x <- toparm(x,iset=iset)
    if( any(! x > rngs[iset,1]) | any(! x < rngs[iset,2]) ) return(-1e15)
    ltbp <- x2list(x)
    ans <- AimDynFLD(LLP,ltbp,stoch=FALSE,graph=FALSE,indat = LLP$LLdat) 
    ans <- ans + logprior(ltbp)
    ans
}

##' Calibrate model parameters
##'
##' content to be written
##' @title getMAP - optimize to find point (MAP) estimate for model
##' @param LLP 
##' @return list(par,value) with par a list of parameters and value log-posterior
##' @author Pete Dodd
##' @export
getMAP <- function(LLP){
    infvar <- rownames(rngs)[c(1:7,9:15)]
    x0 <- rowMeans(rngs[infvar,])
    x0['foi'] <- 5e-3
    x0['bet'] <- 3
    infset <- rep(NA,length(x0))
    for(i in 1:length(x0)) infset[i] <- which(rownames(rngs)==names(x0)[i])
    x1 <- fromparm(x0,infset)           #initial value in transorfmed space
    ## optimize
    sola <- optim(par=x1,fn=function(x)-LLfun(x,LLP,infset),
                  method='Nelder-Mead',
                  control=list(maxit=6e2))
    if(!sola$convergence) warning('Optimization routine did not converge')
    x0 <- toparm(sola$par,infset)
    tbpo <- x2list(x0)
    list(par=tbpo,value=sola$value)
}

##' Markov-chain Monte Carlo model parameter inference
##'
##' content to be written
##' @title runMCMC - run MCMC approach to parameter inference
##' @param tbp 
##' @param LLP 
##' @param NI 
##' @param NW 
##' @param previous 
##' @return data.frame of multiple MCMC chains 
##' @author Pete Dodd
##' @export
runMCMC <- function(tbp,LLP,NI=200,NW=20,previous=NULL){
    iss <- 1:nrow(rngs)
    xx0 <- unlist(tbp)
    xx0 <- fromparm(xx0)
    NT <- 1
    XX0 <- matrix(xx0,nrow=nrow(rngs),ncol=NW*NT) *
        rnorm(n=nrow(rngs)*NW*NT,mean=1,sd=1e-2)
    zcs <- which(abs(xx0)<1e-10)
    if(length(zcs)>0) XX0[zcs,] <- matrix(rnorm(n=length(zcs)*NW*NT,mean=1,sd=1e-2),nrow=length(zcs),ncol=NW*NT)
    rownames(XX0) <- names(xx0)
    XX0 <- XX0[iss,]
    ## run mcmc
    if(is.null(previous)){
        testmc <- PTMCMC(x0=XX0,
                         sigprop=1e-1/length(iss),
                         ntemps=NT,nwalkers=NW,niter=NI,swapn=1,tempfac = 2,
                         adaptstart = 1e2/NW,lag=50,bet=.1,
                         LL=function(x) LLfun(x,LLP,iss) )    
    }
    for(i in 1:dim(testmc$X)[2]) testmc$X[,i,] <- t(toparmM(testmc$X[,i,],iss))
    dimnames(testmc$X) <-list(variable=row.names(rngs)[iss],walker=1:NW,
                              iteration=1:dim(testmc$X)[3])
    coldat <- matrix(c(testmc$X),ncol=dim(testmc$X)[1],nrow=NI*NW,byrow=T)
    colnames(coldat) <- row.names(rngs)[iss]
    as.data.frame(coldat)
}

##' A function to generate parameter sample from MCMC data 
##'
##' content to be written 
##' @title prepMCMC - post-process MCMC data from runMCMC
##' @param X 
##' @param nsamps 
##' @param burnin 
##' @param NW 
##' @return a data.frame of sampled parameters
##' @author Pete Dodd
##' @export
prepMCMC <- function(X,nsamps=150,burnin=150,NW=20){
    if(nrow(X)%%NW) stop('nrow(X) must have NW as a divisor!')
    if(NW*burnin>=nrow(X)){
        warning('NW*burnin>=nrow(X)! Using 3/4 of nrow(X)...')
        keep <- round(.75*nrow(X)):nrow(X)
    } else keep <- (NW*burnin):nrow(X)
    X <- X[keep,]
    if(nsamps>nrow(X)) stop('nsamps > availabe rows!')
    keep <- seq(from=1,to=nrow(X),len=nsamps)
    X[keep,]
}

