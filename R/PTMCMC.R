## trace function
tr <- function(x)sum(diag(x))           #trace

## covariance calculation, including shrinkage estimators
getEmpSig <- function(X,shrink='none',verbose=FALSE){
    DIM <- nrow(X)
    n <- ncol(X)
    if(any(is.na(X))){
        sna <- which.min(apply(X,2,function(x)any(!is.na(x)))) #start of NAs
    } else {sna <- ncol(X)+1}
    sna <- sna-1
    XS <- X[,1:sna]                     #X w/o the NAs
    mn <- rowMeans(XS)
    XS <- XS - matrix(mn,ncol = sna,nrow=nrow(XS),byrow=FALSE) #centre
    S <- XS %*% t(XS) / sna                                   #empiric
    if(shrink=='none') return(S)
    F <- diag(tr(S)/DIM,DIM)            #shrink target
    if(shrink=='OAS'){                  #
        ## OAS shrinkage
        rho <- .1
        alpha <- mean(S^2)
        mu <- tr(S)/DIM
        num <- alpha + mu^2
        den <- (sna + 1) * (alpha - mu^2/DIM)
        rho <- 1
        if(den>0) rho <- min(num/den,1)
        SS <- (1-rho) * S + rho * diag(mu,DIM)
        if(verbose) cat('rho = ',rho,'\n')
    }
    return(SS)
}


## parallel tempered MCMC with adaptive metropolis
PTMCMC <- function(x0,sigprop,ntemps=10,nwalkers=10,niter=1e3,swapn=1,tempfac=2,
                   adaptstart = Inf,lag=250,bet=.1, #adaptation
                   LL){
    nchains <- ntemps*nwalkers
    if(is.null(dim(x0))){
        DIM <- length(x0)
        x0 <- matrix(rep(x0,nchains) + rnorm(n=DIM*nchains,sd=mean(sigprop)),
                     nrow=DIM,ncol=nchains)
    } else { DIM <- nrow(x0)}
    X <- array(NA,dim=c(DIM,nchains,niter))
    LLchn <- array(NA,dim=c(nchains,niter))
    stemps <- tempfac^((0:(ntemps-1))/2)
    cat('temperature ladder = ',stemps,'\n')
    temps <- rep(stemps,each = nwalkers)
    X[,,1] <- x0
    Vbef <- apply(x0,2,LL)
    Vbef <- Vbef / temps
    LLchn[,1] <- Vbef
    acc <- rep(0,ntemps)
    acc2 <- 0
    if(is.null(dim(sigprop))){
        if(length(sigprop)==1){
            cat('growing siprop...\n')
            sigprop <- rep(sigprop,ntemps)
            sigprop <- sigprop * stemps^1
        }
        cat('making SIGMA...\n')
        sigmat <- array(NA,dim=c(DIM,DIM,ntemps))
        for(i in 1:ntemps)
            sigmat[,,i] <- diag(sigprop[i],DIM)
    } else {
        if(!all(dim(sigprop)==c(DIM,DIM,ntemps))){
            stop('If sigprop array, dims must be npars,npars,ntemps!')
        }
        sigmat <- sigprop
        cat('using SIGMA supplied...\n')
    }
    x <- array(NA,dim=c(DIM,nchains))
    pflag <- TRUE                       #for displaying adaptation
    for( i in 2:niter){
        if(!i%%ceiling(niter/10)) cat(round(1e2*i/niter),' % complete \n')
        ## normal proposing
        for(j in 1:ntemps)
            x[,(1+(j-1)*nwalkers):(j*nwalkers)] <- t( mvtnorm::rmvnorm(n=nwalkers,mean=rep(0,DIM),sigma=sigmat[,,j]) ) 
        x <- X[,,i-1] + x
        Vnow <- apply(x,2,LL)
        Vnow <- Vnow / temps
        test <- ( log(runif(nchains)) <= (Vnow - Vbef) ) #accept
        X[,,i] <- X[,,i-1]              #previous
        if(any(is.na(test))) print(test)
        X[,test,i] <- x[,test]          #set to new
        for(j in 1:ntemps)
            acc[j] <- acc[j] + sum(test[1:nwalkers + (j-1)*nwalkers])/nwalkers
        Vbef[test] <- Vnow[test]
        ## swapping
        if(!i%%swapn){
            swap <- sample(nchains,2)        #propose swap
            Vnow <- Vbef
            Vnow[rev(swap)] <- Vnow[swap] * temps[swap]/ temps[rev(swap)]  #real LL
            test <- ( log(runif(1)) <= sum(Vnow[swap] - Vbef[swap]) ) #accept
            if(test){
                ## print(X[1,swap,i])
                X[,rev(swap),i] <- X[,swap,i]
                Vbef[swap] <- Vnow[swap]
                acc2 <- acc2 + 1
            }
        }
        LLchn[,i] <- Vbef
        ## adaptation
        if(i > adaptstart){
            if(pflag){
                cat('Starting adaptive phase...\n')
                pflag <- FALSE
            }
            for(j in 1:ntemps){
                if(runif(1)<bet){
                    sigmat[,,j] <- diag(1e-2/DIM,DIM) #safe one
                } else {
                    xtmp <- X[,(1+(j-1)*nwalkers):(j*nwalkers), max(1,i-lag):i]
                    xtmp <- matrix(xtmp,nrow=DIM,ncol=prod(dim(xtmp)[2:3]))
                    sigmat[,,j] <- 2.38^2*getEmpSig(xtmp)/DIM
                    ## if( tr(sigmat[,,j])<1e-7 | is.nan(sum(sigmat[,,j])) |
                    ##     (sum(sigmat[,,j])==Inf) | is.na(sum(sigmat[,,j])) ) sigmat[,,j] <- diag(1e-2/DIM,DIM) #safe one
                    if( tr(sigmat[,,j])<1e-7 ) sigmat[,,j] <- diag(1e-2/DIM,DIM) #safe one
                }
            }
        }
    }    
    ## end
    acc <- round(1e2*acc/niter)
    acc2 <- round(1e2*acc2*swapn/niter)
    print(acc)
    cat('acc2 = ',acc2,'\n')
    list(X = X,LLchn = LLchn,acc=acc,acc2=acc2,sigmat=sigmat,
         ntemps=ntemps,nwalkers=nwalkers,niter=niter,swapn=swapn,tempfac=tempfac,
         lag=lag,bet=bet,adaptstart=adaptstart,LL=LL)
}

