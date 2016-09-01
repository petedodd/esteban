##' A function for running the AIM dynamics with infection (fortran version)
##'
##' Content to be written
##' @title AIM dynamics with LTBI and deterministic approximation
##' @param P 
##' @param graph 
##' @return ans list .... TODO
##' @author Pete Dodd
##' @import scales
##' @import ggplot2
##' @export AimDynFLD
##' @useDynLib esteban
AimDynFLD <- function(P,tbp,stoch=FALSE,coarse=FALSE,graph=FALSE,indat=NULL){
    list2env(P,envir=environment())
    X <- DT <- array(0,dim=dz,dimnames=list(nmsex,nmage,nmhiv,nmart))
    X[1,,1,1] <- UPD$bp[year==syear & sex==1 & age<81,value] #men
    X[2,,1,1] <- UPD$bp[year==syear & sex==2 & age<81,value] #women
    zr <- yrz*0
    nsmll <- length(arttarget)
    tbp <- as.numeric(tbp)
    rec <- array(0,dim=c(ntime,20))
    LLmode <- FALSE
    if(!is.null(indat)){
        if(is.data.frame(indat)) indat <- as.matrix(indat)
        indat[is.na(indat)] <- -1       #for fortran version; no negative data allowed
        if(!( all(dim(indat)==c(nsmll,3)) | all(dim(indat)==c(nsmll,4)) )) stop(paste0('Expecting indat to have 3 or 4 cols,',nsmll,' rows; got ', ncol(indat),' columns and ',nrow(indat),' rows!'))
        if(ncol(indat)==3) indat <- cbind(indat,rep(0,nrow(indat))) #add another row
        indat[indat[,3]<=0,4] <- 0                                  #safety - variance must be zero
        LLmode <- TRUE
    } else
        indat <- rep(-1,4*nsmll)  
    ## fortran call
    ans <- .Fortran('aimdynLD',           #name
                    ntime=as.integer(ntime),nsmll=as.integer(nsmll), #dims
                    nmsex=as.integer(dz[1]),nmage=as.integer(dz[2]),
                    nmhiv=as.integer(dz[3]),nmart=as.integer(dz[4]), 
                    popt=as.double(zr),poph=as.double(zr),pop1549=as.double(zr), #res
                    poph1549=as.double(zr),popa=as.double(zr),yrz=as.double(yrz),
                    U=as.double(X),L=as.double(0*X),tstep=as.double(tstep),     #initial state
                    FINC=as.double(0*X), #final incidence
                    hivtarget=as.double(hivtarget),hivtarget2=as.double(hivtarget2),
                    arttarget=as.double(arttarget),arttarget2=as.double(arttarget2),
                    nbmv=as.double(nbmv), nbfv=as.double(nbfv),srv=as.double(srv),
                    migrates=as.double(migrates),mortrates=as.double(mortrates),
                    ar=as.double(AR),cdinit=as.double(cdinit[,,,1]),pp=as.double(pp),
                    mh=as.double(mh),ma=as.double(ma),
                    tbp=as.double(tbp),record=as.double(rec),DT=as.double(DT),
                    stoch=as.integer(stoch),coarse=as.integer(coarse),
                    indat=as.integer(indat),LL=as.double(0))
    
    ## back to array
    X <- array(ans$U+ans$L,dim=dz,dimnames=list(nmsex,nmage,nmhiv,nmart))
    U <- array(ans$U,dim=dz,dimnames=list(nmsex,nmage,nmhiv,nmart))
    L <- array(ans$L,dim=dz,dimnames=list(nmsex,nmage,nmhiv,nmart))
    FINC <- array(ans$FINC,dim=dz,dimnames=list(nmsex,nmage,nmhiv,nmart))
    DT <- array(ans$DT,dim=dz,dimnames=list(nmsex,nmage,nmhiv,nmart))
    rec <- matrix(ans$record,nrow=ntime,ncol=20)
    if(coarse) rec <- rec[1:nsmll,]
    colnames(rec) <- c('t','cdr','DD','DU','VR','notes','unnotes','deaths','dreg','TBI',
                    'ARI','pop','vDD','vDU','vnotes','vprev100k','vdeaths','vdreg','vcdr','prev100k')
    rec <- as.data.frame(rec)
    rec$t <- 1:nrow(rec)                #shift time
    rec$t <- rec$t - 1
    if(!stoch | (stoch & !coarse)) rec$t <- rec$t/P$npy
    rec$t <- rec$t + P$syear
    rec$prev <- rec$DD + rec$DU
    if(!stoch | (stoch & !coarse)) rec[,c('TBI','notes','deaths','dreg')] <- rec[,c('TBI','notes','deaths','dreg')]*npy #
    ## graphing
    if(graph){
        AIMplotterL(P,U,L,yrz,ans$popt,ans$poph,ans$poph1549,ans$pop1549,ans$popa)
        AIMplotterLD(P,U,L,yrz,ans$popt,ans$poph,ans$poph1549,ans$pop1549,ans$popa,rec,DT,FINC) 
    }
    if(LLmode)
        return(ans$LL)
    else
        return(list(U=U,L=L,rec=rec,DT=DT,FINC=FINC))
}
