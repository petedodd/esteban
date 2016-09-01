##' A function to prepare country specific parameteres for use in demographic/AIM model
##'
##' Content to be written
##' @title Parameter preparation for AIM
##' @param iso 
##' @return a list of parameters. See Details.
##' @author Pete Dodd
##' @import data.table
##' @export
getCountryAimParms <- function(iso){
    ans <- list()
    ## which region are we in: details from episty
    if(! iso %in% unique(episty$iso3)) stop('iso code not find in hash!')
    ans$iso3 <- iso
    ans$Region <- episty[iso3==iso,Region]
    ans$cn <- episty[iso3==iso,country]
    ans$type <- episty[iso3==iso,Epidemic.Type]
    ans$syear <- episty[iso3==iso,Start.year]
    if(is.na(ans$syear)) ans$syear <- 1970
    ans$ison <- episty[iso3==iso,iso_numeric]
    ans$g_whoregion <- episty[iso3==iso,g_whoregion]
    ## progression parameters, hiv and ART mortality
    ans$pp <- ans$mh <- ans$ma <- ans$cdinit <- array(0,dim=dz,dimnames=list(nmsex,nmage,nmhiv,nmart))
    for( sx in nmsex ){
        for(cd4 in nmhiv[-1]){               #CD4
            for(ag in nmagg){               #aggregated ages
                if(cd4!='<50')
                    ans$pp[sx,aglst[[ag]],cd4,1] <- prog[region==ans$Region & sex==sx  & age==ag & CD4==cd4,rate] #hiv=1
                ans$mh[sx,aglst[[ag]],cd4,1] <- muh[region==ans$Region & sex==sx  & age==ag & CD4==cd4,rate] #hiv=1
                if(TRUE) ans$mh[sx,aglst[[ag]],cd4,1] <-  ma(ans$mh[sx,aglst[[ag]],cd4,1],n=6)
                for(aa in nmart[-1]){
                    ans$ma[sx,aglst[[ag]],cd4,aa] <- mua[region==ans$Region & sex==sx  & age==ag & CD4==cd4 & tART==aa,rate] #3 times on ART
                }
            }
        }
    }
    ## CD4 infection distribution
    for(cd4 in rownames(aim_inf))
        for(ag in colnames(aim_inf))
            ans$cdinit[,aglst[[ag]],cd4,1] <- aim_inf[cd4,ag] #NB done for hiv=1
    ## age ratio, sex ratio, todo: these may need reshaping based on use in dynamics
    whom <- which(names(sr)==ans$type)
    ans$SR <- unlist(sr[,whom,with=FALSE])
    ans$AR <- AR[type==ans$type,list(M=M,F=F)][c(rep(1:12,each=5),rep(12,6)),]
    ans$AR <- as.matrix(ans$AR)
    if(any(is.na(ans$AR))) ans$AR <- matrix(1,nrow=66,ncol=2)  #safety
    ans$AR <- rbind(matrix(0,nrow=15,ncol=2),ans$AR)
    ans$elig <- data.table(year=ans$syear:2015,cd4=4) #top category eligible in nmhiv
    ## HIV incidence, prevalence etc. from EAA
    ans$data <- EAA[Code==ans$ison,list(year=year,Country=Country,hivi=hivi,artm=artm,artf=artf,ltfu=ltfu,hivp=hivp)]
    ## UNPD data
    ans$UPD <- UDP[[iso]]
    ## ---- new data pre-computed for use in dynamics -----
    ans$tstep <- .1
    ans$syear <- 1970
    ans$eyear <- 2015
    ans$ntime <- ceiling((ans$eyear-ans$syear)/ans$tstep) + 1
    ans$yrz <- seq(from=ans$syear,to=ans$eyear,by=ans$tstep)
    ans$npy <- ceiling(1/ans$tstep)                            #times per year
    ans$MB <- d$B[iso3==ans$iso3 & sex=='male',list(years,value=1e3*value/5)] #for male births
    ans$FB <- d$B[iso3==ans$iso3 & sex=='female',list(years,value=1e3*value/5)] #female births
    ans$yrbrks <- as.numeric(unlist(strsplit(as.character(ans$MB[,years]),split='-'))[seq(from=1,to=2*nrow(ans$MB),by=2)]) #year breaks
    if(nrow(HA$AP[iso3==ans$iso3])){
        ans$inHAAP <- TRUE
        ans$xyr <- HA$AP[iso3==ans$iso3,max(year)]
        ans$myr <- HA$AP[iso3==ans$iso3,min(year)]
        ans$ady <- HA$AP[iso3==ans$iso3 & sex=='both' & type=='mid' & age=='15+' & year==ans$myr+1,1e-2*acov]
        ans$ady <- ans$ady-HA$AP[iso3==ans$iso3 & sex=='both' & type=='mid' & age=='15+' & year==ans$myr,1e-2*acov] #initial gradient
        ## myr-asy = min/ady
        ans$asy <- ans$myr - HA$AP[iso3==ans$iso3 & sex=='both' & type=='mid' & age=='15+' & year==ans$myr,1e-2*acov]/ans$ady
        ans$asy <- ceiling(ans$asy)                   #art start year
        ## print(asy)
        ans$ady <- HA$AP[iso3==ans$iso3 & sex=='both' & type=='mid' & age=='15+' & year==ans$myr,1e-2*acov]/(ans$myr-ans$asy + 1e-6) #recompute gradient with whole year round for asy
    } else {
        ans$inHAAP <- FALSE
    }
    tmp <- getAIMtargets(ans)           #precomputed targets and other data
    ans <- c(ans,tmp)
    tmp <- getCountryTBData(ans)        #include TB data in this
    ans <- c(ans,tmp)
    ## ======= return ==========
    ans
}


getAIMtargets <- function(L){
    list2env(L,envir=environment())

    ## quantities: HIV/ART; M/F births; mortality(!); HIV SR; migration
    yrzs <- unique(floor(yrz))
    nsmll <- length(yrzs)
    nbm <- nbf <- srv <- rep(0,nsmll)
    hivtarget <- hivtarget2 <- arttarget <- arttarget2 <- rep(0,nsmll)
    migrates <- mortrates <- array(0,dim=c(nsmll,dz[1:2])) #yrzs x M/F x age
    j <- 0
    for(i in 1:ntime){
        if(!(i-1)%%(1/tstep)){           #NB assumes 1/tstep is integer
            j <- j+1
            yr <- floor(yrz[i])
            ## HIV & ART targets
            yr2 <- yr + 1
            if(yr2>eyear)
                yr2 <- yr

            if(inHAAP){
                hivtarget[j] <- data[year==yr,hivp*1e-2]
                hivtarget2[j] <- data[year==yr2,hivp*1e-2]
                if(yr >= myr ){ #period with data
                    if(yr>xyr)  #over end
                        arttarget[j] <- HA$AP[iso3==L$iso3 & sex=='both' & type=='mid' &
                                              age=='15+' & year==xyr,1e-2*acov]
                    else
                        arttarget[j] <- HA$AP[iso3==L$iso3 & sex=='both' & type=='mid' &
                                              age=='15+' & year==yr,1e-2*acov]
                    if( yr+1 > xyr ){
                        tmp <- HA$AP[iso3==L$iso3 & sex=='both' & type=='mid' &
                                     age=='15+' & year==xyr,1e-2*acov]
                        if(length(tmp)!=1)stop(paste0('art prob ',yr,' ',yr2,' ',length(tmp)))
                        arttarget2[j] <- tmp
                        adf <- (HA$AP[iso3==L$iso3 & sex=='both' & type=='mid' & age=='15+' & year==xyr,1e-2*acov] - HA$AP[iso3==L$iso3 & sex=='both' & type=='mid' & age=='15+' & year==xyr-1,1e-2*acov]) 
                        arttarget2[j] <- arttarget2[j] + adf*(yr+1-xyr)
                        arttarget[j] <- arttarget[j] + adf*(yr-xyr)
                    } else{
                        arttarget2[j] <- HA$AP[iso3==L$iso3 & sex=='both' & type=='mid' &
                                               age=='15+' & year==yr+1,1e-2*acov]
                    }
                } else if(yr+1 > asy){      #early extrapolation period
                    arttarget2[j] <- (yr+1-asy)*ady
                    arttarget[j] <- max( (yr-asy)*ady, 0)
                }
            } else {
                
            }
            
            ## births
            who <- sum(yr>yrbrks)       #which base element
            yover <- yr - yrbrks[who]   #how far to next 5 yrs
            p <- yover/5
            nbm[j] <- (1-p)*MB[who,value] + p*MB[who+1,value] #boys
            nbf[j] <- (1-p)*FB[who,value] + p*FB[who+1,value] #girls
            ## migration rate
            migrates[j,1,] <- UPD$migr[year==yr & sex==1,value]
            migrates[j,2,] <- UPD$migr[year==yr & sex==2,value]
            ## mortality
            mortrates[j,1,] <- UPD$lfts[year==yr & sex==1 & age<81,Sx]^tstep #tstep mu
            mortrates[j,2,] <- UPD$lfts[year==yr & sex==2 & age<81,Sx]^tstep #tstep mu

            ## hiv disaggregation
            hsyear <- syear; if(is.na(hsyear)) hsyear <- 1975
            ypost <- floor(yrz[i] - hsyear) #years after HIV start
            sryr <- min(31,max(1,ypost))                #in between
            sr <- 1                                     #default in case not there
            if('SR' %in% names(L))
                sr <- L$SR[sryr]                            #sex ratio HIV
            sr <- sr/(1+sr)                             #fraction male HIV
            srv[j] <- sr
        }
    }

    return(list(hivtarget=hivtarget,hivtarget2=hivtarget2,
                arttarget=arttarget,arttarget2=arttarget2,
                nbmv=nbm,nbfv=nbf,srv=srv,
                migrates=migrates,mortrates=mortrates))

}

ma <- function(x,n=5){y <- filter(x,rep(1/n,n), sides=2); y[is.na(y)] <- x[is.na(y)]; y}
