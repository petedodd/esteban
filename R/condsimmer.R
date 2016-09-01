##' A function for conditional simulation
##'
##' Content to be written
##' @title condsimmer (un)conditional simulation and graphing
##' @param dat data to condition to (as in the dat slot of country parameters); NULL if no conditioning
##' @param PL country parameters as built by getCountryAimParms 
##' @param tbpl a list of model parameters
##' @param n the number of simulations; overwritten by nrow(tbpl) if this is a data.frame
##' @param quants the unobserved quantities of interest
##' @param obs the observed quanttities
##' @param V the s.d. of prevalence observations
##' @return X a long-form data frame of observed and unobserved quantities
##' @author Pete Dodd
##' @import scales
##' @export
condsimmer <- function(dat = NULL, PL, tbpl, n = 150,
                       quants = c("TBI", "deaths", "prev","prev100k"),
                       obs = c("notes", "dreg", "prev100k"),V=0,
                       show=FALSE,truedat=NULL,fname=NULL,stoch=TRUE,coarse=TRUE,west=FALSE){
    ## safety
    if(!is.null(dat)){
        if(length(V)==nrow(dat)){
            if(any(is.na(dat$prev100k) & !is.na(V))) stop('V must be NA where prevalence is not!')
            if(any(!is.na(dat$prev100k) & is.na(V))) stop('V not NA where prevalence is!')
        }
        if(!(length(V)==1 | length(V)==nrow(dat)))
            stop('V must be of length 1  or pad out to nrow(dat) with NAs!')
        if(length(V)==1){ V <- rep(V,nrow(dat)); V[is.na(dat$prev100k)] <- NA;}
        if( any(!is.na(dat$prev100k))){
            cat('Found prevalence data...\n')
        }else{  cat('No prevalence data found...\n') }
        if( any(!is.na(dat$notes))){
            cat('Found notification data...\n')
        }else{ cat('No notification data found...\n') }
        if( any(!is.na(dat$dreg))){
            cat('Found VR data...\n')
        }else{ cat('No VR data found...\n') }
    }
    if(!stoch){
        coarse <- FALSE
        if(n>0) n <- 1
        warning('Since deterministic, any n>0 is set to 1 and coarse is set FALSE.')
    }

    ## handle matrix or df input
    manypar <- FALSE
    if(is.matrix(tbpl))
        if(is.null(colnames(tbpl)))
            stop('If tbpl is a matrix, it must have colnames set!')
    if(is.data.frame(tbpl)){
        n <- nrow(tbpl)                 #use the input data to overwrite number of runs?
        tblist <- split(tbpl,seq(nrow(tbpl)))
        manypar <- TRUE
    }

    nmz <- unique(c(quants,obs))
    
    ## simulations
    if(n>0){
        ans <- list()
        for (i in 1:n) {
            if(manypar) tbpl <- tblist[[i]]
            rec <- AimDynFLD(PL, tbpl, stoch = stoch, coarse = coarse, graph = 0)$rec
            if (!is.null(dat)){
                if (nrow(rec) != nrow(dat)) {
                    print(dim(rec))
                    print(dim(dat))
                    stop("data/simulation mismatch!")
                }
            }
            rec$replicate <- i
            ans[[i]] <- rec
        }
        ans <- do.call("rbind", ans)
    } else {
        ans <- data.frame(t=PL$syear:PL$eyear,replicate=1)
        for(nm in nmz) ans[[nm]] <- NA
        for(nm in nmz) ans[[nm]] <- as.numeric(ans[[nm]])
        bad <- c()
        for(nm in nmz){
            ditch <- TRUE
            if(!is.null(truedat[[nm]]))
                if(any(!is.na(truedat[[nm]]))) ditch <- FALSE
            if(!is.null(dat[[nm]]))
                if(any(!is.na(dat[[nm]]))) ditch <- FALSE
            if(west)
                if(nrow(wHL[iso3==PL$iso3 & variable==nm])) ditch <- FALSE
            if(ditch)
                bad <- c(bad,nm)
        }
        nmz <- nmz[!nmz%in%bad]
        ans <- ans[,c('t','replicate',nmz)]
    }
    
    ## conditioning
    if (is.null(dat)){                  #just reshape if no data
        X <- ans
        X <- reshape2::melt(X,id=c('t','replicate'))
    } else {                            #data
        X <- ans[, c('t','replicate',nmz)]
        ## reshape X
        Xm <- reshape2::melt(X,id=c('t','replicate'))
        X <- matrix(NA,ncol=max(Xm$replicate),nrow=nrow(Xm)/max(Xm$replicate))
        for(i in 1:ncol(X))
            X[,i] <- Xm[Xm$replicate==i,'value']
        nmzi <- as.character(Xm[Xm$replicate==1,'variable'])
        X <- trans(X)
        nmziy <- paste0(nmzi,'_',as.character(Xm[Xm$replicate==1,'t']))
        
        ## make H and d
        ## make d
        d <- c()
        for(qn in obs)
            d <- c(d,dat[,qn])
        dnmz <- rep(obs,each=nrow(dat))
        dnmzy <- paste0(dnmz,'_',1:nrow(dat) + PL$syear-1)
        VL <- rep(V,length(obs))
        for(i in 1:length(obs))
            if( obs[i] != 'prev100k' ) VL[1:nrow(dat) + (i-1)*nrow(dat)] <- 0

        
        dnmz <- dnmz[!is.na(d)]         #names of data in vector
        dnmzy <- dnmzy[!is.na(d)]         #names of data in vector
        VL <- VL[!is.na(d)]
        d <- d[!is.na(d)]               #data
        d <- trans(d)                   #transform
        nobs <- length(d)
        ## make H
        H <- matrix(0,nrow=nobs,ncol=nrow(X))
        for(i in 1:length(dnmzy)){
            who <- which(nmziy==dnmzy[i])
            H[i,who] <- 1
        }
        
        if(n>10){
            Xc2 <- enkf(X=X,H=H,d=d,eps = 1e2 + VL)
        } else {
            if(n>0) warning('n<=10 so conditioning not performed!')
            Xc2 <- X
        }
        Xc2 <- itrans(Xc2)

        ## reshape output
        X <- data.frame(t=1:nrow(dat) - 1 + PL$syear,
                        replicate=rep(1:ncol(Xc2),each=nrow(Xc2)),
                        variable=nmzi,
                        value=c(Xc2))
    }

    ## plotting
    if(show){
        
        ## if true data supplied
        if(!is.null(truedat)){
            umz <- intersect(nmz,names(truedat))
            truez <- truedat[,c('t',umz)]
            truez$t <- 1:nrow(truez) -1 +PL$syear
            truez <- reshape2::melt(truez,id='t')
            truez$replicate <- 1
        }

        ## prevalence surveys
        if(!is.null(dat)){
            if( sum(!is.na(dat$prev100k))>0 ){
                pss <- which(!is.na(dat$prev100k))
                exdat <- data.frame(t=(1:nrow(dat))[pss] - 1 + PL$syear,
                                    variable='prev100k',value=dat$prev100k[pss],replicate=1,
                                    ymax=dat$prev100k[pss] + V[pss]*1.96,
                                    ymin=dat$prev100k[pss] - V[pss]*1.96)
            }
            indat <- dat[,c('notes','dreg','t')]
            indat$t <- 1:nrow(indat) - 1 + PL$syear
            indat <- reshape2::melt(indat,id='t')
            indat$replicate <- 1
        }

        uv <- c('deaths','dreg','notes','TBI','prev','prev100k') #default
        if(max(X$replicate)>=1){
            ## new version with ribbons
            XS <- data.table::as.data.table(X)[,list(value=median(value,na.rm=TRUE),
                                         ub=quantile(value,probs = .975,na.rm=TRUE),
                                         lb=quantile(value,probs = .025,na.rm=TRUE)),
                                   by=list(t,variable)]
            XS$lb <- pmax(XS$lb,0)
            XS$ub <- pmax(XS$ub,0)
            XS <- XS[XS$variable %in% uv,]
        } else {  ## 1  rep
            if(!is.null(dat)){
                XS <- indat[!indat$variable %in% bad,]
                indat <- indat[!indat$variable %in% bad,]
                if( sum(!is.na(dat$prev100k))>0 )
                    XS <- rbind(XS,exdat[,1:4])
            }
            if(west){
                tmp <- wHL[iso3==PL$iso3,list(t=year,variable,value=mid)]
                tmp$replicate <- 1
                XS <- rbind(XS,tmp)
            }
            XS <- XS[XS$variable %in% nmz,]
            XS <- XS[!is.na(XS$value),]
        }
        
        gp <- ggplot(data=XS,aes(x=t,y=value)) + geom_blank() + xlab('Year') + ggtitle(PL$iso3) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
        gp <- gp + facet_wrap(~variable,scales = 'free_y')
        gp <- gp + scale_y_continuous(label=comma) + expand_limits(y=0)
        if(max(X$replicate)>10)
            gp <- gp + geom_ribbon(aes(ymin=lb,ymax=ub),fill='grey')
        if(n>0)
            gp <- gp + geom_line()
        if(!is.null(truedat))           #if true data supplied
            gp <- gp + geom_point(data=truez,col=6,alpha=.5)
        if(!is.null(dat)){           #if prevalence data supplied?
            gp <- gp + geom_errorbar(data=exdat,aes(x=t,y=value,ymin=ymin,ymax=ymax),col=4,width=1) + geom_point(data=exdat,shape=7,size=2.5,col=4)
            gp <- gp  + geom_point(data=indat,col=4)
        }
        if(west)
            gp <- gp + geom_line(data=wHL[iso3==PL$iso3],aes(x=year,y=mid),col=2)+
                geom_line(data=wHL[iso3==PL$iso3],aes(x=year,y=hi),col=2,lty=2)+
                geom_line(data=wHL[iso3==PL$iso3],aes(x=year,y=lo),col=2,lty=2)
        

        ## save/show
        if(!is.null(fname))
            ggsave(gp,file=fname)
        else
            print(gp)

    }

    X
}

trans <- function(x) x
itrans <- function(x) x

