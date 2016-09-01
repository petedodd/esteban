getCountryTBData <- function(P){
    ## needs to make notes, dreg, prev, and V
    dat <- data.frame(t=1:length(P$syear:P$eyear),notes=NA,dreg=NA,prev100k=NA)
    V <- rep(NA,nrow(dat))
    ## notifications
    if(P$iso3 %in% notifs$iso3){
        cat('Notification data found...\n')
        tmp <- notifs[iso3==P$iso3 & year %in% P$syear:P$eyear]
        tmp <- tmp[order(tmp$year)]
        dat$notes[min(tmp$year):max(tmp$year)-P$syear+1] <- tmp$notes
    } else cat('Notification not data found...\n')
    ## deaths
    if(P$iso3 %in% TBMs$iso3){
        cat('VR data found...\n')
        tmp <- TBMs[iso3==P$iso3 & year %in% P$syear:P$eyear]
        tmp <- tmp[order(tmp$year)]
        dat$dreg[min(tmp$year):max(tmp$year)-P$syear+1] <- tmp$deaths
    } else cat('VR data not found...\n')
    ## prevalence
    if(P$iso3 %in% prev$iso3){
        cat('Prevalence data found...\n')
        tmp <- prev[iso3==P$iso3 & year %in% P$syear:P$eyear]
        tmp <- tmp[order(tmp$year)]
        ## NB not handling smr only surveys. Uses adj for preference then bc.
        dat$prev100k[tmp$year-P$syear+1] <- tmp$prev.bc.100k
        V[tmp$year-P$syear+1] <- ((tmp$prev.bc.100k.hi-tmp$prev.bc.100k.lo)/3.92)^1
        tmp2 <- tmp$prev.adj.100k
        tmp3 <- ((tmp$prev.adj.100k.hi-tmp$prev.adj.100k.lo)/3.92)^1
        dat$prev100k[tmp$year-P$syear+1][!is.na(tmp2)] <- tmp2
        V[tmp$year-P$syear+1][!is.na(tmp3)] <- tmp3
    } else cat('Prevalence data not found...\n')
    LLdat <- dat[,-1]
    LLdat$V <- V
    LLdat[is.na(LLdat)] <- -1
    list(dat=dat,V=V,LLdat=LLdat)
}
