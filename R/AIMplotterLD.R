AIMplotterLD <- function(P,U,L,yrz,popt,poph,poph1549,pop1549,popa,record,DT,FINC){
    ## ##other plots
    AIMplotterL(P,U,L,yrz,popt,poph,poph1549,pop1549,popa)
    
    ## ARI plot
    plt <- ggplot(data=record,aes(P$syear+t,ARI)) + geom_line() + scale_y_continuous(labels = percent) + xlab('Time since start') + ylab('ARI') + expand_limits(y=0)
    ggsave(plt,file=paste0(P$iso3,'_tbchk_ARI.pdf'),width=7,height=7)        #

    ## TB prevalence
    plt <- ggplot(data=record,aes(t,1e5*(DU+DD)/pop)) + geom_line() + xlab('Time since start') + ylab('TB prevalence per 100,000') + expand_limits(y=0)
    ggsave(plt,file=paste0(P$iso3,'_tbchk_TBP.pdf'),width=7,height=7)        #

    ## TB incidence
    plt <- ggplot(data=record,aes(t,1e5*TBI/pop)) + geom_line()+ xlab('Time since start') + ylab('TB incidence per 100,000 per year') + expand_limits(y=0)
    ggsave(plt,file=paste0(P$iso3,'_tbchk_TBI.pdf'),width=7,height=7)        #

    ## TB notifications
    plt <- ggplot(data=record,aes(t,1e5*notes/pop)) + geom_line()+ xlab('Time since start') + ylab('TB notifications per 100,000 per year')  + expand_limits(y=0)
    ggsave(plt,file=paste0(P$iso3,'_tbchk_TBN.pdf'),width=7,height=7)        #

    ## TB deaths
    plt <- ggplot(data=record,aes(t,1e5*deaths/pop)) + geom_line()+ xlab('Time since start') + ylab('TB mortality per 100,000 per year') + expand_limits(y=0)
    ggsave(plt,file=paste0(P$iso3,'_tbchk_TBM.pdf'),width=7,height=7)        #

    ## in addition, by age!
    ## ## -----------------
    tmpU <- reshape2::melt(U+L)
    tmpL <- reshape2::melt(DT)
    names(tmpU) <- names(tmpL) <- c('sex','age','hiv','art','value')
    tmpUL <- tmpU
    tmpUL$TB <- tmpL$value
    names(tmpUL)[5] <- 'X'
    tmpUL <- as.data.table(tmpUL)
    tmpUL <- tmpUL[age!='[80,Inf)',list(X=sum(X),TB=sum(TB)),by=list(age)] 
    plt <- ggplot(data=tmpUL,aes(x=age,y=1e5*TB/X))+ geom_bar(stat='identity')+ylab('TB prevalence per 100,000') + xlab('Age')+ theme(axis.text.x = element_text(size = rel(0.5),angle=90))  + ggtitle(P$iso3)
    ggsave(plt,file=paste0(P$iso3,'_tbchk_TBPage.pdf'),width=7,height=7)        #age LTBI
    ## return(tmplat)

    ## TBP by age and HIV
    tmp <- reshape2::melt(DT)
    names(tmp) <- c('sex','age','hiv','art','value')
    tmp <- as.data.table(tmp)
    tmp[,hivart:='HIV-ve']
    tmp[hiv!='hiv-ve' & art!='none',hivart:='HIV+ve/ART+ve']
    tmp[hiv!='hiv-ve' & art=='none',hivart:='HIV+ve/ART-ve']
    tmp <- tmp[age!='[80,Inf)',list(value=sum(value)),by=list(age,hivart)] 
    plt <- ggplot(data=tmp,aes(x=age,y=value,fill=hivart))+ geom_bar(stat='identity')+ylab('TB prevalence') + xlab('Age')+ theme(axis.text.x = element_text(size = rel(0.5),angle=90)) + ggtitle(P$iso3)
    ggsave(plt,file=paste0(P$iso3,'_tbchk_TBPageHIV.pdf'),width=10,height=7)        #age LTBI

    ## TBI
    ## TBP by age and HIV
    tmp <- reshape2::melt(FINC)
    names(tmp) <- c('sex','age','hiv','art','value')
    tmp <- as.data.table(tmp)
    tmp[,hivart:='HIV-ve']
    tmp[hiv!='hiv-ve' & art!='none',hivart:='HIV+ve/ART+ve']
    tmp[hiv!='hiv-ve' & art=='none',hivart:='HIV+ve/ART-ve']
    tmp <- tmp[age!='[80,Inf)',list(value=sum(value)),by=list(age,hivart)] 
    plt <- ggplot(data=tmp,aes(x=age,y=value,fill=hivart))+ geom_bar(stat='identity')+ylab('TB incidence per year') + xlab('Age')+ theme(axis.text.x = element_text(size = rel(0.5),angle=90)) + ggtitle(P$iso3)
    ggsave(plt,file=paste0(P$iso3,'_tbchk_TBIageHIV.pdf'),width=10,height=7)        #age LTBI
    
}
