AIMplotterL <- function(P,U,L,yrz,popt,poph,poph1549,pop1549,popa){
    ## todo: as with HIV plot -- correct the order
    ## do the other plots
    AIMplotter(P=P,X=U+L,yrz=yrz,popt=popt,poph=poph,poph1549=poph1549,pop1549=pop1549,popa=popa)
    ## ## -----------------
    tmpU <- reshape2::melt(U)
    tmpL <- reshape2::melt(L)
    names(tmpU) <- names(tmpL) <- c('sex','age','hiv','art','value')
    tmpUL <- tmpU
    tmpUL$L <- tmpL$value
    names(tmpUL)[5] <- 'U'
    tmpUL <- as.data.table(tmpUL)
    tmpUL <- tmpUL[,list(U=sum(U),L=sum(L)),by=list(age)]
    tmpU$LTBI <- FALSE; tmpL$LTBI <- TRUE;
    tmp <- rbind(tmpU,tmpL)
    tmp <- as.data.table(tmp)
    tmplat <- tmp[,list(value=sum(value)),by=list(sex,age,LTBI)]
    tmplat <- tmplat[rev(order(LTBI))]
    tmplat$LTBI <- factor(tmplat$LTBI,levels=c('TRUE','FALSE'),ordered=TRUE)      #new
    undmg <- d$N[iso3==P$iso3 & year== max(yrz) &
                     !age %in% c('100-','95-99','90-94','85-89','80-84') ,
                 list(sex2=sex,age2=age,value=1e3*value/5)]
    undmg$sex <- 'M'
    undmg$sex[undmg$sex2=='female'] <- 'F'
    undmg$HIV <- FALSE
    undmg$LTBI <- FALSE
    undmg$age <- as.numeric(unlist(strsplit(as.character(undmg$age2),split='-'))[seq(from=1,to=nrow(undmg),by=2)])
    ## demographic + LTBI
    plt <- ggplot(data=tmplat,aes(x=age,fill=LTBI))+ geom_bar(data=tmplat[sex=='M'],aes(y=value),stat='identity')+ geom_bar(data=tmplat[sex=='F'],aes(y=-value),stat='identity') + coord_flip()+geom_abline(intercept=0,slope=0) + ylab('Number') + xlab('Age')+ theme(axis.text.y = element_text(size = rel(0.5)))+ geom_step(data=undmg[sex=='F'],aes(y=-value)) + geom_step(data=undmg[sex=='M'],aes(y=value)) + scale_y_continuous(label=comma) + ggtitle(P$iso3)
    ggsave(plt,file=paste0(P$iso3,'_chk_ageLTBIdemo.pdf'),width=7,height=7)        #age LTBI
    ## LTBI proportion
    plt <- ggplot(data=tmpUL,aes(x=age,y=L/(U+L)))+ geom_bar(stat='identity')+ylab('LTBI') + xlab('Age')+ theme(axis.text.x = element_text(size = rel(0.5),angle=90)) + scale_y_continuous(label=percent,limits=c(0,1)) + ggtitle(P$iso3)
    ggsave(plt,file=paste0(P$iso3,'_chk_ageLTBI.pdf'),width=7,height=7)        #age LTBI
    ## return(tmplat)
}
