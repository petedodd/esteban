AIMplotter <- function(P,X,yrz,popt,poph,poph1549,pop1549,popa){
    tmp <- d$N[iso3==P$iso3 & year %in% yrz &
                   !age %in% c('100-','95-99','90-94','85-89','80-84') ,
               list(value=1e3*sum(value)),by=year]
    modpop <- data.frame(yrz,popt,poph,poph1549,pop1549,popa)
    ggplot(data=modpop,aes(x=yrz,y=popt)) + geom_line(color=2) + scale_y_continuous(label=comma)+xlab('year')  + ylab('Total population') +geom_point(data=tmp,aes(x=year,y=value))+ ggtitle(P$iso3)
    ggsave(paste0(P$iso3,'_chk_totpop.pdf'),width=7,height=7)
    ## -----------------        
    hatmp <- HA$HP[iso3==P$iso3 & year %in% yrz & type=='mid',list(hp1549=mean(hp1549)),by=year]
    ggplot(data=modpop,aes(x=yrz,y=poph1549/pop1549)) + geom_line(color=2)+ geom_point(data=P$data,aes(x=year,y=hivp*1e-2)) + scale_y_continuous(label=percent) + xlab('year') + ylab('HIV prevalence in 15-49 year olds') + ggtitle(P$iso3) #+ geom_point(data=hatmp,aes(x=year,y=hp1549*1e-2),col=2)
    ggsave(paste0(P$iso3,'_chk_hivp.pdf'),width=7,height=7)
    ## -----------------
    aatmp <- HA$AP[iso3==P$iso3 & year %in% yrz & type=='mid' & age=='15+' & sex=='both',list(acov,year)]
    ggplot(data=modpop,aes(x=yrz,y=popa/(poph+1e-6))) + geom_line(color=2) + scale_y_continuous(label=percent) + xlim(c(2000,max(modpop$yrz)))+ xlab('year') + ylab('ART coverage in HIV') +  geom_point(data=aatmp,aes(x=year,y=acov*1e-2))
    ggsave(paste0(P$iso3,'_chk_artp.pdf'),width=7,height=7)
    ## -----------------
    tmp <- reshape2::melt(X)
    names(tmp) <- c('sex','age','hiv','art','value')
    tmp <- as.data.table(tmp)
    tmpart <- tmp[,list(value=sum(value)),by=list(age,hiv,art)]
    tmp <- tmp[,list(value=sum(value)),by=list(sex,age,hiv)]
    tmp$value <- pmax(tmp$value,0)
    tmp[,HIV:=!(hiv=='hiv-ve')]
    tmpa <- tmp[,list(value=sum(value)),by=list(sex,age)]
    tmpah <- tmp[,list(value=sum(value)),by=list(sex,age,HIV)]
    tmpah <- tmpah[(order(HIV))]        
    tmpah$HIV <- factor(tmpah$HIV,levels=c('TRUE','FALSE'),ordered=TRUE)      #new
    undmg <- d$N[iso3==P$iso3 & year== max(yrz) &
                     !age %in% c('100-','95-99','90-94','85-89','80-84') ,
                 list(sex2=sex,age2=age,value=1e3*value/5)]
    undmg$sex <- 'M'
    undmg$sex[undmg$sex2=='female'] <- 'F'
    undmg$HIV <- FALSE
    undmg$age <- as.numeric(unlist(strsplit(as.character(undmg$age2),split='-'))[seq(from=1,to=nrow(undmg),by=2)])
    ggplot(data=tmpa,aes(x=age,fill=sex))+ geom_bar(data=tmpa[sex=='M'],aes(y=value),stat='identity')+ geom_bar(data=tmpa[sex=='F'],aes(y=-value),stat='identity') + coord_flip()+geom_abline(intercept=0,slope=0) + ylab('Number') + xlab('Age') + theme(axis.text.y = element_text(size = rel(0.5))) + geom_step(data=undmg[sex=='F'],aes(y=-value)) + geom_step(data=undmg[sex=='M'],aes(y=value)) + scale_y_continuous(label=comma)+ ggtitle(P$iso3)
    ggsave(paste0(P$iso3,'_chk_demog.pdf'),width=7,height=7)         #demography
    ## -----------------plyr::desc
    ggplot(data=tmpah,aes(x=age,fill=HIV,order=(HIV)))+ geom_bar(data=tmpah[sex=='M'],aes(y=value),stat='identity')+ geom_bar(data=tmpah[sex=='F'],aes(y=-value),stat='identity') + coord_flip()+geom_abline(intercept=0,slope=0) + ylab('Number') + xlab('Age')+ theme(axis.text.y = element_text(size = rel(0.5)))+ geom_step(data=undmg[sex=='F'],aes(y=-value)) + geom_step(data=undmg[sex=='M'],aes(y=value)) + scale_y_continuous(label=comma) + ggtitle(P$iso3)
    ggsave(paste0(P$iso3,'_chk_agehiv.pdf'),width=7,height=7)        #age HIV: 
    ## -----------------
    tmpart$hiv <- factor(tmpart$hiv,levels=rev(c(levels(tmpart$hiv)[-1],levels(tmpart$hiv)[1])),ordered=TRUE)
    tmpart$hiv <- factor(tmpart$hiv,levels=rev(levels(tmpart$hiv)),ordered=TRUE)
    ggplot(data=tmpart[hiv!='hiv-ve',],aes(x=age,fill=hiv))+ geom_bar(aes(y=value),stat='identity')+ geom_abline(intercept=0,slope=0) + ylab('Number') + xlab('Age') + facet_wrap(~art,scales='free_y') + scale_y_continuous(label=comma)+ theme(axis.text.x = element_text(size = rel(0.5),angle = 90)) + labs(fill='CD4 category')+ ggtitle(P$iso3)
    ggsave(paste0(P$iso3,'_chk_agecd4.pdf'),width=14,height=7)        #age HIV
}
