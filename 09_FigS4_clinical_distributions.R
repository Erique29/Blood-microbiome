#!/usr/bin/env Rscript
# Compare clinical characteristics across age and ER groups.
# Fig. S4; Data S4. Input: script 01.

local({
  .files <- Filter(Negate(is.null), lapply(sys.frames(), function(f) f$ofile))
  .file <- if (length(.files)) .files[[length(.files)]] else NULL
  .cli <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  .dir <- if (!is.null(.file)) dirname(normalizePath(.file)) else if (length(.cli)) {
    dirname(normalizePath(sub("^--file=", "", .cli[1])))
  } else getwd()
  assign("SKCCC_SCRIPTS_DIR", .dir, envir = .GlobalEnv)
  sys.source(file.path(.dir, "00_config_and_analysis_framework.R"), envir = .GlobalEnv)
})

OUT_DIR <- file.path(PROJECT_DIR,'Revision2/fig1-3/preview_ER_subtype/output')
suppressPackageStartupMessages({library(phyloseq);library(vegan);library(ggplot2);library(cowplot);library(openxlsx)})
for(d in c('','Tables','Audit','QA'))dir.create(file.path(OUT_DIR,d),recursive=TRUE,showWarnings=FALSE)
FONT<-11;INK<-'#000000';BLUE<-'#4D76A6';PINK<-'#C77896';set.seed(666)
th<-function()theme_classic(base_size=FONT,base_family='sans')+theme(text=element_text(size=FONT,colour=INK),axis.title=element_text(size=FONT,colour=INK),axis.text=element_text(size=FONT,colour=INK),legend.title=element_text(size=FONT,colour=INK),legend.text=element_text(size=FONT,colour=INK),plot.title=element_blank(),plot.subtitle=element_blank(),plot.caption=element_blank(),axis.line=element_line(colour=INK,linewidth=.4),plot.margin=margin(8,10,8,8))
save_pdf<-function(p,name,w,h){dev<-if(capabilities('aqua'))function(filename,width,height,bg='white',...)grDevices::quartz(type='pdf',file=filename,width=width,height=height,bg=bg,family='Arial') else grDevices::pdf;ggsave(file.path(OUT_DIR,paste0(name,'.pdf')),p,width=w,height=h,device=dev,bg='white')}
wt<-function(z,n,d='Tables')write.csv(z,file.path(OUT_DIR,d,paste0(n,'.csv')),row.names=FALSE,na='')

stat_annotation<-function(first,p,hjust=-.03){list(
 annotate('text',x=-Inf,y=Inf,label=first,hjust=hjust,vjust=1.15,size=FONT/(72.27/25.4),colour=INK),
 annotate('text',x=-Inf,y=Inf,label=sprintf('italic(P) == %.3f',p),parse=TRUE,hjust=hjust,vjust=2.8,size=FONT/(72.27/25.4),colour=INK))}
input<-file.path(PROJECT_DIR,'Output/physeq.rds');ps<-readRDS(input);m<-data.frame(sample_data(ps),check.names=FALSE);a<-as(otu_table(ps),'matrix');if(taxa_are_rows(ps))a<-t(a)
stopifnot(nrow(a)>1,identical(rownames(m),rownames(a)),!anyDuplicated(m$Patient_id),all(a>=0),all(is.finite(a)),max(abs(rowSums(a)-1))<1e-6,all(m$ER%in%c('N','Y')))
m$ER_group<-factor(m$ER,levels=c('N','Y'),labels=c('ER negative','ER positive'))
cols<-c('ER negative'=BLUE,'ER positive'=PINK)

fields<-c('Race: Black'='race','HER2 positive'='HER2','Surgery: Yes'='Pre.enrollment_surgery','Chemotherapy: Yes'='Pre.enrollment_Chemotherapy','Radiotherapy: Yes'='Pre.enrollment_Radiation','Hormonal therapy: Yes'='Pre.enrollment_Hormonal.Therapy','Targeted therapy: Yes'='Pre.enrollment_Targeted.Therapy')
proportions<-list();age_summary<-list()
for(er in levels(m$ER_group)){
 md<-m[m$ER_group==er,];ag<-md$age_MBC_dx
 age_summary[[er]]<-data.frame(ER=er,n_total=nrow(md),n_age=sum(is.finite(ag)),median=median(ag,na.rm=TRUE),Q1=unname(quantile(ag,.25,na.rm=TRUE)),Q3=unname(quantile(ag,.75,na.rm=TRUE)),min=min(ag,na.rm=TRUE),max=max(ag,na.rm=TRUE))
 for(j in seq_along(fields)){
  x<-as.character(md[[fields[j]]]);missing<-is.na(x)|trimws(x)=='';value<-if(j==1)'Black' else 'Y';valid<-!missing
  stopifnot(all(x[valid]%in%if(j==1)c('Black','White') else c('N','Y')))
  yes<-sum(x[valid]==value);denom<-sum(valid)
  proportions[[length(proportions)+1]]<-data.frame(variable=names(fields)[j],field=unname(fields[j]),ER=er,n_total=nrow(md),n_observed=denom,n_missing=sum(missing),n_category=yes,percent=100*yes/denom,row=length(fields)-j+1)
 }
}
proportions<-do.call(rbind,proportions);age_summary<-do.call(rbind,age_summary)
proportions$ER<-factor(proportions$ER,levels=levels(m$ER_group));proportions$y<-proportions$row+ifelse(proportions$ER=='ER negative',.14,-.14)
wt(proportions,'ER_categorical_distributions');wt(age_summary,'ER_age_summary');wt(m[,c('Sample_id','Patient_id','age_MBC_dx','ER','race','HER2',unname(fields[-c(1,2)]))],'clinical_plot_input','Audit')

age_test<-wilcox.test(m$age_MBC_dx~m$ER_group,exact=FALSE,correct=TRUE)
er_tests<-list(data.frame(variable='Age',field='age_MBC_dx',test='Wilcoxon rank-sum, asymptotic tie/continuity corrected',n_negative=sum(is.finite(m$age_MBC_dx)&m$ER=='N'),n_positive=sum(is.finite(m$age_MBC_dx)&m$ER=='Y'),missing_negative=sum(!is.finite(m$age_MBC_dx)&m$ER=='N'),missing_positive=sum(!is.finite(m$age_MBC_dx)&m$ER=='Y'),statistic=unname(age_test$statistic),P=age_test$p.value,OR_positive_vs_negative=NA_real_,CI_low=NA_real_,CI_high=NA_real_,row=NA_real_))
for(j in seq_along(fields)){
 x<-as.character(m[[fields[j]]]);ok<-!is.na(x)&trimws(x)!='';target<-if(j==1)'Black'else'Y'
 yy<-factor(ifelse(x[ok]==target,'Yes','No'),levels=c('No','Yes'));tab<-table(m$ER_group[ok],yy)
 ft<-fisher.test(tab,alternative='two.sided')
 er_tests[[length(er_tests)+1]]<-data.frame(variable=names(fields)[j],field=unname(fields[j]),test='Fisher exact two-sided',n_negative=sum(ok&m$ER=='N'),n_positive=sum(ok&m$ER=='Y'),missing_negative=sum(!ok&m$ER=='N'),missing_positive=sum(!ok&m$ER=='Y'),statistic=NA_real_,P=ft$p.value,OR_positive_vs_negative=unname(ft$estimate),CI_low=ft$conf.int[1],CI_high=ft$conf.int[2],row=length(fields)-j+1)
 wt(as.data.frame.matrix(tab),paste0('ER_test_contingency_',j),'Audit')
}
er_tests<-do.call(rbind,er_tests);er_tests$q_BH8<-p.adjust(er_tests$P,'BH')
er_tests$P_label<-ifelse(er_tests$P<.001,'P < 0.001',sprintf('P = %.3f',er_tests$P))
er_tests$stars_rawP<-ifelse(er_tests$P<.001,'***',ifelse(er_tests$P<.01,'**',ifelse(er_tests$P<.05,'*','')))
wt(er_tests,'ER_group_tests')
pAge<-ggplot(m,aes(ER_group,age_MBC_dx,fill=ER_group))+geom_boxplot(width=.4,outlier.shape=NA,alpha=.25,colour=INK,linewidth=.45)+geom_point(aes(colour=ER_group),position=position_jitter(width=.09,seed=666),size=2.2,alpha=.85)+scale_fill_manual(values=cols,guide='none')+scale_colour_manual(values=cols,guide='none')+labs(x=NULL,y='Age at MBC diagnosis (years)')+th()+theme(aspect.ratio=1)+scale_y_continuous(expand=expansion(mult=c(.03,.16)))+annotate('text',x=-Inf,y=Inf,label=sprintf('italic(P) == %.3f',er_tests$P[1]),parse=TRUE,hjust=-.08,vjust=1.15,size=FONT/(72.27/25.4),colour=INK)
pProp<-ggplot(proportions,aes(percent,y,colour=ER))+geom_segment(aes(x=0,xend=percent,yend=y),linewidth=.6,alpha=.35)+geom_point(size=2.8)+geom_text(aes(label=paste0(n_category,'/',n_observed)),nudge_x=2.3,hjust=0,size=FONT/(72.27/25.4),colour=INK)+scale_colour_manual(values=cols,name=NULL)+geom_text(data=er_tests[!is.na(er_tests$row)&er_tests$P<.05,],aes(x=108,y=row,label=stars_rawP),inherit.aes=FALSE,hjust=0,size=FONT/(72.27/25.4),colour=INK)+scale_x_continuous(limits=c(0,118),breaks=c(0,25,50,75,100),labels=function(x)paste0(x,'%'),expand=expansion(mult=c(.015,0)))+scale_y_continuous(breaks=1:7,labels=rev(names(fields)),limits=c(.5,7.5),expand=c(0,0))+labs(x='Participants with the indicated characteristic (%)',y=NULL)+th()+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),legend.position='bottom')
pClinical<-plot_grid(pAge,pProp,nrow=1,labels=c('A','B'),label_size=FONT,label_colour=INK,rel_widths=c(1,1.9),align='h',axis='tb')
save_pdf(pClinical,'ER_clinical_covariates',12.8,5.3)

m$age_group<-factor(ifelse(m$age_MBC_dx<50,'<50 years','>=50 years'),levels=c('<50 years','>=50 years'))
stopifnot(!anyNA(m$age_group),all(table(m$age_group)>0))
age_fields<-c('ER positive'='ER',fields)
age_proportions<-list();age_tests<-list()
for(j in seq_along(age_fields)){
 field<-unname(age_fields[j]);label<-names(age_fields)[j];x<-as.character(m[[field]])
 ok<-!is.na(x)&trimws(x)!='';target<-if(field=='race')'Black'else'Y'
 stopifnot(all(x[ok]%in%if(field=='race')c('Black','White')else c('N','Y')))
 binary<-factor(ifelse(x[ok]==target,'Yes','No'),levels=c('No','Yes'));tab<-table(m$age_group[ok],binary)
 ft<-fisher.test(tab,alternative='two.sided')
 age_tests[[j]]<-data.frame(variable=label,field=field,test='Fisher exact two-sided',n_under50=sum(ok&m$age_group=='<50 years'),n_50plus=sum(ok&m$age_group=='>=50 years'),missing_under50=sum(!ok&m$age_group=='<50 years'),missing_50plus=sum(!ok&m$age_group=='>=50 years'),P=ft$p.value,OR_50plus_vs_under50=unname(ft$estimate),CI_low=ft$conf.int[1],CI_high=ft$conf.int[2],row=length(age_fields)-j+1)
 wt(as.data.frame.matrix(tab),paste0('Age_test_contingency_',j),'Audit')
 for(g in levels(m$age_group)){
  ii<-m$age_group==g;valid<-ok&ii;den<-sum(valid);ncat<-sum(x[valid]==target)
  age_proportions[[length(age_proportions)+1]]<-data.frame(variable=label,field=field,age_group=g,n_total=sum(ii),n_observed=den,n_missing=sum(ii&!ok),n_category=ncat,percent=100*ncat/den,row=length(age_fields)-j+1)
 }
}
age_tests<-do.call(rbind,age_tests);age_tests$q_BH8<-p.adjust(age_tests$P,'BH');age_tests$stars_rawP<-ifelse(age_tests$P<.001,'***',ifelse(age_tests$P<.01,'**',ifelse(age_tests$P<.05,'*','')))
age_proportions<-do.call(rbind,age_proportions);age_proportions$age_group<-factor(age_proportions$age_group,levels=levels(m$age_group));age_proportions$y<-age_proportions$row+ifelse(age_proportions$age_group=='<50 years',.14,-.14)

all_clinical_tests<-rbind(data.frame(grouping='ER',variable=er_tests$variable,P=er_tests$P,q_BH8=er_tests$q_BH8),data.frame(grouping='Age <50 vs >=50',variable=age_tests$variable,P=age_tests$P,q_BH8=age_tests$q_BH8))
all_clinical_tests$q_BH16<-p.adjust(all_clinical_tests$P,'BH')
wt(age_proportions,'Age_group_distributions');wt(age_tests,'Age_group_tests');wt(all_clinical_tests,'All_clinical_group_tests')
wt(m[,c('Sample_id','Patient_id','age_MBC_dx','age_group','ER','race','HER2',unname(fields[-c(1,2)]))],'clinical_plot_input','Audit')
agecols<-setNames(c(BLUE,PINK),levels(m$age_group))
pAgeGroups<-ggplot(age_proportions,aes(percent,y,colour=age_group))+geom_segment(aes(x=0,xend=percent,yend=y),linewidth=.6,alpha=.35)+geom_point(size=2.8)+geom_text(aes(label=paste0(n_category,'/',n_observed)),nudge_x=2.3,hjust=0,size=FONT/(72.27/25.4),colour=INK)+scale_colour_manual(values=agecols,name=NULL,labels=c('<50 years','≥50 years'))+geom_text(data=age_tests[age_tests$P<.05,],aes(x=108,y=row,label=stars_rawP),inherit.aes=FALSE,hjust=0,size=FONT/(72.27/25.4),colour=INK)+scale_x_continuous(limits=c(0,118),breaks=c(0,25,50,75,100),labels=function(x)paste0(x,'%'),expand=expansion(mult=c(.015,0)))+scale_y_continuous(breaks=1:8,labels=rev(names(age_fields)),limits=c(.5,8.5),expand=c(0,0))+labs(x='Participants with the indicated characteristic (%)',y=NULL)+th()+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),legend.position='bottom')

pAge<-pAge+theme(aspect.ratio=1)
pAgeGroups<-pAgeGroups+theme(aspect.ratio=1)
pProp<-pProp+theme(aspect.ratio=1)
pExpanded<-plot_grid(pAge,pAgeGroups,pProp,nrow=1,labels=c('A','B','C'),label_size=FONT,label_colour=INK,align='hv',axis='tblr')
save_pdf(pExpanded,'FigS4_clinical_distributions',20.5,6.6)

save_session(file.path(OUT_DIR,"Audit"))
