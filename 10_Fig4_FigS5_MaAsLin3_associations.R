#!/usr/bin/env Rscript
# Plot univariate MaAsLin3 results and fit the specified age/ER-adjusted models.
# Fig. 4 and Fig. S5; Data S5. Inputs: scripts 01 and 05.

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

suppressPackageStartupMessages({library(phyloseq); library(maaslin3); library(dplyr); library(tidyr); library(ggplot2); library(cowplot)})
set.seed(666)
module <- file.path(project,'Revision2/Fig4')
out <- file.path(module,'output')
for(d in c('Main','Supplementary','Tables','Audit','Models')) dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
phy <- readRDS(file.path(project,'Output/physeq.rds'))
abund <- as(otu_table(phy),'matrix'); if(taxa_are_rows(phy)) abund <- t(abund)
meta <- data.frame(sample_data(phy),check.names=FALSE)[rownames(abund),,drop=FALSE]
stopifnot(nrow(abund)>1L, identical(rownames(abund),rownames(meta)), !anyDuplicated(rownames(abund)),all(is.finite(abund)),all(abund>=0))
meta$ctc2 <- ifelse(meta$CTC_count_cat %in% c('0','0-5'),'0-5',ifelse(meta$CTC_count_cat=='>5','>5',NA))
meta$age_group <- factor(meta$age_group,levels=c('<50','>=50'))
meta$ER <- factor(meta$ER,levels=c('N','Y'))
meta$race <- factor(meta$race,levels=c('Black','White'))
meta$ctc2 <- factor(meta$ctc2,levels=c('0-5','>5'))
spec <- list(
 age=list(var='age_group',folder='age_group_only_0.5',levels=c('<50','>=50'),labels=c('<50','≥50'),cols=c('#D6A0B6','#8C2857')),
 race=list(var='race',folder='race_only_0.5',levels=c('Black','White'),labels=c('Black','White'),cols=c('#95BCE0','#225E9C')),
 ER=list(var='ER',folder='ER_only_0.5',levels=c('N','Y'),labels=c('Negative','Positive'),cols=c('#96C9BE','#147D70')),
 CTC=list(var='ctc2',folder='ctc2_only_0.5',levels=c('0-5','>5'),labels=c('0–5','>5'),cols=c('#B6AAD8','#65419B')))
read_results <- function(path) read.delim(path,check.names=FALSE,stringsAsFactors=FALSE)
legacy <- lapply(spec,function(s)read_results(file.path(project,'Output/Maaslin3/Paper used analysis',s$folder,'all_results.tsv')))

audit <- lapply(names(spec),function(k){s<-spec[[k]]; z<-read.delim(file.path(project,'Output/Maaslin3/Paper used analysis',s$folder,'features/filtered_data.tsv'),row.names=1,check.names=FALSE); b<-abund[rownames(z),colnames(z),drop=FALSE]; stopifnot(max(abs(as.matrix(z)-b),na.rm=TRUE)<1e-12, all(b[is.na(z)]==0)); data.frame(panel=k,N=nrow(z),taxa=ncol(z),max_abs_difference=max(abs(as.matrix(z)-b),na.rm=TRUE),legacy_NA_are_current_zero=TRUE)})
write.csv(bind_rows(audit),file.path(out,'Audit/legacy_input_verification.csv'),row.names=FALSE)
source_files<-c(file.path(project,'Output/physeq.rds'),vapply(spec,function(s)file.path(project,'Output/Maaslin3/Paper used analysis',s$folder,'all_results.tsv'),character(1)))
write.csv(data.frame(path=source_files,md5=unname(tools::md5sum(source_files))),file.path(out,'Audit/input_manifest.csv'),row.names=FALSE)

write.csv(cbind(sample_id=rownames(meta),meta[,c('age_group','ER','race','CTC_count_cat','ctc2')]),file.path(out,'Audit/sample_inclusion.csv'),row.names=FALSE)
fitspec <- list(age_ER=list(formula='~ age_group + ER',vars=c('age_group','ER'),n=nrow(meta)),race_age_ER=list(formula='~ race + age_group + ER',vars=c('race','age_group','ER'),n=nrow(meta)),CTC_age_ER=list(formula='~ ctc2 + age_group + ER',vars=c('ctc2','age_group','ER'),n=sum(!is.na(meta$ctc2))))
fitone <- function(k){f<-fitspec[[k]]; m<-droplevels(meta[complete.cases(meta[,f$vars,drop=FALSE]),f$vars,drop=FALSE]); a<-abund[rownames(m),,drop=FALSE]; X<-model.matrix(as.formula(f$formula),m); stopifnot(nrow(m)==f$n,qr(X)$rank==ncol(X)); dest<-file.path(out,'Models',k)
 cfg<-list(formula=f$formula,metadata=m,input_md5=unname(tools::md5sum(file.path(project,'Output/physeq.rds'))),package=as.character(packageVersion('maaslin3')),seed=666,parameters=list(min_prevalence=.5,normalization='TSS',transform='LOG',augment=TRUE,standardize=TRUE,median_comparison_abundance=TRUE,median_comparison_prevalence=FALSE,max_significance=.25,correction='BH',cores=4))
 cache<-file.path(dest,'configuration.rds'); result<-file.path(dest,'all_results.tsv')
 if(file.exists(result)){stopifnot(file.exists(cache),identical(readRDS(cache),cfg)); message('Using matching cache: ',k)}else{dir.create(dest,recursive=TRUE,showWarnings=FALSE);saveRDS(cfg,cache);set.seed(666); do.call(maaslin3,c(list(input_data=a,input_metadata=m,output=dest,formula=f$formula,plot_summary_plot=FALSE,plot_associations=FALSE,max_pngs=0),cfg$parameters))}
 read_results(result)
}
if(Sys.getenv('FIG4_PLOT_ONLY','false')=='true') {
 adjusted_all<-lapply(names(fitspec),function(k)read_results(file.path(out,'Models',k,'all_results.tsv'))); names(adjusted_all)<-names(fitspec)
}else adjusted_all <- setNames(lapply(names(fitspec),fitone),names(fitspec))
adjusted <- list(age=subset(adjusted_all$age_ER,metadata=='age_group'),race=subset(adjusted_all$race_age_ER,metadata=='race'),ER=subset(adjusted_all$age_ER,metadata=='ER'),CTC=subset(adjusted_all$CTC_age_ER,metadata=='ctc2'))
alltable <- bind_rows(lapply(names(spec),function(k)bind_rows(base::transform(legacy[[k]],analysis='Main univariate',panel=k),base::transform(adjusted[[k]],analysis='Supplement adjusted',panel=k))))
write.csv(alltable,file.path(out,'Tables/Fig4_complete_target_results.csv'),row.names=FALSE)
write.csv(bind_rows(lapply(names(adjusted_all),function(k)base::transform(adjusted_all[[k]],model_formula=fitspec[[k]]$formula))),file.path(out,'Tables/Adjusted_all_predictors_results.csv'),row.names=FALSE)
counts <- alltable %>% group_by(analysis,panel,model) %>% summarise(taxa=n(),joint_q_lt_025=sum(qval_joint<.25,na.rm=TRUE),individual_q_lt_025=sum(qval_individual<.25,na.rm=TRUE),errors=sum(!is.na(error)&nzchar(error)),.groups='drop')
write.csv(counts,file.path(out,'Tables/Evidence_summary.csv'),row.names=FALSE)

fdr_audit <- bind_rows(lapply(names(adjusted_all),function(k){z<-adjusted_all[[k]]; u<-z[!duplicated(z[,c('feature','metadata','value','name')]),]; data.frame(model=k,joint_tests=nrow(u),max_q_difference=max(abs(p.adjust(u$pval_joint,method='BH')-u$qval_joint),na.rm=TRUE))}))
stopifnot(all(fdr_audit$max_q_difference<1e-12))
write.csv(fdr_audit,file.path(out,'Audit/FDR_families.csv'),row.names=FALSE)

prevalence <- bind_rows(lapply(names(spec),function(k){s<-spec[[k]]; features<-unique(legacy[[k]]$feature);bind_rows(lapply(seq_along(s$levels),function(i){ids<-which(!is.na(meta[[s$var]])&as.character(meta[[s$var]])==s$levels[i]);data.frame(panel=k,feature=features,group=s$levels[i],label=s$labels[i],detected=colSums(abund[ids,features,drop=FALSE]>0),N=length(ids),prevalence=100*colMeans(abund[ids,features,drop=FALSE]>0))}))}))
write.csv(prevalence,file.path(out,'Tables/Prevalence_named_groups.csv'),row.names=FALSE)
font_pt<-11; text_mm<-font_pt/ggplot2::.pt; star_mm<-16/ggplot2::.pt
phylum_cols<-c(Proteobacteria='#D7B348',Firmicutes='#8C5A2B',Actinobacteriota='#D45740')
base_theme<-theme_bw(base_size=font_pt,base_family='Arial')+theme(text=element_text(size=font_pt,colour='black'),axis.text=element_text(size=font_pt,colour='black'),axis.title=element_text(size=font_pt,colour='black'),strip.text=element_text(size=font_pt,colour='black'),legend.text=element_text(size=font_pt,colour='black'),legend.title=element_text(size=font_pt,colour='black'),panel.grid=element_blank(),plot.margin=margin(4,3,4,3))

maxcoef<-max(abs(alltable$coef[grepl('abundance',alltable$model)]),na.rm=TRUE)
lim<-2.5
stopifnot(maxcoef<=lim)
panel <- function(k,df){s<-spec[[k]]; d<-df[grepl('abundance',df$model),];d<-d[order(d$coef),];d$species<-gsub('_',' ',vapply(strsplit(d$feature,';',fixed=TRUE),tail,character(1),1));d$phylum<-vapply(strsplit(d$feature,';',fixed=TRUE),function(x)x[2],character(1)); d$species<-factor(d$species,levels=unique(d$species)); d$direction<-ifelse(d$coef>=0,'pos','neg');d$star_x<-d$coef+ifelse(d$coef<0,-1,1)*.22
 pc<-ggplot(d,aes(y=species))+geom_vline(xintercept=0,linetype='dashed',colour='grey45',linewidth=.4)+geom_segment(aes(x=0,xend=coef,yend=species,colour=direction),linewidth=.8)+geom_point(aes(x=coef,colour=direction),size=2.8)+geom_text(data=d[!is.na(d$qval_joint)&d$qval_joint<.25,],aes(x=star_x,label='*'),size=star_mm,colour='black',vjust=.55)+scale_colour_manual(values=c(neg=s$cols[1],pos=s$cols[2]),guide='none')+scale_x_continuous(limits=c(-lim,lim),breaks=c(-2.5,0,2.5),labels=c('-2.5','0','2.5'),expand=expansion(mult=c(.01,.01)))+scale_y_discrete(limits=levels(d$species),expand=expansion(add=.6))+labs(x='Coefficient',y=NULL)+base_theme+theme(axis.text.y=element_text(face='italic'),plot.margin=margin(4,3,4,14))
 pp<-ggplot(d,aes(x=1,y=species,fill=phylum))+geom_tile(width=.8,colour='white')+scale_fill_manual(values=phylum_cols,drop=FALSE,guide='none')+scale_y_discrete(limits=levels(d$species),expand=expansion(add=.6))+theme_void(base_size=font_pt,base_family='Arial')+theme(plot.margin=margin(4,0,4,0))
 pr<-prevalence[prevalence$panel==k,];pr$species<-factor(as.character(d$species)[match(pr$feature,d$feature)],levels=levels(d$species));pr$label<-factor(pr$label,levels=s$labels)
 bottom<-setNames(paste0('Prevalence in\n',s$labels,' (%)'),s$labels)
 pv<-ggplot(pr,aes(x=prevalence,y=species,fill=label))+geom_col(width=.7)+facet_grid(.~label,switch='x',labeller=labeller(label=bottom))+scale_fill_manual(values=setNames(s$cols,s$labels),guide='none')+scale_x_continuous(limits=c(0,100),breaks=c(0,50,100),expand=expansion(mult=c(0,.08)))+scale_y_discrete(limits=levels(d$species),expand=expansion(add=.6))+labs(x=NULL,y=NULL)+base_theme+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),strip.placement='outside',strip.background=element_blank(),panel.spacing=grid::unit(14,'pt'),plot.margin=margin(4,10,4,3))
 plot_grid(pc,pp,pv,nrow=1,align='h',axis='tb',rel_widths=c(1.5,.055,1.0))
}
legend_plot<-function(vals,title){d<-data.frame(x=names(vals),y=1);get_legend(ggplot(d,aes(x=x,y=y,fill=x))+geom_tile()+scale_fill_manual(values=vals,breaks=names(vals),name=title)+guides(fill=guide_legend(ncol=1,title.position='top',title.hjust=0))+theme_void(base_size=font_pt,base_family='Arial')+theme(legend.position='top',legend.title=element_text(size=font_pt,colour='black'),legend.text=element_text(size=font_pt,colour='black'),legend.key.size=grid::unit(9,'pt'),legend.margin=margin(0,0,0,0),legend.spacing.x=grid::unit(3,'pt')))}
legtop<-legend_plot(phylum_cols,'Phylum')
legs<-lapply(names(spec),function(k){s<-spec[[k]];legend_plot(setNames(s$cols,s$labels),switch(k,age='Age',race='Race',ER='ER',CTC='CTC'))})
legendrow<-plot_grid(plotlist=c(list(legtop),legs),nrow=1,rel_widths=c(1.5,1.3,1,1,1))
makefigure<-function(results,path){ps<-lapply(names(spec),function(k)panel(k,results[[k]]));body<-plot_grid(plotlist=ps,nrow=2,labels=c('A','B','C','D'),label_size=font_pt,label_fontface='bold',label_fontfamily='Arial',hjust=0,vjust=1); fig<-plot_grid(legendrow,body,ncol=1,rel_heights=c(.10,.90)); dev<-function(filename,width,height,bg='white',...)grDevices::quartz(type='pdf',file=filename,width=width,height=height,bg=bg,family='Arial');ggsave(path,fig,width=20,height=13,device=dev,bg='white',limitsize=FALSE)}
makefigure(legacy,file.path(out,'Main/Fig4.pdf'))
makefigure(adjusted,file.path(out,'Supplementary/FigS_Fig4_adjusted.pdf'))
comparison<-merge(subset(alltable,analysis=='Main univariate'),subset(alltable,analysis=='Supplement adjusted'),by=c('panel','feature','metadata','value','model'),suffixes=c('_univariate','_adjusted'),all=TRUE)
write.csv(comparison,file.path(out,'Tables/Univariate_vs_adjusted.csv'),row.names=FALSE)
openxlsx::write.xlsx(list(Evidence_summary=as.data.frame(counts),Complete_target_results=alltable,Adjusted_all_predictors=bind_rows(adjusted_all,.id='fit'),Univariate_vs_adjusted=comparison,Prevalence=prevalence,FDR_families=fdr_audit),file.path(out,'Tables/Fig4_results.xlsx'),overwrite=TRUE)

writeLines(capture.output(sessionInfo()),file.path(out,'Audit/sessionInfo.txt'))
writeLines(c(paste('Shared coefficient limits:',-lim,lim),'All coefficients fit within -2.5 to 2.5; no boundary arrows.','Main stars: exact legacy qval_joint < 0.25.','Supplement stars: package qval_joint < 0.25 from each full adjusted model.','Joint support must not be described as abundance-specific support.','Prevalence labels joined using group names, not column positions.'),file.path(out,'Audit/plot_settings.txt'))
print(counts);message('Finished: ',out)
