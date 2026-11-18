#!/usr/bin/env Rscript
# Analyze clinical subtypes and treatment exposures: descriptive statistics, diversity and MaAsLin3.
# Table 1, Data S3/S5; subtype statistics for Fig. 3. Input: script 01.

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

OUT_DIR <- file.path(PROJECT_DIR,'Revision2/Major4/output')
suppressPackageStartupMessages({library(phyloseq);library(vegan);library(maaslin3);library(openxlsx)})
for(d in c('Tables','Audit','QA'))dir.create(file.path(OUT_DIR,d),recursive=TRUE,showWarnings=FALSE)
write_table<-function(z,name,dir='Tables')write.csv(z,file.path(OUT_DIR,dir,paste0(name,'.csv')),row.names=FALSE,na='')
dir.create(file.path(OUT_DIR,'Audit/MaAsLin3'),recursive=TRUE,showWarnings=FALSE)
SEED<-666L;set.seed(SEED)
input<-file.path(PROJECT_DIR,'Output/physeq.rds');ps<-readRDS(input);m<-data.frame(sample_data(ps),check.names=FALSE);a<-as(otu_table(ps),'matrix');if(taxa_are_rows(ps))a<-t(a);storage.mode(a)<-'double';stopifnot(nrow(m)>1,identical(rownames(m),rownames(a)),!anyDuplicated(m$Patient_id),all(a>=0),all(is.finite(a)),max(abs(rowSums(a)-1))<1e-6)
m$subtype<-factor(ifelse(m$HER2=='Y','HER2-positive',ifelse(m$ER=='Y'|m$PR=='Y','Luminal','TNBC')),levels=c('Luminal','HER2-positive','TNBC'));stopifnot(all((m$subtype=='TNBC')==(m$TNBC=='Y')))
fields<-c(Surgery='Pre.enrollment_surgery',Chemotherapy='Pre.enrollment_Chemotherapy',Radiation='Pre.enrollment_Radiation',Hormonal_therapy='Pre.enrollment_Hormonal.Therapy',Targeted_therapy='Pre.enrollment_Targeted.Therapy')
for(v in names(fields))m[[v]]<-factor(m[[fields[[v]]]],levels=c('N','Y'))
oktx<-complete.cases(m[,fields,drop=FALSE]);m$treatment_count<-ifelse(oktx,rowSums(m[,fields,drop=FALSE]=='Y'),NA);m$treatment_group<-ifelse(is.na(m$treatment_count),'Missing',ifelse(m$treatment_count==0,'None',ifelse(m$treatment_count==1,'Single','Multiple')))
m$Shannon<-apply(a,1,function(z)-sum(z[z>0]*log(z[z>0])))
write_table(m[,c('Sample_id','Patient_id','ER','PR','HER2','TNBC','subtype','CTC_count','CTC_count_cat',names(fields),'treatment_count','treatment_group','Shannon')],'clinical_input','Audit')
n<-nrow(m);summaryrows<-list();add<-function(section,label,count){summaryrows[[length(summaryrows)+1]]<<-data.frame(section=section,label=label,n=count,percent=100*count/n,display=sprintf('%d (%.1f)',count,100*count/n))}
for(s in levels(m$subtype))add('subtype',s,sum(m$subtype==s))
for(v in names(fields))for(g in c('Y','N',NA)){nn<-if(is.na(g))sum(is.na(m[[v]]))else sum(m[[v]]==g,na.rm=TRUE);add(v,if(is.na(g))'Missing'else if(g=='Y')'Yes'else 'No',nn)}
for(s in c('None','Single','Multiple','Missing'))add('treatment_group',s,sum(m$treatment_group==s))
clinical<-do.call(rbind,summaryrows);write_table(clinical,'Table1_added_data')

params<-list(min_prevalence=.5,normalization='TSS',transform='LOG',augment=TRUE,standardize=TRUE,median_comparison_abundance=TRUE,median_comparison_prevalence=FALSE,max_significance=.25,cores=4L,correction='BH',seed=SEED)
writeLines(c('Approved scope: subtype Luminal/HER2-positive/TNBC plus five independent Yes/No treatment comparisons.','No adjustment covariates. No Single/Multiple modeling. Corrected physeq.rds input.','Subtype Shannon: Kruskal-Wallis. Treatment Shannon: Wilcoxon, using the original rstatix default exact policy (automatic when possible).','Bray-Curtis PERMANOVA: 999 permutations, original setting. Dispersion diagnostic: 999 permutations.','MaAsLin3 settings match the original 6_Maaslin.rmd; plots disabled only to prevent unwanted raster/PDF output. All abundance/prevalence/joint rows and warnings retained.','Package BH q values retained. Joint q is counted once per taxon/contrast, and individual component evidence is reported separately. No alternate-method or parameter search.'),file.path(OUT_DIR,'Analysis_plan.md'))
vars<-c('subtype',names(fields));alpha<-list();beta<-list();dispersion<-list();fits<-list();sample_manifest<-list();support<-list()
for(i in seq_along(vars)){
 v<-vars[i];ok<-!is.na(m[[v]]);md<-droplevels(m[ok,]);aa<-a[ok,,drop=FALSE];g<-droplevels(md[[v]])
 t<-if(v=='subtype')kruskal.test(md$Shannon,g) else wilcox.test(md$Shannon~g)
 alpha[[v]]<-data.frame(variable=v,n=nrow(md),group_sizes=paste(names(table(g)),as.integer(table(g)),sep='=',collapse='; '),test=if(v=='subtype')'Kruskal-Wallis'else 'Wilcoxon rank-sum',statistic=unname(t$statistic),p=t$p.value)
 ds<-vegdist(aa,method='bray');set.seed(SEED);b<-adonis2(ds~g,permutations=999);beta[[v]]<-data.frame(variable=v,n=nrow(md),df=b$Df[1],R2=b$R2[1],F=b$F[1],p=b$`Pr(>F)`[1],permutations=999)
 set.seed(SEED);bd<-betadisper(ds,g,bias.adjust=TRUE);pb<-permutest(bd,permutations=999);dispersion[[v]]<-data.frame(variable=v,n=nrow(md),F=pb$tab$F[1],p=pb$tab$`Pr(>F)`[1])
 sample_manifest[[v]]<-data.frame(analysis=v,sample_id=md$Sample_id,group=as.character(g))
 out<-file.path(OUT_DIR,'Audit/MaAsLin3',v);cfg<-list(abundance=aa,metadata=md[,v,drop=FALSE],formula=paste0('~',v),parameters=params,package=as.character(packageVersion('maaslin3')))
 cache<-file.path(out,'configuration.rds');results<-file.path(out,'all_results.tsv')
 if(file.exists(results)&&(!file.exists(cache)||!identical(readRDS(cache),cfg)))stop('Input/settings changed; archive existing MaAsLin3 output first: ',v)
 if(!file.exists(results)){
  set.seed(SEED);maaslin3(input_data=as.data.frame(aa),input_metadata=md[,v,drop=FALSE],output=out,formula=paste0('~',v),min_prevalence=.5,normalization='TSS',transform='LOG',augment=TRUE,standardize=TRUE,median_comparison_abundance=TRUE,median_comparison_prevalence=FALSE,max_significance=.25,cores=4,plot_summary_plot=FALSE,plot_associations=FALSE,max_pngs=0,verbosity='WARN');saveRDS(cfg,cache)
 }
 z<-read.delim(results,check.names=FALSE);z$analysis<-v;z$formula<-paste0('~',v);z$Species<-gsub('_',' ',sub('.*;','',z$feature));z$valid<-is.na(z$error)|z$error=='';z$individual_support<-z$valid&!is.na(z$qval_individual)&z$qval_individual<.25;z$joint_support<-z$valid&!is.na(z$qval_joint)&z$qval_joint<.25;fits[[v]]<-z
 keep<-colSums(aa>0)>.5*nrow(aa);write_table(data.frame(feature=colnames(aa),nonzero_n=colSums(aa>0),n=nrow(aa),retained=keep),paste0(v,'_filter'),'Audit')
 for(lev in levels(g)){ii<-g==lev;support[[paste(v,lev)]]<-data.frame(analysis=v,level=lev,feature=colnames(aa)[keep],n=sum(ii),nonzero_n=colSums(aa[ii,keep,drop=FALSE]>0))}
 cat('Completed ',v,'; n=',nrow(md),'\n')
}
alpha<-do.call(rbind,alpha);beta<-do.call(rbind,beta);dispersion<-do.call(rbind,dispersion);all<-do.call(rbind,fits)

alpha$q_BH_five_treatments<-NA_real_;beta$q_BH_five_treatments<-NA_real_;ii<-alpha$variable!='subtype';alpha$q_BH_five_treatments[ii]<-p.adjust(alpha$p[ii],'BH');ii<-beta$variable!='subtype';beta$q_BH_five_treatments[ii]<-p.adjust(beta$p[ii],'BH')
write_table(alpha,'Shannon_tests');write_table(beta,'PERMANOVA_tests');write_table(dispersion,'Dispersion_tests');write_table(all,'MaAsLin3_all_results');write_table(all[all$individual_support|all$joint_support,],'MaAsLin3_supported_rows');write_table(do.call(rbind,sample_manifest),'sample_inclusion','Audit');write_table(do.call(rbind,support),'nonzero_support_by_group','Audit')
counts<-do.call(rbind,lapply(vars,function(v){z<-fits[[v]];do.call(rbind,lapply(unique(z$value),function(val){zz<-z[z$value==val,];data.frame(analysis=v,contrast=val,n=unique(alpha$n[alpha$variable==v]),features=length(unique(zz$feature)),abundance_individual_n=sum(zz$individual_support&zz$model=='abundance'),prevalence_individual_n=sum(zz$individual_support&zz$model=='prevalence'),joint_n=sum(zz$joint_support&zz$model=='abundance'),warning_rows=sum(!zz$valid))}))}));write_table(counts,'MaAsLin3_evidence_counts')
w<-createWorkbook();tabs<-list(Table1_additions=clinical,Shannon=alpha,PERMANOVA=beta,Dispersion=dispersion,DA_counts=counts,DA_all=all,DA_supported=all[all$individual_support|all$joint_support,]);for(nm in names(tabs)){addWorksheet(w,nm);writeData(w,nm,tabs[[nm]]);freezePane(w,nm,firstRow=TRUE);setColWidths(w,nm,cols=seq_len(ncol(tabs[[nm]])),widths=pmin(48,pmax(12,nchar(names(tabs[[nm]]))+3)))};saveWorkbook(w,file.path(OUT_DIR,'Tables/Major4_results.xlsx'),overwrite=TRUE)
qc<-data.frame(check=c('38 unique patients','Luminal/HER2/TNBC=15/12/11','TNBC concordance','Each treatment n=36','No PNG or PDF output','All models univariate'),pass=c(nrow(m)==38&&!anyDuplicated(m$Patient_id),identical(as.integer(table(m$subtype)),c(15L,12L,11L)),all((m$subtype=='TNBC')==(m$TNBC=='Y')),all(alpha$n[alpha$variable!='subtype']==36),!any(grepl('\\.(png|pdf|tiff?)$',list.files(OUT_DIR,recursive=TRUE))),all(vapply(fits,function(z)length(unique(z$metadata))==1,logical(1)))));write_table(qc,'analysis_checks','QA');stopifnot(all(qc$pass))
write_table(data.frame(path=input,md5=unname(tools::md5sum(input))),'input_manifest','Audit');writeLines(capture.output(sessionInfo()),file.path(OUT_DIR,'Audit/sessionInfo.txt'));saveRDS(params,file.path(OUT_DIR,'Audit/approved_parameters.rds'))
writeLines(c('# Major4 results: confirmed univariate analyses',capture.output(print(clinical)),capture.output(print(alpha)),capture.output(print(beta)),capture.output(print(dispersion)),capture.output(print(counts)),'Only treatment categories with individual Yes/No records are modeled. None/Single/Multiple/Missing counts are for Table1. Luminal is the agreed label for HR-positive/HER2-negative, not a claim of measured Luminal A/B subtype.','Treatments overlap; models are unadjusted associations and are not causal treatment-effect estimates. Individual abundance and prevalence q values support their respective components; joint q alone does not specify the component. CTC categories remain exactly as the reference Word table. No additional covariate models were fitted.'),file.path(OUT_DIR,'Major4_results_notes.md'))
cat('Major4 analyses complete\n')
