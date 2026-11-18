#!/usr/bin/env Rscript
# Build the four-taxon MBR score and analyze survival, ROC and internal evaluation.
# Fig. 5B-F and Fig. S6; Data S6. Inputs: scripts 01 and 06.
# Fixed-four evaluation is conditional; upstream batch correction remains fixed.

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

OUT_DIR <- file.path(PROJECT_DIR,'Revision2','Fig5','output')
CACHE_DIR <- file.path(OUT_DIR,'Cache')
FONT_SIZE <- 11
N_SELECTION_REPEATS <- 50L
N_BOOTSTRAP <- 500L
N_OUTER_REPEATS <- 1L
ALPHAS <- seq(.1,.9,.1)
BASE_SEED <- 666L
PREVALENCE <- .10
PSEUDOCOUNT <- 1e-6
FREQUENCY_CUTOFF <- .50
FORCE <- tolower(Sys.getenv('FIG5_FORCE_RECOMPUTE','false'))=='true'
REDRAW_ONLY <- tolower(Sys.getenv('FIG5_REDRAW_ONLY','false'))=='true'
required <- c('phyloseq','survival','glmnet','ggplot2','cowplot','timeROC','ggbeeswarm','openxlsx')
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop('Missing packages (no automatic installation): ',paste(missing,collapse=', '))
suppressPackageStartupMessages({library(phyloseq);library(survival);library(glmnet);library(ggplot2);library(cowplot)})
options(stringsAsFactors=FALSE)
for(d in c('Main_Fig5','Supplementary_Fig5','Tables','Analysis','Cache','QA')) dir.create(file.path(OUT_DIR,d),recursive=TRUE,showWarnings=FALSE)
dir.create(CACHE_DIR,recursive=TRUE,showWarnings=FALSE)
ANA<-file.path(OUT_DIR,'Analysis');TBL<-file.path(ANA,'Audit_tables');dir.create(TBL,showWarnings=FALSE)
write_table<-function(x,name)write.csv(x,file.path(TBL,paste0(name,'.csv')),row.names=FALSE,na='')
species4<-c('Sphingopyxis_sp001990265','Bacillus_pumilus','Ralstonia_insidiosa','Paraburkholderia_tropica')
pretty_species<-function(z)gsub('_',' ',sub('.*;','',z))
hash_object<-function(z){f<-tempfile();saveRDS(z,f);h<-unname(tools::md5sum(f));unlink(f);h}
input_path<-file.path(PROJECT_DIR,'Output','physeq.rds')
input_md5<-unname(tools::md5sum(input_path))
ps<-readRDS(input_path);meta<-data.frame(sample_data(ps),check.names=FALSE)
raw<-as(otu_table(ps),'matrix');if(taxa_are_rows(ps))raw<-t(raw)
storage.mode(raw)<-'double'
stopifnot(identical(rownames(meta),rownames(raw)),!anyDuplicated(rownames(raw)),!anyNA(raw),all(raw>=0))

time<-as.numeric(meta$survival_months);event<-as.integer(meta$vital_status)
stopifnot(nrow(raw)>1,sum(event)>0,all(time>0),all(event%in%c(0,1)),!anyDuplicated(meta$Patient_id))
meta$age_binary<-as.numeric(meta$age_MBC_dx>=50)
meta$ER_positive<-as.numeric(meta$ER=='Y');meta$HER2_positive<-as.numeric(meta$HER2=='Y');meta$White<-as.numeric(meta$race=='White')
stopifnot(!anyNA(meta[,c('age_binary','ER_positive','HER2_positive','White')]))
taxa4<-vapply(species4,function(s){v<-grep(paste0(';',s,'$'),colnames(raw),value=TRUE);stopifnot(length(v)==1);v},character(1))
fold_ids<-function(status,k,seed){set.seed(seed);f<-integer(length(status));for(s in unique(status)){ix<-which(status==s);f[ix]<-sample(rep(1:k,length.out=length(ix)))};f}
recipe<-function(x)colnames(x)[colSums(x>0)>=ceiling(PREVALENCE*nrow(x))]
clr_apply<-function(x,features){z<-log(x[,features,drop=FALSE]+PSEUDOCOUNT);z-rowMeans(z)}
training_z<-function(tx,vx){
 if(!ncol(tx))return(list(train=tx,test=vx))
 mn<-colMeans(tx);ss<-apply(tx,2,sd)

 ss[ss==0]<-1
 list(train=sweep(sweep(tx,2,mn,'-'),2,ss,'/'),test=sweep(sweep(vx,2,mn,'-'),2,ss,'/'))
}
concord<-function(t,e,p){if(any(!is.finite(p)))return(NA_real_);survival::concordance(Surv(t,e)~p,reverse=TRUE)$concordance}
fit_cox<-function(x,t,e,newx=x,itermax=50,strict=TRUE){
 x<-as.matrix(x);newx<-as.matrix(newx)
 if(!ncol(x))return(list(fit=NULL,lp=rep(0,length(t)),test=rep(0,nrow(newx)),warning='',error='',unstable=FALSE))
 nm<-paste0('v',seq_len(ncol(x)));d<-data.frame(time=t,status=e,x,check.names=FALSE);names(d)<-c('time','status',nm)
 nd<-as.data.frame(newx);names(nd)<-nm;ww<-character()
 ans<-tryCatch(withCallingHandlers({f<-coxph(Surv(time,status)~.,data=d,x=TRUE,y=TRUE,control=coxph.control(iter.max=itermax));
  list(fit=f,lp=as.numeric(predict(f,type='lp')),test=as.numeric(predict(f,newdata=nd,type='lp')))},
  warning=function(w){ww<<-c(ww,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)list(error=conditionMessage(e)))
 ans$warning<-paste(ww,collapse=' | ');ans$unstable<-any(grepl('infinite|converg',ww,ignore.case=TRUE))
 if(is.null(ans$error))ans$error<-''
 if(!is.null(ans$fit)&&any(!is.finite(coef(ans$fit))))ans$error<-'Nonfinite or aliased coefficient'
 if(strict&&ans$unstable)ans$error<-ans$warning
 ans
}
extract_cox<-function(f,terms,unit){s<-summary(f);data.frame(term=terms,unit=unit,coefficient=coef(f),SE=s$coefficients[,3],HR=s$conf.int[,1],CI_low=s$conf.int[,3],CI_high=s$conf.int[,4],p=s$coefficients[,5],row.names=NULL)}
config<-list(input_md5=input_md5,repeats=N_SELECTION_REPEATS,alphas=ALPHAS,prevalence=PREVALENCE,pseudocount=PSEUDOCOUNT,cutoff=FREQUENCY_CUTOFF,base_seed=BASE_SEED,
 glmnet_version=as.character(packageVersion('glmnet')),survival_version=as.character(packageVersion('survival')),R_version=R.version.string)

discover<-function(rr,tt,ee,label){
 key<-hash_object(list(config=config,raw=rr,time=tt,event=ee));cache<-file.path(CACHE_DIR,paste0(label,'.rds'))
 if(!FORCE&&file.exists(cache)){v<-readRDS(cache);if(identical(v$key,key)){cat('Using verified cache:',label,'\n');return(v)}}
 if(REDRAW_ONLY)stop('Missing or incompatible cache: ',label)
 features<-recipe(rr);x<-clr_apply(rr,features);sel<-array(FALSE,dim=c(length(features),N_SELECTION_REPEATS,length(ALPHAS)),dimnames=list(features,NULL,as.character(ALPHAS)))
 rows<-list();foldrows<-list();warnrows<-list()
 for(r in seq_len(N_SELECTION_REPEATS)){
  fold<-fold_ids(ee,5,BASE_SEED+r*1000);foldrows[[r]]<-data.frame(repetition=r,sample_id=rownames(rr),fold=fold)
  for(j in seq_along(ALPHAS)){
   a<-ALPHAS[j];ww<-character()
   ans<-tryCatch(withCallingHandlers({cv<-cv.glmnet(x,Surv(tt,ee),family='cox',alpha=a,type.measure='C',foldid=fold,standardize=TRUE)
    f<-glmnet(x,Surv(tt,ee),family='cox',alpha=a,lambda=cv$lambda.min,standardize=TRUE);b<-as.matrix(coef(f))[,1]
    list(lambda=cv$lambda.min,C=max(cv$cvm),selected=b!=0)},warning=function(w){ww<<-c(ww,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)list(error=conditionMessage(e)))
   if(!is.null(ans$error))stop('Grid fit failed in ',label,': ',ans$error)
   sel[,r,j]<-ans$selected
   rows[[length(rows)+1]]<-data.frame(repetition=r,alpha=a,lambda=ans$lambda,internal_C=ans$C,n_selected=sum(ans$selected))
   if(length(ww))warnrows[[length(warnrows)+1]]<-data.frame(repetition=r,alpha=a,warning=paste(unique(ww),collapse=' | '))
  }
  if(r%%10==0){cat(label,':',r,'/',N_SELECTION_REPEATS,'repetitions complete\n');flush.console()}
 }
 cv<-do.call(rbind,rows);ag<-aggregate(internal_C~alpha,cv,mean);ag<-ag[order(-ag$internal_C,-ag$alpha),];a<-ag$alpha[1];j<-which.min(abs(ALPHAS-a))
 fr<-data.frame(taxon=features,selection_fraction=rowMeans(sel[,,j]),original_four=features%in%taxa4)
 fr$species<-sub('.*;','',fr$taxon);fr$stable<-fr$selection_fraction>FREQUENCY_CUTOFF
 fr<-fr[order(-fr$selection_fraction,fr$taxon),]
 v<-list(key=key,alpha=a,features=features,stable=fr$taxon[fr$stable],frequencies=fr,grid=cv,alpha_mean=ag,fold_manifest=do.call(rbind,foldrows),warnings=if(length(warnrows))do.call(rbind,warnrows) else data.frame(repetition=integer(),alpha=numeric(),warning=character()))
 saveRDS(v,cache);v
}
cat('Input:',nrow(raw),'patients;',sum(event),'events\n')
sortix<-order(rownames(raw))
discovery<-discover(raw[sortix,,drop=FALSE],time[sortix],event[sortix],'full_cohort_grid')
write_table(discovery$alpha_mean,'S01_alpha_grid_mean_C');write_table(discovery$grid,'S02_all_full_cohort_CV_runs')
write_table(discovery$frequencies,'S03_full_cohort_selection_frequencies');write_table(discovery$fold_manifest,'S04_full_cohort_CV_folds')
write_table(discovery$warnings,'S05_full_cohort_grid_warnings')
if(!setequal(discovery$stable,unname(taxa4)))stop('Discovery did not reproduce the original four taxa. Results saved; do not force the signature or redraw as four taxa.')
x<-clr_apply(raw,discovery$features)[,taxa4,drop=FALSE];colnames(x)<-species4
xs<-scale(x);mu<-attr(xs,'scaled:center');sds<-attr(xs,'scaled:scale')
fit<-fit_cox(x,time,event);zfit<-fit_cox(xs,time,event)
stopifnot(fit$error=='',zfit$error=='',max(abs(fit$lp-zfit$lp))<1e-8)
score<-fit$lp;group<-factor(ifelse(score>=median(score),'High risk','Low risk'),levels=c('Low risk','High risk'))
taxon_coef<-extract_cox(zfit$fit,species4,'Per SD of CLR abundance');taxon_coef$species<-pretty_species(taxon_coef$term)
taxon_coef$selection_fraction<-discovery$frequencies$selection_fraction[match(unname(taxa4),discovery$frequencies$taxon)]
taxon_coef$direction<-ifelse(taxon_coef$coefficient>0,'Higher risk','Lower risk')
write_table(taxon_coef,'S06_four_taxa_unit_scaled_Cox');write_table(extract_cox(fit$fit,species4,'Per CLR unit'),'S07_four_taxa_original_CLR_Cox')
write_table(data.frame(species=species4,CLR_mean=mu,CLR_SD=sds),'S08_scaling_parameters')
write_table(data.frame(sample_id=rownames(meta),patient_id=meta$Patient_id,time=time,event=event,MBR=score,risk_group=group),'S09_sample_MBR_scores')
write_table(data.frame(metric=c('Max LP difference after standardization','Apparent C'),value=c(max(abs(fit$lp-zfit$lp)),concord(time,event,score))),'S10_scaling_checks')
ph<-cox.zph(fit$fit)$table;write_table(data.frame(term=c(species4,'GLOBAL'),ph,check.names=FALSE),'S11_proportional_hazards')
hrrows<-list()
for(u in c('Per one MBR unit','Per SD MBR','Per IQR MBR','High versus low MBR')){
 v<-switch(u,'Per one MBR unit'=score,'Per SD MBR'=as.numeric(scale(score)),'Per IQR MBR'=score/IQR(score),'High versus low MBR'=as.numeric(group=='High risk'))
 f<-fit_cox(matrix(v,ncol=1),time,event);hrrows[[u]]<-extract_cox(f$fit,'MBR',u)
}
hr<-do.call(rbind,hrrows);write_table(hr,'S12_MBR_HR_scales')
lr<-survdiff(Surv(time,event)~group);logrank_p<-pchisq(lr$chisq,df=1,lower.tail=FALSE)
km<-survfit(Surv(time,event)~group)
write_table(data.frame(test='Log-rank',p=logrank_p),'S13_KM_logrank')

firsttax<-as.data.frame(tax_table(ps));otaxa<-rownames(firsttax)[firsttax$Species%in%species4]
oldx<-clr_apply(raw,discovery$features)[,otaxa,drop=FALSE];colnames(oldx)<-gsub('[^A-Za-z0-9]','_',colnames(oldx))
dd<-data.frame(survival_months=time,vital_status=event,oldx,check.names=FALSE)
ff<-as.formula(paste('Surv(survival_months,vital_status) ~',paste(colnames(oldx),collapse=' + ')))
quiet_old<-function(tr,nd=tr){ww<-character();ans<-tryCatch(withCallingHandlers({f<-coxph(ff,data=tr);list(fit=f,train=as.numeric(predict(f,newdata=tr,type='lp')),test=as.numeric(predict(f,newdata=nd,type='lp')))},warning=function(w){ww<<-c(ww,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)list(error=conditionMessage(e)));ans$warning<-paste(ww,collapse=' | ');ans}
app<-concord(time,event,score);set.seed(20260612);boots<-vector('list',N_BOOTSTRAP)
boot_rank<-matrix(NA_real_,nrow(raw),N_BOOTSTRAP)
boot_high<-matrix(NA_real_,nrow(raw),N_BOOTSTRAP)
boot_beta<-matrix(NA_real_,4,N_BOOTSTRAP)
for(b in seq_len(N_BOOTSTRAP)){
 ix<-sample(seq_along(time),replace=TRUE);f<-quiet_old(dd[ix,,drop=FALSE],dd)
 cb<-co<-NA_real_;if(is.null(f$error)){
  cb<-concord(time[ix],event[ix],f$train);co<-concord(time,event,f$test)
  if(all(is.finite(f$test))){boot_rank[,b]<-(rank(f$test,ties.method='average')-.5)/length(time);boot_high[,b]<-as.numeric(f$test>=median(f$test))}
  boot_beta[,b]<-coef(f$fit)[match(unname(taxa4),otaxa)]
 }
 boots[[b]]<-data.frame(iteration=b,C_bootstrap=cb,C_original=co,optimism=cb-co,warning=f$warning,error=if(is.null(f$error))'' else f$error)
}
boots<-do.call(rbind,boots);loo<-vapply(seq_along(time),function(i){f<-quiet_old(dd[-i,,drop=FALSE],dd[i,,drop=FALSE]);if(is.null(f$error))f$test else NA_real_},numeric(1))
validation<-data.frame(analysis=c('Apparent','Fixed-four bootstrap corrected','Fixed-four LOOCV'),C=c(app,app-mean(boots$optimism,na.rm=TRUE),concord(time,event,loo)),scope='Fixed signature; selection not repeated')
write_table(boots,'S14_fixed_four_bootstrap_attempts');write_table(validation,'S15_fixed_four_validation')
write_table(data.frame(sample_id=rownames(meta),time=time,event=event,prediction=loo),'S16_fixed_four_LOOCV_predictions')
original_rank<-(rank(score,ties.method='average')-.5)/length(time)
brho<-apply(boot_rank,2,function(v)if(all(is.finite(v)))cor(original_rank,v,method='spearman') else NA_real_)
bagree<-apply(boot_high,2,function(v)if(all(is.finite(v)))mean(v==as.numeric(group=='High risk')) else NA_real_)
bootstrap_stability<-data.frame(iteration=seq_len(N_BOOTSTRAP),rank_Spearman=brho,median_group_agreement=bagree,warning=boots$warning,error=boots$error,
 scope='Conditional fixed-four refitting in development cohort; discovery not repeated')
valid_rank<-colSums(is.finite(boot_rank))==nrow(raw)
qrank<-t(apply(boot_rank[,valid_rank,drop=FALSE],1,quantile,probs=c(.025,.5,.975),na.rm=TRUE))
subject_stability<-data.frame(sample_id=rownames(meta),original_MBR=score,original_percentile=original_rank,risk_group=group,
 percentile_low=qrank[,1],percentile_median=qrank[,2],percentile_high=qrank[,3],bootstrap_high_group_fraction=rowMeans(boot_high[,valid_rank,drop=FALSE]))
conditional_summary<-data.frame(metric=c('Attempted bootstrap fits','Finite rank predictions','Warned bootstrap fits','Median rank Spearman','Rank Spearman 2.5th percentile','Rank Spearman 97.5th percentile','Median group agreement'),
 value=c(N_BOOTSTRAP,sum(valid_rank),sum(nzchar(boots$warning)),median(brho,na.rm=TRUE),quantile(brho,.025,na.rm=TRUE),quantile(brho,.975,na.rm=TRUE),median(bagree,na.rm=TRUE)),
 scope='Conditional on the four taxa and development cohort; no claim of independent validation')
direction_summary<-data.frame(species=species4,sign_agreement=rowMeans(sign(boot_beta)==sign(coef(fit$fit)),na.rm=TRUE),finite_coefficients=rowSums(is.finite(boot_beta)))
write_table(bootstrap_stability,'S36_conditional_MBR_bootstrap_stability');write_table(subject_stability,'S37_subject_MBR_percentile_intervals');write_table(conditional_summary,'S38_conditional_stability_summary');write_table(direction_summary,'S39_bootstrap_coefficient_signs')
set.seed(20260612);perm<-replicate(1000,{ii<-sample(seq_along(time));concord(time[ii],event[ii],score)})
permp<-(sum(perm>=app)+1)/1001
write_table(data.frame(permutation=seq_along(perm),C=perm),'S17_fixed_score_permutation')
write_table(data.frame(test='Outcome permutation relative to fixed learned score',p=permp,scope='Descriptive conditional test; discovery/fitting not rerun under permutation'),'S18_permutation_scope')

base<-as.matrix(meta[,c('age_binary','ER_positive','White')]);h2<-as.matrix(meta[,c('age_binary','ER_positive','White','HER2_positive')])
specs<-list(Clinical_original=base,Clinical_HER2=h2,MBR=matrix(score,ncol=1,dimnames=list(NULL,'MBR')),Combined_original=cbind(MBR=score,base),Combined_HER2=cbind(MBR=score,h2))
roc<-list();model_coef<-list();model_diag<-list();aucrows<-list()
for(n in names(specs)){
 f<-fit_cox(specs[[n]],time,event,strict=FALSE)
 model_diag[[n]]<-data.frame(model=n,n=length(time),events=sum(event),n_predictors=ncol(specs[[n]]),events_per_predictor=sum(event)/ncol(specs[[n]]),warnings=f$warning,error=f$error,unstable=f$unstable)
 if(f$error!=''||f$unstable)stop('Full-cohort clinical model unstable: ',n)
 model_coef[[n]]<-cbind(model=n,extract_cox(f$fit,colnames(specs[[n]]),'Original clinical coding; MBR per score unit'))
 roc[[n]]<-timeROC::timeROC(T=time,delta=event,marker=f$lp,cause=1,times=c(12,36),iid=FALSE)
 aucrows[[n]]<-data.frame(model=n,time_months=c(12,36),apparent_AUC=roc[[n]]$AUC[match(c(12,36),roc[[n]]$times)])
}
auc<-do.call(rbind,aucrows);write_table(auc,'S19_clinical_HER2_apparent_AUC');write_table(do.call(rbind,model_coef),'S20_clinical_Cox_coefficients');write_table(do.call(rbind,model_diag),'S21_clinical_model_diagnostics')

set.seed(BASE_SEED);te<-unlist(lapply(c(0,1),function(s){ii<-which(event==s);sample(ii,ceiling(.25*length(ii)))}));tr<-setdiff(seq_along(time),te)
hold<-discover(raw[tr,,drop=FALSE],time[tr],event[tr],'holdout_training_grid')
htx<-clr_apply(raw[tr,,drop=FALSE],hold$features)[,hold$stable,drop=FALSE];hvx<-clr_apply(raw[te,,drop=FALSE],hold$features)[,hold$stable,drop=FALSE]
hz<-training_z(htx,hvx);hf<-fit_cox(hz$train,time[tr],event[tr],hz$test)
hc<-if(hf$error=='')concord(time[te],event[te],hf$test) else NA_real_
holdsum<-data.frame(n_train=length(tr),events_train=sum(event[tr]),n_test=length(te),events_test=sum(event[te]),alpha=hold$alpha,n_selected=length(hold$stable),n_original_four=sum(hold$stable%in%taxa4),selected_species=paste(pretty_species(hold$stable),collapse='; '),C=hc,error=hf$error,scope='Downstream train-only selection; full-cohort upstream batch correction fixed')
write_table(holdsum,'S22_holdout_summary');write_table(hold$frequencies,'S23_holdout_selection_frequencies');write_table(hold$alpha_mean,'S24_holdout_alpha_grid')
write_table(hold$fold_manifest,'S25_holdout_training_CV_folds');write_table(data.frame(sample_id=rownames(meta),partition=ifelse(seq_along(time)%in%te,'Test','Train'),time=time,event=event),'S26_holdout_partition')
if(hf$error=='')write_table(data.frame(sample_id=rownames(meta)[te],time=time[te],event=event[te],prediction=hf$test),'S27_holdout_test_predictions')

outerrows<-list();outerdiag<-list();outerfreq<-list();outerman<-list()
for(r in seq_len(N_OUTER_REPEATS)){
 outerfold<-fold_ids(event,5,BASE_SEED+r*1000)
 outerman[[r]]<-data.frame(repetition=r,sample_id=rownames(meta),fold=outerfold)
 for(k in 1:5){
  it<-which(outerfold!=k);iv<-which(outerfold==k);lab<-sprintf('outer_%02d_fold_%02d_grid',r,k)
  ds<-discover(raw[it,,drop=FALSE],time[it],event[it],lab)
  write.csv(ds$grid,file.path(ANA,paste0(lab,'_CV_runs.csv')),row.names=FALSE)
  write.csv(ds$fold_manifest,file.path(ANA,paste0(lab,'_inner_fold_manifest.csv')),row.names=FALSE)
  write.csv(ds$warnings,file.path(ANA,paste0(lab,'_glmnet_warnings.csv')),row.names=FALSE)
  tx<-clr_apply(raw[it,,drop=FALSE],ds$features)[,ds$stable,drop=FALSE];vx<-clr_apply(raw[iv,,drop=FALSE],ds$features)[,ds$stable,drop=FALSE]
  tz<-training_z(tx,vx);f<-fit_cox(tz$train,time[it],event[it],tz$test)
  pp<-if(f$error=='')f$test else rep(NA_real_,length(iv))
  numerical<-if(is.null(f$test))rep(NA_real_,length(iv)) else f$test
  outerrows[[length(outerrows)+1]]<-data.frame(repetition=r,fold=k,sample_id=rownames(meta)[iv],time=time[iv],event=event[iv],prediction=pp,numerical_prediction_even_if_unstable=numerical)
  outerdiag[[length(outerdiag)+1]]<-data.frame(repetition=r,fold=k,alpha=ds$alpha,n_train=length(it),n_test=length(iv),events_train=sum(event[it]),events_test=sum(event[iv]),n_selected=length(ds$stable),n_original_four=sum(ds$stable%in%taxa4),selected_species=paste(pretty_species(ds$stable),collapse='; '),C_within_fold=if(f$error=='')concord(time[iv],event[iv],pp) else NA_real_,C_numerical_even_if_unstable=concord(time[iv],event[iv],numerical),error=f$error,warning=f$warning)
  outerfreq[[length(outerfreq)+1]]<-data.frame(repetition=r,fold=k,ds$frequencies)
  write.csv(do.call(rbind,outerdiag),file.path(ANA,'outer_progress.csv'),row.names=FALSE)
  cat('Completed matching discovery-pipeline outer fold',k,'of 5\n');flush.console()
 }
}
op<-do.call(rbind,outerrows);od<-do.call(rbind,outerdiag)
outer_summary<-do.call(rbind,lapply(split(op,op$repetition),function(p)data.frame(repetition=p$repetition[1],complete=nrow(p)==nrow(raw)&&all(is.finite(p$prediction)),pooled_C=if(nrow(p)==nrow(raw)&&all(is.finite(p$prediction)))concord(p$time,p$event,p$prediction) else NA_real_,scope='Matching grid-frequency-Cox pipeline; upstream batch correction fixed')))

outer_summary$pair_weighted_within_fold_C<-vapply(outer_summary$repetition,function(r){
 pp<-op[op$repetition==r,];if(nrow(pp)!=nrow(raw)||any(!is.finite(pp$prediction)))return(NA_real_)
 cc<-lapply(split(pp,pp$fold),function(z){p<-z$prediction;survival::concordance(Surv(z$time,z$event)~p,reverse=TRUE)$count})
 den<-sum(vapply(cc,function(z)sum(z[1:3]),numeric(1)))
 sum(vapply(cc,function(z)z[1]+.5*z[3],numeric(1)))/den
},numeric(1))
outer_summary$failed_folds<-vapply(outer_summary$repetition,function(r)sum(od$error[od$repetition==r]!=''),integer(1))
outer_summary$numerical_C_including_unstable_fits<-vapply(outer_summary$repetition,function(r){
 pp<-op[op$repetition==r,];if(nrow(pp)!=nrow(raw)||any(!is.finite(pp$numerical_prediction_even_if_unstable)))return(NA_real_)
 cc<-lapply(split(pp,pp$fold),function(z){p<-z$numerical_prediction_even_if_unstable;survival::concordance(Surv(z$time,z$event)~p,reverse=TRUE)$count})
 den<-sum(vapply(cc,function(z)sum(z[1:3]),numeric(1)))
 sum(vapply(cc,function(z)z[1]+.5*z[3],numeric(1)))/den
},numeric(1))
write_table(op,'S28_matching_pipeline_outer_predictions');write_table(od,'S29_matching_pipeline_outer_signatures');write_table(do.call(rbind,outerfreq),'S30_matching_pipeline_outer_frequencies');write_table(do.call(rbind,outerman),'S31_matching_pipeline_outer_folds');write_table(outer_summary,'S32_matching_pipeline_outer_summary')

penalized_refit<-function(tx,vx,tt,ee,seed){
 if(ncol(tx)==0)return(list(prediction=rep(0,nrow(vx)),lambda=NA_real_,warning='',error='',fit=NULL))

 if(ncol(tx)==1){tx<-cbind(tx,zero_padding=0);vx<-cbind(vx,zero_padding=0)}
 ww<-character()
 ans<-tryCatch(withCallingHandlers({
  fd<-fold_ids(ee,5,seed)
  cf<-cv.glmnet(tx,Surv(tt,ee),family='cox',alpha=0,type.measure='deviance',foldid=fd,standardize=TRUE)
  list(prediction=as.numeric(predict(cf,newx=vx,s='lambda.1se',type='link')),lambda=cf$lambda.1se,fit=cf,folds=fd,error='')
 },warning=function(w){ww<<-c(ww,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)list(prediction=rep(NA_real_,nrow(vx)),lambda=NA_real_,fit=NULL,error=conditionMessage(e)))
 ans$warning<-paste(unique(ww),collapse=' | ');ans
}
prrows<-list();prdiag<-list();prmodels<-list()
for(r in seq_len(N_OUTER_REPEATS))for(k in 1:5){
 mf<-outerman[[r]];it<-which(mf$fold!=k);iv<-which(mf$fold==k)
 ds<-readRDS(file.path(CACHE_DIR,sprintf('outer_%02d_fold_%02d_grid.rds',r,k)))
 tx<-clr_apply(raw[it,,drop=FALSE],ds$features)[,ds$stable,drop=FALSE];vx<-clr_apply(raw[iv,,drop=FALSE],ds$features)[,ds$stable,drop=FALSE]
 zz<-training_z(tx,vx);pf<-penalized_refit(zz$train,zz$test,time[it],event[it],BASE_SEED+50000+r*100+k)
 prrows[[length(prrows)+1]]<-data.frame(repetition=r,fold=k,sample_id=rownames(meta)[iv],time=time[iv],event=event[iv],prediction=pf$prediction)
 prdiag[[length(prdiag)+1]]<-data.frame(repetition=r,fold=k,n_selected=length(ds$stable),ridge_alpha=0,lambda_rule='lambda.1se; training-only deviance CV',lambda=pf$lambda,C=concord(time[iv],event[iv],pf$prediction),warning=pf$warning,error=pf$error)
 prmodels[[sprintf('outer_%d_%d',r,k)]]<-pf
 if(!is.null(pf$folds))write.csv(data.frame(sample_id=rownames(meta)[it],fold=pf$folds),file.path(ANA,sprintf('penalized_refit_inner_folds_%d_%d.csv',r,k)),row.names=FALSE)
}
prp<-do.call(rbind,prrows);prd<-do.call(rbind,prdiag)
pair_weighted<-function(pp){
 if(nrow(pp)!=nrow(raw)||any(!is.finite(pp$prediction)))return(NA_real_)
 cc<-lapply(split(pp,pp$fold),function(z){p<-z$prediction;survival::concordance(Surv(z$time,z$event)~p,reverse=TRUE)$count})
 sum(vapply(cc,function(z)z[1]+.5*z[3],numeric(1)))/sum(vapply(cc,function(z)sum(z[1:3]),numeric(1)))
}
prsummary<-data.frame(analysis='Post-selection ridge refit sensitivity',n=nrow(raw),events=sum(event),successful_folds=sum(prd$error==''&is.finite(prd$C)),total_folds=nrow(prd),
 within_fold_C=pair_weighted(prp),pooled_C=concord(prp$time,prp$event,prp$prediction),scope='Modified regularized estimator; upstream correction fixed; not validation of original Cox MBR')
phridge<-penalized_refit(hz$train,hz$test,time[tr],event[tr],BASE_SEED+60000)
prhold<-data.frame(analysis='Holdout post-selection ridge sensitivity',n=length(te),events=sum(event[te]),C=concord(time[te],event[te],phridge$prediction),lambda=phridge$lambda,warning=phridge$warning,error=phridge$error,
 scope='Secondary modified-estimator sensitivity on the previously examined holdout; not an untouched confirmation test')
write_table(prp,'S40_penalized_refit_outer_predictions');write_table(prd,'S41_penalized_refit_fold_results');write_table(prsummary,'S42_penalized_refit_summary');write_table(prhold,'S43_penalized_refit_holdout')
saveRDS(prmodels,file.path(ANA,'penalized_refit_models.rds'))
cat('Targeted ridge sensitivity:\n');print(prsummary);print(prhold)

contrib<-sweep(sweep(x,2,colMeans(x),'-'),2,coef(fit$fit),'*')
stopifnot(max(abs(rowSums(contrib)-(score-mean(score))))<1e-8)
write_table(data.frame(sample_id=rownames(meta),contrib,check.names=FALSE),'S33_exact_linear_SHAP')
write_table(data.frame(sample_id=rownames(meta),xs,check.names=FALSE),'S34_four_taxa_CLR_z_matrix')

NAVY<-'#3B6FB6';CORAL<-'#D95F59';TEAL<-'#299D91';PURPLE<-NAVY;ROSE<-CORAL;INK<-'#222222';GREY<-'#A0A0A0'
riskcols<-c('Low risk'=NAVY,'High risk'=CORAL)
dircols<-c('Higher risk'=CORAL,'Lower risk'=NAVY)
theme_fig<-function(size=FONT_SIZE)theme_classic(base_size=FONT_SIZE,base_family='sans')+theme(text=element_text(colour=INK,size=FONT_SIZE),axis.text=element_text(colour=INK,size=FONT_SIZE),axis.title=element_text(size=FONT_SIZE),strip.text=element_text(size=FONT_SIZE),axis.line=element_line(linewidth=.45),axis.ticks=element_line(linewidth=.4),plot.title=element_blank(),plot.margin=margin(5,6,5,5),legend.title=element_text(size=FONT_SIZE),legend.text=element_text(size=FONT_SIZE))
inside<-function(x,y)theme(legend.position='inside',legend.position.inside=c(x,y),legend.justification=c(1,0),legend.background=element_rect(fill=scales::alpha('white',.9),colour=NA),legend.margin=margin(1,2,1,2),legend.key.height=grid::unit(.3,'cm'),legend.key.width=grid::unit(.45,'cm'))
save_plot<-function(p,folder,name,w,h){

  tmp<-tempfile(fileext='.pdf')
   if(isTRUE(capabilities('aqua'))){

    native_pdf<-function(filename,width,height,bg='white',...)grDevices::quartz(type='pdf',file=filename,width=width,height=height,bg=bg,family='Arial')
    ggsave(tmp,p,width=w,height=h,device=native_pdf,bg='white')
   } else ggsave(tmp,p,width=w,height=h,device=grDevices::pdf,useDingbats=FALSE,bg='white')
  dest<-file.path(OUT_DIR,folder,paste0(name,'.pdf'));part<-paste0(dest,'.writing')
  stopifnot(file.copy(tmp,part,overwrite=TRUE),file.rename(part,dest));unlink(tmp)
}
order_species<-taxon_coef$species[order(-taxon_coef$coefficient)]
taxon_coef$species<-factor(taxon_coef$species,levels=rev(order_species))
ay<-scale_y_discrete(drop=FALSE,expand=expansion(add=c(.95,.55)))
paf<-ggplot(taxon_coef,aes(selection_fraction,species,colour=direction))+geom_segment(aes(x=0,xend=selection_fraction,yend=species),colour='grey80',linewidth=.65)+geom_point(size=2.8)+scale_colour_manual(values=dircols,guide='none')+geom_vline(xintercept=.5,linetype=2,colour='grey55',linewidth=.4)+scale_x_continuous(limits=c(0,1),breaks=c(0,.5,1),labels=c('0','50','100'))+ay+labs(x='Selection frequency (%)',y=NULL)+theme_fig()+theme(axis.text.y=element_text(face='italic',size=FONT_SIZE))
pah<-ggplot(taxon_coef,aes(HR,species,colour=direction))+geom_vline(xintercept=1,linetype=2,colour='grey55',linewidth=.4)+geom_segment(aes(x=CI_low,xend=CI_high,yend=species),linewidth=.7)+geom_point(size=2.6)+scale_x_log10(limits=c(.15,6),breaks=c(.2,.5,1,2,5))+ay+scale_colour_manual(values=dircols,guide='none')+labs(x='HR per SD (95% CI)',y=NULL)+theme_fig()+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
pac<-ggplot(taxon_coef,aes(coefficient,species,fill=direction))+geom_vline(xintercept=0,linetype=2,colour='grey55',linewidth=.4)+geom_col(width=.55,colour=INK,linewidth=.25)+ay+scale_fill_manual(values=dircols)+scale_x_continuous(expand=expansion(mult=c(.08,.12)))+labs(x='Coefficient per SD',y=NULL,fill=NULL)+theme_fig()+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())+inside(.98,.025)
pa<-plot_grid(paf,pah,pac,nrow=1,rel_widths=c(1.55,.95,1.1),align='h',axis='tb')

ks<-summary(km,censored=TRUE)
kd<-data.frame(time=ks$time,surv=ks$surv,strata=sub('group=','',as.character(ks$strata)),n_censor=ks$n.censor)
kd<-rbind(data.frame(time=0,surv=1,strata=levels(group),n_censor=0),kd);kd$strata<-factor(kd$strata,levels=levels(group))
hrg<-hr[hr$unit=='High versus low MBR',]
kmlabel<-sprintf('HR = %.2f (95%% CI %.2f-%.2f)\nLog-rank P = %.3f',hrg$HR,hrg$CI_low,hrg$CI_high,logrank_p)
pb<-ggplot(kd,aes(time,surv,colour=strata))+geom_step(linewidth=.85)+geom_point(data=kd[kd$n_censor>0,],shape=3,size=1.7,stroke=.45)+scale_colour_manual(values=riskcols)+scale_x_continuous(limits=c(0,130),breaks=seq(0,125,25),expand=expansion(mult=c(.01,.01)))+scale_y_continuous(limits=c(0,1.02),breaks=seq(0,1,.25),expand=c(0,0))+annotate('text',x=2,y=.075,label=kmlabel,hjust=0,vjust=0,size=FONT_SIZE/(72.27/25.4),family='sans')+labs(x='Time (months)',y='Overall survival probability',colour=NULL)+theme_fig()+inside(.98,.79)
rt<-seq(0,125,25);rk<-summary(km,times=rt,extend=TRUE)
rd<-data.frame(time=rk$time,n=rk$n.risk,group=sub('group=','',as.character(rk$strata)))
rd$group<-factor(rd$group,levels=rev(levels(group)))
prisk<-ggplot(rd,aes(time,group,label=n,colour=group))+geom_text(size=FONT_SIZE/(72.27/25.4),show.legend=FALSE)+scale_colour_manual(values=riskcols)+scale_x_continuous(limits=c(0,130),breaks=rt,expand=expansion(mult=c(.01,.01)))+labs(x=NULL,y=NULL)+theme_fig()+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(size=FONT_SIZE),plot.margin=margin(0,6,0,5))
pbfull<-plot_grid(pb,prisk,ncol=1,rel_heights=c(1,.18),align='v',axis='lr')
rocmods<-c('Clinical_HER2','MBR','Combined_HER2');rocnames<-c('Age + ER + Race + HER2','MBR','Clinical + MBR');roccols<-setNames(c('#1B9E77','#8C6BB1','#F07818'),rocnames)
rrows<-list();rlabels<-character()
for(i in seq_along(rocmods)){a<-roc[[rocmods[i]]];j<-match(36,a$times);rrows[[i]]<-data.frame(FP=a$FP[,j],TP=a$TP[,j],model=rocnames[i]);rlabels[i]<-sprintf('%s (AUC %.3f)',rocnames[i],a$AUC[j])}
rr<-do.call(rbind,rrows);rr$model<-factor(rr$model,levels=rocnames)
pc<-ggplot(rr,aes(FP,TP,colour=model))+geom_abline(intercept=0,slope=1,linetype=2,colour='grey65',linewidth=.45)+geom_path(linewidth=.85)+coord_equal(xlim=c(0,1),ylim=c(0,1),expand=FALSE)+scale_colour_manual(values=roccols,labels=rlabels)+scale_x_continuous(breaks=seq(0,1,.25))+scale_y_continuous(breaks=seq(0,1,.25))+labs(x='1 - Specificity',y='Sensitivity',colour=NULL)+theme_fig()+inside(.98,.025)+theme(legend.text=element_text(size=FONT_SIZE))
write_table(rr,'S35_Fig5C_ROC_coordinates')

datcon<-data.frame(sample_id=rep(rownames(meta),times=4),species=rep(pretty_species(species4),each=nrow(raw)),contribution=as.vector(contrib),CLR_z=as.vector(xs))
datcon$species<-factor(datcon$species,levels=rev(order_species))
pd1<-ggplot(datcon,aes(contribution,species,fill=CLR_z))+geom_vline(xintercept=0,colour='grey70',linewidth=.4)+ggbeeswarm::geom_quasirandom(groupOnX=FALSE,width=.25,size=1.8,alpha=.95,shape=21,colour='#554F5A',stroke=.2)+scale_fill_gradient2(low=PURPLE,mid='#F5F3F5',high=ROSE,midpoint=0,limits=c(-2,2),oob=scales::squish,name='CLR z-score')+scale_y_discrete(drop=FALSE)+labs(x='Contribution to MBR (SHAP)',y=NULL)+theme_fig()+theme(axis.text.y=element_text(face='italic',size=FONT_SIZE),legend.position='bottom',legend.key.width=grid::unit(.5,'cm'),legend.key.height=grid::unit(.18,'cm'),legend.title=element_text(size=FONT_SIZE),legend.text=element_text(size=FONT_SIZE))
ord<-order(group,-score);heat<-data.frame(patient=rep(seq_along(ord),times=4),species=rep(pretty_species(species4),each=nrow(raw)),value=as.vector(xs[ord,,drop=FALSE]))
heat$species<-factor(heat$species,levels=rev(order_species))
pd2<-ggplot(heat,aes(patient,species,fill=value))+geom_tile(colour='white',linewidth=.18)+scale_fill_gradient2(low=PURPLE,mid='white',high=ROSE,midpoint=0,limits=c(-2,2),oob=scales::squish,name='CLR z-score')+scale_x_continuous(expand=c(0,0),limits=c(.5,nrow(raw)+.5))+scale_y_discrete(expand=c(0,0),drop=FALSE)+labs(x='Patients ordered by risk group and MBR',y=NULL)+theme_fig()+theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.line=element_blank(),legend.position='right',legend.key.height=grid::unit(.3,'cm'),legend.key.width=grid::unit(.18,'cm'),legend.text=element_text(size=FONT_SIZE),legend.title=element_text(size=FONT_SIZE))
strip<-data.frame(patient=seq_along(ord),group=group[ord])
pstrip<-ggplot(strip,aes(patient,1,fill=group))+geom_tile(colour='white',linewidth=.18)+scale_fill_manual(values=riskcols)+scale_x_continuous(expand=c(0,0),limits=c(.5,nrow(raw)+.5))+scale_y_continuous(expand=c(0,0))+labs(fill=NULL)+theme_void(base_family='sans')+theme(legend.position='none',plot.margin=margin(0,6,0,5))
pstrip<-pstrip+geom_text(data=aggregate(patient~group,strip,mean),aes(patient,1,label=group),inherit.aes=FALSE,colour='white',size=FONT_SIZE/(72.27/25.4))

pd2full<-plot_grid(pstrip,pd2,ncol=1,rel_heights=c(.08,1),align='v',axis='lr')
pd1<-pd1+theme(legend.position='none')
pd1full<-plot_grid(NULL,pd1,ncol=1,rel_heights=c(.08,1))
pd<-plot_grid(pd1full,pd2full,nrow=1,rel_widths=c(1.08,1),align='h',axis='tb')
save_plot(pa,'Main_Fig5','Fig5B',12.5,2.8);save_plot(pbfull,'Main_Fig5','Fig5C',6.8,4.4);save_plot(pc,'Main_Fig5','Fig5D',5.5,4.4);save_plot(pd,'Main_Fig5','Fig5E',12.5,3.1)

ag<-discovery$alpha_mean
pag<-ggplot(ag,aes(alpha,internal_C))+geom_line(colour=NAVY,linewidth=.7)+geom_point(size=2,colour=NAVY)+geom_point(data=ag[ag$alpha==discovery$alpha,],colour=CORAL,size=3)+scale_x_continuous(breaks=ALPHAS)+labs(x=expression(alpha),y='Mean internal CV C-index')+theme_fig()+annotate('text',x=.12,y=max(ag$internal_C),label=sprintf('Selected alpha = %.1f',discovery$alpha),hjust=0,vjust=1,size=FONT_SIZE/(72.27/25.4))
top<-head(discovery$frequencies,10);ppf<-ggplot(top,aes(selection_fraction,reorder(pretty_species(species),selection_fraction),fill=original_four))+geom_col(width=.65)+geom_vline(xintercept=.5,linetype=2,colour='grey45')+scale_fill_manual(values=c('TRUE'=TEAL,'FALSE'='#B8D8D3'),guide='none')+scale_x_continuous(limits=c(0,1),labels=function(z)paste0(100*z,'%'))+labs(x='Full-cohort selection frequency',y=NULL)+theme_fig(9)+theme(axis.text.y=element_text(face='italic',size=FONT_SIZE))
pv<-ggplot(validation,aes(C,reorder(analysis,-C),colour=analysis))+geom_segment(aes(x=.5,xend=C,yend=reorder(analysis,-C)),linewidth=.7,alpha=.65)+geom_point(size=3)+scale_colour_manual(values=setNames(c('#7366AD',TEAL,'#F07818'),validation$analysis),guide='none')+geom_text(aes(label=sprintf('%.3f',C)),nudge_y=.22,size=FONT_SIZE/(72.27/25.4))+scale_x_continuous(limits=c(.5,1),breaks=seq(.5,1,.1))+labs(x='C-index: fixed four-taxon model',y=NULL)+theme_fig(9)
popt<-ggplot(boots,aes(optimism))+geom_histogram(bins=25,fill=TEAL,colour='white',linewidth=.2)+geom_vline(xintercept=mean(boots$optimism,na.rm=TRUE),colour=CORAL,linewidth=.7)+labs(x='Bootstrap optimism (fixed four taxa)',y='Resamples')+theme_fig(9)

phf<-ggplot(bootstrap_stability,aes(rank_Spearman))+geom_histogram(bins=25,fill='#7366AD',colour='white',linewidth=.2)+geom_vline(xintercept=median(brho,na.rm=TRUE),colour=INK,linewidth=.7)+labs(x='Rank correlation with original MBR',y='Bootstrap refits')+theme_fig()
ss<-subject_stability[order(subject_stability$original_percentile),];ss$patient_order<-seq_len(nrow(ss))
po<-ggplot(ss,aes(patient_order,percentile_median,colour=risk_group))+geom_hline(yintercept=.5,linetype=2,colour='grey65',linewidth=.4)+geom_errorbar(aes(ymin=percentile_low,ymax=percentile_high),width=0,linewidth=.45,alpha=.65)+geom_line(aes(y=original_percentile),colour='grey45',linetype=2,linewidth=.45)+geom_point(size=1.8)+scale_colour_manual(values=riskcols)+scale_y_continuous(limits=c(0,1),breaks=c(0,.25,.5,.75,1))+scale_x_continuous(breaks=NULL)+labs(x='Patients ordered by original MBR',y='Bootstrap MBR percentile',colour=NULL)+theme_fig()+theme(legend.position='bottom')
hh<-hr;hh$unit<-factor(hh$unit,levels=rev(hh$unit))
phr<-ggplot(hh,aes(HR,unit))+geom_vline(xintercept=1,linetype=2,colour='grey50')+geom_segment(aes(x=CI_low,xend=CI_high,yend=unit),colour=TEAL,linewidth=.7)+geom_point(colour=TEAL,size=2.6)+scale_x_log10()+labs(x='MBR hazard ratio (95% CI)',y=NULL)+theme_fig()
prplot<-ggplot(prd,aes(factor(fold),C))+geom_hline(yintercept=.5,linetype=2,colour='grey65')+geom_point(size=2.8,colour=INK)+scale_y_continuous(limits=c(0,1),breaks=c(0,.25,.5,.75,1))+labs(x='Outer test fold',y='C-index: penalized-refit sensitivity')+theme_fig()

supp_top<-plot_grid(pag,ppf,phr,pv,ncol=2,labels=LETTERS[1:4],label_size=FONT_SIZE,align='hv')
supp_bottom<-plot_grid(popt,phf,po,nrow=1,labels=LETTERS[5:7],label_size=FONT_SIZE,align='h',rel_widths=c(1,1,1.2))
supp_bottom<-ggdraw()+draw_plot(supp_bottom,x=.1525,y=0,width=.8475,height=1)
supp1<-plot_grid(supp_top,supp_bottom,ncol=1,rel_heights=c(2,1))
save_plot(supp1,'Supplementary_Fig5','Supplementary_Fig5_validation',13.5,10.5)

wb<-openxlsx::createWorkbook();csvfiles<-list.files(TBL,pattern='\\.csv$',full.names=TRUE)
for(f in csvfiles){nm<-substr(sub('\\.csv$','',basename(f)),1,31);openxlsx::addWorksheet(wb,nm);d<-read.csv(f,check.names=FALSE);openxlsx::writeData(wb,nm,d);openxlsx::freezePane(wb,nm,firstRow=TRUE);openxlsx::setColWidths(wb,nm,cols=seq_len(ncol(d)),widths=pmin(45,pmax(12,nchar(names(d))+2)))}
openxlsx::saveWorkbook(wb,file.path(ANA,'Fig5_All_Audit_Tables.xlsx'),overwrite=TRUE)

validation_pub<-rbind(data.frame(analysis=validation$analysis,C=validation$C,n=nrow(raw),events=sum(event),scope=validation$scope),
 data.frame(analysis='Matching grid-frequency-Cox outer CV',C=outer_summary$pair_weighted_within_fold_C[1],n=nrow(raw),events=sum(event),scope=sprintf('%d/5 training Cox fits failed stability checks; primary C omitted if any fold failed; upstream correction fixed',outer_summary$failed_folds[1])),
 data.frame(analysis='Training-only discovery holdout',C=hc,n=length(te),events=sum(event[te]),scope='Different training-derived signature; upstream correction fixed'))
freqpub<-data.frame(species=pretty_species(taxa4),full_cohort_frequency=discovery$frequencies$selection_fraction[match(taxa4,discovery$frequencies$taxon)],holdout_training_frequency=hold$frequencies$selection_fraction[match(taxa4,hold$frequencies$taxon)])
pubtables<-list(SuppTable1_four_taxa=taxon_coef,SuppTable2_validation=validation_pub,SuppTable3_MBR_HR_scales=hr,SuppTable4_HER2_AUC=auc,SuppTable5_selection_stability=freqpub,SuppTable6_outer_signatures=od,SuppTable7_conditional_MBR=conditional_summary,SuppTable8_penalized_refit=prd)
pw<-openxlsx::createWorkbook()
for(n in names(pubtables)){
 d<-pubtables[[n]];write.csv(d,file.path(OUT_DIR,'Tables',paste0(n,'.csv')),row.names=FALSE,na='')
 sn<-substr(n,1,31);openxlsx::addWorksheet(pw,sn);openxlsx::writeData(pw,sn,d);openxlsx::freezePane(pw,sn,firstRow=TRUE);openxlsx::setColWidths(pw,sn,cols=seq_len(ncol(d)),widths=18)
}
openxlsx::saveWorkbook(pw,file.path(OUT_DIR,'Tables','Fig5_Supplementary_Tables.xlsx'),overwrite=TRUE)
saveRDS(list(config=config,discovery=discovery,fit=fit$fit,zfit=zfit$fit,KM=km,ROC=roc,holdout=holdsum,outer_summary=outer_summary,SHAP=contrib),file.path(ANA,'Fig5_analysis_objects.rds'))
writeLines(capture.output(sessionInfo()),file.path(ANA,'sessionInfo.txt'))
write.csv(data.frame(input=input_path,md5=input_md5,n=nrow(raw),events=sum(event),config_hash=hash_object(config)),file.path(ANA,'input_manifest.csv'),row.names=FALSE)
qa<-data.frame(check=c('Unique patients and valid events','Grid stable taxa match original four','Unit scaling leaves LP unchanged','Additive SHAP sums to centered MBR','Every outer prediction is held out','Outer fold counts complete','Fixed-four corrected C is finite'),pass=c(nrow(raw)>1&&sum(event)>0,setequal(discovery$stable,unname(taxa4)),max(abs(fit$lp-zfit$lp))<1e-8,max(abs(rowSums(contrib)-(score-mean(score))))<1e-8,all(vapply(seq_len(nrow(op)),function(i){mf<-outerman[[op$repetition[i]]];mf$fold[match(op$sample_id[i],mf$sample_id)]==op$fold[i]},logical(1))),nrow(op)==nrow(raw)*N_OUTER_REPEATS,is.finite(validation$C[2])))
write.csv(qa,file.path(OUT_DIR,'QA','numerical_checks.csv'),row.names=FALSE);stopifnot(all(qa$pass))
isolation<-list(data.frame(check='Holdout patients excluded from discovery CV',pass=!any(rownames(meta)[te]%in%hold$fold_manifest$sample_id)))
for(r in seq_len(N_OUTER_REPEATS))for(k in 1:5){
 mf<-outerman[[r]];ds<-readRDS(file.path(CACHE_DIR,sprintf('outer_%02d_fold_%02d_grid.rds',r,k)))
 isolation[[length(isolation)+1]]<-data.frame(check=sprintf('Outer %d fold %d: test patients excluded from discovery CV',r,k),pass=!any(mf$sample_id[mf$fold==k]%in%ds$fold_manifest$sample_id)&&setequal(mf$sample_id[mf$fold!=k],unique(ds$fold_manifest$sample_id)))
}
isolation<-do.call(rbind,isolation);write.csv(isolation,file.path(OUT_DIR,'QA','partition_isolation_checks.csv'),row.names=FALSE);stopifnot(all(isolation$pass))

suppressPackageStartupMessages({library(ggplot2);library(cowplot);library(ggbeeswarm)})
module <- file.path(project,'Revision2/Fig5')
out <- file.path(project,'Revision2/Fig5/output')
control_out <- file.path(project,'Revision2/Fig1_controls/output')
for(d in c('Main_Fig5','QA'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
tables<-file.path(project,'Revision2/Fig5/output/Analysis/Audit_tables')
shap<-read.csv(file.path(tables,'S33_exact_linear_SHAP.csv'),check.names=FALSE)
z<-read.csv(file.path(tables,'S34_four_taxa_CLR_z_matrix.csv'),check.names=FALSE)
scores<-read.csv(file.path(tables,'S09_sample_MBR_scores.csv'))
cf<-read.csv(file.path(tables,'S06_four_taxa_unit_scaled_Cox.csv'))
controls<-read.csv(file.path(control_out,'Tables/Four_taxa_group_summary.csv'))
sp<-cf$term[order(cf$coefficient,decreasing=TRUE)];labels<-gsub('_',' ',sp)
stopifnot(identical(shap$sample_id,z$sample_id),identical(shap$sample_id,scores$sample_id),all(sp%in%names(shap)),all(sp%in%names(z)))
N<-nrow(shap);stopifnot(N==nrow(scores))
plotdata<-data.frame(species=rep(sp,each=N),y=rep(4:1,each=N),contribution=as.vector(as.matrix(shap[,sp])),CLR_z=as.vector(as.matrix(z[,sp])))
grp<-factor(scores$risk_group,levels=c('Low risk','High risk'));ord<-order(grp,-scores$MBR)
heat<-data.frame(patient=rep(seq_len(N),times=4),y=rep(4:1,each=N),value=as.vector(as.matrix(z[ord,sp])))
strip<-data.frame(patient=seq_len(N),group=grp[ord]);riskcols<-c('Low risk'='#218C7E','High risk'='#D4A137')
fs<-controls;fs$y<-5-match(fs$species,labels);stopifnot(!anyNA(fs$y));fs$group<-factor(fs$group,levels=c('MBC','Healthy','Negative'))
FONT<-11;TM<-FONT/ggplot2::.pt;INK<-'#000000'
old_null_device<-cowplot::set_null_device(function(width,height)grDevices::quartz(type='pdf',file=tempfile(fileext='.pdf'),width=width,height=height,family='Arial'))
th<-theme_classic(base_size=FONT,base_family='Arial')+theme(text=element_text(size=FONT,colour=INK),axis.text=element_text(size=FONT,colour=INK),axis.title=element_text(size=FONT,colour=INK),legend.text=element_text(size=FONT,colour=INK),legend.title=element_text(size=FONT,colour=INK),plot.margin=margin(8,6,5,4),plot.title=element_blank())
yscale<-function(show=TRUE)scale_y_continuous(limits=c(.5,5),breaks=1:4,labels=if(show)rev(labels) else NULL,expand=c(0,0))

E1<-ggplot(plotdata,aes(contribution,y,fill=CLR_z,group=species))+geom_vline(xintercept=0,colour='grey70',linewidth=.4)+ggbeeswarm::geom_quasirandom(groupOnX=FALSE,width=.25,size=1.8,shape=21,colour='#554F5A',stroke=.2)+scale_fill_gradient2(low='#3B6FB6',mid='white',high='#D95F59',midpoint=0,limits=c(-2,2),oob=scales::squish,guide='none')+yscale(TRUE)+labs(x='Contribution to MBR (SHAP)',y=NULL)+th+theme(axis.text.y=element_text(face='italic'))
E2<-ggplot(heat,aes(patient,y,fill=value))+geom_tile(width=1,height=.97,colour='white',linewidth=.18)+scale_fill_gradient2(low='#3B6FB6',mid='white',high='#D95F59',midpoint=0,limits=c(-2,2),oob=scales::squish,name='CLR z-score',guide=guide_colourbar(barheight=grid::unit(.9,'inches'),barwidth=grid::unit(.12,'inches')))+scale_x_continuous(limits=c(.5,N+.5),expand=c(0,0),breaks=NULL)+yscale(FALSE)+labs(x='Patients ordered by\nrisk group and MBR',y=NULL)+th+theme(axis.line=element_blank(),axis.ticks=element_blank(),legend.position='right',legend.key.height=grid::unit(.25,'cm'),legend.key.width=grid::unit(.15,'cm'))
E2<-E2+annotate('rect',xmin=strip$patient-.5,xmax=strip$patient+.5,ymin=4.42,ymax=4.64,fill=unname(riskcols[as.character(strip$group)]),colour=NA)+geom_text(data=aggregate(patient~group,strip,mean),aes(patient,4.84,label=group),inherit.aes=FALSE,size=TM,colour=INK,family='Arial')
control_plot<-function(show_species=FALSE){ggplot(fs,aes(group,y))+geom_point(data=fs[fs$detected>0,],aes(size=prevalence,fill=mean_abundance_pct),shape=21,colour='#815752',stroke=.35)+scale_size_area(limits=c(0,1),breaks=c(0,.5,1),labels=c('0%','50%','100%'),max_size=8,name='Prevalence')+scale_fill_gradient(low='#FFF0E5',high='#B84B4B',trans='sqrt',limits=c(0,max(fs$mean_abundance_pct)),breaks=c(0,.10,.20),labels=c('0.00','0.10','0.20'),name='Mean relative\nabundance (%)')+scale_x_discrete(labels=c(MBC='MBC',Healthy='Healthy',Negative='Negative\ncontrols'),expand=expansion(add=.55))+yscale(show_species)+labs(x=NULL,y=NULL)+th+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.y=element_text(face='italic'),legend.position='right',legend.key.height=grid::unit(.24,'cm'),legend.key.width=grid::unit(.17,'cm'),legend.spacing.y=grid::unit(3,'pt'),legend.margin=margin(0,0,0,0))+guides(size=guide_legend(order=1,override.aes=list(stroke=c(0,.35,.35))),fill=guide_colourbar(order=2,barheight=grid::unit(.85,'inches'),barwidth=grid::unit(.12,'inches')))}
F<-control_plot(FALSE)

pdfdev<-function(filename,width,height,bg='white',...)grDevices::quartz(type='pdf',file=filename,width=width,height=height,bg=bg,family='Arial')
ggsave(file.path(out,'Main_Fig5/Fig5E_contributions.pdf'),E1,width=7.5,height=3.6,device=pdfdev,bg='white')
ggsave(file.path(out,'Main_Fig5/Fig5E_abundance_heatmap.pdf'),E2,width=6,height=3.6,device=pdfdev,bg='white')
ggsave(file.path(out,'Main_Fig5/Fig5F.pdf'),control_plot(TRUE),width=7.2,height=3.6,device=pdfdev,bg='white')
cowplot::set_null_device(old_null_device)
