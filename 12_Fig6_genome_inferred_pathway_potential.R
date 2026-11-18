#!/usr/bin/env Rscript
# Calculate genome-inferred pathway potential weighted by corrected relative abundance.
# Fig. 6; Data S7. Inputs: scripts 01 and 11, plus gapseq prediction tables.

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

suppressPackageStartupMessages({library(data.table);library(dplyr);library(purrr)})
sp_names<-c(Sphingopyxis="Sphingopyxis_sp001990265",Bacillus="Bacillus_pumilus",Ralstonia="Ralstonia_insidiosa",Paraburkholderia="Paraburkholderia_tropica")
path_dir<-file.path(PROJECT_DIR,"gapseq");dir.create(path_dir,recursive=TRUE,showWarnings=FALSE)
lst<-lapply(names(sp_names),function(genus){
  source<-private_input("gapseq","4genome",genus,paste0(genus,"-all-Pathways.tbl"))
  dest<-file.path(path_dir,"4genome",genus,basename(source));dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
  if (!file.copy(source,dest,overwrite=TRUE)) stop("Could not stage the private genome prediction input.")
  data.table::fread(source,data.table=FALSE,check.names=FALSE)
})
names(lst)<-unname(sp_names)
clean_tbl<-function(df){df %>% dplyr::select(Name,ID,Prediction) %>% mutate(Prediction=ifelse(tolower(Prediction)=="true",1,0))}
df_merged<-purrr::reduce(purrr::map(lst,clean_tbl),function(x,y)full_join(x,y,by=c("ID","Name")))
colnames(df_merged)<-c("Pathway","ID",names(lst))
df_merged<-df_merged[,c("Pathway","ID","Sphingopyxis_sp001990265","Bacillus_pumilus","Ralstonia_insidiosa","Paraburkholderia_tropica")]
write.csv(df_merged,file.path(path_dir,"pathway.csv"),row.names=FALSE)

suppressPackageStartupMessages({library(phyloseq);library(VSURF);library(ComplexHeatmap);library(circlize);library(grid);library(ggplot2)})
module <- file.path(project,'Revision2/Fig6')
out <- file.path(module,'output')
for(d in c('Main','Supplementary','Tables','Audit','QA','Cache'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)

grDevices::pdfFonts(Arial=grDevices::pdfFonts('Helvetica')[[1]])
font_size<-11;seed<-666L
sp<-c('Sphingopyxis_sp001990265','Bacillus_pumilus','Ralstonia_insidiosa','Paraburkholderia_tropica')
phyfile<-file.path(project,'Output/physeq.rds');pathfile<-file.path(project,'gapseq/pathway.csv');scorefile<-file.path(project,'Revision2/Fig5/output/Analysis/Audit_tables/S09_sample_MBR_scores.csv')
phy<-readRDS(phyfile);ab<-as(otu_table(phy),'matrix');if(taxa_are_rows(phy))ab<-t(ab)
meta<-data.frame(sample_data(phy),check.names=FALSE)[rownames(ab),,drop=FALSE]
stopifnot(nrow(ab)>1,identical(rownames(ab),rownames(meta)),all(is.finite(ab)),all(ab>=0),!anyDuplicated(meta$Sample_id))
cols<-vapply(sp,function(s){f<-grep(paste0(';',s,'$'),colnames(ab),value=TRUE);stopifnot(length(f)==1);f},character(1))
pw<-read.csv(pathfile,check.names=FALSE);scores<-read.csv(scorefile,check.names=FALSE)
stopifnot(!anyDuplicated(pw$Pathway),!anyDuplicated(pw$ID),setequal(scores$sample_id,meta$Sample_id))
scores<-scores[match(meta$Sample_id,scores$sample_id),,drop=FALSE]
y<-scores$MBR;names(y)<-meta$Sample_id
grp<-factor(scores$risk_group,levels=c('Low risk','High risk'))
stopifnot(all(as.character(grp)==ifelse(y>=median(y),'High risk','Low risk')),!anyDuplicated(meta$Order))
patient<-paste0('P',rank(meta$Order,ties.method='first'))
W<-as.matrix(pw[,sp]);storage.mode(W)<-'double';rownames(W)<-pw$Pathway
stopifnot(all(W%in%c(0,1)),!anyNA(W))

genome_folders<-c('Sphingopyxis','Bacillus','Ralstonia','Paraburkholderia')
genome_checks<-lapply(seq_along(sp),function(i){f<-file.path(project,'gapseq/4genome',genome_folders[i],paste0(genome_folders[i],'-all-Pathways.tbl'));d<-data.table::fread(f,data.table=FALSE,check.names=FALSE);j<-match(pw$ID,d$ID);stopifnot(!anyNA(j),!anyDuplicated(d$ID));v<-as.integer(tolower(as.character(d$Prediction[j]))=='true');bad<-sum(v!=W[,i]);stopifnot(bad==0);data.frame(species=sp[i],records=nrow(d),presence_mismatches=bad,input=f)})
genome_checks<-do.call(rbind,genome_checks)
A4<-ab[,cols,drop=FALSE];colnames(A4)<-sp
potential<-A4%*%t(W);rownames(potential)<-meta$Sample_id
stopifnot(all(potential>=0),identical(colnames(potential),pw$Pathway))
active<-apply(potential,2,sd)>0

profile_code<-apply(W,1,paste0,collapse='');pattern_id<-paste0('G',match(profile_code,sort(unique(profile_code))))
patterns<-data.frame(pathway=pw$Pathway,pathway_ID=pw$ID,pattern_id=pattern_id,profile_code=profile_code,W,check.names=FALSE)
audit<-data.frame(metric=c('Cases','Pathway records','Nonconstant pathways','Distinct presence patterns'),value=c(nrow(ab),nrow(W),sum(active),length(unique(profile_code))))
write.csv(audit,file.path(out,'Audit/input_scale_audit.csv'),row.names=FALSE)
write.csv(genome_checks,file.path(out,'Audit/genome_prediction_checks.csv'),row.names=FALSE)
write.csv(data.frame(sample_id=meta$Sample_id,patient=patient,MBR=y,risk_group=grp,A4,check.names=FALSE),file.path(out,'Audit/sample_MBR_and_relative_abundance.csv'),row.names=FALSE)
write.csv(data.frame(sample_id=meta$Sample_id,potential,check.names=FALSE),file.path(out,'Audit/pathway_potential_matrix.csv'),row.names=FALSE)
write.csv(patterns,file.path(out,'Tables/pathway_species_presence.csv'),row.names=FALSE)

vsurf_method<-getS3method('VSURF_thres','default',envir=asNamespace('VSURF'))
writeLines(capture.output(formals(vsurf_method)),file.path(out,'Audit/VSURF_threshold_formals.txt'))
X<-as.data.frame(potential[,active,drop=FALSE]);mtry<-max(floor(ncol(X)/3),1)
hashobj<-function(z){f<-tempfile();saveRDS(z,f);h<-unname(tools::md5sum(f));unlink(f);h}
cfg<-list(input_hash=hashobj(list(X,y)),seed=seed,ntree.thres=2000L,nfor.thres=20L,mtry=mtry,parallel=FALSE,RFimplem='randomForest',VSURF_version=as.character(packageVersion('VSURF')),randomForest_version=as.character(packageVersion('randomForest')))
cache<-file.path(out,'Cache/VSURF_threshold_relative_potential.rds')
if(file.exists(cache)){v<-readRDS(cache);stopifnot(identical(v$config,cfg));vr<-v$result;fit_warnings<-v$warnings;if(is.null(fit_warnings))fit_warnings<-character();message('Using matching current-module importance cache.')}else{set.seed(seed);message('Fresh VSURF threshold ranking: ',ncol(X),' predictors; ',cfg$nfor.thres,' forests; ',cfg$ntree.thres,' trees per forest; mtry=',mtry);fit_warnings<-character();vr<-withCallingHandlers(VSURF::VSURF_thres(x=X,y=y,ntree.thres=cfg$ntree.thres,nfor.thres=cfg$nfor.thres,mtry=mtry,parallel=FALSE,RFimplem='randomForest',verbose=TRUE),warning=function(w){fit_warnings<<-c(fit_warnings,conditionMessage(w));invokeRestart('muffleWarning')});saveRDS(list(config=cfg,result=vr,warnings=fit_warnings),cache)}
write.csv(data.frame(warning=fit_warnings),file.path(out,'Audit/VSURF_warnings.csv'),row.names=FALSE)
ix<-vr$imp.mean.dec.ind
importance<-data.frame(rank=seq_along(ix),pathway=colnames(X)[ix],importance=vr$imp.mean.dec,importance_sd=vr$imp.sd.dec)
stopifnot(nrow(importance)==sum(active),all(diff(importance$importance)<=0))

getcors<-function(M){lapply(seq_len(ncol(M)),function(j){v<-M[,j];if(sd(v)==0)return(data.frame(pathway=colnames(M)[j],rho=NA_real_,P=NA_real_));r<-suppressWarnings(cor.test(v,y,method='spearman',exact=FALSE));data.frame(pathway=colnames(M)[j],rho=unname(r$estimate),P=r$p.value)})}
cors<-do.call(rbind,getcors(potential));cors$q_BH<-NA_real_;cors$q_BH[active]<-p.adjust(cors$P[active],method='BH');cors$direction<-ifelse(!is.finite(cors$rho),'Constant',ifelse(cors$rho>0,'Positive','Negative'))
cors$mean_low<-colMeans(potential[grp=='Low risk',,drop=FALSE]);cors$mean_high<-colMeans(potential[grp=='High risk',,drop=FALSE]);cors$pattern_id<-pattern_id
allresults<-merge(importance,cors,by='pathway',all.y=TRUE,sort=FALSE);allresults<-allresults[order(allresults$rank,na.last=TRUE),]
top<-head(importance$pathway,30);topresults<-allresults[match(top,allresults$pathway),]
write.csv(allresults,file.path(out,'Tables/all_pathway_results.csv'),row.names=FALSE)
write.csv(topresults,file.path(out,'Tables/Fig6_top30_results.csv'),row.names=FALSE)
write.csv(importance,file.path(out,'Tables/VSURF_complete_importance.csv'),row.names=FALSE)

riskcols<-c('Low risk'='#159A8C','High risk'='#C73652');effectcols<-c('negative'='#3B6FB6','positive'='#D95F59');heatcols<-circlize::colorRamp2(c(-2,-.7,0,2),c('#735269','#BCA0B3','#E9D7B0','#C77B39'))
z<-t(scale(log2(potential[,top,drop=FALSE]+1e-6)));rownames(z)<-top;colnames(z)<-patient
rho<-topresults$rho
grDevices::quartz(type='pdf',file=file.path(out,'Main/Fig6.pdf'),width=13.4,height=9.4,family='Arial',bg='white')
fg<-gpar(fontsize=font_size,fontfamily='Arial',fontface='plain',col='black')
ha<-HeatmapAnnotation(MBR=anno_points(y,pch=16,size=unit(1.2,'mm'),gp=gpar(col='#34495E'),border=FALSE,height=unit(12,'mm'),axis_param=list(at=range(y),labels=sprintf('%.1f',range(y)),gp=fg)),`Risk Group`=grp,col=list(`Risk Group`=riskcols),simple_anno_size=unit(3.5,'mm'),annotation_name_gp=fg,annotation_name_side='right',show_annotation_name=c(MBR=TRUE,`Risk Group`=FALSE),gap=unit(1.5,'mm'),annotation_legend_param=list(`Risk Group`=list(title='Risk Group',title_gp=fg,labels_gp=fg,nrow=1)))
ra<-rowAnnotation(spacer=anno_empty(which='row',border=FALSE,width=unit(3,'mm')),`Spearman rho`=anno_points(rho,which='row',pch=16,size=unit(1.7,'mm'),gp=gpar(col=ifelse(rho>=0,effectcols['positive'],effectcols['negative'])),ylim=c(-1,1),extend=0,axis=TRUE,axis_param=list(at=c(-1,0,1),labels=c('-1','0','1'),side='bottom',labels_rot=0,gp=fg),border=FALSE,width=unit(25,'mm')),show_annotation_name=c(spacer=FALSE,`Spearman rho`=TRUE),gap=unit(0,'mm'),annotation_name_gp=fg,annotation_name_side='top',annotation_name_rot=0)
ht<-Heatmap(z,name='Pathway potential\n(log2, row Z)',col=heatcols,top_annotation=ha,right_annotation=ra,column_split=grp,cluster_column_slices=FALSE,cluster_columns=TRUE,clustering_distance_columns='euclidean',clustering_method_columns='average',column_dend_reorder=FALSE,column_title=NULL,column_gap=unit(3,'mm'),cluster_rows=TRUE,clustering_distance_rows='euclidean',clustering_method_rows='average',row_dend_reorder=FALSE,row_dend_side='left',row_dend_width=unit(7,'mm'),column_dend_height=unit(7,'mm'),row_names_side='left',row_names_gp=fg,row_names_max_width=max_text_width(top,gp=fg)+unit(2,'mm'),column_names_gp=fg,column_names_rot=90,rect_gp=gpar(col='white',lwd=.25),border=FALSE,width=unit(152,'mm'),height=unit(177,'mm'),heatmap_legend_param=list(title_gp=fg,labels_gp=fg,at=c(-2,0,2),direction='horizontal',legend_width=unit(45,'mm')))
drawn<-draw(ht,show_heatmap_legend=FALSE,show_annotation_legend=FALSE,padding=unit(c(10,6,10,6),'mm'))
decorate_annotation('Spearman rho', {grid.lines(x=unit(c(.5,.5),'npc'),y=unit(c(0,1),'npc'),gp=gpar(col='grey70',lty=2,lwd=.5))})

risk_lgd<-Legend(title='Risk Group',labels=names(riskcols),legend_gp=gpar(fill=riskcols,col=NA),nrow=1,title_position='topleft',title_gp=fg,labels_gp=fg,grid_height=unit(3,'mm'),grid_width=unit(3,'mm'))
heat_lgd<-Legend(title='Pathway potential\n(log2, row Z)',col_fun=heatcols,at=c(-2,0,2),direction='horizontal',legend_width=unit(40,'mm'),title_position='topleft',title_gp=fg,labels_gp=fg)
lgd<-packLegend(risk_lgd,heat_lgd,direction='horizontal',gap=unit(6,'mm'))
upViewport(0)
draw(lgd,x=unit(12,'mm'),y=unit(1,'npc')-unit(13,'mm'),just=c('left','top'))
ro<-row_order(drawn);co<-unlist(column_order(drawn),use.names=FALSE)
write.csv(data.frame(display_row=seq_along(ro),pathway=rownames(z)[ro],rho=rho[ro]),file.path(out,'Audit/main_row_order.csv'),row.names=FALSE)
write.csv(data.frame(display_column=seq_along(co),patient=patient[co],sample_id=meta$Sample_id[co],MBR=y[co],risk_group=grp[co]),file.path(out,'Audit/main_column_order.csv'),row.names=FALSE)
dev.off()

saveRDS(list(potential=potential,top30=top,display_z=z,MBR=y,risk_group=grp,patient=patient,main_row_order=ro,main_column_order=co,VSURF_config=cfg),file.path(out,'Audit/Fig6_plot_inputs.rds'))
sources<-c(phyfile,pathfile,scorefile,genome_checks$input)
write.csv(data.frame(path=sources,md5=unname(tools::md5sum(sources))),file.path(out,'Audit/input_manifest.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'Audit/sessionInfo.txt'))
writeLines(c('Formula: potential[i,p] = sum_s(corrected_relative_abundance[i,s] * predicted_presence[p,s]); four selected taxa only.','No closure over pathways; no renormalization to the four-taxon total.','All direction statistics use unscaled potential and continuous unchanged MBR.','For display only, apply log2(potential + 1e-6), then row-standardize; colors saturate outside -2 and +2. All directions and RF use unlogged, unscaled potential.','Risk Group annotation and legend titles use plain font, not bold.','Rows use Euclidean/average clustering of row-standardized profiles; columns cluster within each risk group.','RF ranking: VSURF threshold mean permutation importance; 20 forests of 2000 trees, seed666, serial, mtry=floor(p/3).','P/q values in tables are descriptive because pathways and MBR share the same four taxa and many pathway profiles are identical.','Do not claim measured pathway abundance/activity, independent prognosis or robust classification.'),file.path(out,'Audit/analysis_and_plot_settings.txt'))
