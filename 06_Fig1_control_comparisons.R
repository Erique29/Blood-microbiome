#!/usr/bin/env Rscript
# Compare case, healthy-control and negative-control profiles.
# Fig. 1; control summaries for Fig. 5F. Inputs: scripts 01-03.

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

options(warn=1)
suppressPackageStartupMessages({library(phyloseq);library(vegan);library(dplyr);library(tidyr);library(ggplot2);library(cowplot);library(openxlsx)})
module <- file.path(project,'Revision2/Fig1_controls')
out <- file.path(module,'output')
layout_choice<-Sys.getenv('FIG1_CONTROLS_LAYOUT','three')
for(d in c('Main','Tables','Audit'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)

if(identical(tolower(Sys.getenv('CONTROLS_REDRAW_ONLY','false')),'true')){
 cache_path<-Sys.getenv('CONTROLS_CACHE',file.path(out,'Audit/control_analysis_objects.rds'))
 cached<-readRDS(cache_path)
 paths<-cached$paths;coords<-cached$coords;vr<-cached$variance
 global<-cached$global;pair_stats<-cached$pairwise;core_summary<-cached$core_summary
}else{
paths<-c(MBC=file.path(project,'Output/physeq.rds'),Healthy=file.path(project,'Output/healthy_controls/phyloseq_health_decontam_clean.rds'),Negative=file.path(project,'Output/physeq_raw.RDS'))
flags_path<-file.path(project,'revision_analysis/02_contamination_control/decontam_flagged_taxa_details.csv')
write_table<-function(x,n)write.csv(x,file.path(out,'Tables',paste0(n,'.csv')),row.names=FALSE,na='')
mat<-function(ps){a<-as(otu_table(ps),'matrix');if(taxa_are_rows(ps))a<-t(a);stopifnot(all(is.finite(a)),all(a>=0),all(rowSums(a)>0));a}
ps<-lapply(paths,readRDS);a<-lapply(ps,mat);raw_all<-a$Negative;raw_meta<-data.frame(sample_data(ps$Negative),check.names=FALSE)
negative_ids<-rownames(raw_meta)[raw_meta$Sample_or_Control=='Control Sample'];a$Negative<-a$Negative[negative_ids,,drop=FALSE]
stopifnot(all(vapply(a,nrow,integer(1)) > 0))
a<-lapply(a,function(x)x/rowSums(x));taxa<-Reduce(union,lapply(a,colnames))
pad<-function(x){z<-matrix(0,nrow(x),length(taxa),dimnames=list(rownames(x),taxa));z[,colnames(x)]<-x;z}
u<-lapply(a,pad);md<-bind_rows(lapply(names(u),function(g)data.frame(source_id=rownames(u[[g]]),sample_id=paste(g,rownames(u[[g]]),sep='::'),group=g,source_path=unname(paths[g]))))
x<-do.call(rbind,u);rownames(x)<-md$sample_id;md$group<-factor(md$group,levels=c('MBC','Healthy','Negative'));rownames(md)<-md$sample_id
stopifnot(identical(rownames(x),rownames(md)),!anyDuplicated(rownames(x)))
write.csv(md,file.path(out,'Audit/sample_manifest.csv'),row.names=FALSE)
write.csv(data.frame(path=c(paths,flags_path),md5=unname(tools::md5sum(c(paths,flags_path)))),file.path(out,'Audit/input_manifest.csv'),row.names=FALSE)
groups<-levels(md$group);pairs<-combn(groups,2,simplify=FALSE)

D<-vegdist(x,'bray');set.seed(20260625);ad<-adonis2(D~group,data=md,permutations=999)
bd<-betadisper(D,md$group,bias.adjust=FALSE);dp<-permutest(bd,permutations=999)
global<-data.frame(comparison='MBC / Healthy / Negative',N=nrow(md),PERMANOVA_R2=ad$R2[1],PERMANOVA_F=ad$F[1],PERMANOVA_P=ad$`Pr(>F)`[1],PERMDISP_F=dp$tab$F[1],PERMDISP_P=dp$tab$`Pr(>F)`[1],permutations=999,seed=20260625)
pair_stats<-bind_rows(lapply(seq_along(pairs),function(i){p<-pairs[[i]];m<-droplevels(md[md$group%in%p,,drop=FALSE]);d<-vegdist(x[rownames(m),,drop=FALSE],'bray');set.seed(666+i);aa<-adonis2(d~group,data=m,permutations=999);bb<-betadisper(d,m$group,bias.adjust=FALSE);dd<-permutest(bb,permutations=999);data.frame(reference=p[1],comparison=p[2],N=nrow(m),N_reference=sum(m$group==p[1]),N_comparison=sum(m$group==p[2]),PERMANOVA_R2=aa$R2[1],PERMANOVA_F=aa$F[1],PERMANOVA_P=aa$`Pr(>F)`[1],PERMDISP_F=dd$tab$F[1],PERMDISP_P=dd$tab$`Pr(>F)`[1],permutations=999,seed=666+i)}))
pair_stats$PERMANOVA_q_BH3<-p.adjust(pair_stats$PERMANOVA_P,'BH');pair_stats$PERMDISP_q_BH3<-p.adjust(pair_stats$PERMDISP_P,'BH')
write_table(global,'Community_global');write_table(pair_stats,'Community_pairwise')
ord<-cmdscale(D,k=2,eig=TRUE);vr<-ord$eig/sum(ord$eig[ord$eig>0])*100
coords<-base::transform(md,PCoA1=ord$points[,1],PCoA2=ord$points[,2],distance_to_centroid=bd$distances)
write_table(coords,'PCoA_coordinates_and_dispersion')

preval<-lapply(u,function(z)colMeans(z>0));meanpct<-lapply(u,function(z)100*colMeans(z))
core<-bind_rows(lapply(c(.2,.5),function(cut){f<-taxa[preval$MBC>=cut];data.frame(threshold=cut,feature=f,MBC_detected_n=colSums(u$MBC[,f,drop=FALSE]>0),Healthy_detected_n=colSums(u$Healthy[,f,drop=FALSE]>0),Negative_detected_n=colSums(u$Negative[,f,drop=FALSE]>0))}))
core$category<-ifelse(core$Negative_detected_n>0,'Negative-control detected',ifelse(core$Healthy_detected_n>0,'MBC + healthy; no negative control','MBC only'))
core_summary<-core%>%count(threshold,category,name='taxa_n')%>%group_by(threshold)%>%mutate(total=sum(taxa_n),percent=100*taxa_n/total)%>%ungroup()
stopifnot(all(core$threshold %in% c(.2,.5)))
write_table(core,'Core_overlap_details');write_table(core_summary,'Core_overlap_summary')

sp<-c('Sphingopyxis_sp001990265','Bacillus_pumilus','Ralstonia_insidiosa','Paraburkholderia_tropica')
four<-vapply(sp,function(s){f<-grep(paste0(';',s,'$'),taxa,value=TRUE);stopifnot(length(f)==1);f},character(1))
fs<-bind_rows(lapply(seq_along(four),function(j)bind_rows(lapply(groups,function(g){v<-u[[g]][,four[j]];ci<-binom.test(sum(v>0),length(v))$conf.int;data.frame(feature=four[j],species=gsub('_',' ',sp[j]),group=g,N=length(v),detected=sum(v>0),prevalence=mean(v>0),prevalence_CI_low=ci[1],prevalence_CI_high=ci[2],mean_abundance_pct=100*mean(v),median_abundance_pct=100*median(v))}))))
fglobal<-bind_rows(lapply(seq_along(four),function(j){v<-x[,four[j]];tt<-rbind(detected=sapply(u,function(z)sum(z[,four[j]]>0)),not_detected=sapply(u,function(z)sum(z[,four[j]]==0)));kk<-kruskal.test(v,md$group);data.frame(species=gsub('_',' ',sp[j]),prevalence_Fisher_P=fisher.test(tt)$p.value,abundance_KW_H=unname(kk$statistic),abundance_KW_P=kk$p.value)}))
fglobal$prevalence_q_BH4<-p.adjust(fglobal$prevalence_Fisher_P,'BH');fglobal$abundance_q_BH4<-p.adjust(fglobal$abundance_KW_P,'BH')
fpair<-bind_rows(lapply(seq_along(four),function(j)bind_rows(lapply(pairs,function(p){v<-u[[p[1]]][,four[j]];w<-u[[p[2]]][,four[j]];tt<-rbind(c(sum(v>0),sum(v==0)),c(sum(w>0),sum(w==0)));wp<-if(length(unique(c(v,w)))==1)1 else wilcox.test(w,v,exact=FALSE,correct=TRUE)$p.value;data.frame(species=gsub('_',' ',sp[j]),reference=p[1],comparison=p[2],N_reference=length(v),N_comparison=length(w),prevalence_difference=mean(w>0)-mean(v>0),mean_abundance_difference_pp=100*(mean(w)-mean(v)),prevalence_Fisher_P=fisher.test(tt)$p.value,abundance_Wilcoxon_P=wp)}))))
fpair$prevalence_q_BH12<-p.adjust(fpair$prevalence_Fisher_P,'BH');fpair$abundance_q_BH12<-p.adjust(fpair$abundance_Wilcoxon_P,'BH')
write_table(fs,'Four_taxa_group_summary');write_table(fglobal,'Four_taxa_global_tests');write_table(fpair,'Four_taxa_pairwise_tests')

flags<-read.csv(flags_path,check.names=FALSE)$taxon;case_before<-raw_all[rownames(a$MBC),,drop=FALSE]
before_tax<-colnames(case_before)[colSums(case_before>0)>0];neg_tax<-colnames(a$Negative)[colSums(a$Negative>0)>0]
counts<-data.frame(metric=c('MBC samples','Healthy samples','Negative samples','MBC taxa before decontamination','All flagged taxa','Flagged taxa also detected in MBC','Flagged taxa only detected in negative controls','Final MBC taxa'),value=c(nrow(u$MBC),nrow(u$Healthy),nrow(u$Negative),length(before_tax),length(flags),sum(flags%in%before_tax),sum(flags%in%neg_tax & !flags%in%before_tax),sum(colSums(u$MBC>0)>0)))
stopifnot(all(is.finite(counts$value)),all(counts$value>=0))
write_table(counts,'Processing_counts')

openxlsx::write.xlsx(list(Community_global=global,Community_pairwise=pair_stats,Core_summary=as.data.frame(core_summary),Core_details=core,Four_taxa_summary=fs,Four_taxa_global=fglobal,Four_taxa_pairwise=fpair,Processing_counts=counts,Sample_manifest=md),file.path(out,'Tables/Control_comparison_results.xlsx'),overwrite=TRUE)
saveRDS(list(paths=paths,x=x,metadata=md,coords=coords,variance=vr[1:2],global=global,pairwise=pair_stats,core=core,core_summary=core_summary,four_summary=fs,four_global=fglobal,four_pairwise=fpair,counts=counts),file.path(out,'Audit/control_analysis_objects.rds'))
writeLines(capture.output(sessionInfo()),file.path(out,'Audit/sessionInfo.txt'))
}

old_null_device<-cowplot::set_null_device(function(width,height)grDevices::quartz(type='pdf',file=tempfile(fileext='.pdf'),width=width,height=height,family='Arial'))
FONT<-11;INK<-'#000000';TM<-FONT/ggplot2::.pt
cols<-c(MBC='#D95F59',Healthy='#3B6FB6',Negative='#BE9638')
glab<-c(MBC='MBC',Healthy='Healthy',Negative='Negative controls')
theme_fig<-function()theme_classic(base_size=FONT,base_family='Arial')+theme(text=element_text(size=FONT,colour=INK),axis.text=element_text(size=FONT,colour=INK),axis.title=element_text(size=FONT,colour=INK),legend.text=element_text(size=FONT,colour=INK),legend.title=element_text(size=FONT,colour=INK),plot.margin=margin(6,8,6,6))
pcoa<-ggplot(coords,aes(PCoA1,PCoA2,colour=group,shape=group))+stat_ellipse(type='norm',level=.95,linetype=2,linewidth=.55,show.legend=FALSE)+geom_point(size=2.5,alpha=.9)+scale_colour_manual(values=cols,labels=glab,name=NULL)+scale_shape_manual(values=c(MBC=16,Healthy=17,Negative=15),labels=glab,name=NULL)+labs(x=sprintf('PCoA1 (%.1f%%)',vr[1]),y=sprintf('PCoA2 (%.1f%%)',vr[2]))+coord_equal()+theme_fig()+theme(legend.position='inside',legend.position.inside=c(.02,.98),legend.justification=c(0,1),legend.background=element_rect(fill=scales::alpha('white',.8),colour=NA),legend.key.height=grid::unit(11,'pt'))+annotate('label',x=Inf,y=-Inf,label=sprintf('PERMANOVA~~R^2=="%.3f"',global$PERMANOVA_R2),parse=TRUE,hjust=1.01,vjust=-1.35,size=TM,family='Arial',fill='white',linewidth=0,label.padding=grid::unit(1,'pt'))+annotate('label',x=Inf,y=-Inf,label=sprintf('italic(P)=="%.3f"',global$PERMANOVA_P),parse=TRUE,hjust=1.01,vjust=-.08,size=TM,family='Arial',fill='white',linewidth=0,label.padding=grid::unit(1,'pt'))

pc_build<-ggplot_build(pcoa)
pc_range<-pc_build$layout$panel_params[[1]]
pc_span<-max(diff(pc_range$x.range),diff(pc_range$y.range))
pc_x<-mean(pc_range$x.range)+c(-.5,.5)*pc_span
pc_y<-mean(pc_range$y.range)+c(-.5,.5)*pc_span
pcoa<-pcoa+coord_equal(xlim=pc_x,ylim=pc_y,expand=FALSE)
cs<-core_summary;cs$threshold<-factor(cs$threshold,levels=c(.2,.5),labels=c('≥20%','≥50%'));cs$category<-factor(cs$category,levels=c('MBC only','MBC + healthy; no negative control','Negative-control detected'))
cc<-c('MBC only'='#D895A2','MBC + healthy; no negative control'='#84ADD0','Negative-control detected'='#D7B45D')
cs$x<-ifelse(cs$threshold=='≥20%',1,3)
pcore<-ggplot(cs,aes(x,taxa_n,fill=category))+geom_col(width=.80,key_glyph=function(data,params,size)grid::rectGrob(width=grid::unit(9,'pt'),height=grid::unit(9,'pt'),gp=grid::gpar(fill=data$fill,col=NA)))+geom_text(aes(label=sprintf('%d (%.1f%%)',taxa_n,percent)),position=position_stack(vjust=.5),size=TM,colour=INK)+scale_fill_manual(values=cc,name=NULL,labels=c('MBC only'='Case core only;\nnegative-control absent','MBC + healthy; no negative control'='Shared with healthy;\nnegative-control absent','Negative-control detected'='Negative-control detected'))+scale_x_continuous(limits=c(.35,3.8),breaks=c(1,3),labels=c('≥20%','≥50%'),expand=c(0,0))+scale_y_continuous(expand=expansion(mult=c(0,.06)))+labs(x='Prevalence threshold in MBC',y='Number of MBC core species')+theme_fig()+theme(legend.position='inside',legend.position.inside=c(.99,.99),legend.justification=c(1,1),legend.background=element_rect(fill='white',colour=NA),legend.margin=margin(2,2,2,2),legend.key.height=grid::unit(10,'pt'),legend.key.width=grid::unit(10,'pt'),legend.key.spacing.y=grid::unit(7,'pt'))+guides(fill=guide_legend(ncol=1))

pz<-pair_stats;pz$contrast<-paste(pz$comparison,'vs',pz$reference);pz$contrast<-factor(pz$contrast,levels=rev(pz$contrast));pz$stars<-ifelse(pz$PERMANOVA_P<.001,'***',ifelse(pz$PERMANOVA_P<.01,'**',ifelse(pz$PERMANOVA_P<.05,'*','')))
ppairs<-ggplot(pz,aes(PERMANOVA_R2,contrast))+geom_col(fill='#80B6B7',width=.58)+geom_text(aes(label=stars),hjust=-.5,size=TM)+scale_x_continuous(expand=expansion(mult=c(0,.18)))+labs(x=expression('Bray–Curtis PERMANOVA '*R^2),y=NULL)+theme_fig()

save_square_row<-function(plots,letters,filename,panel_inches=3.9){
 metric_pdf<-tempfile(fileext='.pdf')
 grDevices::quartz(type='pdf',file=metric_pdf,width=18,height=6,family='Arial')
 grobs<-lapply(plots,ggplotGrob)
 stopifnot(length(unique(vapply(grobs,function(g)length(g$heights),integer(1))))==1L)
 common_heights<-do.call(grid::unit.pmax,lapply(grobs,`[[`,'heights'))
 for(i in seq_along(grobs)){
  g<-grobs[[i]];g$heights<-common_heights
  k<-which(g$layout$name=='panel');stopifnot(length(k)==1L)
  g$widths[g$layout$l[k]]<-grid::unit(panel_inches,'inches')
  g$heights[g$layout$t[k]]<-grid::unit(panel_inches,'inches')
  g<-gtable::gtable_add_rows(g,grid::unit(16,'pt'),pos=0)
  g<-gtable::gtable_add_grob(g,grid::textGrob(letters[i],x=grid::unit(2,'pt'),y=.7,hjust=0,gp=grid::gpar(fontsize=FONT,fontface='bold',fontfamily='Arial',col=INK)),t=1,l=1,r=ncol(g),clip='off')
  grobs[[i]]<-g
 }
 widths<-vapply(grobs,function(g)grid::convertWidth(sum(g$widths),'inches',valueOnly=TRUE),numeric(1))
 heights<-vapply(grobs,function(g)grid::convertHeight(sum(g$heights),'inches',valueOnly=TRUE),numeric(1))
 grDevices::dev.off();unlink(metric_pdf)
 grDevices::quartz(type='pdf',file=filename,width=sum(widths),height=max(heights),family='Arial',bg='white')
 grid::grid.newpage();grid::pushViewport(grid::viewport(layout=grid::grid.layout(1,length(grobs),widths=grid::unit(widths,'inches'),heights=grid::unit(max(heights),'inches'))))
 for(i in seq_along(grobs)){
  grid::pushViewport(grid::viewport(layout.pos.row=1,layout.pos.col=i))
  grid::grid.draw(grobs[[i]]);grid::popViewport()
 }
 grDevices::dev.off()
 write.csv(data.frame(panel=letters,panel_width_inches=panel_inches,panel_height_inches=panel_inches,outer_width_inches=widths,outer_height_inches=heights),file.path(out,'Audit/Fig1_panel_dimensions.csv'),row.names=FALSE)
}
if(layout_choice=='three')save_square_row(list(pcoa,ppairs,pcore),c('A','B','C'),file.path(out,'Main/Fig1.pdf')) else save_square_row(list(pcoa,pcore),c('A','B'),file.path(out,'Main/Fig1.pdf'))
writeLines(c('Inputs: current batch-corrected MBC, clean healthy, raw negative controls; within-sample closure only. No extra batch correction.','Community PERMANOVA and PERMDISP: three groups overall and all three pairwise comparisons, 999 permutations.','Pairwise community P values are BH-adjusted across three contrasts, separately for PERMANOVA and PERMDISP.','Core prevalence thresholds are >=0.20 and >=0.50; core overlap is descriptive.','Four-taxon prevalence: exact Fisher tests; relative abundance including zeros: Kruskal-Wallis overall and Wilcoxon pairwise.','Four-taxon global BH family: four taxa per outcome; pairwise BH family: four taxa x three contrasts =12 per outcome.','No full-taxon screen, survival refit, or new feature selection.','Four-taxon descriptive panel belongs to Fig5F, not Fig1.'),file.path(out,'Audit/methods_contract.txt'))
cowplot::set_null_device(old_null_device)
print(global);print(pair_stats);message('Finished controls analysis: ',out)
