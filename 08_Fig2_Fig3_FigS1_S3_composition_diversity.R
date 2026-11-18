#!/usr/bin/env Rscript
# Analyze species composition and diversity across clinical groups.
# Figs. 2-3 and S1-S3. Inputs: scripts 01, 04 and 07.

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

OUT_DIR <- file.path(PROJECT_DIR,'Revision2/fig1-3/output')
suppressPackageStartupMessages({library(phyloseq);library(vegan);library(ggplot2);library(cowplot);library(openxlsx)})
for(d in c('Main','Supplementary','Tables','Audit','QA'))dir.create(file.path(OUT_DIR,d),recursive=TRUE,showWarnings=FALSE)
FONT<-11;SEED<-666L;PC<-1e-6;BLUE<-'#4D76A6';PINK<-'#C77896';TEAL<-'#518B86';INK<-'#000000';set.seed(SEED)
wt<-function(z,n,dir='Tables')write.csv(z,file.path(OUT_DIR,dir,paste0(n,'.csv')),row.names=FALSE,na='')
pretty_species<-function(x)gsub('_',' ',sub('.*;','',x))
theme_fig<-function()theme_classic(base_size=FONT,base_family='sans')+theme(text=element_text(size=FONT,colour=INK),axis.text=element_text(size=FONT,colour=INK),axis.title=element_text(size=FONT,colour=INK),legend.title=element_text(size=FONT,colour=INK),legend.text=element_text(size=FONT,colour=INK),strip.text=element_text(size=FONT,colour=INK),plot.title=element_blank(),plot.margin=margin(6,8,6,6),axis.line=element_line(colour=INK,linewidth=.4),axis.ticks=element_line(colour=INK,linewidth=.4))
save_pdf<-function(p,folder,name,w,h){tmp<-tempfile(fileext='.pdf');dev<-if(capabilities('aqua'))function(filename,width,height,bg='white',...)grDevices::quartz(type='pdf',file=filename,width=width,height=height,bg=bg,family='Arial') else grDevices::pdf;ggsave(tmp,p,width=w,height=h,device=dev,bg='white');dest<-file.path(OUT_DIR,folder,paste0(name,'.pdf'));st<-paste0(dest,'.writing');stopifnot(file.copy(tmp,st,overwrite=TRUE),file.rename(st,dest));unlink(tmp)}

stat_annotation<-function(first,p,hjust=-.03){list(
 annotate('text',x=-Inf,y=Inf,label=first,hjust=hjust,vjust=1.15,size=FONT/(72.27/25.4),colour=INK),
 annotate('text',x=-Inf,y=Inf,label=sprintf('italic(P) == %.3f',p),parse=TRUE,hjust=hjust,vjust=2.8,size=FONT/(72.27/25.4),colour=INK))}
ps<-readRDS(file.path(PROJECT_DIR,'Output/physeq.rds'));m<-data.frame(sample_data(ps),check.names=FALSE);a<-as(otu_table(ps),'matrix');if(taxa_are_rows(ps))a<-t(a);tax<-data.frame(tax_table(ps));tax<-tax[colnames(a),];stopifnot(identical(rownames(m),rownames(a)),nrow(a)>1,!anyDuplicated(m$Patient_id),!anyDuplicated(tax$Species),all(a>=0),all(is.finite(a)),max(abs(rowSums(a)-1))<1e-6)
m$age_group<-factor(ifelse(m$age_MBC_dx<50,'<50','>=50'),levels=c('<50','>=50'));m$race<-factor(m$race,levels=c('Black','White'));m$ER<-factor(m$ER,levels=c('N','Y'));m$ctc2<-factor(ifelse(is.na(m$CTC_count),NA,ifelse(m$CTC_count<=5,'0-5','>5')),levels=c('0-5','>5'));m$Shannon<-apply(a,1,function(z)-sum(z[z>0]*log(z[z>0])))
vars<-c('age_group','race','ER','ctc2');labels<-c(age_group='Age',race='Race',ER='ER',ctc2='CTC');contrast_names<-c(age_group='Age: >=50 vs <50',race='Race: White vs Black',ER='ER: positive vs negative',ctc2='CTC: >5 vs 0-5')

top<-names(sort(colMeans(a),decreasing=TRUE))[1:20];topnames<-pretty_species(top);heat<-list();alpha<-list();beta<-list();disp<-list();inclusion<-list();oldbars<-list();ordination<-list()
for(i in seq_along(vars)){
 v<-vars[i];ok<-!is.na(m[[v]]);md<-droplevels(m[ok,]);aa<-a[ok,,drop=FALSE];g<-md[[v]];lev<-levels(g)
 for(f in top){x<-aa[,f];r<-x[g==lev[1]];c<-x[g==lev[2]];p<-if(length(unique(x))==1)1 else wilcox.test(c,r,alternative='two.sided',exact=FALSE,correct=TRUE)$p.value;heat[[length(heat)+1]]<-data.frame(variable=v,contrast=contrast_names[[v]],feature=f,species=pretty_species(f),n_reference=length(r),n_comparison=length(c),mean_reference=mean(r),mean_comparison=mean(c),median_reference=median(r),median_comparison=median(c),nonzero_reference=sum(r>0),nonzero_comparison=sum(c>0),difference_mean_pp=100*(mean(c)-mean(r)),log2_mean_ratio=log2((mean(c)+PC)/(mean(r)+PC)),P=p)}
 w<-wilcox.test(md$Shannon~g);alpha[[v]]<-data.frame(variable=v,n=nrow(md),P=w$p.value,test=w$method)
 ds<-vegdist(aa,'bray');set.seed(123);b<-adonis2(ds~g,permutations=999);beta[[v]]<-data.frame(variable=v,n=nrow(md),R2=b$R2[1],F=b$F[1],P=b$`Pr(>F)`[1],seed=123,permutations=999)
 set.seed(SEED);bd<-permutest(betadisper(ds,g,bias.adjust=TRUE),permutations=999);disp[[v]]<-data.frame(variable=v,n=nrow(md),P=bd$tab$`Pr(>F)`[1],F=bd$tab$F[1])
 inclusion[[v]]<-data.frame(variable=v,sample_id=md$Sample_id,group=as.character(g))

 gm<-rbind(colMeans(aa[g==lev[1],,drop=FALSE]),colMeans(aa[g==lev[2],,drop=FALSE]));rownames(gm)<-lev;tt<-names(sort(colMeans(gm),decreasing=TRUE))[1:20];gm2<-cbind(gm[,tt,drop=FALSE],Others=1-rowSums(gm[,tt,drop=FALSE]));oldbars[[v]]<-data.frame(variable=v,group=rep(lev,each=21),species=rep(c(pretty_species(tt),'Others'),2),percent=as.vector(t(gm2))*100)

 dd<-as.dist(as.matrix(ds));ord<-phyloseq::ordinate(prune_samples(rownames(md),ps),method='PCoA',distance=dd);ev<-ord$values$Eigenvalues;pts<-as.data.frame(ord$vectors[,1:2]);ordination[[v]]<-data.frame(pts,sample_id=rownames(pts),group=g[match(rownames(pts),rownames(md))]);attr(ordination[[v]],'variance')<-100*ev[1:2]/sum(ev)
}
heat<-do.call(rbind,heat)

heat$q_BH80<-p.adjust(heat$P,'BH')
m$subtype<-factor(ifelse(m$HER2=='Y','HER2-positive',ifelse(m$ER=='Y'|m$PR=='Y','Luminal','TNBC')),levels=c('Luminal','HER2-positive','TNBC'))
stopifnot(all(table(m$subtype)>0))
pairs<-list(sub_HER2_vs_Luminal=c('Luminal','HER2-positive'),sub_TNBC_vs_Luminal=c('Luminal','TNBC'),sub_TNBC_vs_HER2=c('HER2-positive','TNBC'))
new_names<-c(sub_HER2_vs_Luminal='Subtype: HER2+ vs Luminal',sub_TNBC_vs_Luminal='Subtype: TNBC vs Luminal',sub_TNBC_vs_HER2='Subtype: TNBC vs HER2+')
contrast_names<-c(contrast_names,new_names);heat_vars<-c(vars,names(pairs))
for(v in names(pairs)){
 lev<-pairs[[v]]
 for(f in top){
  r<-a[m$subtype==lev[1],f];cc<-a[m$subtype==lev[2],f]
  pp<-if(length(unique(c(r,cc)))==1)1 else wilcox.test(cc,r,alternative='two.sided',exact=FALSE,correct=TRUE)$p.value
  heat<-rbind(heat,data.frame(variable=v,contrast=contrast_names[[v]],feature=f,species=pretty_species(f),n_reference=length(r),n_comparison=length(cc),mean_reference=mean(r),mean_comparison=mean(cc),median_reference=median(r),median_comparison=median(cc),nonzero_reference=sum(r>0),nonzero_comparison=sum(cc>0),difference_mean_pp=100*(mean(cc)-mean(r)),log2_mean_ratio=log2((mean(cc)+PC)/(mean(r)+PC)),P=pp,q_BH80=NA_real_))
 }
}

heat$q_BH140<-p.adjust(heat$P,'BH')
heat$stars_q_BH140<-ifelse(heat$q_BH140<.001,'***',ifelse(heat$q_BH140<.01,'**',ifelse(heat$q_BH140<.05,'*','')))

heat$stars<-ifelse(heat$P<.001,'***',ifelse(heat$P<.01,'**',ifelse(heat$P<.05,'*','')))

omnibus<-do.call(rbind,lapply(top,function(f){kw<-kruskal.test(a[,f],m$subtype);data.frame(feature=f,species=pretty_species(f),n=nrow(m),H=unname(kw$statistic),df=unname(kw$parameter),P=kw$p.value)}))
omnibus$q_BH20<-p.adjust(omnibus$P,'BH');wt(omnibus,'Fig1_subtype_omnibus')
wt(heat[heat$variable%in%names(pairs),],'Fig1_subtype_pairwise')
alpha<-do.call(rbind,alpha);beta<-do.call(rbind,beta);disp<-do.call(rbind,disp);bars<-do.call(rbind,oldbars)
corage<-cor.test(m$age_MBC_dx,m$Shannon,method='spearman',exact=FALSE);age<-data.frame(n=nrow(m),rho=unname(corage$estimate),P=corage$p.value)
wt(heat,'Fig1_shared_Top20_effects_and_tests');wt(data.frame(rank=seq_along(top),feature=top,species=topnames,overall_mean=colMeans(a)[top]),'Top20_selection');wt(beta,'PERMANOVA_original_comparisons');wt(disp,'PERMDISP_diagnostics');wt(alpha,'Shannon_group_comparisons');wt(age,'Shannon_continuous_age');wt(do.call(rbind,inclusion),'sample_inclusion','Audit');wt(bars,'Supplementary_composition_values');wt(data.frame(sample_id=m$Sample_id,Shannon=m$Shannon,age=m$age_MBC_dx),'Shannon_values','Audit')

heat$contrast<-factor(heat$contrast,levels=rev(contrast_names[heat_vars]))
cluster_matrix<-sapply(heat_vars,function(v)heat$log2_mean_ratio[match(top,heat$feature[heat$variable==v])+which(heat$variable==v)[1]-1])
rownames(cluster_matrix)<-topnames;colnames(cluster_matrix)<-heat_vars
stopifnot(all(is.finite(cluster_matrix)))
hc<-hclust(dist(cluster_matrix,method='euclidean'),method='ward.D2')
species_order<-hc$labels[hc$order]
wt(data.frame(order=1:20,species=species_order),'Fig1_species_cluster_order')
wt(data.frame(species=rownames(cluster_matrix),cluster_matrix,check.names=FALSE),'Fig1_clustering_matrix')
saveRDS(hc,file.path(OUT_DIR,'Audit/Fig1_species_hclust.rds'))

node_x<-numeric(19);node_y<-hc$height;segments<-list()
for(i in 1:19){
 ch<-hc$merge[i,];xx<-sapply(ch,function(z)if(z<0)match(-z,hc$order) else node_x[z]);yy<-sapply(ch,function(z)if(z<0)0 else node_y[z]);node_x[i]<-mean(xx)
 segments[[length(segments)+1]]<-data.frame(x=c(xx[1],xx[2],xx[1]),xend=c(xx[1],xx[2],xx[2]),y=c(yy[1],yy[2],node_y[i]),yend=rep(node_y[i],3))
}
segments<-do.call(rbind,segments)
dend<-ggplot(segments)+geom_segment(aes(x,y,xend=xend,yend=yend),colour=INK,linewidth=.35)+scale_x_continuous(limits=c(.5,20.5),expand=c(0,0))+scale_y_continuous(expand=expansion(mult=c(0,.02)))+theme_void()+theme(plot.margin=margin(2,8,0,6))
heat$species<-factor(heat$species,levels=species_order)
cluster_heatmap<-function(data,star_column){
 data$species<-factor(data$species,levels=species_order)
 h<-ggplot(data,aes(species,contrast,fill=log2_mean_ratio))+geom_tile(colour='white',linewidth=.6)+geom_text(aes(label=.data[[star_column]]),colour=INK,size=FONT/(72.27/25.4))+scale_fill_gradient2(low=BLUE,mid='white',high=PINK,midpoint=0,limits=c(-3,3),breaks=c(-3,0,3),labels=c('<= -3','0','>= 3'),oob=scales::squish,name='Log2 mean\nabundance ratio')+scale_x_discrete(expand=c(0,0))+scale_y_discrete(labels=setNames(c('Age: >=50 vs <50','Race: White vs Black','ER: ER+ vs ER-','CTC: >5 vs 0-5','HER2+ vs Luminal','TNBC vs Luminal','TNBC vs HER2+'),contrast_names[heat_vars]),expand=c(0,0))+labs(x=NULL,y=NULL)+theme_fig()+theme(axis.text.x=element_text(face='italic',angle=45,hjust=1,vjust=1,size=FONT,colour=INK),axis.ticks=element_blank(),axis.line=element_blank(),legend.position='right',legend.key.height=grid::unit(.7,'cm'),plot.margin=margin(0,8,6,6))

 legend<-get_legend(h)
 h<-h+theme(legend.position='none')
 body<-plot_grid(dend,h,ncol=1,rel_heights=c(.65,4.6),align='v',axis='lr')
 plot_grid(body,legend,nrow=1,rel_widths=c(1,.12))
}
p1<-cluster_heatmap(heat,'stars')
save_pdf(p1,'Main','Fig2',15.5,8.6)

B <- 19999L
set.seed(SEED)
sensitivity <- list()
for (v in vars) {
 ok <- !is.na(m[[v]]); xx <- a[ok,top,drop=FALSE]; g <- droplevels(m[[v]][ok])
 n <- nrow(xx); k <- sum(g==levels(g)[2]); nr <- n-k
 G <- t(replicate(B, {z<-numeric(n);z[sample.int(n,k)]<-1;z}))
 obsG <- as.numeric(g==levels(g)[2])
 rankx <- apply(xx,2,rank,ties.method='average')
 rankobs <- colSums(rankx*obsG)/k-colSums(rankx*(1-obsG))/nr
 rankperm <- G %*% rankx/k-(matrix(colSums(rankx),B,ncol(xx),byrow=TRUE)-G%*%rankx)/nr
 rankP <- (1+colSums(abs(rankperm)>=matrix(abs(rankobs)-1e-12,B,ncol(xx),byrow=TRUE)))/(B+1)
 cm <- G%*%xx/k; rm <- (matrix(colSums(xx),B,ncol(xx),byrow=TRUE)-G%*%xx)/nr
 cv <- pmax(0,(G%*%(xx^2)-k*cm^2)/(k-1))
 rv <- pmax(0,((matrix(colSums(xx^2),B,ncol(xx),byrow=TRUE)-G%*%(xx^2))-nr*rm^2)/(nr-1))
 ts <- (cm-rm)/sqrt(cv/k+rv/nr)
 for (j in seq_along(top)) {
  xc<-xx[obsG==1,j];xr<-xx[obsG==0,j]
  w<-tryCatch(t.test(xc,xr,var.equal=FALSE),error=function(e)NULL)
  if (is.null(w)) {wp<-NA_real_;to<-NA_real_;pp<-NA_real_;note<-'Welch not estimable'} else {
   wp<-w$p.value;to<-unname(w$statistic)
   pp<-(1+sum(abs(ts[,j])>=abs(to)-1e-12,na.rm=TRUE))/(B+1);note<-'ok'
  }
  h<-heat[heat$variable==v & heat$feature==top[j],]
  sensitivity[[length(sensitivity)+1]]<-data.frame(variable=v,species=pretty_species(top[j]),method=c('Wilcoxon_asymptotic_primary','Permutation_rank','Welch_mean','Permutation_Welch'),P=c(h$P,rankP[j],wp,pp),note=c('primary','exchangeable-label null',note,'exchangeable-label null'),permutations=c(NA,B,NA,B))
 }
}
sensitivity<-do.call(rbind,sensitivity)
sensitivity$q_BH80<-ave(sensitivity$P,sensitivity$method,FUN=function(z)p.adjust(z,'BH'))
summary_methods<-do.call(rbind,lapply(split(sensitivity,sensitivity$method),function(d)data.frame(method=d$method[1],tests=nrow(d),estimable=sum(is.finite(d$P)),nominal_P_lt_05=sum(d$P<.05,na.rm=TRUE),BH80_q_lt_05=sum(d$q_BH80<.05,na.rm=TRUE),unique_species_nominal=length(unique(d$species[which(d$P<.05)])),unique_species_FDR=length(unique(d$species[which(d$q_BH80<.05)])))))
wt(sensitivity,'Fig1_all_method_sensitivity');wt(summary_methods,'Fig1_method_summary')
wt(heat[heat$P<.05,],'Fig1_nominal_P_candidates')

nominal<-heat;nominal$nominal_stars<-ifelse(nominal$P<.001,'***',ifelse(nominal$P<.01,'**',ifelse(nominal$P<.05,'*','')))
pnominal<-cluster_heatmap(nominal,'nominal_stars')
save_pdf(pnominal,'Supplementary','Fig1_nominal_P_exploratory',15.5,8.6)

legacy_alpha<-read.csv(file.path(PROJECT_DIR,'Output/physeq_alpha_div.csv'))
alpha_match<-match(m$Sample_id,legacy_alpha$Sample_id)
stopifnot(!anyNA(alpha_match),!anyDuplicated(legacy_alpha$Sample_id))
m$Observed<-legacy_alpha$observed[alpha_match]
stopifnot(all(m$Observed==rowSums(a>0)))
wt(data.frame(sample_id=m$Sample_id,Observed=m$Observed,source='existing physeq_alpha_div.csv; corrected-profile nonzero species'), 'Observed_source','Audit')
alpha_overview<-list();beta_overview<-beta
bv<-vegdist(a,'bray');set.seed(123);bb<-adonis2(bv~age_MBC_dx,data=m,permutations=999)
beta_overview[beta_overview$variable=='age_group',c('R2','F','P')]<-c(bb$R2[1],bb$F[1],bb$`Pr(>F)`[1])
beta_overview$model<-ifelse(beta_overview$variable=='age_group','continuous age','original binary contrast')
for (v in vars) for (metric in c('Shannon','Observed')) {
 x<-if(v=='age_group') m$age_MBC_dx else as.numeric(m[[v]])-1
 ok<-is.finite(x)&is.finite(m[[metric]]);co<-cor.test(x[ok],m[[metric]][ok],method='spearman',exact=FALSE)
 alpha_overview[[length(alpha_overview)+1]]<-data.frame(variable=v,metric=metric,n=sum(ok),rho=unname(co$estimate),P=co$p.value,coding=if(v=='age_group')'continuous age (years)' else contrast_names[[v]])
}
alpha_overview<-do.call(rbind,alpha_overview)

alpha_overview$q_BH4_per_metric<-ave(alpha_overview$P,alpha_overview$metric,FUN=function(z)p.adjust(z,'BH'))
beta_overview$q_BH4<-p.adjust(beta_overview$P,'BH')
alpha_overview$effect<-alpha_overview$rho
alpha_overview$method<-'Spearman rho'
alpha_overview$statistic<-NA_real_

m$subtype<-factor(ifelse(m$HER2=='Y','HER2-positive',ifelse(m$ER=='Y'|m$PR=='Y','Luminal','TNBC')),levels=c('Luminal','HER2-positive','TNBC'))
stopifnot(all(table(m$subtype)>0))
for(metric in c('Shannon','Observed')) {
 kw<-kruskal.test(m[[metric]],m$subtype)

 eps<-min(1,max(0,unname(kw$statistic)/(nrow(m)-1)))
 alpha_overview<-rbind(alpha_overview,data.frame(variable='subtype',metric=metric,n=nrow(m),rho=NA_real_,P=kw$p.value,coding='Luminal / HER2-positive / TNBC (unordered)',q_BH4_per_metric=NA_real_,effect=eps,method='Kruskal-Wallis rank epsilon squared: H/(n-1)',statistic=unname(kw$statistic)))
}
set.seed(666);bs<-adonis2(bv~subtype,data=m,permutations=999)
bs_ref_path<-file.path(PROJECT_DIR,'Revision2/Major4/output/Tables/PERMANOVA_tests.csv')
bsref<-read.csv(bs_ref_path);bsref<-bsref[bsref$variable=='subtype',]
stopifnot(abs(bs$R2[1]-bsref$R2)<1e-12,bs$`Pr(>F)`[1]==bsref$p)
beta_overview<-rbind(beta_overview,data.frame(variable='subtype',n=nrow(m),R2=bs$R2[1],F=bs$F[1],P=bs$`Pr(>F)`[1],seed=666,permutations=999,model='three-group clinical subtype; original Major4',q_BH4=NA_real_))
set.seed(666);subdisp<-permutest(betadisper(bv,m$subtype,bias.adjust=TRUE),permutations=999)
wt(data.frame(variable='subtype',n=nrow(m),F=subdisp$tab$F[1],P=subdisp$tab$`Pr(>F)`[1]),'Fig2_subtype_dispersion')
alpha_overview$q_BH5_per_metric<-ave(alpha_overview$P,alpha_overview$metric,FUN=function(z)p.adjust(z,'BH'))
beta_overview$q_BH5<-p.adjust(beta_overview$P,'BH')
star<-function(p)ifelse(p<.001,'***',ifelse(p<.01,'**',ifelse(p<.05,'*','')))
order_v<-beta_overview$variable[order(beta_overview$R2)]
row_labels<-c(age_group='Age (years)',race='Race',ER='ER status',ctc2='CTC category',subtype='Clinical subtype')
alpha_overview$row<-factor(alpha_overview$variable,levels=order_v)
alpha_overview$metric<-factor(alpha_overview$metric,levels=c('Shannon','Observed'))
alpha_overview$cell<-paste0(sprintf('%.2f',alpha_overview$effect),star(alpha_overview$P))

alpha_overview$fill_hex<-vapply(seq_len(nrow(alpha_overview)),function(i){e<-alpha_overview$effect[i];if(alpha_overview$variable[i]=='subtype')scales::col_numeric(c('#FAFAFA','#518B86'),domain=c(0,1))(e) else scales::col_numeric(c(BLUE,'#FAFAFA',PINK),domain=c(-.5,.5))(max(-.5,min(.5,e)))},character(1))
beta_overview$row<-factor(beta_overview$variable,levels=order_v)
beta_overview$stars<-star(beta_overview$P)
wt(alpha_overview,'Fig2A_alpha_overview');wt(beta_overview,'Fig2A_PERMANOVA_overview')

ph<-ggplot(alpha_overview,aes(metric,row,fill=fill_hex))+geom_tile(colour=INK,linewidth=.35,width=.94,height=.84)+geom_text(aes(label=sprintf('%.2f',effect)),size=FONT/(72.27/25.4),colour=INK)+geom_text(aes(label=star(P)),nudge_x=.24,nudge_y=.08,size=15/(72.27/25.4),colour=INK)+scale_fill_identity()+scale_x_discrete(position='top')+scale_y_discrete(labels=row_labels,expand=expansion(add=.55))+labs(x=NULL,y=NULL)+theme_fig()+theme(axis.text.x=element_text(angle=90,hjust=0,vjust=.5),axis.line=element_blank(),axis.ticks=element_blank(),plot.margin=margin(6,2,6,6))
pbar<-ggplot(beta_overview,aes(R2,row))+geom_col(width=.70,fill='#78BDB6',colour=INK,linewidth=.4)+geom_text(aes(x=R2+.002,label=stars),hjust=0,size=15/(72.27/25.4),colour=INK)+scale_x_continuous(position='top',labels=scales::label_percent(accuracy=1),limits=c(0,max(beta_overview$R2)*1.18),expand=c(0,0))+scale_y_discrete(expand=expansion(add=.55))+labs(x='Bray-Curtis PERMANOVA R²',y=NULL)+theme_fig()+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),plot.margin=margin(6,12,6,2))

gh<-ggplotGrob(ph);gb<-ggplotGrob(pbar)
for(nm in c('axis-t','xlab-t')) {ih<-gh$layout[gh$layout$name==nm,'t'];ib<-gb$layout[gb$layout$name==nm,'t'];if(length(ih)&&length(ib)){hh<-grid::unit.pmax(gh$heights[ih],gb$heights[ib]);gh$heights[ih]<-hh;gb$heights[ib]<-hh}}
pOverview<-plot_grid(gh,gb,nrow=1,rel_widths=c(1.2,1.4),align='h',axis='tb')

pa<-ggplot(m,aes(age_MBC_dx,Shannon))+geom_smooth(method='lm',se=TRUE,colour=TEAL,fill='#DBE8E5',linewidth=.85)+geom_point(shape=21,fill=BLUE,colour='white',stroke=.3,size=2.5,alpha=.9)+labs(x='Age at MBC diagnosis (years)',y='Shannon diversity')+theme_fig()+theme(aspect.ratio=1)+stat_annotation(sprintf('Spearman rho = %.3f',age$rho),age$P,hjust=-.05)
erdd<-ordination[['ER']];ervv<-attr(erdd,'variance');erdd$group<-factor(erdd$group,levels=c('N','Y'))
pc<-ggplot(erdd,aes(Axis.1,Axis.2,colour=group))+geom_point(size=2.5,alpha=.9)+stat_ellipse(type='norm',linetype=2,linewidth=.6)+scale_colour_manual(values=c(N=BLUE,Y=PINK),labels=c('ER negative','ER positive'),name=NULL)+labs(x=sprintf('PCoA1 (%.1f%%)',ervv[1]),y=sprintf('PCoA2 (%.1f%%)',ervv[2]))+theme_fig()+theme(legend.position='inside',legend.position.inside=c(.98,.02),legend.justification.inside=c(1,0),legend.direction='vertical',legend.background=element_rect(fill='white',colour=NA),legend.margin=margin(3,3,3,3),aspect.ratio=1)+stat_annotation(sprintf('PERMANOVA R² = %.3f',beta$R2[beta$variable=='ER']),beta$P[beta$variable=='ER'])+guides(colour=guide_legend(ncol=1))

subdd<-erdd;subdd$subtype<-m$subtype[match(subdd$sample_id,rownames(m))]
stopifnot(!anyNA(subdd$subtype))
subcols<-c(Luminal='#0072B2','HER2-positive'='#E69F00',TNBC='#CC3377')
pd<-ggplot(subdd,aes(Axis.1,Axis.2,colour=subtype))+geom_point(size=2.5,alpha=.9)+stat_ellipse(type='norm',linetype=2,linewidth=.6)+scale_colour_manual(values=subcols,name=NULL)+guides(colour=guide_legend(ncol=1))+labs(x=sprintf('PCoA1 (%.1f%%)',ervv[1]),y=sprintf('PCoA2 (%.1f%%)',ervv[2]))+theme_fig()+theme(legend.position='inside',legend.position.inside=c(.98,.02),legend.justification.inside=c(1,0),legend.direction='vertical',legend.background=element_rect(fill='white',colour=NA),legend.margin=margin(3,3,3,3),aspect.ratio=1)+stat_annotation(sprintf('PERMANOVA R² = %.3f',bs$R2[1]),bs$`Pr(>F)`[1])
wt(subdd,'Fig2D_subtype_coordinates')

main2<-plot_grid(pOverview,pa,pc,pd,nrow=2,ncol=2,labels=c('A','B','C','D'),label_size=FONT,label_colour=INK,align='hv',axis='tblr',rel_widths=c(1.15,1),rel_heights=c(1,1))
save_pdf(main2,'Main','Fig3',13.8,11.8)

boxplot_alpha<-function(v){ok<-!is.na(m[[v]]);d<-droplevels(m[ok,]);d$g<-d[[v]];lev<-levels(d$g);lab<-switch(v,age_group='Age at MBC diagnosis',race='Race',ER='ER status',ctc2='CTCs / 7.5 mL');ll<-switch(v,age_group=c('<50','>=50'),race=c('Black','White'),ER=c('Negative','Positive'),ctc2=c('0-5','>5'));tt<-alpha[alpha$variable==v,];set.seed(SEED);ggplot(d,aes(g,Shannon,fill=g,colour=g))+geom_violin(trim=TRUE,alpha=.18,linewidth=.4)+geom_boxplot(width=.18,fill='white',outlier.shape=NA,colour=INK,linewidth=.4)+geom_jitter(width=.1,size=1.7,alpha=.85)+scale_fill_manual(values=setNames(c(BLUE,PINK),lev),guide='none')+scale_colour_manual(values=setNames(c(BLUE,PINK),lev),guide='none')+scale_x_discrete(labels=setNames(paste0(ll,'\n(n=',as.integer(table(d$g)),')'),lev))+labs(x=lab,y='Shannon diversity')+theme_fig()+annotate('text',x=1,y=max(d$Shannon)+.35,label=sprintf("paste('Wilcoxon ', italic(P), ' = %.3f')",tt$P),parse=TRUE,hjust=0,size=FONT/(72.27/25.4),colour=INK)+scale_y_continuous(expand=expansion(mult=c(.02,.15)))}
sa<-plot_grid(plot_grid(pa+theme(aspect.ratio=NULL),labels='A',label_size=FONT,label_colour=INK),plot_grid(plotlist=lapply(vars,boxplot_alpha),ncol=2,labels=LETTERS[2:5],label_size=FONT,label_colour=INK),ncol=1,rel_heights=c(1,2));save_pdf(sa,'Supplementary','FigS2_alpha_diversity',11.5,11)

pcs<-lapply(vars,function(v){dd<-ordination[[v]];vv<-attr(dd,'variance');lev<-levels(m[[v]]);dd$group<-factor(dd$group,levels=lev);ll<-switch(v,age_group=c('<50','>=50'),race=c('Black','White'),ER=c('Negative','Positive'),ctc2=c('0-5','>5'));b<-beta[beta$variable==v,];ggplot(dd,aes(Axis.1,Axis.2,colour=group))+geom_point(size=2.1,alpha=.9)+stat_ellipse(type='norm',linetype=2,linewidth=.55)+scale_colour_manual(values=setNames(c(BLUE,PINK),lev),labels=ll,name=switch(v,age_group='Age group',race='Race',ER='ER status',ctc2='CTC category'))+labs(x=sprintf('PCoA1 (%.1f%%)',vv[1]),y=sprintf('PCoA2 (%.1f%%)',vv[2]))+theme_fig()+theme(legend.position='bottom')+annotate('text',x=-Inf,y=Inf,label=sprintf("paste('R2 = %.3f; ', italic(P), ' = %.3f')",b$R2,b$P),parse=TRUE,hjust=-.02,vjust=1.05,size=FONT/(72.27/25.4),colour=INK)})
save_pdf(plot_grid(plotlist=pcs,ncol=2,labels=LETTERS[1:4],label_size=FONT,label_colour=INK),'Supplementary','FigS3_beta_diversity',11.5,8.5)

allspecies<-unique(bars$species);named<-setdiff(allspecies,'Others');cols<-setNames(grDevices::hcl(h=seq(15,375,length.out=length(named)+1)[-1],c=45,l=65),named);cols<-c(cols,Others='#D9DDE1');legorder<-c(sort(named),'Others')
stacks<-lapply(vars,function(v){dd<-bars[bars$variable==v,];dd$species<-factor(dd$species,levels=legorder);dd$group<-factor(dd$group,levels=levels(m[[v]]));pp<-ggplot(dd,aes(group,percent,fill=species))+geom_col(width=.58,colour='white',linewidth=.15,show.legend=TRUE)+scale_fill_manual(values=cols,drop=FALSE)+scale_y_continuous(breaks=c(0,25,50,75,100),expand=c(0,0))+coord_cartesian(ylim=c(0,100))+labs(x=switch(v,age_group='Age group',race='Race',ER='ER status',ctc2='CTC category'),y='Relative abundance (%)',fill=NULL)+theme_fig()+theme(legend.position='none');pp})
sl<-get_legend(stacks[[1]]+theme(legend.position='bottom',legend.text=element_text(face='italic',size=FONT,colour=INK))+guides(fill=guide_legend(ncol=3)));sb<-plot_grid(plot_grid(plotlist=stacks,ncol=2,labels=LETTERS[1:4],label_size=FONT,label_colour=INK),sl,ncol=1,rel_heights=c(1,.48));save_pdf(sb,'Supplementary','FigS1_composition',12,11.2)

w<-createWorkbook();for(nm in c('Fig1','Top20','PERMANOVA','Dispersion','Shannon','Age','Method_summary','All_methods','Fig2A_alpha','Fig2A_beta','Subtype_omnibus')){z<-switch(nm,Fig1=heat,Top20=data.frame(rank=1:20,species=topnames,mean=colMeans(a)[top]),PERMANOVA=beta,Dispersion=disp,Shannon=alpha,Age=age,Method_summary=summary_methods,All_methods=sensitivity,Fig2A_alpha=alpha_overview,Fig2A_beta=beta_overview,Subtype_omnibus=omnibus);addWorksheet(w,nm);writeData(w,nm,z);freezePane(w,nm,firstRow=TRUE);setColWidths(w,nm,cols=seq_len(ncol(z)),widths=pmin(48,pmax(13,nchar(names(z))+3)))};saveWorkbook(w,file.path(OUT_DIR,'Tables/Fig1-3_results.xlsx'),overwrite=TRUE)
old<-read.csv(file.path(PROJECT_DIR,'Output/physeq_alpha_div.csv'));delta<-max(abs(m$Shannon-old$diversity_shannon[match(m$Sample_id,old$Sample_id)]));qc<-data.frame(check=c('Corrected case input','Shared20species x7 comparisons','Shannon unchanged','CTC complete-case inclusion','PDF only','No MaAsLin3 fit'),pass=c(nrow(a)>1,nrow(heat)==140&&length(unique(heat$feature))==20,delta<1e-8,beta$n[beta$variable=='ctc2']==sum(!is.na(m$ctc2)),!any(grepl('\\.(png|tiff?)$',list.files(OUT_DIR,recursive=TRUE))),!exists('maaslin_fit')));wt(qc,'numerical_checks','QA');stopifnot(all(qc$pass))
sources<-c(bs_ref_path,file.path(PROJECT_DIR,'Output/physeq.rds'),file.path(PROJECT_DIR,'Output/physeq_alpha_div.csv'))
wt(data.frame(path=sources,md5=unname(tools::md5sum(sources))),'input_manifest','Audit')
writeLines(capture.output(sessionInfo()),file.path(OUT_DIR,'Audit/sessionInfo.txt'))
saveRDS(list(heatmap=p1,main3=main2,main3_panels=list(A=pOverview,B=pa,C=pc,D=pd),alpha_supplement=sa,beta_supplement=pcs,composition_supplement=sb),file.path(OUT_DIR,'Audit/plot_objects.rds'))
