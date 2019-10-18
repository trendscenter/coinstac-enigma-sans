#R
#system('mkdir -p /Users/hgazula/Desktop/PROJECTS/Publications/papers/TheoG.M.vanErp/enigma_cortical/results_cortical/EffectSizes_SZ_only_SANSTOT_withAge_thickness_asis_105ROIs')
args = commandArgs(trailingOnly=TRUE)
baseDir = args[1]
transferDir = args[2]
outputDir = args[3]
setwd(outputDir)

pdf(paste(outputDir, 'EffectSizes_SZ_only_SANSTOT_withAge_thickness_asis_105ROIs.pdf', sep='/'))

#need the metafor library
library(metafor)

## Added by Harsh
#system('mkdir -p /Users/hgazula/Desktop/PROJECTS/Publications/papers/TheoG.M.vanErp/enigma_cortical/data3c/ASRB/ASRB_CortRegs3c')
#system('mkdir -p /Users/hgazula/Desktop/PROJECTS/Publications/papers/TheoG.M.vanErp/enigma_cortical/data3c/Dublin/Dublin_CortRegs3c')
## Added by Harsh
#Read in effect size tables for each group

ASRB <- new.env()
load(paste(baseDir, 'local0', "EffectSizes_SZ_only_SANSTOT_CONVERTEASY_withAge_thickness_asis.Rdata", sep='/'),ASRB)
ASRB_r.noicv <- get('r.cort',ASRB)
ASRBm <- new.env()
load(paste(baseDir, 'local0', "Models_SZ_only_SANSTOT_CONVERTEASY_withAge_thickness_asis.Rdata", sep='/'),ASRBm)
ASRBm_models.cort <- get('models.cort',ASRBm)
# need to write a for loop: for i in 1:108 that fills up ASRB_npat
# length(ASRBm_models.cort[[i]]$model[,2])
ASRB_npat=rep(0,108)
for (i in 1:108) {
  ASRB_npat[i]<-length(ASRBm_models.cort[[i]]$model[,2])
}

Dublin <- new.env()
load(paste(baseDir, 'local1', "EffectSizes_SZ_only_SANSTOT_withAge_thickness_asis.Rdata", sep='/'),Dublin)
Dublin_r.noicv <- get('r.cort',Dublin)
Dublinm <- new.env()
load(paste(baseDir, 'local1', "Models_SZ_only_SANSTOT_withAge_thickness_asis.Rdata", sep='/'),Dublinm)
Dublinm_models.cort <- get('models.cort',Dublinm)
# need to write a for loop: for i in 1:108 that fills up Dublin_npat
# length(Dublinm_models.cort[[i]]$model[,2])
Dublin_npat=rep(0,108)
for (i in 1:108) {
  Dublin_npat[i]<-length(Dublinm_models.cort[[i]]$model[,2])
}

#Adjusted Cohen's d MA without ICV covariate

#d.noicv=cbind(AMC_d.noicv, ASRB_d.noicv, CAMH_d.noicv, CIAM_d.noicv, CLING_d.noicv, COBRE_d.noicv, Dublin_d.noicv, ESO_d.noicv, EdinburghEHRS_d.noicv, EdinburghFunc_d.noicv, EdinburghSFMH_d.noicv, FBIRN_d.noicv, FIDMAG_d.noicv, Frankfurt_d.noicv, GAP_d.noicv, Galway_d.noicv, HMS_d.noicv, HUBIN_d.noicv, Huilong1_d.noicv, Huilong2_d.noicv, KASP_d.noicv, MCIC_d.noicv, MPRC1_d.noicv, MPRC2_d.noicv, NU_d.noicv, OLIN_d.noicv, Osaka_d.noicv, PAFIP1.5T_d.noicv, PAFIP3T_d.noicv, RSCZ_d.noicv, RomeSL_d.noicv, SCORE_d.noicv, SNUH_d.noicv, SaoPaulo_d.noicv, TOP_d.noicv, UMCU_d.noicv, UMCUS_d.noicv, UNIBA_d.noicv, UPENN_d.noicv)
#r.noicv=cbind(AMC_r.noicv, CAMH_r.noicv)
r.noicv=cbind(ASRB_r.noicv, Dublin_r.noicv)

# exclude surface area regions
r.noicv<-r.noicv[-c(71,72,108),]

# no SE needed for r meta-analyses, only n
#se.noicv=cbind(AMC_se.noicv, ASRB_se.noicv, CAMH_se.noicv, CIAM_se.noicv, CLING_se.noicv, COBRE_se.noicv, Dublin_se.noicv, ESO_se.noicv, EdinburghEHRS_se.noicv, EdinburghFunc_se.noicv, EdinburghSFMH_se.noicv, FBIRN_se.noicv, FIDMAG_se.noicv, Frankfurt_se.noicv, GAP_se.noicv, Galway_se.noicv, HMS_se.noicv, HUBIN_se.noicv, Huilong1_se.noicv, Huilong2_se.noicv, KASP_se.noicv, MCIC_se.noicv, MPRC1_se.noicv, MPRC2_se.noicv, NU_se.noicv, OLIN_se.noicv, Osaka_se.noicv, PAFIP1.5T_se.noicv, PAFIP3T_se.noicv, RSCZ_se.noicv, RomeSL_se.noicv, SCORE_se.noicv, SNUH_se.noicv, SaoPaulo_se.noicv, TOP_se.noicv, UMCU_se.noicv, UMCUS_se.noicv, UNIBA_se.noicv, UPENN_se.noicv)

# exclude surface area regions
#se.noicv<-se.noicv[-c(71,72,108),]

# no controls
#n.ctl=cbind(AMC_nctl, ASRB_nctl, CAMH_nctl, CIAM_nctl, CLING_nctl, COBRE_nctl, Dublin_nctl, ESO_nctl, EdinburghEHRS_nctl, EdinburghFunc_nctl, EdinburghSFMH_nctl, FBIRN_nctl, FIDMAG_nctl, Frankfurt_nctl, GAP_nctl, Galway_nctl, HMS_nctl, HUBIN_nctl, Huilong1_nctl, Huilong2_nctl, KASP_nctl, MCIC_nctl, MPRC1_nctl, MPRC2_nctl, NU_nctl, OLIN_nctl, Osaka_nctl, PAFIP1.5T_nctl, PAFIP3T_nctl, RSCZ_nctl, RomeSL_nctl, SCORE_nctl, SNUH_nctl, SaoPaulo_nctl, TOP_nctl, UMCU_nctl, UMCUS_nctl, UNIBA_nctl, UPENN_nctl)

# exclude surface area regions
#n.ctl<-n.ctl[-c(71,72,108),]

#n.pat=cbind(AMC_npat, ASRB_npat, CAMH_npat, CIAM_npat, CLING_npat, COBRE_npat, Dublin_npat, ESO_npat, EdinburghEHRS_npat, EdinburghFunc_npat, EdinburghSFMH_npat, FBIRN_npat, FIDMAG_npat, Frankfurt_npat, GAP_npat, Galway_npat, HMS_npat, HUBIN_npat, Huilong1_npat, Huilong2_npat, KASP_npat, MCIC_npat, MPRC1_npat, MPRC2_npat, NU_npat, OLIN_npat, Osaka_npat, PAFIP1.5T_npat, PAFIP3T_npat, RSCZ_npat, RomeSL_npat, SCORE_npat, SNUH_npat, SaoPaulo_npat, TOP_npat, UMCU_npat, UMCUS_npat, UNIBA_npat, UPENN_npat)
#n.pat=cbind(AMC_npat,CAMH_npat);
n.pat=cbind(ASRB_npat, Dublin_npat)

# exclude surface area regions
n.pat<-n.pat[-c(71,72,108),]

# r;se;zval;pval;ci.lb;ci.ub;tau2;tause;i2;h2;npat;fdr.p;fdr.p.mean

#meta.noicv.d=rep(0,nrow(d.noicv))
meta.noicv.r=rep(0,nrow(r.noicv))
meta.noicv.se=rep(0,nrow(r.noicv))
meta.noicv.zval=rep(0,nrow(r.noicv))
meta.noicv.pval=rep(0,nrow(r.noicv))
meta.noicv.ci.lb=rep(0,nrow(r.noicv))
meta.noicv.ci.ub=rep(0,nrow(r.noicv))
meta.noicv.tau2=rep(0,nrow(r.noicv))
meta.noicv.tause=rep(0,nrow(r.noicv))
meta.noicv.i2=rep(0,nrow(r.noicv))
meta.noicv.h2=rep(0,nrow(r.noicv))
#meta.noicv.nctl=rep(0,nrow(d.noicv))
meta.noicv.npat=rep(0,nrow(r.noicv))

#Run meta-analysis, random-effects
 for(x in 1:nrow(r.noicv)){

   meta.noicv.r[x] = rma.uni(ri=r.noicv[x,which(!is.na(r.noicv[x,]))],ni=n.pat[x,which(!is.na(r.noicv[x,]))],measure="COR",control=list(maxiter=10000,stepadj=0.00000001),method="REML")$b
   meta.noicv.se[x] = rma.uni(ri=r.noicv[x,which(!is.na(r.noicv[x,]))],ni=n.pat[x,which(!is.na(r.noicv[x,]))],measure="COR",control=list(maxiter=10000,stepadj=0.00000001),method="REML")$se
   meta.noicv.zval[x] = rma.uni(ri=r.noicv[x,which(!is.na(r.noicv[x,]))],ni=n.pat[x,which(!is.na(r.noicv[x,]))],measure="COR",control=list(maxiter=10000,stepadj=0.00000001),method="REML")$zval
   meta.noicv.pval[x] = rma.uni(ri=r.noicv[x,which(!is.na(r.noicv[x,]))],ni=n.pat[x,which(!is.na(r.noicv[x,]))],measure="COR",control=list(maxiter=10000,stepadj=0.00000001),method="REML")$pval
   meta.noicv.ci.lb[x] = rma.uni(ri=r.noicv[x,which(!is.na(r.noicv[x,]))],ni=n.pat[x,which(!is.na(r.noicv[x,]))],measure="COR",control=list(maxiter=10000,stepadj=0.00000001),method="REML")$ci.lb
   meta.noicv.ci.ub[x] = rma.uni(ri=r.noicv[x,which(!is.na(r.noicv[x,]))],ni=n.pat[x,which(!is.na(r.noicv[x,]))],measure="COR",control=list(maxiter=10000,stepadj=0.00000001),method="REML")$ci.ub
   meta.noicv.tau2[x] = rma.uni(ri=r.noicv[x,which(!is.na(r.noicv[x,]))],ni=n.pat[x,which(!is.na(r.noicv[x,]))],measure="COR",control=list(maxiter=10000,stepadj=0.00000001),method="REML")$tau2
   meta.noicv.tause[x] = rma.uni(ri=r.noicv[x,which(!is.na(r.noicv[x,]))],ni=n.pat[x,which(!is.na(r.noicv[x,]))],measure="COR",control=list(maxiter=10000,stepadj=0.00000001),method="REML")$se.tau2
   meta.noicv.i2[x] = rma.uni(ri=r.noicv[x,which(!is.na(r.noicv[x,]))],ni=n.pat[x,which(!is.na(r.noicv[x,]))],measure="COR",control=list(maxiter=10000,stepadj=0.00000001),method="REML")$I2
   meta.noicv.h2[x] = rma.uni(ri=r.noicv[x,which(!is.na(r.noicv[x,]))],ni=n.pat[x,which(!is.na(r.noicv[x,]))],measure="COR",control=list(maxiter=10000,stepadj=0.00000001),method="REML")$H2
   meta.noicv.npat[x] = sum(n.pat[x,which(!is.na(r.noicv[x,]))])
   
      #meta.noicv.d[x] = rma.uni(yi=d.noicv[x,which(!is.na(d.noicv[x,]))],sei=se.noicv[x,which(!is.na(d.noicv[x,]))],control=list(maxiter=10000,stepadj=0.00000001),method="REML")$b
   #              rma.uni(ri=r.noicv[which(!is.na(r.noicv[x,])),x],ni=n.noicv[which(!is.na(n.noicv[x,])),x],measure="COR",control=list(maxiter=10000,stepadj=0.00000001),method="REML")
   # }

 		#meta.noicv.d[x] = rma.uni(yi=d.noicv[x,which(!is.na(d.noicv[x,]))],sei=se.noicv[x,which(!is.na(d.noicv[x,]))],control=list(maxiter=10000,stepadj=0.00000001),method="REML")$b
 		#meta.noicv.se[x] = rma.uni(yi=d.noicv[x,which(!is.na(d.noicv[x,]))],sei=se.noicv[x,which(!is.na(d.noicv[x,]))],control=list(maxiter=10000,stepadj=0.00000001),method="REML")$se
 		#meta.noicv.zval[x] = rma.uni(yi=d.noicv[x,which(!is.na(d.noicv[x,]))],sei=se.noicv[x,which(!is.na(d.noicv[x,]))],control=list(maxiter=10000,stepadj=0.00000001),method="REML")$zval
 		#meta.noicv.pval[x] = rma.uni(yi=d.noicv[x,which(!is.na(d.noicv[x,]))],sei=se.noicv[x,which(!is.na(d.noicv[x,]))],control=list(maxiter=10000,stepadj=0.00000001),method="REML")$pval
 		#meta.noicv.ci.lb[x] = rma.uni(yi=d.noicv[x,which(!is.na(d.noicv[x,]))],sei=se.noicv[x,which(!is.na(d.noicv[x,]))],control=list(maxiter=10000,stepadj=0.00000001),method="REML")$ci.lb
 		#meta.noicv.ci.ub[x] = rma.uni(yi=d.noicv[x,which(!is.na(d.noicv[x,]))],sei=se.noicv[x,which(!is.na(d.noicv[x,]))],control=list(maxiter=10000,stepadj=0.00000001),method="REML")$ci.ub
		#meta.noicv.tau2[x] = rma.uni(yi=d.noicv[x,which(!is.na(d.noicv[x,]))],sei=se.noicv[x,which(!is.na(d.noicv[x,]))],control=list(maxiter=10000,stepadj=0.00000001),method="REML")$tau2
		#meta.noicv.tause[x] = rma.uni(yi=d.noicv[x,which(!is.na(d.noicv[x,]))],sei=se.noicv[x,which(!is.na(d.noicv[x,]))],control=list(maxiter=10000,stepadj=0.00000001),method="REML")$se.tau2
		#meta.noicv.i2[x] = rma.uni(yi=d.noicv[x,which(!is.na(d.noicv[x,]))],sei=se.noicv[x,which(!is.na(d.noicv[x,]))],control=list(maxiter=10000,stepadj=0.00000001),method="REML")$I2
		#meta.noicv.h2[x] = rma.uni(yi=d.noicv[x,which(!is.na(d.noicv[x,]))],sei=se.noicv[x,which(!is.na(d.noicv[x,]))],control=list(maxiter=10000,stepadj=0.00000001),method="REML")$H2
		#meta.noicv.nctl[x] = sum(n.ctl[x,which(!is.na(d.noicv[x,]))])
		#meta.noicv.npat[x] = sum(n.pat[x,which(!is.na(d.noicv[x,]))])
		
    # forest plot
		 if(x > 70){
		   tmp = rma.uni(ri=r.noicv[x,which(!is.na(r.noicv[x,]))],ni=n.pat[x,which(!is.na(r.noicv[x,]))],measure="COR",control=list(maxiter=10000,stepadj=0.00000001),method="REML")
			#tmp = rma.uni(yi=r.noicv[x,which(!is.na(.noicv[x,]))],sei=se.noicv[x,which(!is.na(d.noicv[x,]))],control=list(maxiter=10000,stepadj=0.00000001),method="REML")
			#nams = c("AMC", "ASRB", "CAMH", "CIAM", "CLING", "COBRE", "Dublin", "ESO", "EdinburghEHRS", "EdinburghFunc", "EdinburghSFMH", "FBIRN", "FIDMAG", "Frankfurt", "GAP", "Galway", "HMS", "HUBIN", "Huilong1", "Huilong2", "KASP", "MCIC", "MPRC1", "MPRC2", "NU", "OLIN", "Osaka", "PAFIP1.5T", "PAFIP3T", "RSCZ", "RomeSL", "SCORE", "SNUH", "SaoPaulo", "TOP", "UMCU", "UMCUS", "UNIBA", "UPENN")
      #nams = c("AMC", "CAMH"); 
      nams = c("ASRB", "Dublin")

		  plot.rois<-c("left bankssts","left caudalanteriorcingulate","left caudalmiddlefrontal","left cuneus","left entorhinal","left fusiform","left inferiorparietal","left inferiortemporal","left isthmuscingulate","left lateraloccipital","left lateralorbitofrontal","left lingual","left medialorbitofrontal","left middletemporal","left parahippocampal","left paracentral","left parsopercularis","left parsorbitalis","left parstriangularis","left pericalcarine","left postcentral","left posteriorcingulate","left precentral","left precuneus","left rostralanteriorcingulate","left rostralmiddlefrontal","left superiorfrontal","left superiorparietal","left superiortemporal","left supramarginal","left frontalpole","left temporalpole","left transversetemporal","left insula","right bankssts","right caudalanteriorcingulate","right caudalmiddlefrontal","right cuneus","right entorhinal","right fusiform","right inferiorparietal","right inferiortemporal","right isthmuscingulate","right lateraloccipital","right lateralorbitofrontal","right lingual","right medialorbitofrontal","right middletemporal","right parahippocampal","right paracentral","right parsopercularis","right parsorbitalis","right parstriangularis","right pericalcarine","right postcentral","right posteriorcingulate","right precentral","right precuneus","right rostralanteriorcingulate","right rostralmiddlefrontal","right superiorfrontal","right superiorparietal","right superiortemporal","right supramarginal","right frontalpole","right temporalpole","right transversetemporal","right insula","left thickness","right thickness","mean bankssts","mean caudalanteriorcingulate","mean caudalmiddlefrontal","mean cuneus","mean entorhinal","mean fusiform","mean inferiorparietal","mean inferiortemporal","mean isthmuscingulate","mean lateraloccipital","mean lateralorbitofrontal","mean lingual","mean medialorbitofrontal","mean middletemporal","mean parahippocampal","mean paracentral","mean parsopercularis","mean parsorbitalis","mean parstriangularis","mean pericalcarine","mean postcentral","mean posteriorcingulate","mean precentral","mean precuneus","mean rostralanteriorcingulate","mean rostralmiddlefrontal","mean superiorfrontal","mean superiorparietal","mean superiortemporal","mean supramarginal","mean frontalpole","mean temporalpole","mean transversetemporal","mean insula","mean thickness")
			#png(paste0("meta_icv_forest_thickness_asis", x, ".png"))
		        #pdf(paste0("meta_icv_forest_thickness_asis",".pdf"))
 			forest(tmp, slab=nams[which(!is.na(r.noicv[x,]))], xlim=c(-4,4),xlab="Partial Correlation")
			text(-.4,41,"Sample",5)
			text(5,41,"Partial Correlation [95% CI]",1.6)
			title(plot.rois[x])
 		#	#dev.off()
		 }
 		
 }

## update these lines
tmp.fdr.p=p.adjust(meta.noicv.pval[1:70],method="fdr")
tmp.fdr.p.mean=p.adjust(meta.noicv.pval[71:length(meta.noicv.pval)],method="fdr") 
fdr.p=rep(1,nrow(r.noicv))
for(x in 1:70){
	fdr.p[x]=tmp.fdr.p[x]
}
fdr.p.mean=rep(1,nrow(r.noicv))
for(x in 71:105){
	fdr.p.mean[x]=tmp.fdr.p.mean[x-70]
}
outmat=cbind(meta.noicv.r,meta.noicv.se,meta.noicv.zval,meta.noicv.pval,meta.noicv.ci.lb,meta.noicv.ci.ub,meta.noicv.tau2,meta.noicv.tause,meta.noicv.i2,meta.noicv.h2,meta.noicv.npat,fdr.p,fdr.p.mean)
##

write.table(outmat,file=paste(outputDir, "MetaAnalysis_ENIGMA_SZ_Cortex_thickness_SZ_only_SANSTOT_withAge_thickness_asis_105ROIs.txt", sep='/'),quote=F,sep="\t")
dev.off() # close pdf
