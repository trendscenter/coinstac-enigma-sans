#R

# clear environment
rm(list=ls())

library(tidyverse)
library(rjson)
args = commandArgs(trailingOnly=TRUE)
baseDir = args[1]
outputDir = args[2]
siteMeta <- fromJSON(args[3])

# create an analysis log file
sink(paste(outputDir,"MetaAnalysisLog.txt", sep='/'), split=TRUE)

# load the metafor library
library(metafor)

#
# for every site, for every analysis
# Read in effect size tables for each group

# for each phenotype (thickness, surf, subcort)
# for each SANS measure (SANSSum, each Global and all 5 Global measures, each factor and all 5 Factors)
# run Age + Sex + SANS measure
# run Age + Sex + SANS measure + (mean thick/area/ICV)
# run Age + Sex + SANS measure + CPZ/IQ/AO/AP

# do another loop for effects of sex and symptoms (no brain!)

#############################################
### Meta-analysis of symptom associations ###
#############################################

for (phenoName in c("Cort", "Surf", "SubCort")) {  # brain measure type loop
#for (phenoName in c("SubCort")) {
   for (predictor in c("SANSSum", "SANSMAP", "SANSEXP",
                        "SANSGlobal1", "SANSGlobal2", "SANSGlobal3", "SANSGlobal4",
                        "SANSFac1", "SANSFac2", "SANSFac3", "SANSFac4", "SANSFac5",
                        "SANSMAPwSANSSum", "SANSEXPwSANSSum",
                        "SANSGlobal1wSANSSum", "SANSGlobal2wSANSSum", "SANSGlobal3wSANSSum", "SANSGlobal4wSANSSum",
                        "SANSFac1wSANSSum", "SANSFac2wSANSSum", "SANSFac3wSANSSum", "SANSFac4wSANSSum", "SANSFac5wSANSSum",
                        "SANSGlobal1wSANSEXP", "SANSGlobal2wSANSEXP", "SANSGlobal3wSANSEXP", "SANSGlobal4wSANSEXP",
                        "SANSFac1wSANSEXP", "SANSFac2wSANSEXP", "SANSFac3wSANSEXP", "SANSFac4wSANSEXP", "SANSFac5wSANSEXP",
                        "SANSGlobal1wSANSMAP", "SANSGlobal2wSANSMAP", "SANSGlobal3wSANSMAP", "SANSGlobal4wSANSMAP",
                        "SANSFac1wSANSMAP", "SANSFac2wSANSMAP", "SANSFac3wSANSMAP", "SANSFac4wSANSMAP", "SANSFac5wSANSMAP"
                        )) {  # All the possible predictors!
        for ( WCov in c("NoCovs","WG", "WSum", "WIQ", "WCPZ", "WAO", "WAP")) {  # loop for other covariates: Global brain, IQ, CPZ equiv, AO, AP group

      # generate the input file names
      #cat(paste(predictor, " ", WCov,"\n"))

      Effectsf= paste0("EffectSizes_SZ_only_",predictor,"_withSexAge_",WCov,phenoName,".Rdata") # this depends on the previous output
      Modelf = paste0("Models_SZ_",predictor,"_withSexAge_",WCov,"_", phenoName,".Rdata")

      #cat(paste(Effectsf,"\n"))
      #cat(paste(Modelf, "\n"))

     # output file
      Outfile=paste0("/MetaAnalysis_SZ_",predictor,"_withSexAge_",WCov,phenoName,".txt")
     # out pdf for forest plots
      Outpdf=paste0(outputDir,"/MetaAnalysis_SZ_",predictor,"_withSexAge_",WCov,phenoName,".pdf")
      pdf(Outpdf)

      cat(paste("Creating: ", Outfile,"\n"))

      nsites = 0 # how many sites, for each analysis

      # create empty plot.sites vector for site labels in forest plots
      plot.sites=c()
      for (siteName in names(siteMeta)) {
        siteDir=paste(baseDir, siteName, siteMeta[siteName][[1]]$cohortDirectory, sep='/')


        if (!dir.exists(siteDir)) {
          # site failed to respond
          cat(siteDir,"does not exist \n")
          cat(paste(siteDir,"\n"))
          }
        else { # site exists, keep going
          # open the files and read them in for the effect sizes and the # subjects included
          # if successful, add to the count of nsites for that analysis
          if(!file.exists(paste(siteDir, Effectsf, sep='/')) || !file.exists(paste(siteDir, Modelf, sep='/'))) {
            cat(siteName, "does not have", paste(siteDir, Effectsf, sep='/'),"\n")
            }
          else {
            # read the files--puts n.patients, models.cort, and r.cort into memory
            # I added n.patients to the EffectsF so it is no longer needed to load the Modelf
            load(paste(siteDir, Effectsf, sep='/'))
            load(paste(siteDir, Modelf, sep='/'))

            # for the first site that has data, save the variables
            # for remaining sites, cbind it

           if(!(exists("npat_all"))) {
             npat_all = n.patients
             r_eff=r.cort
             models_all = models.cort
             plot.sites = siteName # siteNames[siteNumber]
           } else {
             npat_all=cbind(npat_all, n.patients)
             r_eff=cbind(r_eff,r.cort)
             models_all = cbind(models_all, models.cort)
             #plot.sites = cbind(plot.sites, siteNames[siteNumber])
             plot.sites = cbind(plot.sites, siteName)
           }

            nsites=nsites+1

          }

        } # end if site exists loop

      } # end  site loop


      if(nsites > 1) { # do the meta analysis and save the results

        #initialize meta analysis structures
        meta.noicv.r=rep(0,nrow(r_eff))
        meta.noicv.se=rep(0,nrow(r_eff))
        meta.noicv.zval=rep(0,nrow(r_eff))
        meta.noicv.pval=rep(0,nrow(r_eff))
        meta.noicv.ci.lb=rep(0,nrow(r_eff))
        meta.noicv.ci.ub=rep(0,nrow(r_eff))
        meta.noicv.tau2=rep(0,nrow(r_eff))
        meta.noicv.tause=rep(0,nrow(r_eff))
        meta.noicv.i2=rep(0,nrow(r_eff))
        meta.noicv.h2=rep(0,nrow(r_eff))
        meta.noicv.npat=rep(0,nrow(r_eff))

        for(x in 1:nrow(models_all)) {
          curr_model = rma.uni(ri=r_eff[x,which(!is.na(r_eff[x,]))],ni=npat_all[x,which(!is.na(r_eff[x,]))],
                               measure="COR",control=list(maxiter=10000,stepadj=0.00000001),method="REML")

          meta.noicv.r[x] = curr_model$b
          meta.noicv.se[x] = curr_model$se
          meta.noicv.zval[x] = curr_model$zval
          meta.noicv.pval[x] = curr_model$pval
          meta.noicv.ci.lb[x] = curr_model$ci.lb
          meta.noicv.ci.ub[x] = curr_model$ci.ub
          meta.noicv.tau2[x] = curr_model$tau2
          meta.noicv.tause[x] = curr_model$se.tau2
          meta.noicv.i2[x] = curr_model$I2
          meta.noicv.h2[x] = curr_model$H2
          meta.noicv.npat[x] = sum(npat_all[x,which(!is.na(r_eff[x,]))])

          # create forest plots
          #plot.sites<-c("1","2","3","4","5")
          if (phenoName == "SubCort") {
          plot.rois<-c("Left Lateral Ventricle", "Right Lateral Ventricle", "Left Thalamus", "Left Caudate", "Left Putamen", "Left Pallidum", "Left Hippocampus", "Left Amygdala", "Left Accumbens", "Right Thalamus", "Right Caudate", "Right Putamen", "Right Palidum", "Right Hippocampus", "Right Amygdala", "Right Accumbens", "Mean Lateral Ventricle", "Mean Thalamus", "Mean Caudate", "Mean Putamen", "Mean Palidum", "Mean Hippocampus", "Mean Amygdala", "Mean Accumbens")
          }
          if (phenoName == "Cort") { # Cort is thickness
          plot.rois<-c("L_bankssts_thickavg","L_caudalanteriorcingulate_thickavg","L_caudalmiddlefrontal_thickavg","L_cuneus_thickavg","L_entorhinal_thickavg","L_fusiform_thickavg","L_inferiorparietal_thickavg","L_inferiortemporal_thickavg","L_isthmuscingulate_thickavg","L_lateraloccipital_thickavg","L_lateralorbitofrontal_thickavg","L_lingual_thickavg","L_medialorbitofrontal_thickavg","L_middletemporal_thickavg","L_parahippocampal_thickavg","L_paracentral_thickavg","L_parsopercularis_thickavg","L_parsorbitalis_thickavg","L_parstriangularis_thickavg","L_pericalcarine_thickavg","L_postcentral_thickavg","L_posteriorcingulate_thickavg","L_precentral_thickavg","L_precuneus_thickavg","L_rostralanteriorcingulate_thickavg","L_rostralmiddlefrontal_thickavg","L_superiorfrontal_thickavg","L_superiorparietal_thickavg","L_superiortemporal_thickavg","L_supramarginal_thickavg","L_frontalpole_thickavg","L_temporalpole_thickavg","L_transversetemporal_thickavg","L_insula_thickavg","R_bankssts_thickavg","R_caudalanteriorcingulate_thickavg","R_caudalmiddlefrontal_thickavg","R_cuneus_thickavg","R_entorhinal_thickavg","R_fusiform_thickavg","R_inferiorparietal_thickavg","R_inferiortemporal_thickavg","R_isthmuscingulate_thickavg","R_lateraloccipital_thickavg","R_lateralorbitofrontal_thickavg","R_lingual_thickavg","R_medialorbitofrontal_thickavg","R_middletemporal_thickavg","R_parahippocampal_thickavg","R_paracentral_thickavg","R_parsopercularis_thickavg","R_parsorbitalis_thickavg","R_parstriangularis_thickavg","R_pericalcarine_thickavg","R_postcentral_thickavg","R_posteriorcingulate_thickavg","R_precentral_thickavg","R_precuneus_thickavg","R_rostralanteriorcingulate_thickavg","R_rostralmiddlefrontal_thickavg","R_superiorfrontal_thickavg","R_superiorparietal_thickavg","R_superiortemporal_thickavg","R_supramarginal_thickavg","R_frontalpole_thickavg","R_temporalpole_thickavg","R_transversetemporal_thickavg","R_insula_thickavg","LThickness","RThickness","M_bankssts_thickavg","M_caudalanteriorcingulate_thickavg","M_caudalmiddlefrontal_thickavg","M_cuneus_thickavg","M_entorhinal_thickavg","M_fusiform_thickavg","M_inferiorparietal_thickavg","M_inferiortemporal_thickavg","M_isthmuscingulate_thickavg","M_lateraloccipital_thickavg","M_lateralorbitofrontal_thickavg","M_lingual_thickavg","M_medialorbitofrontal_thickavg","M_middletemporal_thickavg","M_parahippocampal_thickavg","M_paracentral_thickavg","M_parsopercularis_thickavg","M_parsorbitalis_thickavg","M_parstriangularis_thickavg","M_pericalcarine_thickavg","M_postcentral_thickavg","M_posteriorcingulate_thickavg","M_precentral_thickavg","M_precuneus_thickavg","M_rostralanteriorcingulate_thickavg","M_rostralmiddlefrontal_thickavg","M_superiorfrontal_thickavg","M_superiorparietal_thickavg","M_superiortemporal_thickavg","M_supramarginal_thickavg","M_frontalpole_thickavg","M_temporalpole_thickavg","M_transversetemporal_thickavg","M_insula_thickavg","MThickness","FullSurfArea")
          }
          if (phenoName == "Surf") { # Surf is surface area
          plot.rois<-c("L_bankssts_thickavg","L_caudalanteriorcingulate_thickavg","L_caudalmiddlefrontal_thickavg","L_cuneus_thickavg","L_entorhinal_thickavg","L_fusiform_thickavg","L_inferiorparietal_thickavg","L_inferiortemporal_thickavg","L_isthmuscingulate_thickavg","L_lateraloccipital_thickavg","L_lateralorbitofrontal_thickavg","L_lingual_thickavg","L_medialorbitofrontal_thickavg","L_middletemporal_thickavg","L_parahippocampal_thickavg","L_paracentral_thickavg","L_parsopercularis_thickavg","L_parsorbitalis_thickavg","L_parstriangularis_thickavg","L_pericalcarine_thickavg","L_postcentral_thickavg","L_posteriorcingulate_thickavg","L_precentral_thickavg","L_precuneus_thickavg","L_rostralanteriorcingulate_thickavg","L_rostralmiddlefrontal_thickavg","L_superiorfrontal_thickavg","L_superiorparietal_thickavg","L_superiortemporal_thickavg","L_supramarginal_thickavg","L_frontalpole_thickavg","L_temporalpole_thickavg","L_transversetemporal_thickavg","L_insula_thickavg","R_bankssts_thickavg","R_caudalanteriorcingulate_thickavg","R_caudalmiddlefrontal_thickavg","R_cuneus_thickavg","R_entorhinal_thickavg","R_fusiform_thickavg","R_inferiorparietal_thickavg","R_inferiortemporal_thickavg","R_isthmuscingulate_thickavg","R_lateraloccipital_thickavg","R_lateralorbitofrontal_thickavg","R_lingual_thickavg","R_medialorbitofrontal_thickavg","R_middletemporal_thickavg","R_parahippocampal_thickavg","R_paracentral_thickavg","R_parsopercularis_thickavg","R_parsorbitalis_thickavg","R_parstriangularis_thickavg","R_pericalcarine_thickavg","R_postcentral_thickavg","R_posteriorcingulate_thickavg","R_precentral_thickavg","R_precuneus_thickavg","R_rostralanteriorcingulate_thickavg","R_rostralmiddlefrontal_thickavg","R_superiorfrontal_thickavg","R_superiorparietal_thickavg","R_superiortemporal_thickavg","R_supramarginal_thickavg","R_frontalpole_thickavg","R_temporalpole_thickavg","R_transversetemporal_thickavg","R_insula_thickavg","LThickness","RThickness","M_bankssts_thickavg","M_caudalanteriorcingulate_thickavg","M_caudalmiddlefrontal_thickavg","M_cuneus_thickavg","M_entorhinal_thickavg","M_fusiform_thickavg","M_inferiorparietal_thickavg","M_inferiortemporal_thickavg","M_isthmuscingulate_thickavg","M_lateraloccipital_thickavg","M_lateralorbitofrontal_thickavg","M_lingual_thickavg","M_medialorbitofrontal_thickavg","M_middletemporal_thickavg","M_parahippocampal_thickavg","M_paracentral_thickavg","M_parsopercularis_thickavg","M_parsorbitalis_thickavg","M_parstriangularis_thickavg","M_pericalcarine_thickavg","M_postcentral_thickavg","M_posteriorcingulate_thickavg","M_precentral_thickavg","M_precuneus_thickavg","M_rostralanteriorcingulate_thickavg","M_rostralmiddlefrontal_thickavg","M_superiorfrontal_thickavg","M_superiorparietal_thickavg","M_superiortemporal_thickavg","M_supramarginal_thickavg","M_frontalpole_thickavg","M_temporalpole_thickavg","M_transversetemporal_thickavg","M_insula_thickavg","MThickness","FullSurfArea")
          }
          #forest(curr_model, slab = paste(plot.sites),xlab="Pearson Correlation")
          forest(curr_model, slab = paste(plot.sites),xlim= c(-4,4),xlab="Pearson Correlation")
          ##text(-.4,18,plot.rois([x],5)
          #text(-.4,17,"Sample",5)
          text(5,17,"Pearson R [95% CI]",1.6)
          title(plot.rois[x])
          # cat(paste0(plot.sites), "\n")

        } # end looping through regions
        # turn off output to the pdf
        dev.off()

        ### Review the FDR correction section
        # fdr correction and write to file
        if ((phenoName == "Cort") || (phenoName == "Surf")) {
           tmp.fdr.p=p.adjust(meta.noicv.pval[1:70],method="fdr")
           tmp.fdr.p.mean=p.adjust(meta.noicv.pval[71:length(meta.noicv.pval)],method="fdr")

           fdr.p=rep(1,nrow(r_eff))
           for(x in 1:70){
             fdr.p[x]=tmp.fdr.p[x]
           }

           fdr.p.mean=rep(1,nrow(r_eff))
           for(x in 71:105){
             fdr.p.mean[x]=tmp.fdr.p.mean[x-70]
           }
        } # end if phenoName is Cort (thickness) or Surf (surface area)

        if (phenoName == "SubCort") {
           tmp.fdr.p=p.adjust(meta.noicv.pval[1:16],method="fdr")
           tmp.fdr.p.mean=p.adjust(meta.noicv.pval[17:length(meta.noicv.pval)],method="fdr")

           fdr.p=rep(1,nrow(r_eff))
           for(x in 1:16){
             fdr.p[x]=tmp.fdr.p[x]
           }

           fdr.p.mean=rep(1,nrow(r_eff))
           for(x in 17:24){
             fdr.p.mean[x]=tmp.fdr.p.mean[x-16]
           }
        } # end if phenoName is SubCort

        outmat=cbind(meta.noicv.r,meta.noicv.se,meta.noicv.zval,meta.noicv.pval,
                     meta.noicv.ci.lb,meta.noicv.ci.ub,
                     meta.noicv.tau2,meta.noicv.tause,
                     meta.noicv.i2,meta.noicv.h2,meta.noicv.npat,fdr.p,fdr.p.mean)

        write.table(outmat,file=paste(outputDir, Outfile, sep='/'),
                    quote=F,sep="\t", row.names=TRUE, col.names=TRUE)

      } else {
        cat("no sites for this analysis! \n")
        quit(status=1)
      }

      # clear npat_all, r_eff, and models_all for the next analysis
      #rm(list=c("npat_all", "r_eff","models_all"))
      rm(list=c("npat_all", "r_eff","models_all","n.patients","r.cort","models.cort"))

    } # end WCov loop
  } # end SANS predictor loop
} # end pheno loop

###########################################
####  Meta-Analysis of the sex effects ####
###########################################
# site loop--read it all in, then do the analysis and write out

# clear the variables
rm(list=c("npat_all", "r_eff","models_all"))

# Meta-Analysis of Sex Effects
for (dependent_variable in  c("SANSSum", "SANSMAP", "SANSEXP",
                              "SANSGlobal1", "SANSGlobal2", "SANSGlobal3", "SANSGlobal4",
                              "SANSFac1", "SANSFac2", "SANSFac3", "SANSFac4", "SANSFac5")) {  # All the possible phenoNames!

  nsites = 0 # how many sites, for each analysis
  for (siteName in names(siteMeta)) { # read in sex effect on dependent variable for each site

    siteDir=paste(baseDir, siteName, siteMeta[siteName][[1]]$cohortDirectory, sep='/')
    ModelF=paste0("EffectSizes_SZ_only_Gender_withAge_on_", dependent_variable,".Rdata")
    cat(paste0(ModelF," from ", siteName, "\n"))

    if(!file.exists(paste(siteDir, ModelF, sep='/'))) {
      cat(siteName, "does not have the needed files for", ModelF, siteDir)
      }
    else { # read the files--puts the npat and effect size into memory
      load(paste(siteDir,ModelF, sep='/'))  # puts genderMod model in memory
      #cat(paste(siteDir, ModelF, sep='/'))

      # for the first site that has data, save the variables
      # for remaining sites, cbind it

      if(!(exists("n.pat_all"))) {
        n.pat_all = n.patients # number of participants in the model
        Sex_d_all=cohen_d  # Sex Effect as Cohen's d
        se_all = se
        n.mal_all = n.mal
        n.fem_all = n.fem
        models_all = genderMod
        }
      else {
        n.pat_all=cbind(n.pat_all, n.patients)
        Sex_d_all=cbind(Sex_d_all, cohen_d)
        se_all = cbind(se_all, se)
        n.mal_all = cbind(n.mal_all, n.mal)
        n.fem_all = cbind(n.fem_all, n.fem)
        models_all = cbind(models_all, genderMod)
        } # end what to do in a site
      } # done with that site
      nsites=nsites+1
      #setwd("..") # back to base
    # } # end file.exists
  } # end for loop siteNumber - to load data from all sites

  #initialize meta analysis structures
  meta.d=rep(0,1)
  meta.se=rep(0,1)
  meta.zval=rep(0,1)
  meta.pval=rep(0,1)
  meta.ci.lb=rep(0,1)
  meta.ci.ub=rep(0,1)
  meta.tau2=rep(0,1)
  meta.tause=rep(0,1)
  meta.i2=rep(0,1)
  meta.h2=rep(0,1)
  meta.n.mal=rep(0,1)
  meta.n.fem=rep(0,1)
  meta.n.pat=rep(0,1)

  # do the meta analysis
  metaAnalysis = rma.uni(yi= Sex_d_all, se = se_all,
                        control=list(maxiter=10000,stepadj=0.00000001),
                        method="REML")

  # we only have 1 measure in each analysis, so x=1
  x=1
  # from previous version for comparison
  meta.d[x]=rma.uni(yi=Sex_d_all[x,which(!is.na(Sex_d_all[x,]))],sei=se_all[x,which(!is.na(Sex_d_all[x,]))],control=list(maxiter=10000,stepadj=0.00000001),method="REML")$b
  meta.se[x] = rma.uni(yi=Sex_d_all[x,which(!is.na(Sex_d_all[x,]))],sei=se_all[x,which(!is.na(Sex_d_all[x,]))],control=list(maxiter=10000,stepadj=0.00000001),method="REML")$se
  meta.zval[x] = rma.uni(yi=Sex_d_all[x,which(!is.na(Sex_d_all[x,]))],sei=se_all[x,which(!is.na(Sex_d_all[x,]))],control=list(maxiter=10000,stepadj=0.00000001),method="REML")$zval
  meta.pval[x] = rma.uni(yi=Sex_d_all[x,which(!is.na(Sex_d_all[x,]))],sei=se_all[x,which(!is.na(Sex_d_all[x,]))],control=list(maxiter=10000,stepadj=0.00000001),method="REML")$pval
  meta.ci.lb[x] = rma.uni(yi=Sex_d_all[x,which(!is.na(Sex_d_all[x,]))],sei=se_all[x,which(!is.na(Sex_d_all[x,]))],control=list(maxiter=10000,stepadj=0.00000001),method="REML")$ci.lb
  meta.ci.ub[x] = rma.uni(yi=Sex_d_all[x,which(!is.na(Sex_d_all[x,]))],sei=se_all[x,which(!is.na(Sex_d_all[x,]))],control=list(maxiter=10000,stepadj=0.00000001),method="REML")$ci.ub
  meta.tau2[x] = rma.uni(yi=Sex_d_all[x,which(!is.na(Sex_d_all[x,]))],sei=se_all[x,which(!is.na(Sex_d_all[x,]))],control=list(maxiter=10000,stepadj=0.00000001),method="REML")$tau2
  meta.tause[x] = rma.uni(yi=Sex_d_all[x,which(!is.na(Sex_d_all[x,]))],sei=se_all[x,which(!is.na(Sex_d_all[x,]))],control=list(maxiter=10000,stepadj=0.00000001),method="REML")$se.tau2
  meta.i2[x] = rma.uni(yi=Sex_d_all[x,which(!is.na(Sex_d_all[x,]))],sei=se_all[x,which(!is.na(Sex_d_all[x,]))],control=list(maxiter=10000,stepadj=0.00000001),method="REML")$I2
  meta.h2[x] = rma.uni(yi=Sex_d_all[x,which(!is.na(Sex_d_all[x,]))],sei=se_all[x,which(!is.na(Sex_d_all[x,]))],control=list(maxiter=10000,stepadj=0.00000001),method="REML")$H2
  meta.n.mal[x] = sum(n.mal_all[x,which(!is.na(Sex_d_all[x,]))])
  meta.n.fem[x] = sum(n.fem_all[x,which(!is.na(Sex_d_all[x,]))])
  meta.n.pat[x] = sum(n.pat_all[x,which(!is.na(Sex_d_all[x,]))])

  outmat = cbind(meta.d,meta.se,meta.zval,meta.pval,meta.ci.lb,meta.ci.ub,meta.tau2,meta.tause,meta.i2,meta.h2,meta.n.mal,meta.n.fem,meta.n.pat)

  # and save it!
  cat(paste0("Creating: MetaAnalysis_Sex_",dependent_variable,".txt\n"))

  write.table(outmat,file=paste0(outputDir,"/MetaAnalysis_Sex_",dependent_variable,".txt"),
              quote=F,sep="\t", row.names=TRUE, col.names=TRUE)

  # remove combined data before staring the next clinical measure
  rm(n.pat_all, Sex_d_all, se_all, n.mal_all, n.fem_all, models_all)

} # end dependent variable loop

# turn sink off while debugging
sink()
