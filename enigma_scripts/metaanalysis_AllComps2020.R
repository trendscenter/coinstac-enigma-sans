#R

# clear environment
rm(list=ls())

args = commandArgs(trailingOnly=TRUE)
baseDir = args[1]
outputDir = args[2]
siteList = as.list(strsplit(args[3], ",")[[1]])
#curDir="./"
curDir=getwd()

# get the list of input directories for meta-analysis - these are the output directories from each of the sites
inputDirs<-list.files(path=".", pattern="^output")
# get the number of input directories for the meta-analysis
# NumDir<-NROW(inputDirs)

# create the output directory for the meta-anlaysis results
# outputDir="./results_meta_analysis"  # for the moment!

if (file.exists(outputDir)){
  #setwd(file.path(curDir, outputDir))
} else {
  #dir.create(file.path(curDir, outputDir))
  #setwd(file.path(curDir, outputDir))
}

# now defined above
#NumDir=2  # this should one of the args

# create an analysis log file
sink(paste(outputDir,"MetaAnalysisLog.txt", sep='/'), split=TRUE)

# load the metafor library
library(metafor)

#
# for every site, for every analysis
# Read in effect size tables for each group

# for each phenotype (thickness, surf, subcort)
# for each SANS measure (SANSSum, each Global , each factor etc.)
# run Age + Sex + SANS measure
# run Age + Sex + SANS measure + (mean thick/area/ICV)
# run Age + Sex + SANS measure + CPZ/IQ/AO/AP

# do another loop for effects of sex and symptoms (no brain!)

#############################################
### Meta-analysis of symptom associations ###
#############################################

for (phenoName in c("Cort", "Surf", "SubCort")) {  # brain measure type loop
   for (predictor in c("SANSSum", "SANSMAP", "SANSEXP",
                       "SANSGlobal1", "SANSGlobal2", "SANSGlobal3", "SANSGlobal4",
                       "SANSFac1", "SANSFac2", "SANSFac3", "SANSFac4", "SANSFac5",
                       "SANSMAPwSANSSum", "SANSEXPwSANSSum",
                       "SANSGlobal1wSANSSum", "SANSGlobal2wSANSSum", "SANSGlobal3wSANSSum", "SANSGlobal4wSANSSum",
                       "SANSFac1wSANSSum", "SANSFac2wSANSSum", "SANSFac3wSANSSum", "SANSFac4wSANSSum", "SANSFac5wSANSSum",
                       "SANSGlobal1wSANSEXP", "SANSGlobal2wSANSEXP", "SANSGlobal3wSANSEXP", "SANSGlobal4wSANSEXP",
                       "SANSFac1wSANSEXP", "SANSFac2wSANSEXP", "SANSFac3wSANSEXP", "SANSFac4wSANSEXP", "SANSFac5wSANSEXP",
                       "SANSGlobal1wSANSMAP", "SANSGlobal2wSANSMAP", "SANSGlobal3wSANSMAP", "SANSGlobal4wSANSMAP",
                       "SANSFac1wSANSMAP", "SANSFac2wSANSMAP", "SANSFac3wSANSMAP", "SANSFac4wSANSMAP", "SANSFac5wSANSMAP")) {  # All the possible predictors!
        for ( WCov in c("NoCovs","WG", "WSum", "WIQ", "WCPZ", "WAO", "WAP")) {  # loop for other covariates: Global brain, IQ, CPZ equiv, AO, AP group

      # generate the input file names
      #cat(paste(predictor, " ", WCov,"\n"))

      Effectsf= paste0("EffectSizes_SZ_only_",predictor,"_withSexAge_",WCov,phenoName,".Rdata") # this depends on the previous output
      Modelf = paste0("Models_SZ_",predictor,"_withSexAge_",WCov,"_", phenoName,".Rdata")

      #cat(paste(Effectsf,"\n"))
      #cat(paste(Modelf, "\n"))

     # output file
      Outfile=paste0("/MetaAnalysis_SZ_",predictor,"_withSexAge_",WCov,phenoName,".txt")

      cat(paste("Creating: ", Outfile,"\n"))

      nsites = 0 # how many sites, for each analysis

      for (siteName in siteList) {
        siteDir=paste(baseDir, siteName, sep='/')

        if (!dir.exists(siteDir)) {
          # site failed to respond
          cat(siteDir,"does not exist \n")
          cat(paste(siteDir,"\n"))
          }
        else { # site exists, keep going
          # open the files and read them in for the effect sizes and the # subjects included
          # if successful, add to the count of nsites for that analysis
          if(!file.exists(paste(siteDir, Effectsf, sep='/')) || !file.exists(paste(siteDir, Modelf, sep='/'))) {
            cat(siteName, "does not have", Effectsf,"\n")
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
           } else {
             npat_all=cbind(npat_all, n.patients)
             r_eff=cbind(r_eff,r.cort)
             models_all = cbind(models_all, models.cort)
           }

            nsites=nsites+1

          }

        } # end if site exists loop

      } # end  site loop


      if(nsites > 0) { # do the meta analysis and save the results

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

        } # end looping through regions

        # fdr correction and write to file
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
        outmat=cbind(meta.noicv.r,meta.noicv.se,meta.noicv.zval,meta.noicv.pval,
                     meta.noicv.ci.lb,meta.noicv.ci.ub,
                     meta.noicv.tau2,meta.noicv.tause,
                     meta.noicv.i2,meta.noicv.h2,meta.noicv.npat,fdr.p,fdr.p.mean)

        write.table(outmat,file=paste(outputDir, Outfile, sep='/'),
                    quote=F,sep="\t", row.names=TRUE, col.names=TRUE)

      } else {
        cat("no sites for this analysis! \n")
      }

      # clear npat_all, r_eff, and models_all for the next analysis
      rm(list=c("npat_all", "r_eff","models_all"))

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
  for (siteName in siteList) { # read in sex effect on dependent variable for each site

    siteDir=paste(baseDir, siteName, sep='/')
    ModelF=paste0("EffectSizes_SZ_only_Gender_withAge_on_", dependent_variable,".Rdata")
    cat(paste0(ModelF," from ", siteName, "\n"))

    if(!file.exists(paste(siteDir, ModelF, sep='/'))) {
      cat(siteName, "does not have the needed files for", ModelF, siteDir)
      }
    else { # read the files--puts the npat and effect size into memory
      load(paste(siteDir, ModelF, sep='/'))  # puts genderMod model in memory
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
        n.pat_all = cbind(n.pat_all, n.patients)
        Sex_d_all = cbind(Sex_d_all, cohen_d)
        se_all = cbind(se_all, se)
        n.mal_all = cbind(n.mal_all, n.mal)
        n.fem_all = cbind(n.fem_all, n.fem)
        models_all = cbind(models_all, genderMod)
        } # end what to do in a site
      } # done with that site
      nsites=nsites+1
      #setwd("..") # back to base
    # } # end file.exists
  } # end for loop siteNames - to load data from all sites

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
  meta.d[x] = rma.uni(yi=Sex_d_all[x,which(!is.na(Sex_d_all[x,]))],sei=se_all[x,which(!is.na(Sex_d_all[x,]))],control=list(maxiter=10000,stepadj=0.00000001),method="REML")$b
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
  cat(paste0("Creating: MetaAnalysis_Sex_", dependent_variable,".txt\n"))

  write.table(outmat, file=paste0(outputDir, "/MetaAnalysis_Sex_", dependent_variable, ".txt"),
              quote=F, sep="\t", row.names=TRUE, col.names=TRUE)

  # remove combined data before staring the next clinical measure
  rm(n.pat_all, Sex_d_all, se_all, n.mal_all, n.fem_all, models_all)

} # end dependent variable loop

# turn sink off while debugging
sink()
