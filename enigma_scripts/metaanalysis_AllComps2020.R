#R
options(warn=-1)

args = commandArgs(trailingOnly=TRUE)
baseDir = args[1]
outputDir = args[2]
NumDir = strtoi(args[3])


#need the metafor library
suppressMessages(library(metafor))

pdf(paste(outputDir, 'SANSoutputs.pdf', sep='/'))

# do we need a sink file here too? 
sink(paste(outputDir, "Meta.txt", sep='/'), split=TRUE)

# for every site, for every analysis
# Read in effect size tables for each group

# for each phenotype (thickness, surf, subcort) 
# for each SANS measure (SANSSum, each Global and all 5 Global measures, each factor and all 5 Factors) 
# run Age + Sex + SANS measure
# run Age + Sex + SANS measure + (mean thick/area/ICV) 
# run Age + Sex + SANS measure + CPZ/IQ/AO

# do another loop for effects of sex and symptoms (no brain!)

for (phenoName in c("Cort", "Surf", "SubCort")) {  # brain measure type loop
  for (predictor in c("SANSSum", "SANSMAP", "SANSEXP",
                      "SANSGlobal1", "SANSGlobal2", "SANSGlobal3", "SANSGlobal4", "SANSGlobal5", 
                      "SANSFac1", "SANSFac2", "SANSFac3", "SANSFac4", "SANSFac5")) {  # All the possible predictors!
    for ( WCov in c("NoCovs","WG", "WSum", "WIQ", "WCPZ", "WAO", "WAP")) {  # loop for other covariates: Global brain, IQ, CPZ equiv, AO, AP group
      
      # generate the input file names
      
      Effectsf= paste0("EffectSizes_SZ_only_",predictor,"_withSexAge_",WCov,phenoName,".Rdata") # this depends on the previous output
      Modelf = paste0("Models_SZ_",predictor,"_withSexAge_",WCov,"_", phenoName,".Rdata")
      
     # output file
      Outfile=paste0(outputDir, '/', "MetaAnalysis_SZ_",predictor,"_withSexAge_",WCov,phenoName,".txt")
      
      nsites = 0 # how many sites, for each analysis
      
      for (site in 0:(NumDir-1)) {
        
        sitedir=paste0("local",site)
        
        if (!dir.exists(paste0(baseDir, '/', sitedir))) {
          # site failed to respond
          cat(sitedir,"does not exist \n")
          
        } else { # site exists, keep going
          
          #setwd(paste0(baseDir, '/', sitedir))  # this may not be the way to do it? 
          
          # open the files and read them in for the effect sizes and the # subjects included
          # if successful, add to the count of nsites for that analysis
          
          if(!file.exists(paste0(baseDir, '/', sitedir, '/', Effectsf)) || !file.exists(paste0(baseDir, '/', sitedir, '/', Modelf))) {
            cat(sitedir, "does not have the needed files for", Effectsf)
          } else {
            # read the files--puts n.patients, models.cort, and r.cort into memory
            load(paste0(baseDir, '/', sitedir, '/', Effectsf))
            load(paste0(baseDir, '/', sitedir, '/', Modelf))
            
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
          
#          setwd("..")  # back to the original directory
          
        } # end if site exists loop
        
      } # end  site loop
      
      
      if(nsites > 0) {# do the meta analysis and save the results
        
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
        
        #fdr correction and write to file
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
        outmat=cbind(meta.noicv.r,meta.noicv.se,meta.noicv.zval,meta.noicv.pval,meta.noicv.ci.lb,meta.noicv.ci.ub,meta.noicv.tau2,meta.noicv.tause,meta.noicv.i2,meta.noicv.h2,meta.noicv.npat,fdr.p,fdr.p.mean)
        
        
        write.table(outmat,file=Outfile, quote=F,sep="\t", row.names=TRUE, col.names=TRUE)
        
        
      } else {
        cat("no sites for this analysis! \n")
      }
      
      # clear npat_all, r_eff, and models_all for the enxt analysis
      rm(list=c("npat_all", "r_eff","models_all"))
      
    } # end WCov loop
  } # end SANS predictor loop
} # end pheno loop

sink()
invisible(dev.off()) # close pdf
