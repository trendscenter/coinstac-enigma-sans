# R

# Based on SZ_CortReg from the original cortical analysis Written by TVE, LS, and DPH for the ENIGMA SZ Working Group

# clear existing workspace
rm(list=ls())

# stuff that COINSTAC needs
args = commandArgs(trailingOnly=TRUE)
baseDir = args[1]
transferDir = args[2]

# specify directory the input_dir that contains the follwing sample .csv files:
# CohortInfo.csv
# Covariates.csv
# PANSS.csv
# metr_FA.csv
# metr_AD.csv
# metr_RD.csv
# metr_MD.csv
inputdir<-paste0(getwd(),"/")
# Optional: specify input directory if the files above are not in the code directory
inputdir="../csv/FBIRN2_data_with_PANSS/"

# First get cohort information. It not present, stop.
CohortInfoFile = paste0(inputdir, "CohortInfo.csv")
if (!file.exists(CohortInfoFile)) {
  stop("CohortInfo.csv file is missing, please make sure it located in your specified inputdir:", inputdir)
}
CohortInfo <- read.csv(CohortInfoFile, header=F) # Read in the Cohort Information, this csv file should contain: "COHORT Acronym, analyst name, analyst e-mail address".
CohortName=CohortInfo[1,1]
analyst=CohortInfo[1,2]
analyst_email=CohortInfo[1,3]
# Replace possible spaces in the cohort name by underscores
cohort <- gsub(" ", "_", CohortName)
#save(CohortInfo, file="CohortInfo_output.Rdata")
#write.table(CohortInfo, file="CohortInfo_output.csv", sep=",",  row.names=FALSE, col.names=FALSE, quote = FALSE)
# Present a message that the script is running
cat(paste0("Your Cohort Name is: ", cohort," - Starting analysis...\n")) 

# Make a new directory relative to the working directory to which output is written
outputdir=paste0(getwd(),"/../results/output_sz_panss_factors_",cohort,"/")
dir.create(outputdir, recursive=TRUE, showWarnings = F)

# save CohortInfo to outputdir
save(CohortInfo, file=paste0(outputdir,"CohortInfo_output.Rdata"))
write.table(CohortInfo, file=paste0(outputdir,"CohortInfo_output.csv"), sep=",",  row.names=FALSE, col.names=FALSE, quote = FALSE)

# Create .log file and redirect all message output to this file
messages <- file(paste0(outputdir,"PANSSscriptOut_RD.txt"), open="wt")
sink(messages, type="message")
sink(messages, type="output", split=TRUE)

cat("ENIGMA-SCHIZOPHRENIA PANSS FACTORS ANALYSES LOG FILE\n")
cat("Version 1.0 April 12, 2022")
cat(paste0("Name of the cohort: ", cohort, "\n"))
cat(paste0("Name of the analyst: ", analyst, "\n"))
cat(paste0("E-mail of the analyst: ", analyst_email, "\n"))
cat(paste0(R.version.string,"\n"))
cat("Started the analysis on ", format(Sys.time(),usetz = TRUE), ".\n\n", sep = "")

# will need to source files here with external functions

source("./FileHeaders.R")
source("./Symptom_Functions.R")  
source("./RegressFunc.R")
source('./WriteRaw.R')  
library(emmeans)

# file names that have to exist this script and the input files
#Rscript = "SZ_PANSSReg_AllComps_RD_04122022.R"  # This file
dtiFile = paste0(inputdir, "metr_RD.csv")
PANSSFile = paste0(inputdir, "PANSS.csv")
CovarFile = paste0(inputdir, "Covariates.csv")

# make sure your working directory is where the .csv files are and this file is in that directory--this is not a complete check!
if (!file.exists(PANSSFile)) {
    stop("PANSS.csv file is missing  please make sure this R script is in the directory with the csv files, and your R working directory is set to that directory")
}

# read and check all files before starting on regressions etc.
RD <- read.csv(dtiFile, header = T)
Covs <- read.csv(CovarFile, header = T)
PANSS <- read.csv(PANSSFile, header = T)

# Check that all of the required columns are present in the metr_RD.csv files

cortcolind = match(DTI_ROI, names(RD))
if (length(which(is.na(cortcolind))) > 0) {
    stop(
        "At least one of the required columns in your metr_AD.csv file is missing. Make sure that the column names are spelled exactly as listed in the protocol\n"
    )
}

# Check that all of the required columns are present (Covs)

colind = match(Covcols, names(Covs))
if (length(which(is.na(colind))) > 0) {
    stop(
        "At least one of the required columns in your Covariates.csv file is missing. Make sure that the column names are spelled exactly as listed:\nIt is possible that the problem column(s) is: ",
        Covcols[which(is.na(colind))]
    )
}

# Check that all of the required columns are present (PANSS)

colind = match(PANSScols, names(PANSS))
if (length(which(is.na(colind))) > 0) {
    stop(
        "At least one of the required columns in your PANSS.csv file is missing. Make sure that the column names are spelled exactly as listed:\nIt is possible that the problem column(s) is: ",
        PANSScols[which(is.na(colind))]
    )
}

## Check for duplicated SubjIDs that may cause issues with merging data sets.
if (anyDuplicated(RD[, c("SubjID")]) != 0) {
    stop(paste0("You have duplicate SubjIDs in your RD  file.\nMake sure there are no repeat SubjIDs.\n"))
}

if (anyDuplicated(Covs[, c("SubjID")]) != 0) {
    stop("You have duplicate SubjIDs in your Covariates.csv file.\nMake sure there are no repeat SubjIDs.")
}

if (anyDuplicated(PANSS[, c("SubjID")]) != 0) {
    stop("You have duplicate SubjIDs in your PANSS.csv file.\nMake sure there are no repeat SubjIDs.")
}

# check that PANSS subjIDs are same as cases in Covs, assuming controls Dx = 0 and cases = 1
# this might be needed but I don't think it is 
#pat = which(Covs$Dx == 1)
#cases = Covs[pat, ]$SubjID
#S <- PANSS$SubjID
#if(sort(as.character(S)) != sort(as.character(cases))) { 
#cat(paste0('WARNING: PANSS and Covs. have non-matching patient SubjIDs.','\n')) 
#cat('Please make sure the patients with Covariates and PANSS are equal.','\n') }
#}

# identify the number of sites included in this dataset (if >1)
n.covs <- ncol(Covs) - 1  #Total number of covariates, -1 removes the SubjectID column
n.sites <- n.covs - 15  #Find the number of site variables, subtract the number of predictirs (Dx, Age, Sex etc.) from n.covs

# set gender and AP group to factors
Covs$AP = as.factor(Covs$AP)
Covs$Sex = as.factor(Covs$Sex)

# check gender, stop if error

n.fem=length(which(Covs$Sex==2))   # Women
n.mal=length(which(Covs$Sex==1))   # Men

#Check that Sex was coded properly
if((n.fem + n.mal) != length(Covs$Sex)){
    stop('Did you remember to code the Sex covariate as Males=1 and Females=2?\n')
}

# calculate the PANSS factorizations here and add to the PANSS matrix--external function 

PANSS = CalcPANSS(PANSS)

# Merge Covariates with RD
merged_orderedRD = merge(Covs, RD, by = "SubjID")

# Merge with PANSS
merged_orderedRD = merge(merged_orderedRD, PANSS, by = "SubjID", all.x = TRUE)

cat(paste0("The number of SubjIDs in the merged_ordered RD data set is: ", nrow(merged_orderedRD), "\n"))


# Get overall raw means for each of the structures
WriteRawFiles(merged_orderedRD,ncol(Covs),outputdir,"RawMeans_RD_asis")

# Check if sex is related to negative symptom total severity
### We want to examine gender effects for PANSS: TOT, NEG, POS, GEN, EXP, and MAP
# generate the site covariates as needed, then an ANCOVA
site = NULL  # These are just string variables that modify the model in R if there are Site variables in the Covariates.csv file
if (n.sites != 0) {
  site = " + "  #Starts the string to add-on Site variables
  for (i in 1:n.sites) {
    if (i != n.sites) {
      site = paste(site, "Site", i, " + ", sep = "")
    } else {
      site = paste(site, "Site", i, sep = "")
    }
  }
}  # end site construction


# # Put a for loop here for PANSS: TOT, NEG, POS, GEN, EXP, and MAP
cat(paste0("Total patients in PANSS ", length(PANSS[,1]), "\n"))
for (phenoName in c("TOT", "NEG", "POS", "GEN", "EXP", "MAP" )) {  # All the possible phenoNames! 

  phenoName_string = PredictCovs(phenoName)  # gets the string for the model
  #phenoName_string = phenoName

  cat(paste0("Running: Sex Differences in ", phenoName, "\n"))

  # make the formula and run it
  sexform = as.formula(paste0(phenoName_string, " ~ factor(Sex) + Age ", site))

  genderMod = lm(sexform, data=merged_orderedRD)

  # save the model, the number of men/women included
  # "d.cort"      "low.ci.cort" "n.controls"  "n.patients"  "se.cort"     "up.ci.cort" 
  # TVE: we should probably added the estimated means / least square means to the table also

  # pull marginal means
  em<-emmeans(genderMod, "Sex")
  emd<-data.frame(em)
  m.adjusted.mean<-emd[1,2]
  f.adjusted.mean<-emd[2,2]
  m.se<-emd[1,3]
  f.se<-emd[2,3]

  n.patients = nobs(genderMod)
  n.mal = length(which(genderMod$model$`factor(Sex)`==1))
  n.fem = n.patients - n.mal
  # Gender Effect
  Sex_eff = coef(summary(genderMod))[2]  # regression coeff
  tvalue=coef(summary(genderMod))[2, 't value']
  pvalue=coef(summary(genderMod))[2, 'Pr(>|t|)']
  cohen_d = d.t.unpaired(tvalue, n.mal, n.fem)
  se = se.d2(tvalue, n.mal, n.fem)
  CI=CI1(ES = cohen_d,se = se)
  low.ci=CI[1]
  up.ci=CI[2]
  # Age Effect
  # get the t, turn it into an r (slightly biased if not normally distributed )
  Age_eff = coef(summary(genderMod))[3]  # regression coeff
  tvalue_age = coef(summary(genderMod))[3,'t value']
  pvalue_age = coef(summary(genderMod))[3,'Pr(>|t|)']
  df = df.residual(genderMod)
  r = tvalue_age/(sqrt(tvalue_age*tvalue_age + df))

  # list the number of men and women for phenoName
  cat(paste0("Number of Men: ", n.mal, "\n"))
  cat(paste0("Number of Women: ", n.fem, "\n"))

  # Exclude components of the lm() object that from which individual-level data can be recovered
  genderMod[["fitted.values"]] <- NULL
  genderMod[["residuals"]] <- NULL
  genderMod[["model"]] <- NULL

  # save Gender effects
  save(genderMod, sexform, emd, m.adjusted.mean, f.adjusted.mean, m.se, f.se, n.patients, n.mal, n.fem, Sex_eff, tvalue, pvalue,
       cohen_d, se, low.ci, up.ci, file = paste0(outputdir,"EffectSizes_SZ_only_Gender_withAge_on_",
                                                 phenoName, ".Rdata"))

  # save Age effects
  tvalue=tvalue_age
  pvalue=pvalue_age
  save(genderMod, sexform, emd, m.adjusted.mean, f.adjusted.mean, m.se, f.se, n.patients, n.mal, n.fem, Age_eff, tvalue, pvalue,
       r, file = paste0(outputdir,"EffectSizes_SZ_only_Age_withSex_on_",
                                                 phenoName, ".Rdata"))

  }

# cross correlation of all PANSS data and factors here, save to file!
CrossCorr = cor(PANSS[,2:dim(PANSS)[2]], use="pairwise.complete.obs")
save(CrossCorr, file=paste0(outputdir,"PANSSCorrs.Rdata"))

cc = "asis"  # for the moment! No complete cases only

cat("Starting pheno loop\n")

for (phenoName in c("RD")) {  # brain measure type loop
  
    # extract data subsets using the phenotypic-specific dataset from above
    # Use the total PANSS we created  
    
    dataset = paste0("merged_ordered", phenoName)  # get the right one: RD
    merged_ordered0 = get(dataset)
    
    merged_ordered_PANSSSum <- merged_ordered0[((merged_ordered0$Dx == 1) & !(is.na(merged_ordered0$TOT))), ]
    
    # Check nsub with data
    L = merged_ordered_PANSSSum[!is.na(merged_ordered_PANSSSum[,18]),18]
    #L = merged_ordered_PANSSSum[!is.na(merged_ordered_PANSSSum$TOT)]
    cat(paste0("Total subjects in PANSS: ", length(L), "\n"))
    
    
    # code for the globalMeasure to include
    # GlobalMeasure= switch(phenoName, "Cort" = "MThickness", "Surf" = "FullSurfArea", "SubCort" = "ICV")
    # The average for each of the 4 DTI measure (FA, AD, RD, and MD) was always named AverageFA
    GlobalMeasure= switch(phenoName, "RD" = "AverageFA")
    
    # analysis for the PANSS: TOT, NEG, POS, GEN, EXP, and MAP
    # analysis for EXP covarying for TOT, GEN, and NEG
    # analsyis for MAP covarying for TOT, GEN, and NEG
    
    
    # We have 2 duplicate analyses: EXPwTOT and MAPwTOT results should be the same as those in the EXP and MAP with WSum as a Covariate files
    for (predictor in c("TOT","NEG","POS","GEN","EXP","MAP",
                        "EXPwTOT", "EXPwGEN", "EXPwNEG",
                        "MAPwTOT", "MAPwGEN", "MAPwNEG" )) {  # All the possible predictors!
            
        #predictor_string = PredictCovs(predictor)  # gets the string for the model
    
        predictor_string = PredictCovs(predictor)

        # This code is outdated, NumPredict is always 1
        #NumPredict = switch (predictor,
        #    "SANSGlobals" = 4,
        #    "SANSFactors" = 5,
        #    1
        #)  # this sets up how many factors' effect sizes to save later
        NumPredict=1
        
        # load data for this pheno 
        merged_ordered = merged_ordered_PANSSSum
        
        #cat(paste0("dataset and pheno & predictor & global off/on ", cc," ", phenoName, " ", predictor_string, " ", GlobalMeasure, "\n"))
        cat(paste0("dataset and pheno & predictor & global off/on ", cc," ", phenoName, " ", predictor_string, " \n"))
        
        n.patients = NULL
        
        site = NULL  # These are just string variables that modify the model in R if there are Site variables in the Covariates.csv file
        if (n.sites != 0) {
            site = " + "  #Starts the string to add-on Site variables
            for (i in 1:n.sites) {
                if (i != n.sites) {
                    site = paste(site, "Site", i, " + ", sep = "")
                } else {
                    site = paste(site, "Site", i, sep = "")
                }
            }
        }  # end site construction
        
        
        if (nrow(merged_ordered) > (n.sites + 5)) {  # do or not do, there is no try (if enough subjects, run)
            # this was 3 but it should be 5
            
            # Create list of covariates with available data - more than 20 subjects
            WCovs<-c("NoCovs","WG","WSum") # initialize array WCov with once that are definitely present: NoCovs and WG 

            # Need to make sure the length checks patients' only data
            # length(Covs$IQ[Covs$Dx == 1 & !is.na(Covs$IQ)])
            for ( Cov in c("IQ", "CPZ", "AO", "AP") ) {
               if ( (Cov == "IQ") && (length(Covs$IQ[Covs$Dx == 1 & !is.na(Covs$IQ)]) ) > 19 ) {
                 Cov_tmp<-paste0("W",Cov)
                 WCovs<-c(WCovs,Cov_tmp)
                 #cat(WCovs,"\n")
                 }
               else if( (Cov == "CPZ") && (length(Covs$CPZ[Covs$Dx == 1 & !is.na(Covs$CPZ)]) ) > 19 ) {
                 Cov_tmp<-paste0("W",Cov)
                 WCovs<-c(WCovs,Cov_tmp)
                 #cat(WCovs,"\n")
               }
               else if ( (Cov == "AO") && (length(Covs$AO[Covs$Dx == 1 & !is.na(Covs$AO)]) ) > 19 ) {
                 Cov_tmp=paste0("W",Cov)
                 WCovs=c(WCovs,Cov_tmp)
                 #cat(WCovs,"\n")
               }
               else if ( (Cov == "AP") && (length(Covs$AP[Covs$Dx == 1 & !is.na(Covs$AP)]) ) > 19 ) {
                 Cov_tmp<-paste0("W",Cov)
                 WCovs<-c(WCovs,Cov_tmp)
                 #cat(WCovs,"\n")
               }
               else {
                 cat(paste0("Less than 20 subjects with Covariate: ", Cov, "\n"))
               }
               #cat(WCovs,"\n")
             }

            #for ( WCov in c("NoCovs","WG", "WSum", "WIQ", "WCPZ", "WAO", "WAP")) {  # loop for other covariates: Global brain, IQ, CPZ equiv, AO, AP group
            for ( WCov in WCovs ) {
            cat(paste0("Running: Regression predictor ", predictor, " against ", phenoName, ", ", cc, ", in SZ patients covary for ", WCov, ", Age and Sex \n"))
            
            # set up model variables that don't change in this loop
                if(WCov == "WG") {
                    glob=paste0("+", GlobalMeasure, sep="")
                } else if (WCov == "WIQ") { 
                    glob="+ IQ "
                } else if (WCov == "WCPZ") {
                    glob="+ CPZ "
                } else if (WCov == "WAO") {
                    glob = "+ AO "
                } else if (WCov == "WAP") {
                    glob = "+ AP "
                } else if (WCov == "WSum") {
                    glob = "+ TOT"
                } else {  # default -- no covariates
                    glob=""
                }
                
                
                # number of brain measures is the total columns subtracting out Covs and PANSS and global means
                # Npheno = ncol(merged_ordered) - ncol(Covs) - ncol(PANSS) - 2
                # Why was -2 removed here?
                Npheno = ncol(merged_ordered) - ncol(Covs) - ncol(PANSS) # - 2
                
                #if(phenoName == "SubCort")
                #    Npheno = Npheno+1  # only one global covariate in the file, not two
                
                Startcol = ncol(Covs) + 1
                Endcol = Startcol + Npheno # this should vary between cortical and subcortical
                 
                
                # Store models for troubleshooting
                models.cort = NULL  # This will become a list where we store all of the models made by lm
                
                # allocate empty vectors to store adjust effect sizes, se, ci (noicv)
                r.cort = matrix(NA, nrow= Npheno+1, ncol = NumPredict)
                n.patients = matrix(NA, nrow = Npheno+1, ncol=1)
                
                # se.cort = NULL
                # low.ci.cort = NULL
                # up.ci.cort = NULL
                
                # Loop through and perform each regression
               
            for (x in (Startcol:Endcol)) {  
                # do the model for each brain measure (column)
                
                pheno = merged_ordered[!is.na(merged_ordered[, x]), x]  #Check to make sure there are observations for a given structure
                
                # check if the phenotype is singular after NA removal
                if (length(pheno) == 0) {
                  next  # Skip to the next measure if there are no observations
                }
                
                # Run the modelfor each brain measure x
                DepVar = merged_ordered[, x]
                
                thisformula = as.formula(paste("DepVar ~ ", predictor_string, " + Age + Sex", glob, site, sep = ""))
                
                tmp=lm(thisformula, data=merged_ordered)
         
                models.cort[[x - ncol(Covs)]] = tmp  #Store the model fit for future reference
                
                # subjects can be dropped if they are missing so we can get the precise number of controls/patients for each region tested
                # n.controls[x - ncol(Covs)] = length(which(tmp$model[, 2] == 0))
                # n.patients[x - ncol(Covs)] = length(which(tmp$model[, 2] == 1))
                
                n.controls = 0 # none used in this analysis
                n.patients[x-ncol(Covs)] = length(tmp$model[,2]) # this works!
                
                
                # extract effect sizes for the predictors
                # partcor.i <- pcor.test(tmp$model[1],tmp$model[,2],tmp$model[,c(3:ncol(tmp$model))])
                
                # cohens = cohens_f(tmp) # from sjstats
                
                # get the t, turn it into an r (slightly biased if not normally distributed )
                tvalue = coef(summary(tmp))[2:(NumPredict+1),'t value']
                df = df.residual(tmp)  
                r = tvalue/(sqrt(tvalue*tvalue + df))
                
                r.cort[x-ncol(Covs),]= r # this should save the effects for the first predictor

                # Exclude components of the lm() object that from which individual-level data can be recovered
                # model
                models.cort[[x - ncol(Covs)]]["model"] <- NULL
                # fitted.values
                models.cort[[x - ncol(Covs)]]["fitted.values"] <- NULL
                # residuals
                models.cort[[x - ncol(Covs)]]["residuals"] <- NULL                
                
            } # end for loop on x (brain regions)
                # save results CHECK THIS
                save(r.cort, n.patients, file = paste0(outputdir, "EffectSizes_SZ_only_", predictor, "_withSexAge_", WCov, "_", phenoName, ".Rdata"))
                
                save(models.cort, file = paste0(outputdir,"Models_SZ_", predictor, "_withSexAge_", WCov,"_", phenoName, ".Rdata"))
        } # end With Global loop   
            
        
            
        } #end do or not do based on sample size
            
        else  # not running that regression
            {
            cat("NOT Running: Regression predictor ", predictor, ",", cc, glob, " in SZ patients covary for age + Sex \n")
        }  

    }  #  end predictor loop, Regressions & AP models
    
} # end phenoname loop for debugging
        #### TVE ##########
        
cat("Done working \n")
sink()  # close the log file

#### TVE/JAT ##########
