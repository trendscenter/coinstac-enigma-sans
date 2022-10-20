# R

# Based on SZ_CortReg from the original cortical analysis Written by TVE, LS, and DPH for the ENIGMA SZ Working Group

# clear existing workspace
rm(list=ls())

# stuff that COINSTAC needs
args = commandArgs(trailingOnly=TRUE)
baseDir = args[1]
transferDir = args[2]

# First get cohort information. It not present, stop.
CortThickFile = file.path(baseDir, args[3])
SurfFile = file.path(baseDir, args[4])
SubCortFile = file.path(baseDir, args[5])
SANSFile = file.path(baseDir, args[6])
CovarFile = file.path(baseDir, args[7])
CohortInfoFile = file.path(baseDir, args[8])

if (!file.exists(CohortInfoFile)) {
  print(CohortInfoFile)
  message("CohortInfo.csv file is missing  please make sure this R script is in the directory with the csv files, and your R working directory is set to that directory")
  quit(status=1)
}

CohortInfo <- read.csv(CohortInfoFile, header=F) # Read in the Cohort Information, this csv file should contain: "COHORT Acronym, analyst name, analyst e-mail address".
CohortName=CohortInfo[1,1]
analyst=CohortInfo[1,2]
analyst_email=CohortInfo[1,3]

# Replace possible spaces in the cohort name by underscores
cohort <- gsub(" ", "_", CohortName)

# Present a message that the script is running
cat(paste0("Your Cohort Name is: ", cohort," - Starting analysis...\n"))

# Make a new directory inside the working directory to which output is written
outputdir=paste0(transferDir, "/output_sz_sans_factors_", cohort, "/")
dir.create(outputdir, showWarnings = F)

# save CohortInfo to outputdir
save(CohortInfo, file=paste(outputdir, "CohortInfo_output.Rdata", sep = "/"))
write.table(CohortInfo, file=paste(outputdir, "CohortInfo_output.csv", sep = "/"), sep=",",  row.names=FALSE, col.names=FALSE, quote = FALSE)

# this still needs a log file
#sink(paste0(outputdir,"SANSscriptOut.txt"), split=TRUE)

# Create .log file and redirect all message output to this file
cat("ENIGMA-SCHIZOPHRENIA SANS FACTORS ANALYSES LOG FILE\n")
cat("Test Phase Feb 2020\n")
cat(paste0("Name of the cohort: ", cohort, "\n"))
cat(paste0("Name of the analyst: ", analyst, "\n"))
cat(paste0("E-mail of the analyst: ", analyst_email, "\n"))
cat(paste0(R.version.string,"\n"))
cat("Started the analysis on ", format(Sys.time(),usetz = TRUE), ".\n\n", sep = "")

# will need to source files here with external functions

source("./enigma_scripts/FileHeaders.R")
source("./enigma_scripts/SANSCalcs.R")
source("./enigma_scripts/RegressFunc.R")
source('./enigma_scripts/WriteRaw.R')
library(emmeans)

# file names that have to exist this script and the input files



# make sure your working directory is where the .csv files are and this file is in that directory--this is not a complete check!
if (!file.exists(SANSFile)) {
    message("SANS.csv file is missing  please make sure this R script is in the directory with the csv files, and your R working directory is set to that directory")
    quit(status=1)
}

# read and check all files before starting on regressions etc.

Cort <- read.csv(CortThickFile, header = T)  #Read in the phenotypes file
Surf = read.csv(SurfFile, header = T)
SubCort <- read.csv(SubCortFile, header = T)
Covs = read.csv(CovarFile, header = T)
SANS = read.csv(SANSFile, header = T)



cortcolind = match(CortCols, names(Cort))
if (length(which(is.na(cortcolind))) > 0) {
    message(
        "At least one of the required columns in your Cortical thickness measures file is missing. Make sure that the column names are spelled exactly as listed in the protocol\n"
    )
    quit(status=1)
}

cortcolind = match(SurfCols, names(Surf))
if (length(which(is.na(cortcolind))) > 0) {
    message(
        "At least one of the required columns in your Cortical surfare area measures file is missing. Make sure that the column names are spelled exactly as listed in the protocol\n"
    )
    quit(status=1)
}

cortcolind = match(SubCortCols, names(SubCort))
if (length(which(is.na(cortcolind))) > 0) {
    message(
        "At least one of the required columns in your LandRvolumesn.csv file is missing. Make sure that the column names are spelled exactly as listed in the protocol\n"
    )
    quit(status=1)
}

# Check that all of the required columns are present (Covs)

colind = match(Covcols, names(Covs))
if (length(which(is.na(colind))) > 0) {
    message(
        "At least one of the required columns in your Covariates.csv file is missing. Make sure that the column names are spelled exactly as listed:\nIt is possible that the problem column(s) is: ",
        Covcols[which(is.na(colind))]
    )
    quit(status=1)
}

# Check that all of the required columns are present (SANS)

colind = match(SANScols, names(SANS))
if (length(which(is.na(colind))) > 0) {
    message(
        "At least one of the required columns in your SANS.csv file is missing. Make sure that the column names are spelled exactly as listed:\nIt is possible that the problem column(s) is: ",
        SANScols[which(is.na(colind))]
    )
    quit(status=1)
}

# Check for duplicated SubjIDs that may cause issues with merging data sets.
if (anyDuplicated(Cort[, c("SubjID")]) != 0) {
    message(paste0("You have duplicate SubjIDs in your cortical measures file.\nMake sure there are no repeat SubjIDs.\n"))
    quit(status=1)
}
if (anyDuplicated(Surf[, c("SubjID")]) != 0) {
    message(paste0("You have duplicate SubjIDs in your cortical measures file.\nMake sure there are no repeat SubjIDs.\n"))
    quit(status=1)
}

if (anyDuplicated(SubCort[, c("SubjID")]) != 0) {
    message(paste0("You have duplicate SubjIDs in your subcortical volumes file.\nMake sure there are no repeat SubjIDs.\n"))
    quit(status=1)
}

if (anyDuplicated(Covs[, c("SubjID")]) != 0) {
    message("You have duplicate SubjIDs in your Covariates.csv file.\nMake sure there are no repeat SubjIDs.")
    quit(status=1)
}

if (anyDuplicated(SANS[, c("SubjID")]) != 0) {
    message("You have duplicate SubjIDs in your SANS.csv file.\nMake sure there are no repeat SubjIDs.")
    quit(status=1)
}

# check that SANS subjIDs are same as cases in Covs, assuming controls Dx = 0 and cases = 1
# this might be needed but I don't think it is
#pat = which(Covs$Dx == 1)
#cases = Covs[pat, ]$SubjID
#S <- SANS$SubjID
#if(sort(as.character(S)) != sort(as.character(cases))) {
#cat(paste0('WARNING: SANS and Covs. have non-matching patient SubjIDs.','\n'))
#cat('Please make sure the patients with Covariates and SANS are equal.','\n') }
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
    message('Did you remember to code the Sex covariate as Males=1 and Females=2?\n')
    quit(status=1)
}

# calculate the SANS factorizations here and add to the SANS matrix--external function

SANS = CalcSans(SANS)

#### calculate mean values --This could be looped...

meanCort = NULL
meanSurf = NULL

# calculate means across hemispheres
for (x in 2:35) {  # This depends on the FS file format and cortical ROI labels

    meanCort = c(meanCort, ((Cort[, x] + Cort[x + 34])/2))
    meanSurf = c(meanSurf, ((Surf[, x] + Surf[x + 34])/2))
}

meanCort = c(meanCort, ((Cort[, 70] + Cort[71])/2))
meanCort = c(meanCort, ((Cort[, 72] + Cort[73])))
meanSurf = c(meanSurf, ((Surf[, 70] + Surf[71])/2))
meanSurf = c(meanSurf, ((Surf[, 72] + Surf[73])))

for (x in 1:34) {
    tmp = strsplit(names(meanCort)[x], "_")
    names(meanCort)[x] = paste0("M_", tmp[[1]][2], "_", tmp[[1]][3])
    tmp2 = strsplit(names(meanSurf)[x], "_")
    names(meanSurf)[x] = paste0("M_", tmp2[[1]][2], "_", tmp2[[1]][3])
}

names(meanCort)[35] = "MThickness"
names(meanCort)[36] = "FullSurfArea"
meanCort = as.data.frame(meanCort)

names(meanSurf)[35] = "MThickness"  # this is not what it was before
names(meanSurf)[36] = "FullSurfArea"
meanSurf = as.data.frame(meanSurf)

# drop ICV from Cort and Surf file
Cort = Cort[, -ncol(Cort)]
Surf = Surf[, -ncol(Surf)]

# combine Cort and Surf file with the newly calculated means
Cort = cbind(Cort, meanCort)
Surf = cbind(Surf, meanSurf)

# now do the same for the subcortical: means and combine, but don't drop the ICV
SubCort$Mvent <- rowMeans(SubCort[, c("LLatVent", "RLatVent")])  #calculate mean Ventricle
SubCort$Mthal <- rowMeans(SubCort[, c("Lthal", "Rthal")])  #calculate mean Thalamus
SubCort$Mcaud <- rowMeans(SubCort[, c("Lcaud", "Rcaud")])  #calculate mean Caudate
SubCort$Mput <- rowMeans(SubCort[, c("Lput", "Rput")])  #calculate mean Putamen
SubCort$Mpal <- rowMeans(SubCort[, c("Lpal", "Rpal")])  #calculate mean Pallidum
SubCort$Mhippo <- rowMeans(SubCort[, c("Lhippo", "Rhippo")])  #calculate mean Hippocampus
SubCort$Mamyg <- rowMeans(SubCort[, c("Lamyg", "Ramyg")])  #calculate mean Amygdala
SubCort$Maccumb <- rowMeans(SubCort[, c("Laccumb", "Raccumb")])  #calculate mean Accumbens

#move ICV to the end column--important for later analyses!
SubCort = SubCort[,c(1:17,19:26,18)]

# combine the files into one dataframe per phenotype
merged_orderedCort = merge(Covs, Cort, by = "SubjID")
merged_orderedSurf = merge(Covs, Surf, by = "SubjID")
merged_orderedSubCort = merge(Covs, SubCort, by = "SubjID")

# Check that the number of rows for brains and covars after merging is the same
#I'm not sure this is necessary, since we need all with SANS and no others
# if (nrow(Cort) != nrow(merged_ordered)) {
    # # cat(paste0('WARNING: ', fsfile, ' and Covariates.csv have non-matching SubjIDs.','\n')) cat('Please make sure the number of
    # subjects in your merged data set are as expected.','\n') cat(paste0('The number of SubjIDs in ', fsfile, ' is: ',nrow(Cort),'\n'))
    # cat("The number of SubjIDs in the merged_ordered Cortical Thickness data set is: ", nrow(merged_orderedCort), "\n")
#}

# merge in the SANS data at the end but don't throw out the controls until we know we aren't doing case/control analyses

merged_orderedCort = merge(merged_orderedCort, SANS, by = "SubjID", all.x = TRUE)
merged_orderedSurf = merge(merged_orderedSurf, SANS, by = "SubjID", all.x = TRUE)
merged_orderedSubCort = merge(merged_orderedSubCort, SANS, by = "SubjID", all.x = TRUE)

cat(paste0("The number of SubjIDs in the merged_ordered Cortical Thickness data set is: ", nrow(merged_orderedCort), "\n"))


# Get overall raw means for each of the structures

WriteRawFiles(merged_orderedCort,ncol(Covs),outputdir,"RawMeans_Cort_asis")
WriteRawFiles(merged_orderedSurf,ncol(Covs),outputdir,"RawMeans_Surf_asis")
WriteRawFiles(merged_orderedSubCort,ncol(Covs),outputdir,"RawMeans_SubCort_asis")

# Check if sex is related to negative symptom total severity
### We want to examine gender effects for SANSSum, SANS$GlobalTotal, Factors 1-5, and MAP and EXP
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

# Put a for loop here for SANSSum "SANSSum" "MAP" "EXP" "AnhedoniaFac" "AsocialityFac" "AvolitionFac" "BluntedAffectFac" "AlogiaFac"
cat(paste0("Total patients in SANS ", length(SANS[,1]), "\n"))
for (phenoName in c("SANSSum", "SANSMAP", "SANSEXP",
                    "SANSGlobal1", "SANSGlobal2", "SANSGlobal3", "SANSGlobal4",
                    "SANSFac1", "SANSFac2", "SANSFac3", "SANSFac4", "SANSFac5")) {  # All the possible phenoNames!

  phenoName_string = PredictCovs(phenoName)  # gets the string for the model

  cat(paste0("Running: Sex Differences in ", phenoName, "\n"))

  # make the formula and run it
  sexform = as.formula(paste0(phenoName_string, " ~ factor(Sex) + Age ", site))

  genderMod = lm(sexform, data=merged_orderedSubCort)

  # save the model, the number of men/women included
  #"d.cort"      "low.ci.cort" "n.controls"  "n.patients"  "se.cort"     "up.ci.cort"
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

# cross correlation of all SANS data and factors here, save to file!
CrossCorr = cor(SANS[,2:dim(SANS)[2]], use="pairwise.complete.obs")
save(CrossCorr, file=paste0(outputdir,"SANSCorrs.Rdata"))

### Regression Analyses create file for resetting to original file after modifications per each regression

swap.merged_orderedCort = merged_orderedCort
swap.merged_orderedSurf = merged_orderedSurf
swap.merged_orderedSubCort = merged_orderedSubCort

cc = "asis"  # for the moment! No complete cases only

# Check nsub with data
L = swap.merged_orderedCort[!is.na(swap.merged_orderedCort[,18]),18]
cat(paste0("Total subjects in Cortical ", length(L), "\n"))

# for each model:
# set the original data frame
# select the subjects
# loop over the brain measures to:
# build the model
# run the model, and  save the results for that model

# for each phenotype (thickness, surf, subcort)
    # for each SANS measure (SANSSum, each Global, each of the 5 and 2 factors separately)
        # run Age + Sex + SANS measure
            # for each brain measure run the model and save it
        # run Age + Sex + SANS measure + (mean thick/area/ICV)
        # run Age + Sex + SANS measure + CPZ/IQ/AO


cat("Starting pheno loop\n")
for (phenoName in c("Cort", "Surf", "SubCort")) {  # brain measure type loop

    # extract data subsets using the phenotypic-specific dataset from above
    # Use the total SANS we created

    dataset = paste0("swap.merged_ordered", phenoName)  # get the right one: Cort,Surf, SubCort
    merged_ordered0 = get(dataset)

    merged_ordered_SANSSum <- merged_ordered0[((merged_ordered0$Dx == 1) & !(is.na(merged_ordered0$SANSSum))), ]

    # Check nsub with data
    L = merged_ordered_SANSSum[!is.na(merged_ordered_SANSSum[,18]),18]
    cat(paste0("Total subjects in SANS: ", length(L), "\n"))


    # code for the globalMeasure to include
    GlobalMeasure= switch(phenoName, "Cort" = "MThickness", "Surf" = "FullSurfArea", "SubCort" = "ICV")

    # one analysis for the total SANS,
    # one analysis for each of the five factors separately, and
    # one analysis for each of the FOUR global scores separately (No Global 5)
    # one analysis for the MAP and EXP separately
    # then doing all sorts of checking for covariation...

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

        predictor_string = PredictCovs(predictor)  # gets the string for the model

        NumPredict = switch (predictor,
            "SANSGlobals" = 4,
            "SANSFactors" = 5,
            1
        )  # this sets up how many factors' effect sizes to save later

        # load data for this pheno
        merged_ordered = merged_ordered_SANSSum

        cat(paste0("dataset and pheno & predictor & global off/on ", cc," ", phenoName, " ", predictor_string, " ", GlobalMeasure, "\n"))

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
                    glob=paste0("+",GlobalMeasure, sep="")
                } else if (WCov == "WIQ") {
                    glob="+ IQ "
                } else if (WCov == "WCPZ") {
                    glob="+ CPZ "
                } else if (WCov == "WAO") {
                    glob = "+ AO "
                } else if (WCov == "WAP") {
                    glob = "+ AP "
                } else if (WCov == "WSum") {
                    glob = "+ SANSSum "
                } else {  # default -- no covariates
                    glob=""
                }


                # number of brain measures is the total columns subtracting out Covs and SANS and global means
                Npheno = ncol(merged_ordered) - ncol(Covs) - ncol(SANS) -2
                if(phenoName == "SubCort")
                    Npheno = Npheno+1  # only one global covariate in the file, not two

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


            } # end for loop on x (brain regions)
                # save results CHECK THIS
                save(r.cort, n.patients, file = paste0(outputdir, "EffectSizes_SZ_only_", predictor,
                  "_withSexAge_", WCov, phenoName, ".Rdata"))

                # Exclude components of the lm() object that from which individual-level data can be recovered
                models.cort[["fitted.values"]] <- NULL
                models.cort[["residuals"]] <- NULL
                models.cort[["model"]] <- NULL

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
