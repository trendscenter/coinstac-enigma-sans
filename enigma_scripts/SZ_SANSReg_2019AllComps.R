# R
options(warn=-1)

# Based on SZ_CortReg from the original cortical analysis Written by TVE, LS, and DPH for the ENIGMA SZ Working Group
args = commandArgs(trailingOnly=TRUE)
baseDir = args[1]
transferDir = args[2]

# this still needs a log file
sink("SANscriptOut.txt", split=TRUE)

# will need to source files here with external functions
source("/computation/enigma_scripts/PackagesNeeded.R")
source("/computation/enigma_scripts/FileHeaders.R")
source("/computation/enigma_scripts/SANSCalcs.R")
source("/computation/enigma_scripts/RegressFunc.R")
source('/computation/enigma_scripts/WriteRaw.R')

# file names that have to exist this script and the input files
Rscript = "/computation/enigma_scripts/SZ_SANSReg_2019AllComps.R"
CortThickFile = file.path(baseDir, args[3])
SurfFile = file.path(baseDir, args[4])
SubCortFile = file.path(baseDir, args[5])
SANSFile = file.path(baseDir, args[6])
CovarFile = file.path(baseDir, args[7])

# make sure your working directory is where the .csv files are and this file is in that directory--this is not a complete check!
if (!file.exists(SANSFile)) {
    stop("SANS.csv file is missing  please make sure this R script is in the directory with the csv files, and your R working directory is set to that directory")
}

# read and check all files before starting on regressions etc.

Cort <- read.csv(CortThickFile, header = T)  #Read in the phenotypes file
Surf = read.csv(SurfFile, header = T)
SubCort <- read.csv(SubCortFile, header = T)
Covs = read.csv(CovarFile, header = T)
SANS = read.csv(SANSFile, header = T)



cortcolind = match(CortCols, names(Cort))
if (length(which(is.na(cortcolind))) > 0) {
    stop(
        "At least one of the required columns in your Cortical thickness measures file is missing. Make sure that the column names are spelled exactly as listed in the protocol\n"
    )
}

cortcolind = match(SurfCols, names(Surf))
if (length(which(is.na(cortcolind))) > 0) {
    stop(
        "At least one of the required columns in your Cortical surfare area measures file is missing. Make sure that the column names are spelled exactly as listed in the protocol\n"
    )
}

cortcolind = match(SubCortCols, names(SubCort))
if (length(which(is.na(cortcolind))) > 0) {
    stop(
        "At least one of the required columns in your LandRvolumesn.csv file is missing. Make sure that the column names are spelled exactly as listed in the protocol\n"
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

# Check that all of the required columns are present (SANS)

colind = match(SANScols, names(SANS))
if (length(which(is.na(colind))) > 0) {
    stop(
        "At least one of the required columns in your SANS.csv file is missing. Make sure that the column names are spelled exactly as listed:\nIt is possible that the problem column(s) is: ",
        SANScols[which(is.na(colind))]
    )
}

# Check for duplicated SubjIDs that may cause issues with merging data sets.
if (anyDuplicated(Cort[, c("SubjID")]) != 0) {
    stop(paste0("You have duplicate SubjIDs in your cortical measures file.\nMake sure there are no repeat SubjIDs.\n"))
}
if (anyDuplicated(Surf[, c("SubjID")]) != 0) {
    stop(paste0("You have duplicate SubjIDs in your cortical measures file.\nMake sure there are no repeat SubjIDs.\n"))
}

if (anyDuplicated(SubCort[, c("SubjID")]) != 0) {
    stop(paste0("You have duplicate SubjIDs in your subcortical volumes file.\nMake sure there are no repeat SubjIDs.\n"))
}

if (anyDuplicated(Covs[, c("SubjID")]) != 0) {
    stop("You have duplicate SubjIDs in your Covariates.csv file.\nMake sure there are no repeat SubjIDs.")
}

if (anyDuplicated(SANS[, c("SubjID")]) != 0) {
    stop("You have duplicate SubjIDs in your SANS.csv file.\nMake sure there are no repeat SubjIDs.")
}

# check that SANS subjIDs are same as cases in Covs, assuming controls Dx = 0 and cases = 1
# this might be needed but I don't think it is 
# pat = which(Covs$Dx == 1)
# cases = Covs[pat, ]$SubjID
# S <- SANS$SubjID
# if(sort(as.character(S)) != sort(as.character(cases))) { 
# cat(paste0('WARNING: SANS and Covs. have non-matching patient SubjIDs.','\n')) 
# cat('Please make sure the patients with Covariates and SANS are equal.','\n') }

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

# calculate the SANS factorizations here and add to the SANS matrix--external function 

SANS = CalcSans(SANS)

#### calculate mean values --This could be looped...

meanCort = NULL
meanSurf = NULL
# calculate means across hemispheres
for (x in 2:35) {  # This depends on the FS file format and cortical ROI labels
    
    meanCort = c(meanCort, ((Cort[, x] + Cort[x + 34])/2))
    meanSurf = c(meanSurf, ((Surf[, x] + Cort[x + 34])/2))
}
meanCort = c(meanCort, ((Cort[, 70] + Cort[71])/2))
meanCort = c(meanCort, ((Cort[, 72] + Cort[73])))
meanSurf = c(meanSurf, ((Surf[, 70] + Surf[71])/2))
meanSurf = c(meanSurf, ((Surf[, 72] + Surf[73])))

for (x in 1:34) {
    tmp = strsplit(names(meanCort)[x], "_")
    names(meanCort)[x] = paste0("M_", tmp[[1]][2], "_", tmp[[1]][3])
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

cat("The number of SubjIDs in the merged_ordered Cortical Thickness data set is: ", nrow(merged_orderedCort), "\n")


# Get overall raw means for each of the structures

WriteRawFiles(merged_orderedCort,ncol(Covs), 'RawMeans_Cort_asis.Rdata')
WriteRawFiles(merged_orderedSurf,ncol(Covs),'RawMeans_Surf_asis.Rdata')
WriteRawFiles(merged_orderedSubCort,ncol(Covs),'RawMeans_SubCort_asis.Rdata')

# Check if sex is related to negative symptom total severity

genderMod = lm(SANSSum ~ Age + as.factor(Sex), data=merged_orderedSubCort)
save(genderMod, file=paste0(transferDir, '/', "GenderSymptoms.Rdata"))

# cross correlation of all SANS data and factors here, save to file!
CrossCorr = cor(SANS[,2:dim(SANS)[2]], use="pairwise.complete.obs")
save(CrossCorr, file=paste0(transferDir, '/', "SANSCorrs.Rdata"))

### Regression Analyses create file for resetting to original file after modifications per each regression

swap.merged_orderedCort = merged_orderedCort
swap.merged_orderedSurf = merged_orderedSurf
swap.merged_orderedSubCort = merged_orderedSubCort

cc = "asis"  # for the moment! No complete cases only

# Check nsub with data
L = swap.merged_orderedCort[!is.na(swap.merged_orderedCort[,18]),18]
cat("Total subjects in Cortical", length(L), "\n")

# for each model: 
# set the original data frame 
# select the subjects 
#loop over the brain measures to: 
# build the model 
# run the model, and  save the results for that model

# for each phenotype (thickness, surf, subcort) 
    # for each SANS measure (SANSSum, each Global and all 5 Global measures, each factor and all 5 Factors) 
        # run Age + Sex + SANS measure
            # for each brain measure run the model and save it 
        # run Age + Sex + SANS measure + (mean thick/area/ICV) 
        # run Age + Sex + SANS measure + CPZ/IQ/AO
        

# need to debug all this to see how to make it parallel across all phenotypes

cat("Starting pheno loop\n")
for (phenoName in c("Cort", "Surf", "SubCort")) {  # brain measure type loop
    
    # extract data subsets using the phenotypic-specific dataset from above
    # Use the total SANS we created 
    
    dataset = paste0("swap.merged_ordered", phenoName)  # get the right one: Cort,Surf, SubCort
    merged_ordered0 = get(dataset)
    
    merged_ordered_SANSSum <- merged_ordered0[((merged_ordered0$Dx == 1) & !(is.na(merged_ordered0$SANSSum))), ]
    
    # Check nsub with data
    L = merged_ordered_SANSSum[!is.na(merged_ordered_SANSSum[,18]),18]
    cat("Total subjects in SANS", length(L), "\n")
    
    
    # code for the globalMeasure to include
    GlobalMeasure= switch(phenoName, "Cort" = "MThickness", "Surf" = "FullSurfArea", "SubCort" = "ICV")
    
    # one analysis for the total SANS, one analysis with all five global scores together, 
    # one for all five factors together,
    # one analysis for each of the five factors separately, and
    # one analysis for each of the global scores separately
    # add in the two factor solution here
    
    for (predictor in c("SANSSum", "SANSMAP", "SANSEXP",
                        "SANSGlobal1", "SANSGlobal2", "SANSGlobal3", "SANSGlobal4", "SANSGlobal5", 
                        "SANSFac1", "SANSFac2", "SANSFac3", "SANSFac4", "SANSFac5")) {  # All the possible predictors!
        
        predictor_string = PredictCovs(predictor)  # gets the string for the model
    
        NumPredict = switch (predictor,
            "SANSGlobals" = 5,
            "SANSFactors" = 5,
            1
        )  # this sets up how many factors' effect sizes to save later

        # load data for this pheno 
        merged_ordered = merged_ordered_SANSSum
        
        cat("dataset and pheno & predictor & global off/on ", cc, phenoName, predictor_string, GlobalMeasure, "\n")
        
        
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
            
            for ( WCov in c("NoCovs","WG", "WSum", "WIQ", "WCPZ", "WAO", "WAP")) {  # loop for other covariates: Global brain, IQ, CPZ equiv, AO, AP group
                
            cat("Running: Regression predictor", predictor, " against", phenoName, ",", cc, WCov,  
                "in SZ patients covary for Age and Sex \n")
            
            # set up model variables that don't change in this loop
                if(WCov == "WG") {
                    glob=paste0("+",GlobalMeasure, sep="")
                } else if (WCov =="WIQ") { 
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
                save(r.cort, n.patients, file = paste0(transferDir, '/',"EffectSizes_SZ_only_", predictor, 
                  "_withSexAge_", WCov, phenoName, ".Rdata"))
                
                save(models.cort, file = paste0(transferDir, '/',"Models_SZ_", predictor, "_withSexAge_", WCov,"_", phenoName, ".Rdata"))
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
