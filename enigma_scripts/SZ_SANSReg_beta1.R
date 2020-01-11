# R
options(warn=-1)

# Based on SZ_CortReg from the original cortical analysis Written by TVE, LS, and DPH for the ENIGMA SZ Working Group
args = commandArgs(trailingOnly=TRUE)
baseDir = args[1]
transferDir = args[2]

# will need to source files here with external functions
source("/computation/enigma_scripts/PackagesNeeded.R")
source("/computation/enigma_scripts/FileHeaders.R")
source("/computation/enigma_scripts/SANSCalcs.R")  # this works!
source("/computation/enigma_scripts/RegressFunc.R")
# source('./WriteRaw.R')

# file names that have to exist this script and the input files
Rscript = "/computation/enigma_scripts/SZ_SANSReg_beta1.R"
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


# calculate the SANS factorizations here and add to the SANS matrix--external function 5 negative domain scores from Strauss 2018
# Figure 1 5 global ratings (8,13,17,22, and 25) SANSTOT = sum(SANS items 1-7, 9-12, 14-16, 18-21, and 23-24) SANSGlobal = sum(SANS
# items 8, 13, 17, 22, and 25)

SANS = CalcSans(SANS)

#### calculate mean values --This could be looped...

meanCort = NULL
meanSurf = NULL
# calculate means across hemispheres
for (x in 2:35) {
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


# Get overall raw means for each of the structures --this needs to be polished in a separate function!  Uncomment when it works

# WriteRawFiles(mergedorderedCort,paste0('RawMeans_Cort',cc,'.Rdata'))
# WriteRawFiles(mergedorderedSurf,paste0('RawMeans_Surf',cc,'.Rdata'))
# WriteRawFiles(mergedorderedSubCort,paste0('RawMeans_Cort',cc,'.Rdata'))

### Regression Analyses create file for resetting to original file after modifications per each regression

swap.merged_orderedCort = merged_orderedCort
swap.merged_orderedSurf = merged_orderedSurf
swap.merged_orderedSubCort = merged_orderedSubCort

# running both complete and asis

# for(cc in c('complete','asis')){ # find the closing loop for this

cc = "asis"  # for the moment!

# this might need to refer to merged_ordered...?
if (cc == "complete") {
    # only keep subjects with all volumes
    Cort = Cort[complete.cases(Cort), ]
    SubCort = SubCort[complete.cases(SubCort), ]
    Surf = Surf[complete.cases(Surf), ]
}

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
    # for each SANS measure (SANSSum, Global measures, Factors)
        # run Age + Sex + SANS measure
            # for each brain measure run the model and save it
        # run Age + Sex + SANS measure + (mean thick/area/ICV)
        # run Age + Sex + SANS measure + IQ --this will be different group of subjects, new data selection!

# need to debug all this to see how to make it parallel across all phenotypes

cat("Starting pheno loop\n")
for (phenoName in c("Cort", "Surf", "SubCort")) {  #

    # extract data subsets using the phenotypic-specific dataset from above
    # Use the total SANS we created

    dataset = paste0("swap.merged_ordered", phenoName)  #Cort,Surf, SubCort
    merged_ordered0 = get(dataset)

    merged_ordered_SANSSum <- merged_ordered0[((merged_ordered0$Dx == 1) & !(is.na(merged_ordered0$SANSSum))), ]

    merged_ordered_SANSGlobals = merged_ordered_SANSSum  # should be all the same subjects who have a SANS Sum
    merged_ordered_SANSFactors = merged_ordered_SANSSum

    # Check nsub with data
    L = merged_ordered_SANSSum[!is.na(merged_ordered_SANSSum[,18]),18]
    cat("Total subjects in SANSSum", length(L), "\n")


    # code for the globalMeasure to include
    GlobalMeasure= switch(phenoName, "Cort" = "MThickness", "Surf" = "FullSurfArea", "SubCort" = "ICV")


    for (predictor in c("SANSSum", "SANSGlobals", "SANSFactors")) {  #

        predictor_string = PredictCovs(predictor)

        data_subset_name = paste0("merged_ordered", sep = "_", predictor)

        # load data for this pheno and SANS measure
        merged_ordered = get(data_subset_name)

        # attach(merged_ordered)  # this is critical for building the model from strings, below

        cat("dataset and pheno & predictor & global off/on ", cc, phenoName, predictor_string, GlobalMeasure, "\n")

        # 1. Regressions predicting phenoName for SZ patients against each SANS measure only covary for Sex and Age
        # (and site if needed)

        # Store models for troubleshooting
        models.cort = NULL  # This will become a list where we store all of the models made by lm

        # allocate empty vectors to store adjust effect sizes, se, ci (noicv)
        r.cort = rep(NA, (ncol(merged_ordered) - ncol(Covs)))
        se.cort = rep(NA, (ncol(merged_ordered) - ncol(Covs)))
        low.ci.cort = rep(NA, (ncol(merged_ordered) - ncol(Covs)))
        up.ci.cort = rep(NA, (ncol(merged_ordered) - ncol(Covs)))

        # n.controls=rep(NA,(ncol(merged_ordered)-ncol(Covs)))

        n.patients = rep(NA, (ncol(merged_ordered) - ncol(Covs)))

        if (nrow(merged_ordered) > (n.sites + 5)) {  # do or not do, there is no try
            # this was 3 but it should be 5

            for ( WGlobalBrain in c("WG","NotWG")) {  # loop both with and without the global measure covariate

            cat("Running: Regression predictor", predictor, " against", phenoName, ",", cc, WGlobalBrain,  "in SZ patients covary for Age and Sex \n")

            # set up model variables that don't change in this loop
            if(WGlobalBrain == "WG") {
                glob=paste0("+",GlobalMeasure, sep="")
            } else {
                glob=""
            }

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

            # Loop through and perform each regression
            # number of brain measures is the total columns subtracting out Covs and SANS
            Npheno = ncol(merged_ordered) - ncol(Covs) - ncol(SANS)
            Startcol = ncol(Covs) + 1
            Endcol = Startcol + Npheno -1 # this should vary between cortical and subcortical

            for (x in (Startcol:Endcol)) { # make this startcol to end col!
                # do the model for each brain measure (column)

                pheno = merged_ordered[!is.na(merged_ordered[, x]), x]  #Check to make sure there are observations for a given structure

                # check if the phenotype is singular after NA removal
                if (length(pheno) == 0) {
                  next  # Skip to the next measure if there are no observations
                }

                # Run the modelfor each x
                IndVar = merged_ordered[, x]

                thisformula = as.formula(paste("IndVar ~ ", predictor_string, " + Age + Sex", glob, site, sep = ""))
                tmp=lm(thisformula, data=merged_ordered)

                models.cort[[x - ncol(Covs)]] = tmp  #Store the model fit for future reference

                # subjects can be dropped if they are missing so we can get the precise number of controls/patients for each region tested
                # n.controls[x - ncol(Covs)] = length(which(tmp$model[, 2] == 0))
                # n.patients[x - ncol(Covs)] = length(which(tmp$model[, 2] == 1))
                n.controls = 0 # none used in this analysis
                n.patients = length(tmp$model[,2]) # Does this always work?? It *should* be the num of subjects included in the model

                # not entirely sure what this is??
                partcor.i <- pcor.test(tmp$model[1], tmp$model[, 2], tmp$model[, c(3:ncol(tmp$model))])
                r.cort[x - ncol(Covs)] = partcor.i[, 1]

            } # end for loop on x (brain regions)
                # save results
                save(r.cort, se.cort, low.ci.cort, up.ci.cort, n.controls, n.patients,
                     file = paste0(transferDir, '/', "EffectSizes_SZ_only_", predictor,
                  "_withSexAge_", WGlobalBrain, phenoName, ".Rdata"))
                # save(models.cort, file=paste0('Models_SZ_only_',predictor,'_withAge_',filetype,'.Rdata'))
                save(models.cort, file = paste0(transferDir, '/', "Models_SZ_only_", predictor, "_withSexAge_", WGlobalBrain, phenoName, ".Rdata"))
        } # end With Global loop

            # detach(merged_ordered)

        } #end do or not do based on sample size

        else  # not running that regression
            {
            cat("NOT Running: Regression predictor ", predictor, ",", cc, glob, " in SZ patients covary for age + Sex \n")
        }

        }  #  end predictor loop, regression 1


    } # end phenoName loop for debugging


#### TVE ##########
