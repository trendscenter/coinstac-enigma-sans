# R

# Written by TVE, LS, and DPH for the ENIGMA SZ Working Group

# Before running the code the first time, 
# you need to download & install the "lsmeans" and "pcorr" packages,
# using the install.packages command. The R "library" or "require" commands
# can then be used to make the package available within R.
#
# If you have write access to the directory where R is installed, you
# can just open "R" by typing "R" on the command line, and once within the R command
# line, type 'install.packages("lsmeans")' and install.packages("pcorr") 
# if you do not you will want to specify a library location (lib.loc), see below.
#
# install.packages("lsmeans", lib="/YourRlibDir")
# library("lsmeans",lib.loc="/YourRlibDir")
# install.packages("ppcor", lib="/YourRlibDir")
# library("pcorr",lib.loc="/YourRlibDir")
#
# Example custom library location:
# install.packages("lsmeans", lib="/master_raid/users/tvanerp/R")
# library("lsmeans",lib.loc="/master_raid/users/tvanerp/R")
#  
# library(lsmeans)
# require(lsmeans)
#
# When you type 'search()', you should see "package:lsmeans" among the installed packages.


library(ppcor)
library(emmeans)

args = commandArgs(trailingOnly=TRUE)
baseDir = args[1]
transferDir = args[2]

#file1 = file.path(baseDir, "CorticalMeasuresENIGMA_ThickAvg.csv")
#file2 = file.path(baseDir, "CorticalMeasuresENIGMA_SurfAvg.csv")
#covarFile = file.path(baseDir, "Covariates.csv")

setwd(baseDir)

file1 = "CorticalMeasuresENIGMA_ThickAvg.csv"
file2 = "CorticalMeasuresENIGMA_SurfAvg.csv"
covarFile = "Covariates.csv"


for(cc in c("complete","asis")){

for(fsfile in c(file1,file2)){

if(fsfile == file1 && cc == "complete"){
	filetype="thickness_complete"
} else if(fsfile == file1 && cc == "asis"){
	filetype="thickness_asis"
} else if(fsfile == file2 && cc == "asis"){
	filetype="surfarea_asis"
} else {
	filetype="surfarea_complete"	
}

#create log file
messages=file(paste0(transferDir, "/", fsfile,paste0("_",cc,".log")), open="wt")
#rest=file("rest.Rout", open="wt")
sink(messages, type="message")
sink(messages, type="output")

cat('Working on ',fsfile,' a(n) ', cc, ' analysis\n')

Cort <- read.csv(fsfile); #Read in the phenotypes file

#check to make sure all of the necessary columns are present
if(fsfile == file1){
	CortCols=c("SubjID", "L_bankssts_thickavg", "L_caudalanteriorcingulate_thickavg", "L_caudalmiddlefrontal_thickavg", "L_cuneus_thickavg", "L_entorhinal_thickavg", "L_fusiform_thickavg", "L_inferiorparietal_thickavg", "L_inferiortemporal_thickavg", "L_isthmuscingulate_thickavg", "L_lateraloccipital_thickavg", "L_lateralorbitofrontal_thickavg", "L_lingual_thickavg", "L_medialorbitofrontal_thickavg", "L_middletemporal_thickavg", "L_parahippocampal_thickavg", "L_paracentral_thickavg", "L_parsopercularis_thickavg", "L_parsorbitalis_thickavg", "L_parstriangularis_thickavg", "L_pericalcarine_thickavg", "L_postcentral_thickavg", "L_posteriorcingulate_thickavg", "L_precentral_thickavg", "L_precuneus_thickavg", "L_rostralanteriorcingulate_thickavg", "L_rostralmiddlefrontal_thickavg", "L_superiorfrontal_thickavg", "L_superiorparietal_thickavg", "L_superiortemporal_thickavg", "L_supramarginal_thickavg", "L_frontalpole_thickavg", "L_temporalpole_thickavg", "L_transversetemporal_thickavg", "L_insula_thickavg", "R_bankssts_thickavg", "R_caudalanteriorcingulate_thickavg", "R_caudalmiddlefrontal_thickavg", "R_cuneus_thickavg", "R_entorhinal_thickavg", "R_fusiform_thickavg", "R_inferiorparietal_thickavg", "R_inferiortemporal_thickavg", "R_isthmuscingulate_thickavg", "R_lateraloccipital_thickavg", "R_lateralorbitofrontal_thickavg", "R_lingual_thickavg", "R_medialorbitofrontal_thickavg", "R_middletemporal_thickavg", "R_parahippocampal_thickavg", "R_paracentral_thickavg", "R_parsopercularis_thickavg", "R_parsorbitalis_thickavg", "R_parstriangularis_thickavg", "R_pericalcarine_thickavg", "R_postcentral_thickavg", "R_posteriorcingulate_thickavg", "R_precentral_thickavg", "R_precuneus_thickavg", "R_rostralanteriorcingulate_thickavg", "R_rostralmiddlefrontal_thickavg", "R_superiorfrontal_thickavg", "R_superiorparietal_thickavg", "R_superiortemporal_thickavg", "R_supramarginal_thickavg", "R_frontalpole_thickavg", "R_temporalpole_thickavg", "R_transversetemporal_thickavg", "R_insula_thickavg", "LThickness", "RThickness", "LSurfArea", "RSurfArea", "ICV")
	cortcolind=match(CortCols,names(Cort))
	if(length(which(is.na(cortcolind))) > 0){
		stop('At least one of the required columns in your ', fsfile, ' file is missing. Make sure that the column names are spelled exactly as listed in the protocol\n')
	}
} else {
	CortCols=c("SubjID", "L_bankssts_surfavg", "L_caudalanteriorcingulate_surfavg", "L_caudalmiddlefrontal_surfavg", "L_cuneus_surfavg", "L_entorhinal_surfavg", "L_fusiform_surfavg", "L_inferiorparietal_surfavg", "L_inferiortemporal_surfavg", "L_isthmuscingulate_surfavg", "L_lateraloccipital_surfavg", "L_lateralorbitofrontal_surfavg", "L_lingual_surfavg", "L_medialorbitofrontal_surfavg", "L_middletemporal_surfavg", "L_parahippocampal_surfavg", "L_paracentral_surfavg", "L_parsopercularis_surfavg", "L_parsorbitalis_surfavg", "L_parstriangularis_surfavg", "L_pericalcarine_surfavg", "L_postcentral_surfavg", "L_posteriorcingulate_surfavg", "L_precentral_surfavg", "L_precuneus_surfavg", "L_rostralanteriorcingulate_surfavg", "L_rostralmiddlefrontal_surfavg", "L_superiorfrontal_surfavg", "L_superiorparietal_surfavg", "L_superiortemporal_surfavg", "L_supramarginal_surfavg", "L_frontalpole_surfavg", "L_temporalpole_surfavg", "L_transversetemporal_surfavg", "L_insula_surfavg", "R_bankssts_surfavg", "R_caudalanteriorcingulate_surfavg", "R_caudalmiddlefrontal_surfavg", "R_cuneus_surfavg", "R_entorhinal_surfavg", "R_fusiform_surfavg", "R_inferiorparietal_surfavg", "R_inferiortemporal_surfavg", "R_isthmuscingulate_surfavg", "R_lateraloccipital_surfavg", "R_lateralorbitofrontal_surfavg", "R_lingual_surfavg", "R_medialorbitofrontal_surfavg", "R_middletemporal_surfavg", "R_parahippocampal_surfavg", "R_paracentral_surfavg", "R_parsopercularis_surfavg", "R_parsorbitalis_surfavg", "R_parstriangularis_surfavg", "R_pericalcarine_surfavg", "R_postcentral_surfavg", "R_posteriorcingulate_surfavg", "R_precentral_surfavg", "R_precuneus_surfavg", "R_rostralanteriorcingulate_surfavg", "R_rostralmiddlefrontal_surfavg", "R_superiorfrontal_surfavg", "R_superiorparietal_surfavg", "R_superiortemporal_surfavg", "R_supramarginal_surfavg", "R_frontalpole_surfavg", "R_temporalpole_surfavg", "R_transversetemporal_surfavg", "R_insula_surfavg", "LThickness", "RThickness", "LSurfArea", "RSurfArea", "ICV")
	cortcolind=match(CortCols,names(Cort))
	if(length(which(is.na(cortcolind))) > 0){
		stop('At least one of the required columns in your ', fsfile, ' file is missing. Make sure that the column names are spelled exactly as listed in the protocol\n')
	}
}

if(cc == "complete"){
	Cort = Cort[complete.cases(Cort),]
}

####
####calculate mean values
####

meanCort=NULL
#calculate means
for(x in 2:35){
	meanCort = c(meanCort, ((Cort[,x] + Cort[x+34])/2))
}
meanCort = c(meanCort, ((Cort[,70] + Cort[71])/2))
meanCort = c(meanCort, ((Cort[,72] + Cort[73])))

for(x in 1:34){
	tmp=strsplit(names(meanCort)[x],"_")
	names(meanCort)[x]=paste0("M_",tmp[[1]][2],"_",tmp[[1]][3])
}
names(meanCort)[35]="MThickness"
names(meanCort)[36]="FullSurfArea"
meanCort=as.data.frame(meanCort)

#drop ICV from Cort file
Cort = Cort[,-ncol(Cort)]

#combine Cort file with the newly calculated means
Cort = cbind(Cort, meanCort)

####
####
####

# Check for duplicated SubjIDs that may cause issues with merging data sets.
if(anyDuplicated(Cort[,c("SubjID")]) != 0) { stop(paste0('You have duplicate SubjIDs in your ', fsfile, 'file.\nMake sure there are no repeat SubjIDs.\n')) }


Covs <- read.csv(covarFile); #Read in the covariates file

# Check for duplicated SubjIDs that may cause issues with merging data sets.
if(anyDuplicated(Covs[,c("SubjID")]) != 0) { stop('You have duplicate SubjIDs in your Covariates.csv file.\nMake sure there are no repeat SubjIDs.') }
	
#Check that all of the required columns are present
mcols=c("SubjID","Dx","Age","Sex","AP","CPZ","AO","DURILL","PANSSTOT","PANSSPOS","PANSSNEG","SAPSTOT","SANSTOT","HAND","PARENTSES","IQ")
colind=match(mcols,names(Covs))
if(length(which(is.na(colind))) > 0){
	stop('At least one of the required columns in your Covariates.csv file is missing. Make sure that the column names are spelled exactly as listed:\nIt is possible that the problem column(s) is: ', mcols[which(is.na(colind))])
}

n.covs <- ncol(Covs) - 1; #Total number of covariates, -1 removes the SubjectID column
n.sites <- n.covs - 15; #Find the number of site variables, subtract the number of predictirs (Dx, Age, Sex etc.) from n.covs

#combine the files into one dataframe
merged_ordered = merge(Covs, Cort, by="SubjID");

#create file for resetting to original file after modifications
swap.merged_ordered=merged_ordered

#Check that the number of rows after merging is the same
if(nrow(Cort) != nrow(merged_ordered)){
       cat(paste0('WARNING: ', fsfile, ' and Covariates.csv have non-matching SubjIDs.','\n'))
       cat('Please make sure the number of subjects in your merged data set are as expected.','\n')
       cat(paste0('The number of SubjIDs in ', fsfile, ' is: ',nrow(Cort),'\n'))
       cat('The number of SubjIDs in the merged_ordered data set is: ',nrow(merged_ordered),'\n')
}

cat('Calculating demographics for ', fsfile,'\n')


#Get overall raw means for each of the structures
raw.means=colMeans(merged_ordered[(ncol(Covs)+1):ncol(merged_ordered)], na.rm=T)

#Get raw sd and number of subjects included for each of the structures
sd.raw=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
n.raw=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
min.raw=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
max.raw=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
for(z in (ncol(Covs)+1):ncol(merged_ordered)){
	sd.raw[z-ncol(Covs)]=sd(merged_ordered[,z], na.rm=T)
	n.raw[z-ncol(Covs)]=length(merged_ordered[which(!is.na(merged_ordered[,z])),z])
	min.raw[z-ncol(Covs)]=min(merged_ordered[,z], na.rm=T)
    max.raw[z-ncol(Covs)]=max(merged_ordered[,z], na.rm=T)
}

#Save raw values
save(raw.means, sd.raw, n.raw, min.raw, max.raw, file=paste0(transferDir, "/","RawMeans_",filetype,".Rdata"))


#Get raw means for each of the structures for the Dx groups:
# patients (Dx==1)
# controls (Dx==0)

for(DXvalue in 0:1){

   DXgroup=which(merged_ordered$Dx==DXvalue)
   raw.means=colMeans(merged_ordered[DXgroup,(ncol(Covs)+1):ncol(merged_ordered)], na.rm=T)

   #Get raw sd and number of subjects included for each of the structures for the Dx groups
   sd.raw=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   n.raw=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   min.raw=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   max.raw=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   for(z in (ncol(Covs)+1):ncol(merged_ordered)){
        sd.raw[z-ncol(Covs)]=sd(merged_ordered[DXgroup,z], na.rm=T)
        n.raw[z-ncol(Covs)]=length(merged_ordered[which(!is.na(merged_ordered[DXgroup,z])),z])
        min.raw[z-ncol(Covs)]=min(merged_ordered[DXgroup,z], na.rm=T)
        max.raw[z-ncol(Covs)]=max(merged_ordered[DXgroup,z], na.rm=T)
   }

   #Save raw values
   save(raw.means, sd.raw, n.raw, min.raw, max.raw, file=paste0(transferDir, "/","RawMeans_DX_",DXvalue,"_",filetype,".Rdata"))    
}


#Get raw means for each of the structures for the antipsychotic user groups
#AP groups:
#AP==0: controls
#AP==1: unmedicated patients
#AP==2: patients on typical antipsychotics
#AP==3: patients on atypical antipsychotics
#AP==4: patients on both typical and atypical antipsychotics  

for(APvalue in 0:4){

   APgroup=which(merged_ordered$AP==APvalue)
   # TVE: without the if statement, this code generate warnings when when the number of subjects in an APgroup=0
   # if(length(APgroup)>0) {
   raw.means=colMeans(merged_ordered[APgroup,(ncol(Covs)+1):ncol(merged_ordered)], na.rm=T)

   #Get raw sd and number of subjects included for each of the structures for the APgroup
   sd.raw=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   n.raw=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   min.raw=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   max.raw=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   for(z in (ncol(Covs)+1):ncol(merged_ordered)){
        sd.raw[z-ncol(Covs)]=sd(merged_ordered[APgroup,z], na.rm=T)
        n.raw[z-ncol(Covs)]=length(merged_ordered[which(!is.na(merged_ordered[APgroup,z])),z])
        min.raw[z-ncol(Covs)]=min(merged_ordered[APgroup,z], na.rm=T)
        max.raw[z-ncol(Covs)]=max(merged_ordered[APgroup,z], na.rm=T)
   }
   #Save raw values
   save(raw.means, sd.raw, n.raw, min.raw, max.raw, file=paste0(transferDir, "/","RawMeans_AP_",APvalue,"_",filetype,".Rdata"))  
  #} # end if (length(APgroup)>0){
}


#Get demographics
age.mu=mean(merged_ordered$Age,na.rm=T) #raw age mean
age.sd=sd(merged_ordered$Age,na.rm=T) #raw age sd
age.range=range(merged_ordered$Age,na.rm=T) #raw age range

age.mu.dx0=mean(merged_ordered$Age[which(merged_ordered$Dx==0)],na.rm=T) #age mean for ctls
age.sd.dx0=sd(merged_ordered$Age[which(merged_ordered$Dx==0)],na.rm=T) #age sd for ctls
age.range.dx0=range(merged_ordered$Age[which(merged_ordered$Dx==0)],na.rm=T) #age range for ctls
age.mu.dx1=mean(merged_ordered$Age[which(merged_ordered$Dx==1)],na.rm=T) #age mean for patients
age.sd.dx1=sd(merged_ordered$Age[which(merged_ordered$Dx==1)],na.rm=T) #age sd for patients
age.range.dx1=range(merged_ordered$Age[which(merged_ordered$Dx==1)],na.rm=T) #age range for patients

parentses.mu.dx0=mean(merged_ordered$PARENTSES[which(merged_ordered$Dx==0)],na.rm=T) #parentses mean for ctls
parentses.sd.dx0=sd(merged_ordered$PARENTSES[which(merged_ordered$Dx==0)],na.rm=T) #parentses sd for ctls
parentses.range.dx0=range(merged_ordered$PARENTSES[which(merged_ordered$Dx==0)],na.rm=T) #parentses range for ctls
parentses.n.dx0=length(which(!is.na(merged_ordered$PARENTSES[which(merged_ordered$Dx==0)]))) #parentses n for ctls

iq.mu.dx0=mean(merged_ordered$IQ[which(merged_ordered$Dx==0)],na.rm=T) #iq mean for ctls
iq.sd.dx0=sd(merged_ordered$IQ[which(merged_ordered$Dx==0)],na.rm=T) #iq sd for ctls
iq.range.dx0=range(merged_ordered$IQ[which(merged_ordered$Dx==0)],na.rm=T) #iq range for ctls
iq.n.dx0=length(which(!is.na(merged_ordered$IQ[which(merged_ordered$Dx==0)]))) #iq n for ctls

iq.mu.dx1=mean(merged_ordered$IQ[which(merged_ordered$Dx==1)],na.rm=T) #iq mean for ctls
iq.sd.dx1=sd(merged_ordered$IQ[which(merged_ordered$Dx==1)],na.rm=T) #iq sd for ctls
iq.range.dx1=range(merged_ordered$IQ[which(merged_ordered$Dx==1)],na.rm=T) #iq range for ctls
iq.n.dx1=length(which(!is.na(merged_ordered$IQ[which(merged_ordered$Dx==1)]))) #iq n for ctls

age.mu.AP0=mean(merged_ordered$Age[which(merged_ordered$AP==0)],na.rm=T)   #age mean for controls
age.sd.AP0=sd(merged_ordered$Age[which(merged_ordered$AP==0)],na.rm=T)   #age sd for controls
age.range.AP0=range(merged_ordered$Age[which(merged_ordered$AP==0)],na.rm=T)   #age range for controls
age.mu.AP1=mean(merged_ordered$Age[which(merged_ordered$AP==1)],na.rm=T)   #age mean for unmedicated patients
age.sd.AP1=sd(merged_ordered$Age[which(merged_ordered$AP==1)],na.rm=T)   #age sd for unmedicated patients
age.range.AP1=range(merged_ordered$Age[which(merged_ordered$AP==1)],na.rm=T)   #age range for unmedicated patients
age.mu.AP2=mean(merged_ordered$Age[which(merged_ordered$AP==2)],na.rm=T)   #age mean for typical AP patients
age.sd.AP2=sd(merged_ordered$Age[which(merged_ordered$AP==2)],na.rm=T)   #age sd for typical AP patients
age.range.AP2=range(merged_ordered$Age[which(merged_ordered$AP==2)],na.rm=T)   #age range for typical AP patients
age.mu.AP3=mean(merged_ordered$Age[which(merged_ordered$AP==3)],na.rm=T)   #age mean for atypical AP patients
age.sd.AP3=sd(merged_ordered$Age[which(merged_ordered$AP==3)],na.rm=T)   #age sd for atypical AP patients
age.range.AP3=range(merged_ordered$Age[which(merged_ordered$AP==3)],na.rm=T)   #age range for atypical AP patients
age.mu.AP4=mean(merged_ordered$Age[which(merged_ordered$AP==4)],na.rm=T)   #age mean for both a and t AP patients
age.sd.AP4=sd(merged_ordered$Age[which(merged_ordered$AP==4)],na.rm=T)   #age sd for both a and t AP patients
age.range.AP4=range(merged_ordered$Age[which(merged_ordered$AP==4)],na.rm=T)   #age range for both a and t AP patients

cpz.mu.dx1=mean(merged_ordered$CPZ[which(merged_ordered$Dx==1)],na.rm=T) #cpz mean for patients
cpz.sd.dx1=sd(merged_ordered$CPZ[which(merged_ordered$Dx==1)],na.rm=T) #cpz sd for patients
cpz.range.dx1=range(merged_ordered$CPZ[which(merged_ordered$Dx==1)],na.rm=T) #cpz range for patients
cpz.mu.AP1=mean(merged_ordered$CPZ[which(merged_ordered$AP==1)],na.rm=T)   #cpz mean for unmedicated patients
cpz.sd.AP1=sd(merged_ordered$CPZ[which(merged_ordered$AP==1)],na.rm=T)   #cpz sd for unmedicated patients
cpz.range.AP1=range(merged_ordered$CPZ[which(merged_ordered$AP==1)],na.rm=T)   #cpz range for unmedicated patients
cpz.mu.AP2=mean(merged_ordered$CPZ[which(merged_ordered$AP==2)],na.rm=T)   #cpz mean for typical AP patients
cpz.sd.AP2=sd(merged_ordered$CPZ[which(merged_ordered$AP==2)],na.rm=T)   #cpz sd for typical AP patients
cpz.range.AP2=range(merged_ordered$CPZ[which(merged_ordered$AP==2)],na.rm=T)   #cpz range for typical AP patients
cpz.mu.AP3=mean(merged_ordered$CPZ[which(merged_ordered$AP==3)],na.rm=T)   #cpz mean for atypical AP patients
cpz.sd.AP3=sd(merged_ordered$CPZ[which(merged_ordered$AP==3)],na.rm=T)   #cpz sd for atypical AP patients
cpz.range.AP3=range(merged_ordered$CPZ[which(merged_ordered$AP==3)],na.rm=T)   #cpz range for atypical AP patients
cpz.mu.AP4=mean(merged_ordered$CPZ[which(merged_ordered$AP==4)],na.rm=T)   #cpz mean for both a and t AP patients
cpz.sd.AP4=sd(merged_ordered$CPZ[which(merged_ordered$AP==4)],na.rm=T)   #cpz sd for both a and t AP patients
cpz.range.AP4=range(merged_ordered$CPZ[which(merged_ordered$AP==4)],na.rm=T)   #cpz range for both a and t AP patients

ao.mu.dx1=mean(merged_ordered$AO[which(merged_ordered$Dx==1)],na.rm=T) #ao mean for patients
ao.sd.dx1=sd(merged_ordered$AO[which(merged_ordered$Dx==1)],na.rm=T) #ao sd for patients
ao.range.dx1=range(merged_ordered$AO[which(merged_ordered$Dx==1)],na.rm=T) #ao range for patients
ao.mu.AP1=mean(merged_ordered$AO[which(merged_ordered$AP==1)],na.rm=T)   #ao mean for unmedicated patients
ao.sd.AP1=sd(merged_ordered$AO[which(merged_ordered$AP==1)],na.rm=T)   #ao sd for unmedicated patients
ao.range.AP1=range(merged_ordered$AO[which(merged_ordered$AP==1)],na.rm=T)   #ao range for unmedicated patients
ao.mu.AP2=mean(merged_ordered$AO[which(merged_ordered$AP==2)],na.rm=T)   #ao mean for typical AP patients
ao.sd.AP2=sd(merged_ordered$AO[which(merged_ordered$AP==2)],na.rm=T)   #ao sd for typical AP patients
ao.range.AP2=range(merged_ordered$AO[which(merged_ordered$AP==2)],na.rm=T)   #ao range for typical AP patients
ao.mu.AP3=mean(merged_ordered$AO[which(merged_ordered$AP==3)],na.rm=T)   #ao mean for atypical AP patients
ao.sd.AP3=sd(merged_ordered$AO[which(merged_ordered$AP==3)],na.rm=T)   #ao sd for atypical AP patients
ao.range.AP3=range(merged_ordered$AO[which(merged_ordered$AP==3)],na.rm=T)   #ao range for atypical AP patients
ao.mu.AP4=mean(merged_ordered$AO[which(merged_ordered$AP==4)],na.rm=T)   #ao mean for both a and t AP patients
ao.sd.AP4=sd(merged_ordered$AO[which(merged_ordered$AP==4)],na.rm=T)   #ao sd for both a and t AP patients
ao.range.AP4=range(merged_ordered$AO[which(merged_ordered$AP==4)],na.rm=T)   #ao range for both a and t AP patients

durill.mu.dx1=mean(merged_ordered$DURILL[which(merged_ordered$Dx==1)],na.rm=T) #durill mean for patients
durill.sd.dx1=sd(merged_ordered$DURILL[which(merged_ordered$Dx==1)],na.rm=T) #durill sd for patients
durill.range.dx1=range(merged_ordered$DURILL[which(merged_ordered$Dx==1)],na.rm=T) #durill range for patients
durill.mu.AP1=mean(merged_ordered$DURILL[which(merged_ordered$AP==1)],na.rm=T)   #durill mean for unmedicated patients
durill.sd.AP1=sd(merged_ordered$DURILL[which(merged_ordered$AP==1)],na.rm=T)   #durill sd for unmedicated patients
durill.range.AP1=range(merged_ordered$DURILL[which(merged_ordered$AP==1)],na.rm=T)   #durill range for unmedicated patients
durill.mu.AP2=mean(merged_ordered$DURILL[which(merged_ordered$AP==2)],na.rm=T)   #durill mean for typical AP patients
durill.sd.AP2=sd(merged_ordered$DURILL[which(merged_ordered$AP==2)],na.rm=T)   #durill sd for typical AP patients
durill.range.AP2=range(merged_ordered$DURILL[which(merged_ordered$AP==2)],na.rm=T)   #durill range for typical AP patients
durill.mu.AP3=mean(merged_ordered$DURILL[which(merged_ordered$AP==3)],na.rm=T)   #durill mean for atypical AP patients
durill.sd.AP3=sd(merged_ordered$DURILL[which(merged_ordered$AP==3)],na.rm=T)   #durill sd for atypical AP patients
durill.range.AP3=range(merged_ordered$DURILL[which(merged_ordered$AP==3)],na.rm=T)   #durill range for atypical AP patients
durill.mu.AP4=mean(merged_ordered$DURILL[which(merged_ordered$AP==4)],na.rm=T)   #durill mean for both a and t AP patients
durill.sd.AP4=sd(merged_ordered$DURILL[which(merged_ordered$AP==4)],na.rm=T)   #durill sd for both a and t AP patients
durill.range.AP4=range(merged_ordered$DURILL[which(merged_ordered$AP==4)],na.rm=T)   #durill range for both a and t AP patients

pansstot.mu.dx1=mean(merged_ordered$PANSSTOT[which(merged_ordered$Dx==1)],na.rm=T) #pansstot mean for patients
pansstot.sd.dx1=sd(merged_ordered$PANSSTOT[which(merged_ordered$Dx==1)],na.rm=T) #pansstot sd for patients
pansstot.range.dx1=range(merged_ordered$PANSSTOT[which(merged_ordered$Dx==1)],na.rm=T) #pansstot range for patients
pansstot.mu.AP1=mean(merged_ordered$PANSSTOT[which(merged_ordered$AP==1)],na.rm=T)   #pansstot mean for unmedicated patients
pansstot.sd.AP1=sd(merged_ordered$PANSSTOT[which(merged_ordered$AP==1)],na.rm=T)   #pansstot sd for unmedicated patients
pansstot.range.AP1=range(merged_ordered$PANSSTOT[which(merged_ordered$AP==1)],na.rm=T)   #pansstot range for unmedicated patients
pansstot.mu.AP2=mean(merged_ordered$PANSSTOT[which(merged_ordered$AP==2)],na.rm=T)   #pansstot mean for typical AP patients
pansstot.sd.AP2=sd(merged_ordered$PANSSTOT[which(merged_ordered$AP==2)],na.rm=T)   #pansstot sd for typical AP patients
pansstot.range.AP2=range(merged_ordered$PANSSTOT[which(merged_ordered$AP==2)],na.rm=T)   #pansstot range for typical AP patients
pansstot.mu.AP3=mean(merged_ordered$PANSSTOT[which(merged_ordered$AP==3)],na.rm=T)   #pansstot mean for atypical AP patients
pansstot.sd.AP3=sd(merged_ordered$PANSSTOT[which(merged_ordered$AP==3)],na.rm=T)   #pansstot sd for atypical AP patients
pansstot.range.AP3=range(merged_ordered$PANSSTOT[which(merged_ordered$AP==3)],na.rm=T)   #pansstot range for atypical AP patients
pansstot.mu.AP4=mean(merged_ordered$PANSSTOT[which(merged_ordered$AP==4)],na.rm=T)   #pansstot mean for both a and t AP patients
pansstot.sd.AP4=sd(merged_ordered$PANSSTOT[which(merged_ordered$AP==4)],na.rm=T)   #pansstot sd for both a and t AP patients
pansstot.range.AP4=range(merged_ordered$PANSSTOT[which(merged_ordered$AP==4)],na.rm=T)   #pansstot range for both a and t AP patients

pansspos.mu.dx1=mean(merged_ordered$PANSSPOS[which(merged_ordered$Dx==1)],na.rm=T) #pansspos mean for patients
pansspos.sd.dx1=sd(merged_ordered$PANSSPOS[which(merged_ordered$Dx==1)],na.rm=T) #pansspos sd for patients
pansspos.range.dx1=range(merged_ordered$PANSSPOS[which(merged_ordered$Dx==1)],na.rm=T) #pansspos range for patients
pansspos.mu.AP1=mean(merged_ordered$PANSSPOS[which(merged_ordered$AP==1)],na.rm=T)   #pansspos mean for unmedicated patients
pansspos.sd.AP1=sd(merged_ordered$PANSSPOS[which(merged_ordered$AP==1)],na.rm=T)   #pansspos sd for unmedicated patients
pansspos.range.AP1=range(merged_ordered$PANSSPOS[which(merged_ordered$AP==1)],na.rm=T)   #pansspos range for unmedicated patients
pansspos.mu.AP2=mean(merged_ordered$PANSSPOS[which(merged_ordered$AP==2)],na.rm=T)   #pansspos mean for typical AP patients
pansspos.sd.AP2=sd(merged_ordered$PANSSPOS[which(merged_ordered$AP==2)],na.rm=T)   #pansspos sd for typical AP patients
pansspos.range.AP2=range(merged_ordered$PANSSPOS[which(merged_ordered$AP==2)],na.rm=T)   #pansspos range for typical AP patients
pansspos.mu.AP3=mean(merged_ordered$PANSSPOS[which(merged_ordered$AP==3)],na.rm=T)   #pansspos mean for atypical AP patients
pansspos.sd.AP3=sd(merged_ordered$PANSSPOS[which(merged_ordered$AP==3)],na.rm=T)   #pansspos sd for atypical AP patients
pansspos.range.AP3=range(merged_ordered$PANSSPOS[which(merged_ordered$AP==3)],na.rm=T)   #pansspos range for atypical AP patients
pansspos.mu.AP4=mean(merged_ordered$PANSSPOS[which(merged_ordered$AP==4)],na.rm=T)   #pansspos mean for both a and t AP patients
pansspos.sd.AP4=sd(merged_ordered$PANSSPOS[which(merged_ordered$AP==4)],na.rm=T)   #pansspos sd for both a and t AP patients
pansspos.range.AP4=range(merged_ordered$PANSSPOS[which(merged_ordered$AP==4)],na.rm=T)   #pansspos range for both a and t AP patients

panssneg.mu.dx1=mean(merged_ordered$PANSSNEG[which(merged_ordered$Dx==1)],na.rm=T) #panssneg mean for patients
panssneg.sd.dx1=sd(merged_ordered$PANSSNEG[which(merged_ordered$Dx==1)],na.rm=T) #panssneg sd for patients
panssneg.range.dx1=range(merged_ordered$PANSSNEG[which(merged_ordered$Dx==1)],na.rm=T) #panssneg range for patients
panssneg.mu.AP1=mean(merged_ordered$PANSSNEG[which(merged_ordered$AP==1)],na.rm=T)   #panssneg mean for unmedicated patients
panssneg.sd.AP1=sd(merged_ordered$PANSSNEG[which(merged_ordered$AP==1)],na.rm=T)   #panssneg sd for unmedicated patients
panssneg.range.AP1=range(merged_ordered$PANSSNEG[which(merged_ordered$AP==1)],na.rm=T)   #panssneg range for unmedicated patients
panssneg.mu.AP2=mean(merged_ordered$PANSSNEG[which(merged_ordered$AP==2)],na.rm=T)   #panssneg mean for typical AP patients
panssneg.sd.AP2=sd(merged_ordered$PANSSNEG[which(merged_ordered$AP==2)],na.rm=T)   #panssneg sd for typical AP patients
panssneg.range.AP2=range(merged_ordered$PANSSNEG[which(merged_ordered$AP==2)],na.rm=T)   #panssneg range for typical AP patients
panssneg.mu.AP3=mean(merged_ordered$PANSSNEG[which(merged_ordered$AP==3)],na.rm=T)   #panssneg mean for atypical AP patients
panssneg.sd.AP3=sd(merged_ordered$PANSSNEG[which(merged_ordered$AP==3)],na.rm=T)   #panssneg sd for atypical AP patients
panssneg.range.AP3=range(merged_ordered$PANSSNEG[which(merged_ordered$AP==3)],na.rm=T)   #panssneg range for atypical AP patients
panssneg.mu.AP4=mean(merged_ordered$PANSSNEG[which(merged_ordered$AP==4)],na.rm=T)   #panssneg mean for both a and t AP patients
panssneg.sd.AP4=sd(merged_ordered$PANSSNEG[which(merged_ordered$AP==4)],na.rm=T)   #panssneg sd for both a and t AP patients
panssneg.range.AP4=range(merged_ordered$PANSSNEG[which(merged_ordered$AP==4)],na.rm=T)   #panssneg range for both a and t AP patients

sapstot.mu.dx1=mean(merged_ordered$SAPSTOT[which(merged_ordered$Dx==1)],na.rm=T) #sapstot mean for patients
sapstot.sd.dx1=sd(merged_ordered$SAPSTOT[which(merged_ordered$Dx==1)],na.rm=T) #sapstot sd for patients
sapstot.range.dx1=range(merged_ordered$SAPSTOT[which(merged_ordered$Dx==1)],na.rm=T) #sapstot range for patients
sapstot.mu.AP1=mean(merged_ordered$SAPSTOT[which(merged_ordered$AP==1)],na.rm=T)   #sapstot mean for unmedicated patients
sapstot.sd.AP1=sd(merged_ordered$SAPSTOT[which(merged_ordered$AP==1)],na.rm=T)   #sapstot sd for unmedicated patients
sapstot.range.AP1=range(merged_ordered$SAPSTOT[which(merged_ordered$AP==1)],na.rm=T)   #sapstot range for unmedicated patients
sapstot.mu.AP2=mean(merged_ordered$SAPSTOT[which(merged_ordered$AP==2)],na.rm=T)   #sapstot mean for typical AP patients
sapstot.sd.AP2=sd(merged_ordered$SAPSTOT[which(merged_ordered$AP==2)],na.rm=T)   #sapstot sd for typical AP patients
sapstot.range.AP2=range(merged_ordered$SAPSTOT[which(merged_ordered$AP==2)],na.rm=T)   #sapstot range for typical AP patients
sapstot.mu.AP3=mean(merged_ordered$SAPSTOT[which(merged_ordered$AP==3)],na.rm=T)   #sapstot mean for atypical AP patients
sapstot.sd.AP3=sd(merged_ordered$SAPSTOT[which(merged_ordered$AP==3)],na.rm=T)   #sapstot sd for atypical AP patients
sapstot.range.AP3=range(merged_ordered$SAPSTOT[which(merged_ordered$AP==3)],na.rm=T)   #sapstot range for atypical AP patients
sapstot.mu.AP4=mean(merged_ordered$SAPSTOT[which(merged_ordered$AP==4)],na.rm=T)   #sapstot mean for both a and t AP patients
sapstot.sd.AP4=sd(merged_ordered$SAPSTOT[which(merged_ordered$AP==4)],na.rm=T)   #sapstot sd for both a and t AP patients
sapstot.range.AP4=range(merged_ordered$SAPSTOT[which(merged_ordered$AP==4)],na.rm=T)   #sapstot range for both a and t AP patients

sanstot.mu.dx1=mean(merged_ordered$SANSTOT[which(merged_ordered$Dx==1)],na.rm=T) #sanstot mean for patients
sanstot.sd.dx1=sd(merged_ordered$SANSTOT[which(merged_ordered$Dx==1)],na.rm=T) #sanstot sd for patients
sanstot.range.dx1=range(merged_ordered$SANSTOT[which(merged_ordered$Dx==1)],na.rm=T) #sanstot range for patients
sanstot.mu.AP1=mean(merged_ordered$SANSTOT[which(merged_ordered$AP==1)],na.rm=T)   #sanstot mean for unmedicated patients
sanstot.sd.AP1=sd(merged_ordered$SANSTOT[which(merged_ordered$AP==1)],na.rm=T)   #sanstot sd for unmedicated patients
sanstot.range.AP1=range(merged_ordered$SANSTOT[which(merged_ordered$AP==1)],na.rm=T)   #sanstot range for unmedicated patients
sanstot.mu.AP2=mean(merged_ordered$SANSTOT[which(merged_ordered$AP==2)],na.rm=T)   #sanstot mean for typical AP patients
sanstot.sd.AP2=sd(merged_ordered$SANSTOT[which(merged_ordered$AP==2)],na.rm=T)   #sanstot sd for typical AP patients
sanstot.range.AP2=range(merged_ordered$SANSTOT[which(merged_ordered$AP==2)],na.rm=T)   #sanstot range for typical AP patients
sanstot.mu.AP3=mean(merged_ordered$SANSTOT[which(merged_ordered$AP==3)],na.rm=T)   #sanstot mean for atypical AP patients
sanstot.sd.AP3=sd(merged_ordered$SANSTOT[which(merged_ordered$AP==3)],na.rm=T)   #sanstot sd for atypical AP patients
sanstot.range.AP3=range(merged_ordered$SANSTOT[which(merged_ordered$AP==3)],na.rm=T)   #sanstot range for atypical AP patients
sanstot.mu.AP4=mean(merged_ordered$SANSTOT[which(merged_ordered$AP==4)],na.rm=T)   #sanstot mean for both a and t AP patients
sanstot.sd.AP4=sd(merged_ordered$SANSTOT[which(merged_ordered$AP==4)],na.rm=T)   #sanstot sd for both a and t AP patients
sanstot.range.AP4=range(merged_ordered$SANSTOT[which(merged_ordered$AP==4)],na.rm=T)   #sanstot range for both a and t AP patients

parentses.mu.dx1=mean(merged_ordered$PARENTSES[which(merged_ordered$Dx==1)],na.rm=T) #parentses mean for patients
parentses.sd.dx1=sd(merged_ordered$PARENTSES[which(merged_ordered$Dx==1)],na.rm=T) #parentses sd for patients
parentses.range.dx1=range(merged_ordered$PARENTSES[which(merged_ordered$Dx==1)],na.rm=T) #parentses range for patients
parentses.n.dx1=length(which(!is.na(merged_ordered$PARENTSES[which(merged_ordered$Dx==1)]))) #parentses range for patients

parentses.mu.AP1=mean(merged_ordered$PARENTSES[which(merged_ordered$AP==1)],na.rm=T)   #sanstot mean for unmedicated patients
parentses.sd.AP1=sd(merged_ordered$PARENTSES[which(merged_ordered$AP==1)],na.rm=T)   #sanstot sd for unmedicated patients
parentses.range.AP1=range(merged_ordered$PARENTSES[which(merged_ordered$AP==1)],na.rm=T)   #sanstot range for unmedicated patients
parentses.mu.AP2=mean(merged_ordered$PARENTSES[which(merged_ordered$AP==2)],na.rm=T)   #sanstot mean for typical AP patients
parentses.sd.AP2=sd(merged_ordered$PARENTSES[which(merged_ordered$AP==2)],na.rm=T)   #sanstot sd for typical AP patients
parentses.range.AP2=range(merged_ordered$PARENTSES[which(merged_ordered$AP==2)],na.rm=T)   #sanstot range for typical AP patients
parentses.mu.AP3=mean(merged_ordered$PARENTSES[which(merged_ordered$AP==3)],na.rm=T)   #sanstot mean for atypical AP patients
parentses.sd.AP3=sd(merged_ordered$PARENTSES[which(merged_ordered$AP==3)],na.rm=T)   #sanstot sd for atypical AP patients
parentses.range.AP3=range(merged_ordered$PARENTSES[which(merged_ordered$AP==3)],na.rm=T)   #sanstot range for atypical AP patients
parentses.mu.AP4=mean(merged_ordered$PARENTSES[which(merged_ordered$AP==4)],na.rm=T)   #sanstot mean for both a and t AP patients
parentses.sd.AP4=sd(merged_ordered$PARENTSES[which(merged_ordered$AP==4)],na.rm=T)   #sanstot sd for both a and t AP patients
parentses.range.AP4=range(merged_ordered$PARENTSES[which(merged_ordered$AP==4)],na.rm=T)   #sanstot range for both a and t AP patients

n.dx0=length(which(merged_ordered$Dx==0))  # Total number of ctls
n.dx1=length(which(merged_ordered$Dx==1))  # Total number of patients

n.AP0=length(which(merged_ordered$AP==0))   # Total number of controls
n.AP1=length(which(merged_ordered$AP==1))   # Total number of unmedicated patients
n.AP2=length(which(merged_ordered$AP==2))   # Total number of typical AP patients
n.AP3=length(which(merged_ordered$AP==3))   # Total number of atypical AP patients
n.AP4=length(which(merged_ordered$AP==4))   # Total number of both a and t AP patients

n.fem=length(which(merged_ordered$Sex==2))   # Women
n.mal=length(which(merged_ordered$Sex==1))   # Men

#Check that Sex was coded properly
if((n.fem + n.mal) != length(merged_ordered$Sex)){
	stop('Did you remember to code the Sex covariate as Males=1 and Females=2?\n')
}

n.fem.dx0=length(which(merged_ordered$Sex==2 & merged_ordered$Dx==0))   # Women ctls
n.mal.dx0=length(which(merged_ordered$Sex==1 & merged_ordered$Dx==0))   # Men ctls
n.fem.dx1=length(which(merged_ordered$Sex==2 & merged_ordered$Dx==1))   # Women patients
n.mal.dx1=length(which(merged_ordered$Sex==1 & merged_ordered$Dx==1))   # Men patients

n.fem.AP0=length(which(merged_ordered$Sex==2 & merged_ordered$AP==0))   # Women controls
n.mal.AP0=length(which(merged_ordered$Sex==1 & merged_ordered$AP==0))   # Men controls
n.fem.AP1=length(which(merged_ordered$Sex==2 & merged_ordered$AP==1))   # Women unmedicated patients
n.mal.AP1=length(which(merged_ordered$Sex==1 & merged_ordered$AP==1))   # Men unmedicated patients
n.fem.AP2=length(which(merged_ordered$Sex==2 & merged_ordered$AP==2))   # Women typical AP patients
n.mal.AP2=length(which(merged_ordered$Sex==1 & merged_ordered$AP==2))   # Men typical AP patients
n.fem.AP3=length(which(merged_ordered$Sex==2 & merged_ordered$AP==3))   # Women atypical AP patients
n.mal.AP3=length(which(merged_ordered$Sex==1 & merged_ordered$AP==3))   # Men atypical AP patients
n.fem.AP4=length(which(merged_ordered$Sex==2 & merged_ordered$AP==4))   # Women both a and t AP patients
n.mal.AP4=length(which(merged_ordered$Sex==1 & merged_ordered$AP==4))   # Men both a and t AP patients

n.hand0.dx0=length(which(merged_ordered$HAND==0 & merged_ordered$Dx==0))   # Number Right Handed ctls
n.hand0.dx1=length(which(merged_ordered$HAND==0 & merged_ordered$Dx==1))   # Number Right Handed patients
n.hand1.dx0=length(which(merged_ordered$HAND==1 & merged_ordered$Dx==0))   # Number Left Handed ctls
n.hand1.dx1=length(which(merged_ordered$HAND==1 & merged_ordered$Dx==1))   # Number Left Handed patients
n.hand2.dx0=length(which(merged_ordered$HAND==2 & merged_ordered$Dx==0))   # Number Ambidextrousctls
n.hand2.dx1=length(which(merged_ordered$HAND==2 & merged_ordered$Dx==1))   # Number Ambidextrous patients

#Save demographic info
save(age.mu, age.sd, age.range, age.mu.dx0, age.sd.dx0, age.range.dx0, age.mu.dx1, age.sd.dx1, age.range.dx1, age.mu, age.sd, age.range, age.mu.dx0, age.sd.dx0, age.range.dx0, age.mu.dx1, age.sd.dx1, age.range.dx1, age.mu.AP0, age.sd.AP0, age.range.AP0, age.mu.AP1, age.sd.AP1, age.range.AP1, age.mu.AP2, age.sd.AP2, age.range.AP2, age.mu.AP3, age.sd.AP3, age.range.AP3, age.mu.AP4, age.sd.AP4, age.range.AP4, cpz.mu.dx1, cpz.sd.dx1, cpz.range.dx1, cpz.mu.AP1, cpz.sd.AP1, cpz.range.AP1, cpz.mu.AP2, cpz.sd.AP2, cpz.range.AP2, cpz.mu.AP3, cpz.sd.AP3, cpz.range.AP3, cpz.mu.AP4, cpz.sd.AP4, cpz.range.AP4, ao.mu.dx1, ao.sd.dx1, ao.range.dx1, ao.mu.AP1, ao.sd.AP1, ao.range.AP1, ao.mu.AP2, ao.sd.AP2, ao.range.AP2, ao.mu.AP3, ao.sd.AP3, ao.range.AP3, ao.mu.AP4, ao.sd.AP4, ao.range.AP4, durill.mu.dx1, durill.sd.dx1, durill.range.dx1, durill.mu.AP1, durill.sd.AP1, durill.range.AP1, durill.mu.AP2, durill.sd.AP2, durill.range.AP2, durill.mu.AP3, durill.sd.AP3, durill.range.AP3, durill.mu.AP4, durill.sd.AP4, durill.range.AP4, pansstot.mu.dx1, pansstot.sd.dx1, pansstot.range.dx1, pansstot.mu.AP1, pansstot.sd.AP1, pansstot.range.AP1, pansstot.mu.AP2, pansstot.sd.AP2, pansstot.range.AP2, pansstot.mu.AP3, pansstot.sd.AP3, pansstot.range.AP3, pansstot.mu.AP4, pansstot.sd.AP4, pansstot.range.AP4, pansspos.mu.dx1, pansspos.sd.dx1, pansspos.range.dx1, pansspos.mu.AP1, pansspos.sd.AP1, pansspos.range.AP1, pansspos.mu.AP2, pansspos.sd.AP2, pansspos.range.AP2, pansspos.mu.AP3, pansspos.sd.AP3, pansspos.range.AP3, pansspos.mu.AP4, pansspos.sd.AP4, pansspos.range.AP4, panssneg.mu.dx1, panssneg.sd.dx1, panssneg.range.dx1, panssneg.mu.AP1, panssneg.sd.AP1, panssneg.range.AP1, panssneg.mu.AP2, panssneg.sd.AP2, panssneg.range.AP2, panssneg.mu.AP3, panssneg.sd.AP3, panssneg.range.AP3, panssneg.mu.AP4, panssneg.sd.AP4, panssneg.range.AP4, sapstot.mu.dx1, sapstot.sd.dx1, sapstot.range.dx1, sapstot.mu.AP1, sapstot.sd.AP1, sapstot.range.AP1, sapstot.mu.AP2, sapstot.sd.AP2, sapstot.range.AP2, sapstot.mu.AP3, sapstot.sd.AP3, sapstot.range.AP3, sapstot.mu.AP4, sapstot.sd.AP4, sapstot.range.AP4, sanstot.mu.dx1, sanstot.sd.dx1, sanstot.range.dx1, sanstot.mu.AP1, sanstot.sd.AP1, sanstot.range.AP1, sanstot.mu.AP2, sanstot.sd.AP2, sanstot.range.AP2, sanstot.mu.AP3, sanstot.sd.AP3, sanstot.range.AP3, sanstot.mu.AP4, sanstot.sd.AP4, sanstot.range.AP4, parentses.mu.dx0, parentses.sd.dx0, parentses.range.dx0, parentses.n.dx0, parentses.mu.dx1,parentses.sd.dx1,parentses.range.dx1, parentses.n.dx1, parentses.mu.AP1,parentses.sd.AP1,parentses.range.AP1,parentses.mu.AP2,parentses.sd.AP2,parentses.range.AP2,parentses.mu.AP3,parentses.sd.AP3,parentses.range.AP3,parentses.mu.AP4,parentses.sd.AP4,parentses.range.AP4,n.dx0, n.dx1, n.AP0, n.AP1, n.AP2, n.AP3, n.AP4, n.fem, n.mal, n.fem.dx0, n.mal.dx0, n.fem.dx1, n.mal.dx1, n.fem.AP0, n.mal.AP0, n.fem.AP1, n.mal.AP1, n.fem.AP2, n.mal.AP2, n.fem.AP3, n.mal.AP3, n.fem.AP4, n.mal.AP4, n.hand0.dx0, n.hand0.dx1, n.hand1.dx0, n.hand1.dx1, n.hand2.dx0, n.hand2.dx1, iq.mu.dx0, iq.sd.dx0, iq.range.dx0, iq.n.dx0, iq.mu.dx1, iq.sd.dx1, iq.range.dx1, iq.n.dx1, file=paste0(transferDir, "/","Demographics_",filetype,".Rdata"))

cat('Done calculating demographics for ', fsfile,'\n')


### Regression Analysis

###Functions used in the code for Cohens d
d.t.unpaired<-function(t.val,n1,n2){
  d<-t.val*sqrt((n1+n2)/(n1*n2))
  names(d)<-"effect size d"
  return(d)
}

partial.d<-function(t.val,df,n1,n2){
  d<-t.val*(n1+n2)/(sqrt(n1*n2)*sqrt(df))
  names(d)<-"effect size d"
  return(d)
}

CI1<-function(ES,se){
  ci<-c((ES-(1.96)*se),(ES+(1.96)*se))
  names(ci)<-c("95% CI lower","95% CI upper")
  return(ci)
}

se.d2<-function(d,n1,n2){
  se<-sqrt((n1+n2)/(n1*n2)+(d^2)/(2*(n1+n2-2)))
  names(se)<-"se for d"
  return(se)
}
##########################

#attach main file
attach(merged_ordered)


# 1. SZ patients vs controls

#check that there is enough data to proceed with this step
if(nrow(merged_ordered[!is.na(merged_ordered$Dx),]) > (n.sites + 3)){

cat('Running: 1. SZ patients vs controls\n')
#Store models for troubleshooting
models.cort=NULL; # This will become a list where we store all of the models made by lm

#allocate empty vectors to store adjust effect sizes, se, ci (noicv)
d.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
se.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
low.ci.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
up.ci.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
n.controls=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
n.patients=rep(NA,(ncol(merged_ordered)-ncol(Covs)))

#Loop through and perform each regression
for(x in (ncol(Covs)+1):ncol(merged_ordered)){
	pheno=merged_ordered[!is.na(merged_ordered[,x]),x] #Check to make sure there are observations for a given structure
	#check if the phenotype is singular after NA removal
	if(length(pheno)==0){
		next; # Skip the whole structure if there are no observations
	}
	site=NULL # These are just string variables that modify the model in R if there are Site variables in the Covariates.csv file
	if(n.sites != 0){
		site=" + " #Starts the string to add-on Site variables
		for(i in 1:n.sites){
			if(i != n.sites){
				site=paste(site, "Site", i, " + ", sep='')
			}
			else{
				site=paste(site, "Site", i, sep='')
			}
		}
	}
	#Run the model
	test=merged_ordered[,x]
	eval(parse(text=paste("tmp=lm(test ~ factor(Dx) + Age + factor(Sex)", site, ")", sep='')))
	models.cort[[x-ncol(Covs)]]=tmp #Store the model fit for future reference
        	

	#subjects can be dropped if they are missing so we can get the precise number of controls/patients for each region tested
	n.controls[x-ncol(Covs)] = length(which(tmp$model[,2] == 0))
	n.patients[x-ncol(Covs)] = length(which(tmp$model[,2] == 1))
	
	#Convert the lm model to a summary format so we can extract statistics
	tmp=summary(tmp)
	tstat=tmp$coefficients[2,3] # Get t-statistic from regression to convert to Cohens d
	tstat.df=tmp$df[2]
	
	#collect effect size data
	d.cort[x-ncol(Covs)]=partial.d(tstat,tstat.df,n.controls[x-ncol(Covs)],n.patients[x-ncol(Covs)])
	se.cort[x-ncol(Covs)]=se.d2(d.cort[x-ncol(Covs)],n.controls[x-ncol(Covs)],n.patients[x-ncol(Covs)])
	bound.cort=CI1(d.cort[x-ncol(Covs)],se.cort[x-ncol(Covs)])
	low.ci.cort[x-ncol(Covs)]=bound.cort[1]
	up.ci.cort[x-ncol(Covs)]=bound.cort[2]
}

#save results
save(d.cort,se.cort,low.ci.cort,up.ci.cort,n.controls,n.patients, file=paste0(transferDir, "/","EffectSizes_SZvHV_",filetype,".Rdata"))
save(models.cort, file=paste0(transferDir, "/","Models_SZvHV_",filetype,".Rdata"))
} else {cat('NOT Running: 1. SZ patients vs controls\n')}


# 2.1 Diagnosis by Sex interaction (all SZ patients vs controls)

#check that there is enough data to proceed with this step
if(nrow(merged_ordered[!is.na(merged_ordered$Dx),]) > (n.sites + 3)){
cat('Running: 2.1 Diagnosis by Sex interaction (all SZ patients vs controls)\n')
	
#Store models for troubleshooting
models.cort=NULL; # This will become a list where we store all of the models made by lm

#allocate empty vectors to store adjust effect sizes, se, ci (noicv)
d.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
se.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
low.ci.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
up.ci.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
n.controls=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
n.patients=rep(NA,(ncol(merged_ordered)-ncol(Covs)))

#Loop through and perform each regression
for(x in (ncol(Covs)+1):ncol(merged_ordered)){
	pheno=merged_ordered[!is.na(merged_ordered[,x]),x] #Check to make sure there are observations for a given structure
	#check if the phenotype is singular after NA removal
	if(length(pheno)==0){
		next; # Skip the whole structure if there are no observations
	}
	site=NULL # These are just string variables that modify the model in R if there are Site variables in the Covariates.csv file
	if(n.sites != 0){
		site=" + " #Starts the string to add-on Site variables
		for(i in 1:n.sites){
			if(i != n.sites){
				site=paste(site, "Site", i, " + ", sep='')
			}
			else{
				site=paste(site, "Site", i, sep='')
			}
		}
	}
	#Run the model
	test=merged_ordered[,x]
	eval(parse(text=paste("tmp=lm(test ~ factor(Dx) + Age + factor(Sex) + factor(Dx):factor(Sex)", site, ")", sep='')))
	models.cort[[x-ncol(Covs)]]=tmp #Store the model fit for future reference
	
	#subjects can be dropped if they are missing so we can get the precise number of controls/patients for each region tested
	n.controls[x-ncol(Covs)] = length(which(tmp$model[,2] == 0))
	n.patients[x-ncol(Covs)] = length(which(tmp$model[,2] == 1))
	
	#Convert the lm model to a summary format so we can extract statistics
	tmp=summary(tmp)
	tstat=tmp$coefficients[(5 + n.sites),3] # Get t-statistic from regression to convert to Cohens d
	tstat.df=tmp$df[2]
	
	#collect effect size data
	d.cort[x-ncol(Covs)]=partial.d(tstat,tstat.df,n.controls[x-ncol(Covs)],n.patients[x-ncol(Covs)])
	se.cort[x-ncol(Covs)]=se.d2(d.cort[x-ncol(Covs)],n.controls[x-ncol(Covs)],n.patients[x-ncol(Covs)])
	bound.cort=CI1(d.cort[x-ncol(Covs)],se.cort[x-ncol(Covs)])
	low.ci.cort[x-ncol(Covs)]=bound.cort[1]
	up.ci.cort[x-ncol(Covs)]=bound.cort[2]
}

#save results
save(d.cort,se.cort,low.ci.cort,up.ci.cort,n.controls,n.patients, file=paste0(transferDir, "/","EffectSizes_SZvHV_DxBySex_",filetype,".Rdata"))
save(models.cort, file=paste0(transferDir, "/","Models_SZvHV_DxBySex_",filetype,".Rdata"))
} else {cat('NOT Running: 2.1 Diagnosis by Sex interaction (all SZ patients vs controls)\n')}


# 2.2 Diagnosis by Age interaction (all SZ patients vs controls)

#check that there is enough data to proceed with this step
if(nrow(merged_ordered[!is.na(merged_ordered$Dx),]) > (n.sites + 3)){
cat('Running: 2.2 Diagnosis by Age interaction (all SZ patients vs controls)\n')

#Store models for troubleshooting
models.cort=NULL; # This will become a list where we store all of the models made by lm

#allocate empty vectors to store adjust effect sizes, se, ci (noicv)
d.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
se.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
low.ci.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
up.ci.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
n.controls=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
n.patients=rep(NA,(ncol(merged_ordered)-ncol(Covs)))

#Loop through and perform each regression
for(x in (ncol(Covs)+1):ncol(merged_ordered)){
        pheno=merged_ordered[!is.na(merged_ordered[,x]),x] #Check to make sure there are observations for a given structure
        #check if the phenotype is singular after NA removal
        if(length(pheno)==0){
                next; # Skip the whole structure if there are no observations
        }
        site=NULL # These are just string variables that modify the model in R if there are Site variables in the Covariates.csv file
        if(n.sites != 0){
                site=" + " #Starts the string to add-on Site variables
                for(i in 1:n.sites){
                        if(i != n.sites){
                                site=paste(site, "Site", i, " + ", sep='')
                        }
                        else{
                                site=paste(site, "Site", i, sep='')
                        }
                }
        }
        #Run the model
        test=merged_ordered[,x]
        eval(parse(text=paste("tmp=lm(test ~ factor(Dx) + Age + factor(Sex) + factor(Dx):Age", site, ")", sep='')))
        models.cort[[x-ncol(Covs)]]=tmp #Store the model fit for future reference

        #subjects can be dropped if they are missing so we can get the precise number of controls/patients for each region tested
        n.controls[x-ncol(Covs)] = length(which(tmp$model[,2] == 0))
        n.patients[x-ncol(Covs)] = length(which(tmp$model[,2] == 1))

        #Convert the lm model to a summary format so we can extract statistics
        tmp=summary(tmp)
        tstat=tmp$coefficients[(5 + n.sites),3] # Get t-statistic from regression to convert to Cohens d
        tstat.df=tmp$df[2]

        #collect effect size data
        d.cort[x-ncol(Covs)]=partial.d(tstat,tstat.df,n.controls[x-ncol(Covs)],n.patients[x-ncol(Covs)])
        se.cort[x-ncol(Covs)]=se.d2(d.cort[x-ncol(Covs)],n.controls[x-ncol(Covs)],n.patients[x-ncol(Covs)])
        bound.cort=CI1(d.cort[x-ncol(Covs)],se.cort[x-ncol(Covs)])
        low.ci.cort[x-ncol(Covs)]=bound.cort[1]
        up.ci.cort[x-ncol(Covs)]=bound.cort[2]
}

#save results
save(d.cort,se.cort,low.ci.cort,up.ci.cort,n.controls,n.patients, file=paste0(transferDir, "/","EffectSizes_SZvHV_DxByAge_",filetype,".Rdata"))
save(models.cort, file=paste0(transferDir, "/","Models_SZvHV_DxByAge_",filetype,".Rdata"))
} else {cat('NOT Running: 2.2 Diagnosis by Age interaction (all SZ patients vs controls)\n')}


# 3. AntiPsychotic (AP) use (typical, atypical, both, none) comparisons

#detach(merged_ordered)
#rm(merged_ordered)
merged_ordered=swap.merged_ordered

# create data subsets
merged_ordered_AP1vsAP0 <- merged_ordered[((merged_ordered$AP!=2) & (merged_ordered$AP!=3) & (merged_ordered$AP!=4) & !(is.na(merged_ordered$AP))),]
merged_ordered_AP2vsAP0 <- merged_ordered[((merged_ordered$AP!=1) & (merged_ordered$AP!=3) & (merged_ordered$AP!=4) & !(is.na(merged_ordered$AP))),]
merged_ordered_AP3vsAP0 <- merged_ordered[((merged_ordered$AP!=2) & (merged_ordered$AP!=1) & (merged_ordered$AP!=4) & !(is.na(merged_ordered$AP))),]
merged_ordered_AP4vsAP0 <- merged_ordered[((merged_ordered$AP!=2) & (merged_ordered$AP!=3) & (merged_ordered$AP!=1) & !(is.na(merged_ordered$AP))),]
merged_ordered_AP2vsAP1 <- merged_ordered[((merged_ordered$AP!=0) & (merged_ordered$AP!=3) & (merged_ordered$AP!=4) & !(is.na(merged_ordered$AP))),]
merged_ordered_AP3vsAP1 <- merged_ordered[((merged_ordered$AP!=0) & (merged_ordered$AP!=2) & (merged_ordered$AP!=4) & !(is.na(merged_ordered$AP))),]
merged_ordered_AP4vsAP1 <- merged_ordered[((merged_ordered$AP!=0) & (merged_ordered$AP!=2) & (merged_ordered$AP!=3) & !(is.na(merged_ordered$AP))),]
merged_ordered_AP3vsAP2 <- merged_ordered[((merged_ordered$AP!=0) & (merged_ordered$AP!=1) & (merged_ordered$AP!=4) & !(is.na(merged_ordered$AP))),]
merged_ordered_AP4vsAP2 <- merged_ordered[((merged_ordered$AP!=0) & (merged_ordered$AP!=1) & (merged_ordered$AP!=3) & !(is.na(merged_ordered$AP))),]
merged_ordered_AP4vsAP3 <- merged_ordered[((merged_ordered$AP!=0) & (merged_ordered$AP!=1) & (merged_ordered$AP!=2) & !(is.na(merged_ordered$AP))),]

# AntiPsychotic (AP) groups vs controls and vs. eachother

for(APcomparison in c("AP1vsAP0","AP2vsAP0","AP3vsAP0","AP4vsAP0","AP2vsAP1","AP3vsAP1","AP4vsAP1","AP3vsAP2","AP4vsAP2","AP4vsAP3")){ 
#for(APcomparison in c("AP3vsAP0")){
  
   # Make sure there are 2 levels and each level has at least 2 data points: minimum required to compute variance
   # Can decide later on reasonable number per group
   #min_range=min(range(merged_ordered_AP3vsAP0$AP,na.rm=T));
   #max_range=max(range(merged_ordered_AP3vsAP0$AP,na.rm=T));
   #n_min_range=length(which(merged_ordered_AP3vsAP0$AP==min_range));
   #n_max_range=length(which(merged_ordered_AP3vsAP0$AP==max_range));
   #if( (min_range != max_range) & n_min_range>2 & n_max_range>2 ) {   
 
   #Store models for troubleshooting
   models.cort=NULL; # This will become a list where we store all of the models made by lm

   #allocate empty vectors to store adjust effect sizes, se, ci (noicv)
   d.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   se.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   low.ci.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   up.ci.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))

   #n.controls=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   #n.patients=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   n.APlow=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   n.APhigh=rep(NA,(ncol(merged_ordered)-ncol(Covs)))

   #attach main file
   detach(merged_ordered)
   data_subset_name=paste0("merged_ordered",sep="_",APcomparison)
   merged_ordered=get(data_subset_name)
   #merged_ordered=swap.merged_ordered
   #merged_ordered <- merged_ordered[((merged_ordered$AD!=1) & !(is.na(merged_ordered$AD))),]
   attach(merged_ordered)


   # Make sure there are 2 levels and each level has at least 2 data points: minimum required to compute variance
   # Can decide later on reasonable number per group
   min_range=min(range(merged_ordered$AP,na.rm=T));
   max_range=max(range(merged_ordered$AP,na.rm=T));
   n_min_range=length(which(merged_ordered$AP==min_range));
   n_max_range=length(which(merged_ordered$AP==max_range));
   if( (min_range != max_range) & n_min_range>3 & n_max_range>3 ) {

   if(nrow(merged_ordered) > (n.sites + 3)){
      cat('Running: 3. Antipsychotic (AP) Group Comparison', APcomparison,'\n')
	
   #Loop through and perform each regression
   for(x in (ncol(Covs)+1):ncol(merged_ordered)){
	pheno=merged_ordered[!is.na(merged_ordered[,x]),x] #Check to make sure there are observations for a given structure
	#check if the phenotype is singular after NA removal
	if(length(pheno)==0){
		next; # Skip the whole structure if there are no observations
	}
	site=NULL # These are just string variables that modify the model in R if there are Site variables in the Covariates.csv file
	if(n.sites != 0){
		site=" + " #Starts the string to add-on Site variables
		for(i in 1:n.sites){
			if(i != n.sites){
				site=paste(site, "Site", i, " + ", sep='')
			}
			else{
				site=paste(site, "Site", i, sep='')
			}
		}
	}
	#Run the model
	test=merged_ordered[,x]
	eval(parse(text=paste("tmp=lm(test ~ factor(AP) + Age + Sex", site, ")", sep='')))
	models.cort[[x-ncol(Covs)]]=tmp #Store the model fit for future reference
	
	#subjects can be dropped if they are missing so we can get the precise number of controls/patients for each region tested
	#n.controls[x-ncol(Covs)] = length(which(tmp$model[,2] == 0))
	#n.patients[x-ncol(Covs)] = length(which(tmp$model[,2] == 1))

        n.APlow[x-ncol(Covs)] = length(which(tmp$model[,2] == min_range))
        n.APhigh[x-ncol(Covs)] = length(which(tmp$model[,2] == max_range))
	
	#Convert the lm model to a summary format so we can extract statistics
	tmp=summary(tmp)
	tstat=tmp$coefficients[2,3] # Get t-statistic from regression to convert to Cohens d
	tstat.df=tmp$df[2]
	
	#collect effect size data
	#d.cort[x-ncol(Covs)]=partial.d(tstat,tstat.df,n.controls[x-ncol(Covs)],n.patients[x-ncol(Covs)])
	#se.cort[x-ncol(Covs)]=se.d2(d.cort[x-ncol(Covs)],n.controls[x-ncol(Covs)],n.patients[x-ncol(Covs)])
        # APhigh, as indicated by subjects with AP==max_range is the experimental condition
        # APlow, as indicated by subjects with AP==max_range is the control condition
        d.cort[x-ncol(Covs)]=partial.d(tstat,tstat.df,n.APlow[x-ncol(Covs)],n.APhigh[x-ncol(Covs)])
        se.cort[x-ncol(Covs)]=se.d2(d.cort[x-ncol(Covs)],n.APlow[x-ncol(Covs)],n.APhigh[x-ncol(Covs)])
	bound.cort=CI1(d.cort[x-ncol(Covs)],se.cort[x-ncol(Covs)])
	low.ci.cort[x-ncol(Covs)]=bound.cort[1]
	up.ci.cort[x-ncol(Covs)]=bound.cort[2]
   }

   #save results
   save(d.cort,se.cort,low.ci.cort,up.ci.cort,n.APlow,n.APhigh, file=paste0(transferDir, "/","EffectSizes_",APcomparison,sep="_",filetype,".Rdata"))
   save(models.cort, file=paste0(transferDir, "/","Models_AP_",APcomparison,sep="_",filetype,".Rdata"))

   } else {cat('NOT Running: 3. Antipsychotic (AP) Group Comparison', APcomparison,'\n')}

   } # end if min and max range not the same and n_min_range and n_max_range > 0

} # end for loop APcomparison


# 4. Regressions predicting cortical Thickness and Surface Area with continuous variables available for SZ patients only covary for Age #####ppcor

merged_ordered=swap.merged_ordered

# create data subsets
merged_ordered_CPZ <- merged_ordered[((merged_ordered$Dx==1) & !(is.na(merged_ordered$CPZ))),];
merged_ordered_AO <- merged_ordered[((merged_ordered$Dx==1) & !(is.na(merged_ordered$AO))),];
merged_ordered_DURILL <- merged_ordered[((merged_ordered$Dx==1) & !(is.na(merged_ordered$DURILL))),];
merged_ordered_PANSSTOT <- merged_ordered[((merged_ordered$Dx==1) & !(is.na(merged_ordered$PANSSTOT))),];
merged_ordered_PANSSPOS <- merged_ordered[((merged_ordered$Dx==1) & !(is.na(merged_ordered$PANSSPOS))),];
merged_ordered_PANSSNEG <- merged_ordered[((merged_ordered$Dx==1) & !(is.na(merged_ordered$PANSSNEG))),];
merged_ordered_SAPSTOT <- merged_ordered[((merged_ordered$Dx==1) & !(is.na(merged_ordered$SAPSTOT))),];
merged_ordered_SANSTOT <- merged_ordered[((merged_ordered$Dx==1) & !(is.na(merged_ordered$SANSTOT))),];

for(predictor in c("CPZ","AO","DURILL","PANSSTOT","PANSSPOS","PANSSNEG","SAPSTOT","SANSTOT")){

   #Store models for troubleshooting
   models.cort=NULL; # This will become a list where we store all of the models made by lm

   #allocate empty vectors to store adjust effect sizes, se, ci (noicv)
   r.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   se.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   low.ci.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   up.ci.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   n.controls=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   n.patients=rep(NA,(ncol(merged_ordered)-ncol(Covs)))

   #attach main file
   detach(merged_ordered)
   data_subset_name=paste0("merged_ordered",sep="_",predictor)
   merged_ordered=get(data_subset_name)
   attach(merged_ordered)

   if(nrow(merged_ordered) > (n.sites + 3)){
   cat('Running: 4. Regression predictor', predictor, ' in SZ patients covary for Age\n')
	
   #Loop through and perform each regression
   for(x in (ncol(Covs)+1):ncol(merged_ordered)){
	pheno=merged_ordered[!is.na(merged_ordered[,x]),x] #Check to make sure there are observations for a given structure
	#check if the phenotype is singular after NA removal
	if(length(pheno)==0){
		next; # Skip the whole structure if there are no observations
	}
	site=NULL # These are just string variables that modify the model in R if there are Site variables in the Covariates.csv file
	if(n.sites != 0){
		site=" + " #Starts the string to add-on Site variables
		for(i in 1:n.sites){
			if(i != n.sites){
				site=paste(site, "Site", i, " + ", sep='')
			}
			else{
				site=paste(site, "Site", i, sep='')
			}
		}
	}
	#Run the model
	test=merged_ordered[,x]
	eval(parse(text=paste("tmp=lm(test ~ ", predictor," + Age + Sex", site, ")", sep='')))
	models.cort[[x-ncol(Covs)]]=tmp #Store the model fit for future reference
	
	#subjects can be dropped if they are missing so we can get the precise number of controls/patients for each region tested
	n.controls[x-ncol(Covs)] = length(which(tmp$model[,2] == 0))
	n.patients[x-ncol(Covs)] = length(which(tmp$model[,2] == 1))
	
	partcor.i <- pcor.test(tmp$model[1],tmp$model[,2],tmp$model[,c(3:ncol(tmp$model))])	
	r.cort[x-ncol(Covs)]=partcor.i[,1]
   }

   #save results
   save(r.cort,se.cort,low.ci.cort,up.ci.cort,n.controls,n.patients, file=paste0(transferDir, "/","EffectSizes_SZ_only_",predictor,"_withAge_",filetype,".Rdata"))
   save(models.cort, file=paste0(transferDir, "/","Models_SZ_only_",predictor,"_withAge_",filetype,".Rdata"))
   } else {cat('NOT Running: 4. Regression predictor ', predictor, ' in SZ patients covary for age\n')}

} # end for loop predictor


# 4a. Regressions/Partial Correlation predicting cortical Thickness and Surface Area with CPZ for SZ patients only covary for Age and Negative Symptom Severity#####ppcor
#     This also runs a negative symptom regression based on SANSTOT or PANSSNEG converted SANSTOT based on Van Erp et al. 2014

# create data subset
#detach(merged_ordered)
#detach(merged_ordered_CPZ)
merged_ordered_CPZ <- swap.merged_ordered[((swap.merged_ordered$Dx==1) & !(is.na(swap.merged_ordered$CPZ))),]
merged_ordered_SANSTOT_CONVERTEASY=merged_ordered_CPZ
#attach(merged_ordered_SANSTOT_CONVERTEASY)
# if SANSTOT is available, run SANSTOT_CONVERTEASY and CPZ analysis controlling for Negative Symptom Severity with available SANSTOT
# if SANTOT is not avaiable but PANSS_NEG is available, run ANSTOT_CONVERTEASY and CPZ analysis controlling for Negative Symptom Severity with PANSS_TOT converted to SANSTOT based on Van Erp et al. 2014
n.SANSTOT=length(which(!(is.na(merged_ordered_SANSTOT_CONVERTEASY$SANSTOT))))
n.PANSSNEG=length(which(!(is.na(merged_ordered_SANSTOT_CONVERTEASY$PANSSNEG))))

# if no SANSTOT and no PANSSNEG, set merged_ordered_SANSTOT_CONVERTEASY$SANSTOT_CONVERTEASY to merged_ordered_SANSTOT_CONVERTEASY$SANSTOT, which is all NA
if(n.SANSTOT==0 & n.PANSSNEG==0){merged_ordered_SANSTOT_CONVERTEASY$SANSTOT_CONVERTEASY=merged_ordered_SANSTOT_CONVERTEASY$SANSTOT}
# if SANSTOT is avaible, set merged_ordered_SANSTOT_CONVERTEASY$SANSTOT_CONVERTEASY to merged_ordered_SANSTOT_CONVERTEASY$SANSTOT
if (n.SANSTOT>0){merged_ordered_SANSTOT_CONVERTEASY$SANSTOT_CONVERTEASY=merged_ordered_SANSTOT_CONVERTEASY$SANSTOT}
# if SANSTOT is NOT avaible but PANSSNEG is available, set merged_ordered_SANSTOT_CONVERTEASY$SANSTOT_CONVERTEASY to converted merged_ordered_SANSTOT_CONVERTEASY$PANSSNEG
if (n.SANSTOT==0 & n.PANSSNEG>0){merged_ordered_SANSTOT_CONVERTEASY$SANSTOT_CONVERTEASY=((merged_ordered_SANSTOT_CONVERTEASY$PANSSNEG*2.1002) - 8.376)}

# create data set for CPZ analysis controling for Negative Symptom Severity
merged_ordered_CPZ_SANSTOT_CONVERTEASY=merged_ordered_SANSTOT_CONVERTEASY
# Making sure no NA for SANSTOT_CONVERTEASY to reduced output
#merged_ordered_CPZ_SANSTOT_CONVERTEASY<-merged_ordered_SANSTOT_CONVERTEASY[( !is.na(merged_ordered_SANSTOT_CONVERTEASY$SANSTOT_CONVERTEASY) ),];
#merged_ordered_SANSTOT_CONVERTEASY=merged_ordered_CPZ_SANSTOT_CONVERTEASY;

for(predictor in c("SANSTOT_CONVERTEASY","CPZ_SANSTOT_CONVERTEASY")){
  
  #Store models for troubleshooting
  models.cort=NULL; # This will become a list where we store all of the models made by lm
  
  #allocate empty vectors to store adjust effect sizes, se, ci (noicv)
  #r.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
  #se.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
  #low.ci.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
  #up.ci.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
  #n.controls=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
  #n.patients=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
  # CAREFUL: I HAVE NOW CREATED A NEW VARIABLE SANSTOT_CONVERTEASY at the end, so I need to reduce the size of ncol(merged_ordered) by 1
  r.cort=rep(NA,((ncol(merged_ordered)-1)-ncol(Covs)))
  se.cort=rep(NA,((ncol(merged_ordered)-1)-ncol(Covs)))
  low.ci.cort=rep(NA,((ncol(merged_ordered)-1)-ncol(Covs)))
  up.ci.cort=rep(NA,((ncol(merged_ordered)-1)-ncol(Covs)))
  n.controls=rep(NA,((ncol(merged_ordered)-1)-ncol(Covs)))
  n.patients=rep(NA,((ncol(merged_ordered)-1)-ncol(Covs))) 
    
  #attach main file
  detach(merged_ordered)
  data_subset_name=paste0("merged_ordered",sep="_",predictor)
  merged_ordered=get(data_subset_name)
  attach(merged_ordered)
  
  if(nrow(merged_ordered) > (n.sites + 3)){
    cat('Running: 4a. Regression predictor', predictor, ' in SZ patients covary for Age\n')
    
    #Loop through and perform each regression
    # CAREFUL: I HAVE NOW CREATED A NEW VARIABLE SANSTOT_CONVERTEASY at the end, so I need to reduce the size of the loop by 1
    #for(x in (ncol(Covs)+1):ncol(merged_ordered)){
    for(x in (ncol(Covs)+1):(ncol(merged_ordered)-1)){
      pheno=merged_ordered[!is.na(merged_ordered[,x]),x] #Check to make sure there are observations for a given structure
      #check if the phenotype is singular after NA removal
      if(length(pheno)==0){
        next; # Skip the whole structure if there are no observations
      }
      site=NULL # These are just string variables that modify the model in R if there are Site variables in the Covariates.csv file
      if(n.sites != 0){
        site=" + " #Starts the string to add-on Site variables
        for(i in 1:n.sites){
          if(i != n.sites){
            site=paste(site, "Site", i, " + ", sep='')
          }
          else{
            site=paste(site, "Site", i, sep='')
          }
        }
      }
      
      # TVE, note, important to code SANSTOT_CONVERTEASY after CPZ as CPZ is the variable of interest
      if(predictor=="CPZ_SANSTOT_CONVERTEASY"){ predictors=c(" CPZ + SANSTOT_CONVERTEASY ") }
      else { predictors=c(" SANSTOT_CONVERTEASY ") }
      #predictors=c(" CPZ ")
      #predictors=c(" CPZ + SANSTOT_CONVERTEASY ")
      ##predictors=c(" SANSTOT_CONVERTEASY ") # looks like the variable SANSTOT_CONVERTEASY does not exist?!
      #predictors=c(" SANSTOT ") # looks like the variable SANSTOT does not exist?!
      
      #Run the model
      test=merged_ordered[,x]
      eval(parse(text=paste("tmp=lm(test ~ ", predictors ," + Age + Sex", site, ")", sep='')))
      #eval(parse(text=paste("tmp=lm(test ~ CPZ + Age + Sex", site, ")", sep='')))
      models.cort[[x-ncol(Covs)]]=tmp #Store the model fit for future reference
      
      #subjects can be dropped if they are missing so we can get the precise number of controls/patients for each region tested
      n.controls[x-ncol(Covs)] = length(which(tmp$model[,2] == 0))
      n.patients[x-ncol(Covs)] = length(which(tmp$model[,2] == 1))
      
      partcor.i <- pcor.test(tmp$model[1],tmp$model[,2],tmp$model[,c(3:ncol(tmp$model))])	
      r.cort[x-ncol(Covs)]=partcor.i[,1]
    }
    
    #save results
    save(r.cort,se.cort,low.ci.cort,up.ci.cort,n.controls,n.patients, file=paste0(transferDir, "/","EffectSizes_SZ_only_",predictor,"_withAge_",filetype,".Rdata"))
    save(models.cort, file=paste0(transferDir, "/","Models_SZ_only_",predictor,"_withAge_",filetype,".Rdata"))
  } else {cat('NOT Running: 4a. Regression predictor ', predictor, ' in SZ patients covary for age\n')}
  
} # end for loop predictor 


# 5. Regressions predicting cortical Thickness and Surface Area with continuous variables available for SZ patients only without covary for Age #####ppcor

# Data sets already created in Step 4.
#merged_ordered=swap.merged_ordered

## create data subsets
merged_ordered_CPZ <- merged_ordered[((merged_ordered$Dx==1) & !(is.na(merged_ordered$CPZ))),];
#merged_ordered_AO <- merged_ordered[((merged_ordered$Dx==1) & !(is.na(merged_ordered$AO))),];
#merged_ordered_DURILL <- merged_ordered[((merged_ordered$Dx==1) & !(is.na(merged_ordered$DURILL))),];
#merged_ordered_PANSSTOT <- merged_ordered[((merged_ordered$Dx==1) & !(is.na(merged_ordered$PANSSTOT))),];
#merged_ordered_PANSSPOS <- merged_ordered[((merged_ordered$Dx==1) & !(is.na(merged_ordered$PANSSPOS))),];
#merged_ordered_PANSSNEG <- merged_ordered[((merged_ordered$Dx==1) & !(is.na(merged_ordered$PANSSNEG))),];
#merged_ordered_SAPSTOT <- merged_ordered[((merged_ordered$Dx==1) & !(is.na(merged_ordered$SAPSTOT))),];
#merged_ordered_SANSTOT <- merged_ordered[((merged_ordered$Dx==1) & !(is.na(merged_ordered$SANSTOT))),];

for(predictor in c("CPZ","AO","DURILL","PANSSTOT","PANSSPOS","PANSSNEG","SAPSTOT","SANSTOT")){

   #Store models for troubleshooting
   models.cort=NULL; # This will become a list where we store all of the models made by lm

   #allocate empty vectors to store adjust effect sizes, se, ci (noicv)
   r.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   se.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   low.ci.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   up.ci.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   n.controls=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   n.patients=rep(NA,(ncol(merged_ordered)-ncol(Covs)))

   #attach main file
   detach(merged_ordered)
   data_subset_name=paste0("merged_ordered",sep="_",predictor)
   merged_ordered=get(data_subset_name)
   attach(merged_ordered)

   if(nrow(merged_ordered) > (n.sites + 3)){
   cat('Running: 5. Regression predictor', predictor, ' in SZ patients without covary for Age\n')

   #Loop through and perform each regression
   for(x in (ncol(Covs)+1):ncol(merged_ordered)){
        pheno=merged_ordered[!is.na(merged_ordered[,x]),x] #Check to make sure there are observations for a given structure
        #check if the phenotype is singular after NA removal
        if(length(pheno)==0){
                next; # Skip the whole structure if there are no observations
        }
        site=NULL # These are just string variables that modify the model in R if there are Site variables in the Covariates.csv file
        if(n.sites != 0){
                site=" + " #Starts the string to add-on Site variables
                for(i in 1:n.sites){
                        if(i != n.sites){
                                site=paste(site, "Site", i, " + ", sep='')
                        }
                        else{
                                site=paste(site, "Site", i, sep='')
                        }
                }
        }
        #Run the model
        test=merged_ordered[,x]
        eval(parse(text=paste("tmp=lm(test ~ ", predictor," + Sex", site, ")", sep='')))
        models.cort[[x-ncol(Covs)]]=tmp #Store the model fit for future reference

        #subjects can be dropped if they are missing so we can get the precise number of controls/patients for each region tested
        n.controls[x-ncol(Covs)] = length(which(tmp$model[,2] == 0))
        n.patients[x-ncol(Covs)] = length(which(tmp$model[,2] == 1))

        partcor.i <- pcor.test(tmp$model[1],tmp$model[,2],tmp$model[,c(3:ncol(tmp$model))])
        r.cort[x-ncol(Covs)]=partcor.i[,1]
   }

   #save results
   save(r.cort,se.cort,low.ci.cort,up.ci.cort,n.controls,n.patients, file=paste0(transferDir, "/","EffectSizes_SZ_only_",predictor,"_withoutAge_",filetype,".Rdata"))
   save(models.cort, file=paste0(transferDir, "/","Models_SZ_only_",predictor,"_withoutAge_",filetype,".Rdata"))
   } else {cat('NOT Running: 5. Regression predictor ', predictor, ' in SZ patients without covary for age\n')}

} # end for loop predictor


# 6. Regressions predicting cortical Thickness and Surface Area with continuous variables available for SZ patients and Controls covary for Age #####ppcor

merged_ordered=swap.merged_ordered

# create data subsets
merged_ordered_SZ_only_IQ <- merged_ordered[((merged_ordered$Dx==1) & !(is.na(merged_ordered$IQ))),];
merged_ordered_HV_only_IQ <- merged_ordered[((merged_ordered$Dx==0) & !(is.na(merged_ordered$IQ))),];
merged_ordered_IQ <- merged_ordered[!(is.na(merged_ordered$IQ)),];

for(predictor in c("IQ","SZ_only_IQ","HV_only_IQ")){

   #Store models for troubleshooting
   models.cort=NULL; # This will become a list where we store all of the models made by lm

   #allocate empty vectors to store adjust effect sizes, se, ci (noicv)
   r.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   se.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   low.ci.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   up.ci.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   n.controls=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   n.patients=rep(NA,(ncol(merged_ordered)-ncol(Covs)))

   #attach main file
   detach(merged_ordered)
   data_subset_name=paste0("merged_ordered",sep="_",predictor)
   merged_ordered=get(data_subset_name)
   attach(merged_ordered)

   if(nrow(merged_ordered) > (n.sites + 3)){
   cat('Running: 6. Regression predictor', predictor, ' covary for Age\n')

   #Loop through and perform each regression
   for(x in (ncol(Covs)+1):ncol(merged_ordered)){
        pheno=merged_ordered[!is.na(merged_ordered[,x]),x] #Check to make sure there are observations for a given structure
        #check if the phenotype is singular after NA removal
        if(length(pheno)==0){
                next; # Skip the whole structure if there are no observations
        }
        site=NULL # These are just string variables that modify the model in R if there are Site variables in the Covariates.csv file
        if(n.sites != 0){
                site=" + " #Starts the string to add-on Site variables
                for(i in 1:n.sites){
                        if(i != n.sites){
                                site=paste(site, "Site", i, " + ", sep='')
                        }
                        else{
                                site=paste(site, "Site", i, sep='')
                        }
                }
        }

        # TVE, note, important to code Dx after IQ as IQ is the variable of interest
        if(predictor=="IQ"){predictors=" + IQ + Dx"}
        else { predictors=" + IQ" }

        #Run the model
        test=merged_ordered[,x]
        eval(parse(text=paste("tmp=lm(test ~ ", predictors," + Age + Sex", site, ")", sep='')))
        models.cort[[x-ncol(Covs)]]=tmp #Store the model fit for future reference

        #subjects can be dropped if they are missing so we can get the precise number of controls/patients for each region tested
        n.controls[x-ncol(Covs)] = length(which(tmp$model[,2] == 0))
        n.patients[x-ncol(Covs)] = length(which(tmp$model[,2] == 1))

        partcor.i <- pcor.test(tmp$model[1],tmp$model[,2],tmp$model[,c(3:ncol(tmp$model))])
        r.cort[x-ncol(Covs)]=partcor.i[,1]
   }

   #save results
   save(r.cort,se.cort,low.ci.cort,up.ci.cort,n.controls,n.patients, file=paste0(transferDir, "/","EffectSizes_",predictor,"_withAge_",filetype,".Rdata"))
   save(models.cort, file=paste0(transferDir, "/","Models_",predictor,"_withAge_",filetype,".Rdata"))
   } else {cat('NOT Running: 6. Regression predictor ', predictor, ' covary for age\n')}

} # end for loop predictor


# 7. Regressions predicting cortical Thickness and Surface Area with Age in SZ patients and Controls separately  (to investigate DX by Age interaction) #####

merged_ordered=swap.merged_ordered

# create data subsets
merged_ordered_SZ_only_Age <- merged_ordered[((merged_ordered$Dx==1) & !(is.na(merged_ordered$Age))),];
merged_ordered_HV_only_Age <- merged_ordered[((merged_ordered$Dx==0) & !(is.na(merged_ordered$Age))),];

for(predictor in c("SZ_only_Age","HV_only_Age")){
  
  #Store models for troubleshooting
  models.cort=NULL; # This will become a list where we store all of the models made by lm
  
  #allocate empty vectors to store adjust effect sizes, se, ci (noicv)
  r.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
  se.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
  low.ci.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
  up.ci.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
  n.males=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
  n.females=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
  
  #attach main file
  detach(merged_ordered)
  data_subset_name=paste0("merged_ordered",sep="_",predictor)
  merged_ordered=get(data_subset_name)
  attach(merged_ordered)
  
  if(nrow(merged_ordered) > (n.sites + 3)){
    cat('Running: 7. Regression predictor', predictor, ' Effect of Age Covary for Sex\n')
    
    #Loop through and perform each regression
    for(x in (ncol(Covs)+1):ncol(merged_ordered)){
      pheno=merged_ordered[!is.na(merged_ordered[,x]),x] #Check to make sure there are observations for a given structure
      #check if the phenotype is singular after NA removal
      if(length(pheno)==0){
        next; # Skip the whole structure if there are no observations
      }
      site=NULL # These are just string variables that modify the model in R if there are Site variables in the Covariates.csv file
      if(n.sites != 0){
        site=" + " #Starts the string to add-on Site variables
        for(i in 1:n.sites){
          if(i != n.sites){
            site=paste(site, "Site", i, " + ", sep='')
          }
          else{
            site=paste(site, "Site", i, sep='')
          }
        }
      }
      
      # TVE, note, important to code Age as first predictor as Age the variable of interest
      if(predictor=="SZ_only_Age"){predictors="+ Age"}
      else { predictors="+ Age" }
      
      #Run the model
      test=merged_ordered[,x]
      eval(parse(text=paste("tmp=lm(test ~ ", predictors," + Sex ", site, ")", sep='')))
      models.cort[[x-ncol(Covs)]]=tmp #Store the model fit for future reference
      
      #subjects can be dropped if they are missing so we can get the precise number of controls/patients for each region tested
      n.males[x-ncol(Covs)] = length(which(tmp$model[,2] == 1))
      n.females[x-ncol(Covs)] = length(which(tmp$model[,2] == 2))
      
      partcor.i <- pcor.test(tmp$model[1],tmp$model[,2],tmp$model[,c(3:ncol(tmp$model))])
      r.cort[x-ncol(Covs)]=partcor.i[,1]
    }
    
    #save results
    save(r.cort,se.cort,low.ci.cort,up.ci.cort,n.males,n.females, file=paste0(transferDir, "/","EffectSizes_",predictor,"_withSex_",filetype,".Rdata"))
    save(models.cort, file=paste0(transferDir, "/","Models_",predictor,"_withSex_",filetype,".Rdata"))
  } else {cat('NOT Running: 7. Regression predictor ', predictor, ' covary for sex\n')}
  
} # end for loop predictor



cat('Rerunning models controlling for either average thickness over the cortex or surface area over the cortex\n')

##################
#######################
if(fsfile == file1){
        whole=" + MThickness"
} else {
        whole=" + FullSurfArea"
}

detach(merged_ordered)
merged_ordered=swap.merged_ordered
#attach main file
attach(merged_ordered)


# 1. SZ patients vs controls

#check that there is enough data to proceed with this step
if(nrow(merged_ordered[!is.na(merged_ordered$Dx),]) > (n.sites + 3)){

cat('Running: 1. SZ patients vs controls control for ', whole, '\n')
#Store models for troubleshooting
models.cort=NULL; # This will become a list where we store all of the models made by lm

#allocate empty vectors to store adjust effect sizes, se, ci (noicv)
d.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
se.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
low.ci.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
up.ci.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
n.controls=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
n.patients=rep(NA,(ncol(merged_ordered)-ncol(Covs)))

#Loop through and perform each regression
for(x in (ncol(Covs)+1):ncol(merged_ordered)){
        pheno=merged_ordered[!is.na(merged_ordered[,x]),x] #Check to make sure there are observations for a given structure
        #check if the phenotype is singular after NA removal
        if(length(pheno)==0){
                next; # Skip the whole structure if there are no observations
        }
        site=NULL # These are just string variables that modify the model in R if there are Site variables in the Covariates.csv file
        if(n.sites != 0){
                site=" + " #Starts the string to add-on Site variables
                for(i in 1:n.sites){
                        if(i != n.sites){
                                site=paste(site, "Site", i, " + ", sep='')
                        }
                        else{
                                site=paste(site, "Site", i, sep='')
                        }
                }
        }
        #Run the model
        test=merged_ordered[,x]
        eval(parse(text=paste("tmp=lm(test ~ factor(Dx) + Age + factor(Sex)", whole, site, ")", sep='')))
        models.cort[[x-ncol(Covs)]]=tmp #Store the model fit for future reference

        #subjects can be dropped if they are missing so we can get the precise number of controls/patients for each region tested
        n.controls[x-ncol(Covs)] = length(which(tmp$model[,2] == 0))
        n.patients[x-ncol(Covs)] = length(which(tmp$model[,2] == 1))

        #Convert the lm model to a summary format so we can extract statistics
        tmp=summary(tmp)
        tstat=tmp$coefficients[2,3] # Get t-statistic from regression to convert to Cohens d
        tstat.df=tmp$df[2]

        #collect effect size data
        d.cort[x-ncol(Covs)]=partial.d(tstat,tstat.df,n.controls[x-ncol(Covs)],n.patients[x-ncol(Covs)])
        se.cort[x-ncol(Covs)]=se.d2(d.cort[x-ncol(Covs)],n.controls[x-ncol(Covs)],n.patients[x-ncol(Covs)])
        bound.cort=CI1(d.cort[x-ncol(Covs)],se.cort[x-ncol(Covs)])
        low.ci.cort[x-ncol(Covs)]=bound.cort[1]
        up.ci.cort[x-ncol(Covs)]=bound.cort[2]
}

#save results
save(d.cort,se.cort,low.ci.cort,up.ci.cort,n.controls,n.patients, file=paste0(transferDir, "/","EffectSizes_SZvHV_",filetype,"_ctlWHOLE.Rdata"))
save(models.cort, file=paste0(transferDir, "/","Models_SZvHV_",filetype,"_ctlWHOLE.Rdata"))
} else {cat('NOT Running: 1. SZ patients vs controls control for ', whole, '\n')}


# 2.1 Diagnosis by Sex interaction (all SZ patients vs controls)

#check that there is enough data to proceed with this step
if(nrow(merged_ordered[!is.na(merged_ordered$Dx),]) > (n.sites + 3)){
cat('Running: 2.1 Diagnosis by Sex interaction (all SZ patients vs controls) control for ', whole, '\n')

#Store models for troubleshooting
models.cort=NULL; # This will become a list where we store all of the models made by lm

#allocate empty vectors to store adjust effect sizes, se, ci (noicv)
d.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
se.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
low.ci.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
up.ci.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
n.controls=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
n.patients=rep(NA,(ncol(merged_ordered)-ncol(Covs)))

#Loop through and perform each regression
for(x in (ncol(Covs)+1):ncol(merged_ordered)){
        pheno=merged_ordered[!is.na(merged_ordered[,x]),x] #Check to make sure there are observations for a given structure
        #check if the phenotype is singular after NA removal
        if(length(pheno)==0){
                next; # Skip the whole structure if there are no observations
        }
        site=NULL # These are just string variables that modify the model in R if there are Site variables in the Covariates.csv file
        if(n.sites != 0){
                site=" + " #Starts the string to add-on Site variables
                for(i in 1:n.sites){
                        if(i != n.sites){
                                site=paste(site, "Site", i, " + ", sep='')
                        }
                        else{
                                site=paste(site, "Site", i, sep='')
                        }
                }
        }
        #Run the model
        test=merged_ordered[,x]
        eval(parse(text=paste("tmp=lm(test ~ factor(Dx) + Age + factor(Sex) + factor(Dx):factor(Sex)", whole, site, ")", sep='')))
        models.cort[[x-ncol(Covs)]]=tmp #Store the model fit for future reference

        #subjects can be dropped if they are missing so we can get the precise number of controls/patients for each region tested
        n.controls[x-ncol(Covs)] = length(which(tmp$model[,2] == 0))
        n.patients[x-ncol(Covs)] = length(which(tmp$model[,2] == 1))

        #Convert the lm model to a summary format so we can extract statistics
        tmp=summary(tmp)
        tstat=tmp$coefficients[(6 + n.sites),3] # Get t-statistic from regression to convert to Cohens d
        tstat.df=tmp$df[2]

        #collect effect size data
        d.cort[x-ncol(Covs)]=partial.d(tstat,tstat.df,n.controls[x-ncol(Covs)],n.patients[x-ncol(Covs)])
        se.cort[x-ncol(Covs)]=se.d2(d.cort[x-ncol(Covs)],n.controls[x-ncol(Covs)],n.patients[x-ncol(Covs)])
        bound.cort=CI1(d.cort[x-ncol(Covs)],se.cort[x-ncol(Covs)])
        low.ci.cort[x-ncol(Covs)]=bound.cort[1]
        up.ci.cort[x-ncol(Covs)]=bound.cort[2]
}

#save results
save(d.cort,se.cort,low.ci.cort,up.ci.cort,n.controls,n.patients, file=paste0(transferDir, "/","EffectSizes_SZvHV_DxBySex_",filetype,"_ctlWHOLE.Rdata"))
save(models.cort, file=paste0(transferDir, "/","Models_SZvHV_DxBySex_",filetype,"_ctlWHOLE.Rdata"))
} else {cat('NOT Running: 2.1 Diagnosis by Sex interaction (all SZ patients vs controls) control for ', whole, '\n')}


# 2.2 Diagnosis by Sex interaction (all SZ patients vs controls)

#check that there is enough data to proceed with this step
if(nrow(merged_ordered[!is.na(merged_ordered$Dx),]) > (n.sites + 3)){
cat('Running: 2.2 Diagnosis by Age interaction (all SZ patients vs controls) control for ', whole, '\n')

#Store models for troubleshooting
models.cort=NULL; # This will become a list where we store all of the models made by lm

#allocate empty vectors to store adjust effect sizes, se, ci (noicv)
d.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
se.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
low.ci.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
up.ci.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
n.controls=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
n.patients=rep(NA,(ncol(merged_ordered)-ncol(Covs)))

#Loop through and perform each regression
for(x in (ncol(Covs)+1):ncol(merged_ordered)){
        pheno=merged_ordered[!is.na(merged_ordered[,x]),x] #Check to make sure there are observations for a given structure
        #check if the phenotype is singular after NA removal
        if(length(pheno)==0){
                next; # Skip the whole structure if there are no observations
        }
        site=NULL # These are just string variables that modify the model in R if there are Site variables in the Covariates.csv file
        if(n.sites != 0){
                site=" + " #Starts the string to add-on Site variables
                for(i in 1:n.sites){
                        if(i != n.sites){
                                site=paste(site, "Site", i, " + ", sep='')
                        }
                        else{
                                site=paste(site, "Site", i, sep='')
                        }
                }
        }
        #Run the model
        test=merged_ordered[,x]
        eval(parse(text=paste("tmp=lm(test ~ factor(Dx) + Age + factor(Sex) + factor(Dx):Age", whole, site, ")", sep='')))
        models.cort[[x-ncol(Covs)]]=tmp #Store the model fit for future reference

        #subjects can be dropped if they are missing so we can get the precise number of controls/patients for each region tested
        n.controls[x-ncol(Covs)] = length(which(tmp$model[,2] == 0))
        n.patients[x-ncol(Covs)] = length(which(tmp$model[,2] == 1))

        #Convert the lm model to a summary format so we can extract statistics
        tmp=summary(tmp)
        tstat=tmp$coefficients[(6 + n.sites),3] # Get t-statistic from regression to convert to Cohens d
        tstat.df=tmp$df[2]

        #collect effect size data
        d.cort[x-ncol(Covs)]=partial.d(tstat,tstat.df,n.controls[x-ncol(Covs)],n.patients[x-ncol(Covs)])
        se.cort[x-ncol(Covs)]=se.d2(d.cort[x-ncol(Covs)],n.controls[x-ncol(Covs)],n.patients[x-ncol(Covs)])
        bound.cort=CI1(d.cort[x-ncol(Covs)],se.cort[x-ncol(Covs)])
        low.ci.cort[x-ncol(Covs)]=bound.cort[1]
        up.ci.cort[x-ncol(Covs)]=bound.cort[2]
}

#save results
save(d.cort,se.cort,low.ci.cort,up.ci.cort,n.controls,n.patients, file=paste0(transferDir, "/","EffectSizes_SZvHV_DxByAge_",filetype,"_ctlWHOLE.Rdata"))
save(models.cort, file=paste0(transferDir, "/","Models_SZvHV_DxByAge_",filetype,"_ctlWHOLE.Rdata"))
} else {cat('NOT Running: 2.2 Diagnosis by Age interaction (all SZ patients vs controls) control for ', whole, '\n')}



# 3. AntiPsychotic (AP) use (typical, atypical, both, none) comparisons

#detach(merged_ordered)
#rm(merged_ordered)
merged_ordered=swap.merged_ordered

# create data subsets
merged_ordered_AP1vsAP0 <- merged_ordered[((merged_ordered$AP!=2) & (merged_ordered$AP!=3) & (merged_ordered$AP!=4) & !(is.na(merged_ordered$AP))),]
merged_ordered_AP2vsAP0 <- merged_ordered[((merged_ordered$AP!=1) & (merged_ordered$AP!=3) & (merged_ordered$AP!=4) & !(is.na(merged_ordered$AP))),]
merged_ordered_AP3vsAP0 <- merged_ordered[((merged_ordered$AP!=2) & (merged_ordered$AP!=1) & (merged_ordered$AP!=4) & !(is.na(merged_ordered$AP))),]
merged_ordered_AP4vsAP0 <- merged_ordered[((merged_ordered$AP!=2) & (merged_ordered$AP!=3) & (merged_ordered$AP!=1) & !(is.na(merged_ordered$AP))),]
merged_ordered_AP2vsAP1 <- merged_ordered[((merged_ordered$AP!=0) & (merged_ordered$AP!=3) & (merged_ordered$AP!=4) & !(is.na(merged_ordered$AP))),]
merged_ordered_AP3vsAP1 <- merged_ordered[((merged_ordered$AP!=0) & (merged_ordered$AP!=2) & (merged_ordered$AP!=4) & !(is.na(merged_ordered$AP))),]
merged_ordered_AP4vsAP1 <- merged_ordered[((merged_ordered$AP!=0) & (merged_ordered$AP!=2) & (merged_ordered$AP!=3) & !(is.na(merged_ordered$AP))),]
merged_ordered_AP3vsAP2 <- merged_ordered[((merged_ordered$AP!=0) & (merged_ordered$AP!=1) & (merged_ordered$AP!=4) & !(is.na(merged_ordered$AP))),]
merged_ordered_AP4vsAP2 <- merged_ordered[((merged_ordered$AP!=0) & (merged_ordered$AP!=1) & (merged_ordered$AP!=3) & !(is.na(merged_ordered$AP))),]
merged_ordered_AP4vsAP3 <- merged_ordered[((merged_ordered$AP!=0) & (merged_ordered$AP!=1) & (merged_ordered$AP!=2) & !(is.na(merged_ordered$AP))),]

# AntiPsychotic (AP) groups vs controls and vs eachother

for(APcomparison in c("AP1vsAP0","AP2vsAP0","AP3vsAP0","AP4vsAP0","AP2vsAP1","AP3vsAP1","AP4vsAP1","AP3vsAP2","AP4vsAP2","AP4vsAP3")){
#for(APcomparison in c("AP3vsAP0")){

   # Make sure there are 2 levels and each level has at least 2 data points: minimum required to compute variance
   # Can decide later on reasonable number per group
   #min_range=min(range(merged_ordered_AP3vsAP0$AP,na.rm=T));
   #max_range=max(range(merged_ordered_AP3vsAP0$AP,na.rm=T));
   #n_min_range=length(which(merged_ordered_AP3vsAP0$AP==min_range));
   #n_max_range=length(which(merged_ordered_AP3vsAP0$AP==max_range));
   #if( (min_range != max_range) & n_min_range>2 & n_max_range>2 ) {

   #Store models for troubleshooting
   models.cort=NULL; # This will become a list where we store all of the models made by lm

   #allocate empty vectors to store adjust effect sizes, se, ci (noicv)
   d.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   se.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   low.ci.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   up.ci.cort=rep(NA,(ncol(merged_ordered)-ncol(Covs)))

   #n.controls=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   #n.patients=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   n.APlow=rep(NA,(ncol(merged_ordered)-ncol(Covs)))
   n.APhigh=rep(NA,(ncol(merged_ordered)-ncol(Covs)))

   #attach main file
   detach(merged_ordered)
   data_subset_name=paste0("merged_ordered",sep="_",APcomparison)
   merged_ordered=get(data_subset_name)
   #merged_ordered=swap.merged_ordered
   #merged_ordered <- merged_ordered[((merged_ordered$AD!=1) & !(is.na(merged_ordered$AD))),]
   attach(merged_ordered)


   # Make sure there are 2 levels and each level has at least 2 data points: minimum required to compute variance
   # Can decide later on reasonable number per group
   min_range=min(range(merged_ordered$AP,na.rm=T));
   max_range=max(range(merged_ordered$AP,na.rm=T));
   n_min_range=length(which(merged_ordered$AP==min_range));
   n_max_range=length(which(merged_ordered$AP==max_range));
   if( (min_range != max_range) & n_min_range>3 & n_max_range>3 ) {

   if(nrow(merged_ordered) > (n.sites + 3)){
      cat('Running: 3. Antipsychotic (AP) Group Comparison', APcomparison,' control for ', whole, '\n')

   #Loop through and perform each regression
   for(x in (ncol(Covs)+1):ncol(merged_ordered)){
        pheno=merged_ordered[!is.na(merged_ordered[,x]),x] #Check to make sure there are observations for a given structure
        #check if the phenotype is singular after NA removal
        if(length(pheno)==0){
                next; # Skip the whole structure if there are no observations
        }
        site=NULL # These are just string variables that modify the model in R if there are Site variables in the Covariates.csv file
        if(n.sites != 0){
                site=" + " #Starts the string to add-on Site variables
                for(i in 1:n.sites){
                        if(i != n.sites){
                                site=paste(site, "Site", i, " + ", sep='')
                        }
                        else{
                                site=paste(site, "Site", i, sep='')
                        }
                }
        }
        #Run the model
        test=merged_ordered[,x]
        eval(parse(text=paste("tmp=lm(test ~ factor(AP) + Age + Sex", whole, site, ")", sep='')))
        models.cort[[x-ncol(Covs)]]=tmp #Store the model fit for future reference

        #subjects can be dropped if they are missing so we can get the precise number of controls/patients for each region tested
        #n.controls[x-ncol(Covs)] = length(which(tmp$model[,2] == 0))
        #n.patients[x-ncol(Covs)] = length(which(tmp$model[,2] == 1))

        n.APlow[x-ncol(Covs)] = length(which(tmp$model[,2] == min_range))
        n.APhigh[x-ncol(Covs)] = length(which(tmp$model[,2] == max_range))

        #Convert the lm model to a summary format so we can extract statistics
        tmp=summary(tmp)
        tstat=tmp$coefficients[2,3] # Get t-statistic from regression to convert to Cohens d
        tstat.df=tmp$df[2]

        #collect effect size data
        #d.cort[x-ncol(Covs)]=partial.d(tstat,tstat.df,n.controls[x-ncol(Covs)],n.patients[x-ncol(Covs)])
        #se.cort[x-ncol(Covs)]=se.d2(d.cort[x-ncol(Covs)],n.controls[x-ncol(Covs)],n.patients[x-ncol(Covs)])
        # APhigh, as indicated by subjects with AP==max_range is the experimental condition
        # APlow, as indicated by subjects with AP==max_range is the control condition
        d.cort[x-ncol(Covs)]=partial.d(tstat,tstat.df,n.APlow[x-ncol(Covs)],n.APhigh[x-ncol(Covs)])
        se.cort[x-ncol(Covs)]=se.d2(d.cort[x-ncol(Covs)],n.APlow[x-ncol(Covs)],n.APhigh[x-ncol(Covs)])
        bound.cort=CI1(d.cort[x-ncol(Covs)],se.cort[x-ncol(Covs)])
        low.ci.cort[x-ncol(Covs)]=bound.cort[1]
        up.ci.cort[x-ncol(Covs)]=bound.cort[2]
   }

   #save results
   save(d.cort,se.cort,low.ci.cort,up.ci.cort,n.APlow,n.APhigh, file=paste0(transferDir, "/","EffectSizes_",APcomparison,sep="_",filetype,"_ctlWHOLE.Rdata"))
   save(models.cort, file=paste0(transferDir, "/","Models_AP_",APcomparison,sep="_",filetype,"_ctlWHOLE.Rdata"))

   } else {cat('NOT Running: 3. Antipsychotic (AP) Group Comparison', APcomparison, ' control for ', whole, '\n')}

   } # end if min and max range not the same and n_min_range and n_max_range > 0

} # end for loop APcomparison



#### TVE ##########

cat('Done working on ', fsfile,'\n')
sink()
detach(merged_ordered)
print(readLines(paste0(transferDir, "/", fsfile,paste0("_",cc,".log"))))

} # end fsfile loop

} # end cc loop

quit("n")

#### TVE ##########
