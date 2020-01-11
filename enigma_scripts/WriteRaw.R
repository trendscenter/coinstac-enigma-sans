WriteRawFiles <- function(merged_ordered, Ncol, outputfilename) {

raw.means=colMeans(merged_ordered[(Ncol+1):ncol(merged_ordered)], na.rm=T)

#Get raw sd and number of subjects included for each of the structures
sd.raw=rep(NA,(ncol(merged_ordered)-Ncol))
n.raw=rep(NA,(ncol(merged_ordered)-Ncol))
min.raw=rep(NA,(ncol(merged_ordered)-Ncol))
max.raw=rep(NA,(ncol(merged_ordered)-Ncol))
for(z in (Ncol+1):ncol(merged_ordered)){
  sd.raw[z-Ncol]=sd(merged_ordered[,z], na.rm=T)
  n.raw[z-Ncol]=length(merged_ordered[which(!is.na(merged_ordered[,z])),z])
  min.raw[z-Ncol]=min(merged_ordered[,z], na.rm=T)
  max.raw[z-Ncol]=max(merged_ordered[,z], na.rm=T)
}

#Save raw values
save(raw.means, sd.raw, n.raw, min.raw, max.raw, file=outputfilename)


#Get raw means for each of the structures for the Dx groups:
# patients (Dx==1)
# controls (Dx==0)

for(DXvalue in 0:1){
  
  DXgroup=which(merged_ordered$Dx==DXvalue)
  raw.means=colMeans(merged_ordered[DXgroup,(Ncol+1):ncol(merged_ordered)], na.rm=T)
  
  #Get raw sd and number of subjects included for each of the structures for the Dx groups
  sd.raw=rep(NA,(ncol(merged_ordered)-Ncol))
  n.raw=rep(NA,(ncol(merged_ordered)-Ncol))
  min.raw=rep(NA,(ncol(merged_ordered)-Ncol))
  max.raw=rep(NA,(ncol(merged_ordered)-Ncol))
  for(z in (Ncol+1):ncol(merged_ordered)){
    sd.raw[z-Ncol]=sd(merged_ordered[DXgroup,z], na.rm=T)
    n.raw[z-Ncol]=length(merged_ordered[which(!is.na(merged_ordered[DXgroup,z])),z])
    min.raw[z-Ncol]=min(merged_ordered[DXgroup,z], na.rm=T)
    max.raw[z-Ncol]=max(merged_ordered[DXgroup,z], na.rm=T)
  }
  
  #Save raw values
  save(raw.means, sd.raw, n.raw, min.raw, max.raw, file=paste0(outputfilename,"_",DXvalue,".Rdata"))    
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
  raw.means=colMeans(merged_ordered[APgroup,(Ncol+1):ncol(merged_ordered)], na.rm=T)
  
  #Get raw sd and number of subjects included for each of the structures for the APgroup
  sd.raw=rep(NA,(ncol(merged_ordered)-Ncol))
  n.raw=rep(NA,(ncol(merged_ordered)-Ncol))
  min.raw=rep(NA,(ncol(merged_ordered)-Ncol))
  max.raw=rep(NA,(ncol(merged_ordered)-Ncol))
  for(z in (Ncol+1):ncol(merged_ordered)){
    sd.raw[z-Ncol]=sd(merged_ordered[APgroup,z], na.rm=T)
    n.raw[z-Ncol]=length(merged_ordered[which(!is.na(merged_ordered[APgroup,z])),z])
    min.raw[z-Ncol]=min(merged_ordered[APgroup,z], na.rm=T)
    max.raw[z-Ncol]=max(merged_ordered[APgroup,z], na.rm=T)
  }
  #Save raw values
  save(raw.means, sd.raw, n.raw, min.raw, max.raw, file=paste0(outputfilename,"_",APvalue,".Rdata"))  
  
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

#Check that Sex was coded properly--already done!
# if((n.fem + n.mal) != length(merged_ordered$Sex)){
#   stop('Did you remember to code the Sex covariate as Males=1 and Females=2?\n')
# }

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

save(age.mu, age.sd, age.range, age.mu.dx0, age.sd.dx0, age.range.dx0, age.mu.dx1, age.sd.dx1, 
     age.range.dx1, age.mu, age.sd, age.range, age.mu.dx0, age.sd.dx0, age.range.dx0, age.mu.dx1, age.sd.dx1, 
     age.range.dx1, age.mu.AP0, age.sd.AP0, age.range.AP0, age.mu.AP1, age.sd.AP1, age.range.AP1, age.mu.AP2, age.sd.AP2, 
     age.range.AP2, age.mu.AP3, age.sd.AP3, age.range.AP3, age.mu.AP4, age.sd.AP4, age.range.AP4, 
     cpz.mu.dx1, cpz.sd.dx1, cpz.range.dx1, cpz.mu.AP1, cpz.sd.AP1, cpz.range.AP1, cpz.mu.AP2, cpz.sd.AP2, cpz.range.AP2, 
     cpz.mu.AP3, cpz.sd.AP3, cpz.range.AP3, cpz.mu.AP4, cpz.sd.AP4, cpz.range.AP4, ao.mu.dx1, ao.sd.dx1, ao.range.dx1, 
     ao.mu.AP1, ao.sd.AP1, ao.range.AP1, ao.mu.AP2, ao.sd.AP2, ao.range.AP2, ao.mu.AP3, ao.sd.AP3, ao.range.AP3, ao.mu.AP4, 
     ao.sd.AP4, ao.range.AP4, 
     durill.mu.dx1, durill.sd.dx1, durill.range.dx1, durill.mu.AP1, durill.sd.AP1, durill.range.AP1, durill.mu.AP2, 
     durill.sd.AP2, durill.range.AP2, durill.mu.AP3, durill.sd.AP3, durill.range.AP3, durill.mu.AP4, durill.sd.AP4, 
     durill.range.AP4, 
     pansstot.mu.dx1, pansstot.sd.dx1, pansstot.range.dx1, pansstot.mu.AP1, pansstot.sd.AP1, pansstot.range.AP1, 
     pansstot.mu.AP2, pansstot.sd.AP2, pansstot.range.AP2, pansstot.mu.AP3, pansstot.sd.AP3, pansstot.range.AP3,
     pansstot.mu.AP4, pansstot.sd.AP4, pansstot.range.AP4, pansspos.mu.dx1, pansspos.sd.dx1, pansspos.range.dx1, 
     pansspos.mu.AP1, pansspos.sd.AP1, pansspos.range.AP1, pansspos.mu.AP2, pansspos.sd.AP2, pansspos.range.AP2, 
     pansspos.mu.AP3, pansspos.sd.AP3, pansspos.range.AP3, pansspos.mu.AP4, pansspos.sd.AP4, pansspos.range.AP4, 
     panssneg.mu.dx1, panssneg.sd.dx1, panssneg.range.dx1, panssneg.mu.AP1, panssneg.sd.AP1, panssneg.range.AP1, 
     panssneg.mu.AP2, panssneg.sd.AP2, panssneg.range.AP2, panssneg.mu.AP3, panssneg.sd.AP3, panssneg.range.AP3, 
     panssneg.mu.AP4, panssneg.sd.AP4, panssneg.range.AP4, 
     sapstot.mu.dx1, sapstot.sd.dx1, sapstot.range.dx1, sapstot.mu.AP1, sapstot.sd.AP1, sapstot.range.AP1, 
     sapstot.mu.AP2, sapstot.sd.AP2, sapstot.range.AP2, sapstot.mu.AP3, sapstot.sd.AP3, sapstot.range.AP3, 
     sapstot.mu.AP4, sapstot.sd.AP4, sapstot.range.AP4, sanstot.mu.dx1, sanstot.sd.dx1, sanstot.range.dx1, 
     sanstot.mu.AP1, sanstot.sd.AP1, sanstot.range.AP1, sanstot.mu.AP2, sanstot.sd.AP2, sanstot.range.AP2, 
     sanstot.mu.AP3, sanstot.sd.AP3, sanstot.range.AP3, sanstot.mu.AP4, sanstot.sd.AP4, sanstot.range.AP4, 
     parentses.mu.dx0, parentses.sd.dx0, parentses.range.dx0, parentses.n.dx0, parentses.mu.dx1,parentses.sd.dx1,
     parentses.range.dx1, parentses.n.dx1, parentses.mu.AP1,parentses.sd.AP1,parentses.range.AP1,parentses.mu.AP2,
     parentses.sd.AP2,parentses.range.AP2,parentses.mu.AP3,parentses.sd.AP3,parentses.range.AP3,parentses.mu.AP4,
     parentses.sd.AP4,parentses.range.AP4,
     n.dx0, n.dx1, n.AP0, n.AP1, n.AP2, n.AP3, n.AP4, n.fem, n.mal, n.fem.dx0, n.mal.dx0, n.fem.dx1, n.mal.dx1, 
     n.fem.AP0, n.mal.AP0, n.fem.AP1, n.mal.AP1, n.fem.AP2, n.mal.AP2, n.fem.AP3, n.mal.AP3, n.fem.AP4, n.mal.AP4, 
     n.hand0.dx0, n.hand0.dx1, n.hand1.dx0, n.hand1.dx1, n.hand2.dx0, n.hand2.dx1, 
     iq.mu.dx0, iq.sd.dx0, iq.range.dx0, iq.n.dx0, iq.mu.dx1, iq.sd.dx1, iq.range.dx1, iq.n.dx1, 
     file=paste0("Demographics_",outputfilename,".Rdata"))

}