# calculate the SANS factorizations here and add to the SANS matrix--external function

# 5 negative domain scores from Strauss 2018 Figure 1--but his SANS numbers not the same as usual? 
# From Figure 1 he uses the following--these are not the standard numbering, he's dropped a few and has some that the 
# Andreassen scale doesn't have; I'm checking with him to confirm!
# Anhedonia = 0.89 * SANS15
# this is labeled "Anhedonia" (maybe the global rating?)
# Asociality = 0.91 * SANS14 + 0.41 * SANS16 + 0.63* SANS17
# this is "Asociality" (SANS18?? global rating?), then SANS19 and SANS20
# Avolition = 0.3 * SANS10 + 0.88 * SANS11 + 0.75*SANS12 + .65 * SANS13
# this is actually SANS14, and what he calls 11 and 12 he calls Current Role function level and quality, then SANS13 is physical anergia
# which is actually SANS16
# Affect = 0.9 * SANS1 + 0.84 * SANS2 + .89*SANS3 + 0.56*SANS4 + .86*SANS5 + .82*SANS6
# This is SANS1-SANS5 but SANS6 is Lack of vocal Inflections, which is SANS7 (he dropped inappropriate affect)
# Alogia = .89*SANS7 + .41* SANS8 + .54* SANS9
# this is actually SANS9, SANS11, and SANS12


# 5 global ratings (8,13,17,22, and 25)
# SANSTOT = sum(SANS items 1-7, 9-12, 14-16, 18-21, and 23-24)
# SANSGlobal = sum(SANS items 8, 13, 17, 22, and 25)

CalcSans <- function(SANS) {
  # already checked that the data frame has SubjID followed by SANS1-SANS25
  # SANS needs to be cast as floats
  
  
  Globals=c("SANS8","SANS13","SANS17", "SANS22", "SANS25")
  SANS$GlobalTotal = rowSums(SANS[,Globals])
  
  Tot=c(1,2,3,4,5,6,7,9,10,11,12,14,15,16,18,19,20,21,23,24)
  SANS$SANSSum = rowSums(SANS[,Tot+1])
  
  #Strauss Factors here--fix when have the right answers!! these column numbers are not completely correct
  # get the numeric parts, drop subjID
  SANSnum=as.matrix(SANS[,2:28])

  SANS$AnhedoniaFac = SANSnum[,15]*0.89
  SANS$AsocialityFac = SANSnum[,18]*.91 + SANSnum[,19]*0.41 + SANSnum[,20]*.63
  SANS$AvolitionFac = SANSnum[,14]*.30 + SANSnum[,15]*.88 + SANSnum[,16]*0.65
  SANS$BluntedAffectFac= SANSnum[,1]*.9 + SANSnum[,2]*.84 + SANSnum[,3]*.89 + SANSnum[,4]*.56 + SANSnum[,5]*.86 + SANSnum[,7]*.82
  SANS$AlogiaFac = SANSnum[,9]*.89 + SANSnum[,11]*.41 + SANSnum[,12]*.54
  
  return(SANS)
  
}
