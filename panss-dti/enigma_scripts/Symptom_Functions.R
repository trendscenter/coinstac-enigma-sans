# Calculate the MAP and EXP factors 

# Calculate PANSS Total, Negative, Positive, General and the EXP and MAP Factors

CalcPANSS <- function(PANSS){

  PANSS$TOT <- (PANSS$PANSS1 + PANSS$PANSS2 + PANSS$PANSS3 + PANSS$PANSS4 + PANSS$PANSS5 + PANSS$PANSS6 + PANSS$PANSS7 + PANSS$PANSS8 +
                    PANSS$PANSS9 + PANSS$PANSS10 + PANSS$PANSS11 + PANSS$PANSS12 + PANSS$PANSS13 + PANSS$PANSS14 + PANSS$PANSS15 + PANSS$PANSS16
                  + PANSS$PANSS17 + PANSS$PANSS18 + PANSS$PANSS19 + PANSS$PANSS20 + PANSS$PANSS21 + PANSS$PANSS22 + PANSS$PANSS23
                  + PANSS$PANSS24 + PANSS$PANSS25 + PANSS$PANSS26 + PANSS$PANSS27 + PANSS$PANSS28 + PANSS$PANSS29 + PANSS$PANSS30)
  
  PANSS$NEG <- (PANSS$PANSS8 + PANSS$PANSS9 + PANSS$PANSS10 + PANSS$PANSS11 + PANSS$PANSS12 + PANSS$PANSS13 + PANSS$PANSS14)
  
  PANSS$POS <- (PANSS$PANSS1 + PANSS$PANSS2 + PANSS$PANSS3 + PANSS$PANSS4 + PANSS$PANSS5 + PANSS$PANSS6 + PANSS$PANSS7)
  
  PANSS$GEN <- (PANSS$PANSS15 + PANSS$PANSS16 + PANSS$PANSS17 + PANSS$PANSS18 + PANSS$PANSS19 + PANSS$PANSS20 + PANSS$PANSS21
                         + PANSS$PANSS22 + PANSS$PANSS23 + PANSS$PANSS24 + PANSS$PANSS25 + PANSS$PANSS26 + PANSS$PANSS27 + PANSS$PANSS28
                         + PANSS$PANSS29 + PANSS$PANSS30)
  
  PANSS$EXP <- (PANSS$PANSS8*.749)+(PANSS$PANSS10*.696)+(PANSS$PANSS13*.783)+(PANSS$PANSS21*.744)
  
  PANSS$MAP <- (PANSS$PANSS9*.666)+(PANSS$PANSS11*.799)+(PANSS$PANSS30*.886)
  
  return(as.data.frame(PANSS))
}

# calculate the SANS factorizations here and add to the SANS matrix--external function

# 5 negative domain scores from Strauss (personal communication)

# 5 global ratings (8,13,17,22, and 25)
# SANSTOT = sum(SANS items 1-7, 9-12, 14-16, 18-21)
# SANSGlobal = sum(SANS items 8, 13, 17, 22)

CalcSans <- function(SANS) { 
  # this adds columns for total Globals, total other than globals, the MAP/EXP factors, and the Strauss 5
  # Not using Attention (global 5) or those questions
  # Not including SANS6 Inappropriate affect
  # Not including SANS10, poverty of speech content
  # already checked that the data frame has SubjID followed by SANS1-SANS25
  # SANS needs to be cast as floats
  
  #two types of global sums
  
  Globals=c("SANS8","SANS13","SANS17", "SANS22") # without Attention question
  SANS$GlobalTotal = rowSums(SANS[,Globals])
  
  Tot=c(1,2,3,4,5,7,9,11,12,14,15,16,18,19,20,21)  # question numbers without #6 and without attention
  SANS$SANSSum = rowSums(SANS[,Tot+1])  # allowing for subjID column
  
  
  #Strauss Factors here--from the excel files they sent in January
  # get the numeric parts, drop subjID
  SANSnum=as.matrix(SANS[,2:28])
  
  # MAP and EXP
  # MAP usually equals sums of questions without globals from Anhedonia/Avolition and Asociality
  # EXP usually equals sums of questions without globals from Alogia, Blunted Affect (drop blocking and pov.of content)
  # Using Strauss's 2 factor solution as he calls it MAP/EXP--weighted versions of the above, includes blocking question
  
  SANS$MAP = SANSnum[,18]*.810 + SANSnum[,19]*.641 + SANSnum[,20]*.688 + SANSnum[,21]*.804 + 
    SANSnum[,14]*.314 + SANSnum[,15]*.553 + SANSnum[,16]*.648
  
  SANS$EXP = SANSnum[,1]*.831 + SANSnum[,2]*.835 + SANSnum[,3]*.901 + SANSnum[,4]*.535 +
    SANSnum[,5]*.654 + SANSnum[,7]*.765 + SANSnum[,9]*.775 + SANSnum[,11]*.291 + SANSnum[,12]*.415
  
  SANS$AnhedoniaFac = SANSnum[,18]*0.707
  SANS$AsocialityFac = SANSnum[,19]*.685 + SANSnum[,20]*0.736 + SANSnum[,21]*.868
  SANS$AvolitionFac = SANSnum[,14]*.395 + SANSnum[,15]*.666 + SANSnum[,16]*0.807
  SANS$BluntedAffectFac= SANSnum[,1]*.835 + SANSnum[,2]*.840 + SANSnum[,3]*.904 + SANSnum[,4]*.541 + SANSnum[,5]*.658 + SANSnum[,7]*.770
  SANS$AlogiaFac = SANSnum[,9]*.867 + SANSnum[,11]*.313 + SANSnum[,12]*.440
  
  return(SANS)
  
}

