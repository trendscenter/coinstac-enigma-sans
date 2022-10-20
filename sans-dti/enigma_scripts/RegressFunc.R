# Regression Functions


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

# this is what builds the strings for the models: Must align with the SANSCalcs functions!

PredictCovs <- function(predictor) {
  predictString <-
    switch(
      predictor,
      # For SANS Measures
      "SANSSum" = "SANSSum",
      "SANSGlobals" = "SANS8 + SANS13 + SANS17 + SANS22",
      "SANSFactors" = "AnhedoniaFac + AsocialityFac + AvolitionFac + BluntedAffectFac + AlogiaFac",
      "SANSMAP" = "MAP",
      "SANSEXP" = "EXP",
      "SANSGlobal1"="SANS8", 
      "SANSGlobal2"="SANS13", 
      "SANSGlobal3"="SANS17", 
      "SANSGlobal4"="SANS22", 
      "SANSFac1"="AnhedoniaFac", 
      "SANSFac2"="AsocialityFac", 
      "SANSFac3"="AvolitionFac", 
      "SANSFac4"="BluntedAffectFac", 
      "SANSFac5"="AlogiaFac",
      "SANSMAPwSANSSum" = "MAP + SANSSum",
      "SANSEXPwSANSSum" = "EXP + SANSSum",
      "SANSGlobal1wSANSSum"="SANS8 + SANSSum", 
      "SANSGlobal2wSANSSum"="SANS13 + SANSSum", 
      "SANSGlobal3wSANSSum"="SANS17 + SANSSum", 
      "SANSGlobal4wSANSSum"="SANS22 + SANSSum", 
      "SANSFac1wSANSSum"="AnhedoniaFac + SANSSum", 
      "SANSFac2wSANSSum"="AsocialityFac + SANSSum", 
      "SANSFac3wSANSSum"="AvolitionFac + SANSSum", 
      "SANSFac4wSANSSum"="BluntedAffectFac + SANSSum", 
      "SANSFac5wSANSsum"="AlogiaFac + SANSSum",
      "SANSGlobal1wSANSEXP"="SANS8 + EXP", 
      "SANSGlobal2wSANSEXP"="SANS13 + EXP", 
      "SANSGlobal3wSANSEXP"="SANS17 + EXP", 
      "SANSGlobal4wSANSEXP"="SANS22 + EXP", 
      "SANSFac1wSANSEXP"="AnhedoniaFac + EXP", 
      "SANSFac2wSANSEXP"="AsocialityFac + EXP", 
      "SANSFac3wSANSEXP"="AvolitionFac + EXP", 
      "SANSFac4wSANSEXP"="BluntedAffectFac + EXP", 
      "SANSFac5wSANSEXP"="AlogiaFac + EXP",
      "SANSGlobal1wSANSMAP"="SANS8 + MAP", 
      "SANSGlobal2wSANSMAP"="SANS13 + MAP", 
      "SANSGlobal3wSANSMAP"="SANS17 + MAP", 
      "SANSGlobal4wSANSMAP"="SANS22 + MAP", 
      "SANSFac1wSANSMAP"="AnhedoniaFac + MAP", 
      "SANSFac2wSANSMAP"="AsocialityFac + MAP", 
      "SANSFac3wSANSMAP"="AvolitionFac + MAP", 
      "SANSFac4wSANSMAP"="BluntedAffectFac + MAP", 
      "SANSFac5wSANSMAP"="AlogiaFac + MAP",
      # For PANSS Measures
      "TOT"="TOT", 
      "NEG"="NEG", 
      "POS"="POS", 
      "GEN"="GEN", 
      "EXP"="EXP", 
      "MAP"="MAP",
      "EXPwTOT"="EXP + TOT",
      "EXPwGEN"="EXP + GEN",
      "EXPwNEG"="EXP + NEG",
      "MAPwTOT"="MAP + TOT",
      "MAPwGEN"="MAP + GEN",
      "MAPwNEG"="MAP + NEG"
    )

  return(predictString)
}
