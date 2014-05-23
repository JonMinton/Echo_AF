rm(list=ls())
wd <- "X:/BMJ Echo AF Manuscript/S8/"
#wd <- "C:/Users/Jon Minton/Documents/R/Echo AF/"

setwd(wd)
source("PSA_funDat.R")

# Look at sens & spec
#require(MASS)
#contour(kde2d(sens.PSA, spec.PSA), xlim=c(0.75, 1), ylim=c(0.25,0.5), xlab=c("Sensitivity"), ylab="Specificity", main="Contour map of estimated sensitivity and specificity")


###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
################### M A I N    A N A L Y S E S  ###################################################################################
###################################################################################################################################
###################################################################################################################################
hrProp.C0.Exp <- with  (CHADS.Data, events[score==0] / total[score==0])
hrProp.C1.Exp <- with  (CHADS.Data, events[score==1] / total[score==1])
hrProp.CV0.Exp <- with (CV.Data,    events[score==0] / total[score==0])
hrProp.CV1.Exp <- with (CV.Data,    events[score==1] / total[score==1])


# Dab C=0

# Write a meta-function to reduce risk of human (i.e. my) error
W_50_F.Analyses <- Meta.Function(Data=Data.WarfC0_50_F, hrProp.Exp=hrProp.C0.Exp, hrProp.PSA=hrProp.C0, sens.Exp=sens.Exp, spec.Exp=spec.Exp, sens.PSA=sens.PSA, spec.PSA=spec.PSA)
W_50_M.Analyses <- Meta.Function(Data=Data.WarfC0_50_M, hrProp.Exp=hrProp.C0.Exp, hrProp.PSA=hrProp.C0, sens.Exp=sens.Exp, spec.Exp=spec.Exp, sens.PSA=sens.PSA, spec.PSA=spec.PSA)

W_65_F.Analyses <- Meta.Function(Data=Data.WarfC0_65_F, hrProp.Exp=hrProp.C0.Exp, hrProp.PSA=hrProp.C0, sens.Exp=sens.Exp, spec.Exp=spec.Exp, sens.PSA=sens.PSA, spec.PSA=spec.PSA)
W_65_M.Analyses <- Meta.Function(Data=Data.WarfC0_65_M, hrProp.Exp=hrProp.C0.Exp, hrProp.PSA=hrProp.C0, sens.Exp=sens.Exp, spec.Exp=spec.Exp, sens.PSA=sens.PSA, spec.PSA=spec.PSA)

R_50_F.Analyses <- Meta.Function(Data=Data.RivC0_50_F, hrProp.Exp=hrProp.C0.Exp, hrProp.PSA=hrProp.C0, sens.Exp=sens.Exp, spec.Exp=spec.Exp, sens.PSA=sens.PSA, spec.PSA=spec.PSA)
R_50_M.Analyses <- Meta.Function(Data=Data.RivC0_50_M, hrProp.Exp=hrProp.C0.Exp, hrProp.PSA=hrProp.C0, sens.Exp=sens.Exp, spec.Exp=spec.Exp, sens.PSA=sens.PSA, spec.PSA=spec.PSA)

R_65_F.Analyses <- Meta.Function(Data=Data.RivC0_65_F, hrProp.Exp=hrProp.C0.Exp, hrProp.PSA=hrProp.C0, sens.Exp=sens.Exp, spec.Exp=spec.Exp, sens.PSA=sens.PSA, spec.PSA=spec.PSA)
R_65_M.Analyses <- Meta.Function(Data=Data.RivC0_65_M, hrProp.Exp=hrProp.C0.Exp, hrProp.PSA=hrProp.C0, sens.Exp=sens.Exp, spec.Exp=spec.Exp, sens.PSA=sens.PSA, spec.PSA=spec.PSA)

D_65_F.Analyses <- Meta.Function(Data=Data.DabC0_65_F, hrProp.Exp=hrProp.C0.Exp, hrProp.PSA=hrProp.C0, sens.Exp=sens.Exp, spec.Exp=spec.Exp, sens.PSA=sens.PSA, spec.PSA=spec.PSA)
D_65_M.Analyses <- Meta.Function(Data=Data.DabC0_65_M, hrProp.Exp=hrProp.C0.Exp, hrProp.PSA=hrProp.C0, sens.Exp=sens.Exp, spec.Exp=spec.Exp, sens.PSA=sens.PSA, spec.PSA=spec.PSA)



# # Bottom line: base.cost, base.qaly, comp.cost, comp.qaly, ICER, and credible intervals

icer.results <-  Extract.mainQI(W_50_F.Analyses)

icer.results <- rbind(icer.results,      Extract.mainQI(W_50_M.Analyses)   )
icer.results <- rbind(icer.results,      Extract.mainQI(W_65_F.Analyses)   )
icer.results <- rbind(icer.results,      Extract.mainQI(W_65_M.Analyses)   )


rownames(icer.results) <- c("warf.50.F", "warf.50.M", "warf.65.F", "warf.65.M")



# clinical quantities of interest
otherQI.results <-                                 apply(W_50_F.Analyses[["PSA.otherQI"]],   2, mean)
otherQI.results <- rbind(otherQI.results,          apply(W_50_M.Analyses[["PSA.otherQI"]],  2, mean) )
otherQI.results <- rbind(otherQI.results,          apply(W_65_M.Analyses[["PSA.otherQI"]],  2, mean) )
otherQI.results <- rbind(otherQI.results,          apply(W_65_F.Analyses[["PSA.otherQI"]],  2, mean) )

rownames(otherQI.results) <- c("warf.50.F", "warf.50.M", "warf.65.F", "warf.65.M")

# PSA

Draw.PSA(W_50_F.Analyses[["PSA"]]   , main="", export="scatter_W_50_F.jpeg" )
Draw.PSA(W_50_M.Analyses[["PSA"]]   , main="", export="scatter_W_50_M.jpeg" )
Draw.PSA(W_65_F.Analyses[["PSA"]]   , main="", export="scatter_W_65_F.jpeg" )
Draw.PSA(W_65_M.Analyses[["PSA"]]   , main="", export="scatter_W_65_M.jpeg" )

Draw.PSA(R_50_F.Analyses[["PSA"]]   , main="", export="scatter_R_50_F.jpeg" )
Draw.PSA(R_50_M.Analyses[["PSA"]]   , main="", export="scatter_R_50_M.jpeg" )
Draw.PSA(R_65_F.Analyses[["PSA"]]   , main="", export="scatter_R_65_F.jpeg" )
Draw.PSA(R_65_M.Analyses[["PSA"]]   , main="", export="scatter_R_65_M.jpeg" )

Draw.PSA(D_65_F.Analyses[["PSA"]]   , main="", export="scatter_D_65_F.jpeg" )
Draw.PSA(D_65_M.Analyses[["PSA"]]   , main="", export="scatter_D_65_M.jpeg" )


# CEAC
tmp <- Draw.CEAC.NB(W_50_F.Analyses)
tmp <- Draw.CEAC.NB(W_50_M.Analyses)
tmp <- Draw.CEAC.NB(W_65_F.Analyses)
tmp <- Draw.CEAC.NB(W_65_M.Analyses)




# EVPI
# want to extract values at 20,000 and 30,000
# tmp <- Draw.EVPI(W_50_F.Analyses[["EVPI"]])    ; evpi.results <-                           with(tmp, c(opLoss[maicer==20000], opLoss[maicer==30000])) 
# tmp <- Draw.EVPI(Dab.CV0.Analyses[["EVPI"]])   ; evpi.results <- rbind(evpi.results,       with(tmp, c(opLoss[maicer==20000], opLoss[maicer==30000])) )
# tmp <- Draw.EVPI(Warf.C0.Analyses[["EVPI"]])   ; evpi.results <- rbind(evpi.results,       with(tmp, c(opLoss[maicer==20000], opLoss[maicer==30000])) )
# tmp <- Draw.EVPI(Warf.C1.Analyses[["EVPI"]])   ; evpi.results <- rbind(evpi.results,       with(tmp, c(opLoss[maicer==20000], opLoss[maicer==30000])) )
# tmp <- Draw.EVPI(Warf.CV0.Analyses[["EVPI"]])  ; evpi.results <- rbind(evpi.results,       with(tmp, c(opLoss[maicer==20000], opLoss[maicer==30000])) )
# tmp <- Draw.EVPI(Warf.CV1.Analyses[["EVPI"]])  ; evpi.results <- rbind(evpi.results,       with(tmp, c(opLoss[maicer==20000], opLoss[maicer==30000])) )
# tmp <- Draw.EVPI(Riv.C0.Analyses[["EVPI"]])    ; evpi.results <- rbind(evpi.results,       with(tmp, c(opLoss[maicer==20000], opLoss[maicer==30000])) )
# tmp <- Draw.EVPI(Riv.C1.Analyses[["EVPI"]])    ; evpi.results <- rbind(evpi.results,       with(tmp, c(opLoss[maicer==20000], opLoss[maicer==30000])) )
# tmp <- Draw.EVPI(Riv.CV0.Analyses[["EVPI"]])   ; evpi.results <- rbind(evpi.results,       with(tmp, c(opLoss[maicer==20000], opLoss[maicer==30000])) )
# tmp <- Draw.EVPI(Riv.CV1.Analyses[["EVPI"]])   ; evpi.results <- rbind(evpi.results,       with(tmp, c(opLoss[maicer==20000], opLoss[maicer==30000])) )
# 
# rownames(evpi.results) <- c("dab.c0", "dab.cv0", "warf.c0", "warf.c1", "warf.cv0", "warf.cv1", "riv.c0", "riv.c1", "riv.cv0", "riv.cv1")
# colnames(evpi.results) <- c("20k", "30k")

# CEAF
Draw.CEAF.NB(W_65_F.Analyses[["PSA"]]   , W_65_F.Analyses[["EVPI"]], legend=c("TTE", "No TTE"))
Draw.CEAF.NB(W_65_M.Analyses[["PSA"]]   , W_65_M.Analyses[["EVPI"]], legend=c("TTE", "No TTE"))
Draw.CEAF.NB(W_50_F.Analyses[["PSA"]]   , W_50_F.Analyses[["EVPI"]], legend=c("TTE", "No TTE"))
Draw.CEAF.NB(W_50_M.Analyses[["PSA"]]   , W_50_M.Analyses[["EVPI"]], legend=c("TTE", "No TTE")) 

Draw.CEAF.NB(R_65_F.Analyses[["PSA"]]   , R_65_F.Analyses[["EVPI"]], legend=c("TTE", "No TTE"))
Draw.CEAF.NB(R_65_M.Analyses[["PSA"]]   , R_65_M.Analyses[["EVPI"]], legend=c("TTE", "No TTE")) 
Draw.CEAF.NB(R_50_F.Analyses[["PSA"]]   , R_50_F.Analyses[["EVPI"]], legend=c("TTE", "No TTE"))
Draw.CEAF.NB(R_50_M.Analyses[["PSA"]]   , R_50_M.Analyses[["EVPI"]], legend=c("TTE", "No TTE"))

Draw.CEAF.NB(D_65_F.Analyses[["PSA"]]   , D_65_F.Analyses[["EVPI"]], legend=c("TTE", "No TTE"))
Draw.CEAF.NB(D_65_M.Analyses[["PSA"]]   , D_65_M.Analyses[["EVPI"]], legend=c("TTE", "No TTE"))


tmp <- Draw.CEAF.NB(W_65_F.Analyses[["PSA"]]   , W_65_F.Analyses[["EVPI"]] , legend=c("TTE", "No TTE"), export="ceaf_W_65_F.jpeg" )
tmp <- Draw.CEAF.NB(W_65_M.Analyses[["PSA"]]   , W_65_M.Analyses[["EVPI"]] , legend=c("TTE", "No TTE"), export="ceaf_W_65_M.jpeg" )
tmp <- Draw.CEAF.NB(W_50_F.Analyses[["PSA"]]   , W_50_F.Analyses[["EVPI"]] , legend=c("TTE", "No TTE"), export="ceaf_W_50_F.jpeg" )
tmp <- Draw.CEAF.NB(W_50_M.Analyses[["PSA"]]   , W_50_M.Analyses[["EVPI"]] , legend=c("TTE", "No TTE"), export="ceaf_W_50_M.jpeg" )

tmp <- Draw.CEAF.NB(R_65_F.Analyses[["PSA"]]   , R_65_F.Analyses[["EVPI"]] , legend=c("TTE", "No TTE"), export="ceaf_R_65_F.jpeg" )
tmp <- Draw.CEAF.NB(R_65_M.Analyses[["PSA"]]   , R_65_M.Analyses[["EVPI"]] , legend=c("TTE", "No TTE"), export="ceaf_R_65_M.jpeg" )
tmp <- Draw.CEAF.NB(R_50_F.Analyses[["PSA"]]   , R_50_F.Analyses[["EVPI"]] , legend=c("TTE", "No TTE"), export="ceaf_R_50_F.jpeg" )
tmp <- Draw.CEAF.NB(R_50_M.Analyses[["PSA"]]   , R_50_M.Analyses[["EVPI"]] , legend=c("TTE", "No TTE"), export="ceaf_R_50_M.jpeg" )

tmp <- Draw.CEAF.NB(D_65_F.Analyses[["PSA"]]   , D_65_F.Analyses[["EVPI"]] , legend=c("TTE", "No TTE"), export="ceaf_D_65_F.jpeg" )
tmp <- Draw.CEAF.NB(D_65_M.Analyses[["PSA"]]   , D_65_M.Analyses[["EVPI"]] , legend=c("TTE", "No TTE"), export="ceaf_D_65_M.jpeg" )


write.csv(icer.results, file="MainIcer_QI_results.csv")
write.csv(otherQI.results, file="otherQI_results.csv")

write.csv(W_50_M.Analyses[["sensSpec"]],    file="SensSpecBlock_W_50_M.csv")
write.csv(W_50_F.Analyses[["sensSpec"]],    file="SensSpecBlock_W_50_F.csv")
write.csv(W_65_M.Analyses[["sensSpec"]],    file="SensSpecBlock_W_65_M.csv")
write.csv(W_65_F.Analyses[["sensSpec"]],    file="SensSpecBlock_W_65_F.csv")


icer.results
otherQI.results
evpi.results

################### Sensitivity analysis: effect of TPHR on CE of TTE cf No TTE for each patient group



trueHRprop <- seq(0,0.2, by = 0.0025)
IcerDF <- data.frame(
  hr=trueHRprop,
  w50f=NA,
  w50m=NA,
  w65f=NA,
  w65m=NA,
  r50f=NA,
  r50m=NA,
  r65f=NA,
  r65m=NA,
  d65f=NA,
  d65m=NA
  )

for (i in 1:length(trueHRprop)){
  
  W_50_F.Analyses <- Meta.Function(Data=Data.WarfC0_50_F, hrProp.Exp=trueHRprop[i], hrProp.PSA=rep(trueHRprop[i], 1000), sens.Exp=sens.Exp, spec.Exp=spec.Exp, sens.PSA=sens.PSA, spec.PSA=spec.PSA)
  W_50_M.Analyses <- Meta.Function(Data=Data.WarfC0_50_M, hrProp.Exp=trueHRprop[i], hrProp.PSA=rep(trueHRprop[i], 1000), sens.Exp=sens.Exp, spec.Exp=spec.Exp, sens.PSA=sens.PSA, spec.PSA=spec.PSA)
  
  W_65_F.Analyses <- Meta.Function(Data=Data.WarfC0_65_F, hrProp.Exp=trueHRprop[i], hrProp.PSA=rep(trueHRprop[i], 1000), sens.Exp=sens.Exp, spec.Exp=spec.Exp, sens.PSA=sens.PSA, spec.PSA=spec.PSA)
  W_65_M.Analyses <- Meta.Function(Data=Data.WarfC0_65_M, hrProp.Exp=trueHRprop[i], hrProp.PSA=rep(trueHRprop[i], 1000), sens.Exp=sens.Exp, spec.Exp=spec.Exp, sens.PSA=sens.PSA, spec.PSA=spec.PSA)
  
  R_50_F.Analyses <- Meta.Function(Data=Data.RivC0_50_F, hrProp.Exp=trueHRprop[i], hrProp.PSA=rep(trueHRprop[i], 1000), sens.Exp=sens.Exp, spec.Exp=spec.Exp, sens.PSA=sens.PSA, spec.PSA=spec.PSA)
  R_50_M.Analyses <- Meta.Function(Data=Data.RivC0_50_M, hrProp.Exp=trueHRprop[i], hrProp.PSA=rep(trueHRprop[i], 1000), sens.Exp=sens.Exp, spec.Exp=spec.Exp, sens.PSA=sens.PSA, spec.PSA=spec.PSA)
  
  R_65_F.Analyses <- Meta.Function(Data=Data.RivC0_65_F, hrProp.Exp=trueHRprop[i], hrProp.PSA=rep(trueHRprop[i], 1000), sens.Exp=sens.Exp, spec.Exp=spec.Exp, sens.PSA=sens.PSA, spec.PSA=spec.PSA)
  R_65_M.Analyses <- Meta.Function(Data=Data.RivC0_65_M, hrProp.Exp=trueHRprop[i], hrProp.PSA=rep(trueHRprop[i], 1000), sens.Exp=sens.Exp, spec.Exp=spec.Exp, sens.PSA=sens.PSA, spec.PSA=spec.PSA)
  
  D_65_F.Analyses <- Meta.Function(Data=Data.DabC0_65_F, hrProp.Exp=trueHRprop[i], hrProp.PSA=rep(trueHRprop[i], 1000), sens.Exp=sens.Exp, spec.Exp=spec.Exp, sens.PSA=sens.PSA, spec.PSA=spec.PSA)
  D_65_M.Analyses <- Meta.Function(Data=Data.DabC0_65_M, hrProp.Exp=trueHRprop[i], hrProp.PSA=rep(trueHRprop[i], 1000), sens.Exp=sens.Exp, spec.Exp=spec.Exp, sens.PSA=sens.PSA, spec.PSA=spec.PSA)
  
  tmp <- Extract.mainQI(W_50_F.Analyses)["mean.icer"] ; IcerDF[IcerDF$hr==trueHRprop[i],"w50f"] <- ifelse(tmp > 0, tmp, Inf)
  tmp <- Extract.mainQI(W_50_M.Analyses)["mean.icer"] ; IcerDF[IcerDF$hr==trueHRprop[i],"w50m"] <- ifelse(tmp > 0, tmp, Inf)

  tmp <- Extract.mainQI(W_65_F.Analyses)["mean.icer"] ; IcerDF[IcerDF$hr==trueHRprop[i],"w65f"] <- ifelse(tmp > 0, tmp, Inf)
  tmp <- Extract.mainQI(W_65_M.Analyses)["mean.icer"] ; IcerDF[IcerDF$hr==trueHRprop[i],"w65m"] <- ifelse(tmp > 0, tmp, Inf)

  tmp <- Extract.mainQI(R_50_F.Analyses)["mean.icer"] ; IcerDF[IcerDF$hr==trueHRprop[i],"r50f"] <- ifelse(tmp > 0, tmp, Inf)
  tmp <- Extract.mainQI(R_50_M.Analyses)["mean.icer"] ; IcerDF[IcerDF$hr==trueHRprop[i],"r50m"] <- ifelse(tmp > 0, tmp, Inf)
  
  tmp <- Extract.mainQI(R_65_F.Analyses)["mean.icer"] ; IcerDF[IcerDF$hr==trueHRprop[i],"r65f"] <- ifelse(tmp > 0, tmp, Inf)
  tmp <- Extract.mainQI(R_65_M.Analyses)["mean.icer"] ; IcerDF[IcerDF$hr==trueHRprop[i],"r65m"] <- ifelse(tmp > 0, tmp, Inf)

  tmp <- Extract.mainQI(D_65_F.Analyses)["mean.icer"] ; IcerDF[IcerDF$hr==trueHRprop[i],"d65f"] <- ifelse(tmp > 0, tmp, Inf)
  tmp <- Extract.mainQI(D_65_M.Analyses)["mean.icer"] ; IcerDF[IcerDF$hr==trueHRprop[i],"d65m"] <- ifelse(tmp > 0, tmp, Inf)
  print(i)  
}

options(scipen=6) # make sure scientific notation is not activated on y axis
# 50 year old females
jpeg("tphr50f.jpeg", width=800, height=800)
plot(w50f ~ hr, data=IcerDF, ylim=c(0, 60000), xlim=c(0,0.2), type="l", xlab="True proportion high risk", ylab="Mean ICER (£/QALY)")
lines(r50f ~ hr, data=IcerDF, lty="dashed")
legend("topleft", legend=c("Warfarin", "Rivaroxaban" ), lty=c(1,2))
dev.off()

jpeg("tphr65f.jpeg", width=800, height=800)
plot(w65f ~ hr, data=IcerDF, ylim=c(0, 60000), xlim=c(0, 0.2), type="l", xlab="True proportion high risk", ylab="Mean ICER (£/QALY)")
lines(r65f ~ hr, data=IcerDF, lty="dashed")
lines(d65f ~ hr, data=IcerDF, lwd=2)
legend("topright", legend=c("Warfarin", "Rivaroxaban", "Dabigatran"), lty=c(1,2,1), lwd=c(1,1,2))
dev.off()

# 50 year old males
jpeg("tphr50m.jpeg", width=800, height=800)
plot(w50m ~ hr, data=IcerDF, ylim=c(0, 60000), xlim=c(0,0.2), type="l", xlab="True proportion high risk", ylab="Mean ICER (£/QALY)")
lines(r50m ~ hr, data=IcerDF, lty="dashed")
legend("topleft", legend=c("Warfarin", "Rivaroxaban" ), lty=c(1,2))
dev.off()

jpeg("tphr65m.jpeg", width=800, height=800)
plot(w65m ~ hr, data=IcerDF, ylim=c(0, 60000), xlim=c(0, 0.2), type="l", xlab="True proportion high risk", ylab="Mean ICER (£/QALY)")
lines(r65m ~ hr, data=IcerDF, lty="dashed")
lines(d65m ~ hr, data=IcerDF, lwd=2)
legend("topright", legend=c("Warfarin", "Rivaroxaban", "Dabigatran"), lty=c(1,2,1), lwd=c(1,1,2))
dev.off()




plot(icerVec ~ trueHRprop, type="l", main="Relationship between mean ICER and true high risk proportion\nDabigatran", ylab="ICER (£/QALY)", ylim=c(0, 80000), xlab="True proportion high risk", lwd=2)

abline(v=0, col="grey")
abline(h=0, col="grey")

trueHRprop <- seq(0,0.2, by = 0.025)
icerVec <- vector("numeric", length(trueHRprop))

for (i in 1:length(trueHRprop)){
  X <- Calc.MeanCostQaly(Data.DabCV0)
  this.cost <- X$cost
  this.qaly <- X$qaly
  X <- Calc.Prop(trueHRprop[i], sens=sens.Exp, spec=spec.Exp)
  this.prop.base <- X$prop.base
  this.prop.comp <- X$prop.comp
  DabCV.Exp.icer <- Calc.CostQALY(prop.base=this.prop.base, prop.comp=this.prop.comp, cost=this.cost, qaly=this.qaly)
  icerVec[i] <- DabCV.Exp.icer$icer
}

lines(icerVec ~ trueHRprop, type="l", lwd=2, lty="dashed")


legend("topright", lty=c("solid", "dashed"), bty="n", inset=0.01, col=c("black", "black", "white", "white"), lwd=c(2,2), legend=c(expression(CHADS[2]),
                                                                                                                                  expression(paste(CHA[2], DS[2],-VASc), sep="")) )


##########################
#########################
## EVPPI
##########################################################################################

dta.Dab.C0 <- data.frame(
  TP.c=Data.DabC0$TP.c,
  FN.c=Data.DabC0$FN.c,
  TN.c=Data.DabC0$TN.c,
  FP.c=Data.DabC0$FP.c,

  TP.q=Data.DabC0$TP.q,
  FN.q=Data.DabC0$FN.q,
  TN.q=Data.DabC0$TN.q,
  FP.q=Data.DabC0$FP.q,
  PSA.inputs)

dta.Dab.CV0 <- data.frame(
  TP.c=Data.DabCV0$TP.c,
  FN.c=Data.DabCV0$FN.c,
  TN.c=Data.DabCV0$TN.c,
  FP.c=Data.DabCV0$FP.c,

  TP.q=Data.DabCV0$TP.q,
  FN.q=Data.DabCV0$FN.q,
  TN.q=Data.DabCV0$TN.q,
  FP.q=Data.DabCV0$FP.q,
  PSA.inputs)

dta.Riv.C0 <- data.frame(
  TP.c=Data.RivC0$TP.c,
  FN.c=Data.RivC0$FN.c,
  TN.c=Data.RivC0$TN.c,
  FP.c=Data.RivC0$FP.c,

  TP.q=Data.RivC0$TP.q,
  FN.q=Data.RivC0$FN.q,
  TN.q=Data.RivC0$TN.q,
  FP.q=Data.RivC0$FP.q,
  PSA.inputs)

dta.Riv.C1 <- data.frame(
  TP.c=Data.RivC1$TP.c,
  FN.c=Data.RivC1$FN.c,
  TN.c=Data.RivC1$TN.c,
  FP.c=Data.RivC1$FP.c,

  TP.q=Data.RivC1$TP.q,
  FN.q=Data.RivC1$FN.q,
  TN.q=Data.RivC1$TN.q,
  FP.q=Data.RivC1$FP.q,
  PSA.inputs)


dta.Riv.CV0 <- data.frame(
  TP.c=Data.RivCV0$TP.c,
  FN.c=Data.RivCV0$FN.c,
  TN.c=Data.RivCV0$TN.c,
  FP.c=Data.RivCV0$FP.c,

  TP.q=Data.RivCV0$TP.q,
  FN.q=Data.RivCV0$FN.q,
  TN.q=Data.RivCV0$TN.q,
  FP.q=Data.RivCV0$FP.q,
  PSA.inputs)

dta.Riv.CV1 <- data.frame(
  TP.c=Data.RivCV1$TP.c,
  FN.c=Data.RivCV1$FN.c,
  TN.c=Data.RivCV1$TN.c,
  FP.c=Data.RivCV1$FP.c,

  TP.q=Data.RivCV1$TP.q,
  FN.q=Data.RivCV1$FN.q,
  TN.q=Data.RivCV1$TN.q,
  FP.q=Data.RivCV1$FP.q,
  PSA.inputs)


dta.Warf.C0 <- data.frame(
  TP.c=Data.WarfC0$TP.c,
  FN.c=Data.WarfC0$FN.c,
  TN.c=Data.WarfC0$TN.c,
  FP.c=Data.WarfC0$FP.c,

  TP.q=Data.WarfC0$TP.q,
  FN.q=Data.WarfC0$FN.q,
  TN.q=Data.WarfC0$TN.q,
  FP.q=Data.WarfC0$FP.q,
  PSA.inputs)

dta.Warf.C1 <- data.frame(
  TP.c=Data.WarfC1$TP.c,
  FN.c=Data.WarfC1$FN.c,
  TN.c=Data.WarfC1$TN.c,
  FP.c=Data.WarfC1$FP.c,

  TP.q=Data.WarfC1$TP.q,
  FN.q=Data.WarfC1$FN.q,
  TN.q=Data.WarfC1$TN.q,
  FP.q=Data.WarfC1$FP.q,
  PSA.inputs)


dta.Warf.CV0 <- data.frame(
  TP.c=Data.WarfCV0$TP.c,
  FN.c=Data.WarfCV0$FN.c,
  TN.c=Data.WarfCV0$TN.c,
  FP.c=Data.WarfCV0$FP.c,

  TP.q=Data.WarfCV0$TP.q,
  FN.q=Data.WarfCV0$FN.q,
  TN.q=Data.WarfCV0$TN.q,
  FP.q=Data.WarfCV0$FP.q,
  PSA.inputs)

dta.Warf.CV1 <- data.frame(
  TP.c=Data.WarfCV1$TP.c,
  FN.c=Data.WarfCV1$FN.c,
  TN.c=Data.WarfCV1$TN.c,
  FP.c=Data.WarfCV1$FP.c,

  TP.q=Data.WarfCV1$TP.q,
  FN.q=Data.WarfCV1$FN.q,
  TN.q=Data.WarfCV1$TN.q,
  FP.q=Data.WarfCV1$FP.q,
  PSA.inputs)



MM.Dab.C0 <- lm(cbind(TP.c, FN.c, TN.c, FP.c, TP.q, FN.q, TN.q, FP.q) ~ 
##########MAIN INPUTS#########################
  Hr.Stroke.LAABN + Hr.Stroke.C0 +Hr.Stroke.C1 + Hr.Stroke.C2 + Hr.Stroke.C3 + Hr.Stroke.C4 + Prob.Death.stroke +
  Prob.DepState.stroke + Prob.IndState.stroke + Prob.Gos2.Ich +Prob.Gos3.Ich + Prob.Gos4.Ich +  Prob.Gos5.Ich + Util.DepState.Stroke +
  Util.IndState.Stroke + Util.Gos3 + Util.Gos4 + Util.Gos5 +  Util.Nich + Cost.Death.Stroke + Cost.DepState.Inst + Cost.DepState.Cont +
  Cost.IndState.Inst + Cost.IndState.Cont + Cost.Gos2.Inst  + Cost.Gos2.Cont  + Cost.Gos3.Inst  + Cost.Gos3.Cont  + Cost.Gos4.Inst  +
  Cost.Gos5.Inst  + Cost.Death.Bleed + Cost.Nich  + Sens  + Spec + 
##########################################################
  hrProp.C0 +         
### Dab Inputs##################################
 Hr.Bleed.lt75.Dab + Hr.Bleed.75over.Dab + Hr.Stroke.Dab + Prob.Death.Bleed.Dab + Prob.Ich.Bleed.Dab + Prob.Nich.Bleed.Dab,
                            data=dta.Dab.C0)

MM.Dab.CV0 <- lm(cbind(TP.c, FN.c, TN.c, FP.c, TP.q, FN.q, TN.q, FP.q) ~ 
##########MAIN INPUTS#########################
  Hr.Stroke.LAABN + Hr.Stroke.C0 +Hr.Stroke.C1 + Hr.Stroke.C2 + Hr.Stroke.C3 + Hr.Stroke.C4 + Prob.Death.stroke +
  Prob.DepState.stroke + Prob.IndState.stroke + Prob.Gos2.Ich +Prob.Gos3.Ich + Prob.Gos4.Ich +  Prob.Gos5.Ich + Util.DepState.Stroke +
  Util.IndState.Stroke + Util.Gos3 + Util.Gos4 + Util.Gos5 +  Util.Nich + Cost.Death.Stroke + Cost.DepState.Inst + Cost.DepState.Cont +
  Cost.IndState.Inst + Cost.IndState.Cont + Cost.Gos2.Inst  + Cost.Gos2.Cont  + Cost.Gos3.Inst  + Cost.Gos3.Cont  + Cost.Gos4.Inst  +
  Cost.Gos5.Inst  + Cost.Death.Bleed + Cost.Nich  + Sens  + Spec +
#########################################################
                         hrProp.CV0 +      
### Dab Inputs##################################
 Hr.Bleed.lt75.Dab + Hr.Bleed.75over.Dab + Hr.Stroke.Dab + Prob.Death.Bleed.Dab + Prob.Ich.Bleed.Dab + Prob.Nich.Bleed.Dab,
                            data=dta.Dab.CV0)

MM.Riv.C0 <- lm(cbind(TP.c, FN.c, TN.c, FP.c, TP.q, FN.q, TN.q, FP.q) ~
##########MAIN INPUTS#########################
  Hr.Stroke.LAABN + Hr.Stroke.C0 +Hr.Stroke.C1 + Hr.Stroke.C2 + Hr.Stroke.C3 + Hr.Stroke.C4 + Prob.Death.stroke +
  Prob.DepState.stroke + Prob.IndState.stroke + Prob.Gos2.Ich +Prob.Gos3.Ich + Prob.Gos4.Ich +  Prob.Gos5.Ich + Util.DepState.Stroke +
  Util.IndState.Stroke + Util.Gos3 + Util.Gos4 + Util.Gos5 +  Util.Nich + Cost.Death.Stroke + Cost.DepState.Inst + Cost.DepState.Cont +
  Cost.IndState.Inst + Cost.IndState.Cont + Cost.Gos2.Inst  + Cost.Gos2.Cont  + Cost.Gos3.Inst  + Cost.Gos3.Cont  + Cost.Gos4.Inst  +
  Cost.Gos5.Inst  + Cost.Death.Bleed + Cost.Nich  + Sens  + Spec + 
#########################################################
  hrProp.C0 +              
### Riv Inputs ################################
 Hr.Bleed.lt75.Riv + Hr.Bleed.75over.Riv + Hr.Stroke.Riv + Prob.Death.Bleed.Riv + Prob.Ich.Bleed.Riv + Prob.Nich.Bleed.Riv,
                            data=dta.Riv.C0)

MM.Riv.C1 <- lm(cbind(TP.c, FN.c, TN.c, FP.c, TP.q, FN.q, TN.q, FP.q) ~
##########MAIN INPUTS#########################
  Hr.Stroke.LAABN + Hr.Stroke.C0 +Hr.Stroke.C1 + Hr.Stroke.C2 + Hr.Stroke.C3 + Hr.Stroke.C4 + Prob.Death.stroke +
  Prob.DepState.stroke + Prob.IndState.stroke + Prob.Gos2.Ich +Prob.Gos3.Ich + Prob.Gos4.Ich +  Prob.Gos5.Ich + Util.DepState.Stroke +
  Util.IndState.Stroke + Util.Gos3 + Util.Gos4 + Util.Gos5 +  Util.Nich + Cost.Death.Stroke + Cost.DepState.Inst + Cost.DepState.Cont +
  Cost.IndState.Inst + Cost.IndState.Cont + Cost.Gos2.Inst  + Cost.Gos2.Cont  + Cost.Gos3.Inst  + Cost.Gos3.Cont  + Cost.Gos4.Inst  +
  Cost.Gos5.Inst  + Cost.Death.Bleed + Cost.Nich  + Sens  + Spec +
#########################################################
             hrProp.C1 +                
### Riv Inputs ################################
 Hr.Bleed.lt75.Riv + Hr.Bleed.75over.Riv + Hr.Stroke.Riv + Prob.Death.Bleed.Riv + Prob.Ich.Bleed.Riv + Prob.Nich.Bleed.Riv,
                            data=dta.Riv.C1)

MM.Riv.CV0 <- lm(cbind(TP.c, FN.c, TN.c, FP.c, TP.q, FN.q, TN.q, FP.q) ~
##########MAIN INPUTS#########################
  Hr.Stroke.LAABN + Hr.Stroke.C0 +Hr.Stroke.C1 + Hr.Stroke.C2 + Hr.Stroke.C3 + Hr.Stroke.C4 + Prob.Death.stroke +
  Prob.DepState.stroke + Prob.IndState.stroke + Prob.Gos2.Ich +Prob.Gos3.Ich + Prob.Gos4.Ich +  Prob.Gos5.Ich + Util.DepState.Stroke +
  Util.IndState.Stroke + Util.Gos3 + Util.Gos4 + Util.Gos5 +  Util.Nich + Cost.Death.Stroke + Cost.DepState.Inst + Cost.DepState.Cont +
  Cost.IndState.Inst + Cost.IndState.Cont + Cost.Gos2.Inst  + Cost.Gos2.Cont  + Cost.Gos3.Inst  + Cost.Gos3.Cont  + Cost.Gos4.Inst  +
  Cost.Gos5.Inst  + Cost.Death.Bleed + Cost.Nich  + Sens  + Spec +
##############################################
                      hrProp.CV0 +            
### Riv Inputs ################################
 Hr.Bleed.lt75.Riv + Hr.Bleed.75over.Riv + Hr.Stroke.Riv + Prob.Death.Bleed.Riv + Prob.Ich.Bleed.Riv + Prob.Nich.Bleed.Riv,
                            data=dta.Riv.CV0)

MM.Riv.CV1 <- lm(cbind(TP.c, FN.c, TN.c, FP.c, TP.q, FN.q, TN.q, FP.q) ~
##########MAIN INPUTS#########################
  Hr.Stroke.LAABN + Hr.Stroke.C0 +Hr.Stroke.C1 + Hr.Stroke.C2 + Hr.Stroke.C3 + Hr.Stroke.C4 + Prob.Death.stroke +
  Prob.DepState.stroke + Prob.IndState.stroke + Prob.Gos2.Ich +Prob.Gos3.Ich + Prob.Gos4.Ich +  Prob.Gos5.Ich + Util.DepState.Stroke +
  Util.IndState.Stroke + Util.Gos3 + Util.Gos4 + Util.Gos5 +  Util.Nich + Cost.Death.Stroke + Cost.DepState.Inst + Cost.DepState.Cont +
  Cost.IndState.Inst + Cost.IndState.Cont + Cost.Gos2.Inst  + Cost.Gos2.Cont  + Cost.Gos3.Inst  + Cost.Gos3.Cont  + Cost.Gos4.Inst  +
  Cost.Gos5.Inst  + Cost.Death.Bleed + Cost.Nich  + Sens  + Spec +
##############################################
                                       hrProp.CV1 +                  
### Riv Inputs ################################
 Hr.Bleed.lt75.Riv + Hr.Bleed.75over.Riv + Hr.Stroke.Riv + Prob.Death.Bleed.Riv + Prob.Ich.Bleed.Riv + Prob.Nich.Bleed.Riv,
                            data=dta.Riv.CV1)

MM.Warf.C0 <- lm(cbind(TP.c, FN.c, TN.c, FP.c, TP.q, FN.q, TN.q, FP.q) ~
##########MAIN INPUTS#########################
  Hr.Stroke.LAABN + Hr.Stroke.C0 +Hr.Stroke.C1 + Hr.Stroke.C2 + Hr.Stroke.C3 + Hr.Stroke.C4 + Prob.Death.stroke +
  Prob.DepState.stroke + Prob.IndState.stroke + Prob.Gos2.Ich +Prob.Gos3.Ich + Prob.Gos4.Ich +  Prob.Gos5.Ich + Util.DepState.Stroke +
  Util.IndState.Stroke + Util.Gos3 + Util.Gos4 + Util.Gos5 +  Util.Nich + Cost.Death.Stroke + Cost.DepState.Inst + Cost.DepState.Cont +
  Cost.IndState.Inst + Cost.IndState.Cont + Cost.Gos2.Inst  + Cost.Gos2.Cont  + Cost.Gos3.Inst  + Cost.Gos3.Cont  + Cost.Gos4.Inst  +
  Cost.Gos5.Inst  + Cost.Death.Bleed + Cost.Nich  + Sens  + Spec + 
##############################################
  hrProp.C0 +                  
### Warf Inputs##################################
 Hr.Bleed.lt75.War + Hr.Bleed.75over.War + Hr.Stroke.War + Prob.Death.Bleed.War + Prob.Ich.Bleed.War  + Prob.Nich.Bleed.War + Cost.War,
                            data=dta.Warf.C0)

MM.Warf.C1 <- lm(cbind(TP.c, FN.c, TN.c, FP.c, TP.q, FN.q, TN.q, FP.q) ~
##########MAIN INPUTS#########################
  Hr.Stroke.LAABN + Hr.Stroke.C0 +Hr.Stroke.C1 + Hr.Stroke.C2 + Hr.Stroke.C3 + Hr.Stroke.C4 + Prob.Death.stroke +
  Prob.DepState.stroke + Prob.IndState.stroke + Prob.Gos2.Ich +Prob.Gos3.Ich + Prob.Gos4.Ich +  Prob.Gos5.Ich + Util.DepState.Stroke +
  Util.IndState.Stroke + Util.Gos3 + Util.Gos4 + Util.Gos5 +  Util.Nich + Cost.Death.Stroke + Cost.DepState.Inst + Cost.DepState.Cont +
  Cost.IndState.Inst + Cost.IndState.Cont + Cost.Gos2.Inst  + Cost.Gos2.Cont  + Cost.Gos3.Inst  + Cost.Gos3.Cont  + Cost.Gos4.Inst  +
  Cost.Gos5.Inst  + Cost.Death.Bleed + Cost.Nich  + Sens  + Spec  + 
###################################################################                 
                hrProp.C1 +        
### Warf Inputs##################################
 Hr.Bleed.lt75.War + Hr.Bleed.75over.War + Hr.Stroke.War + Prob.Death.Bleed.War + Prob.Ich.Bleed.War  + Prob.Nich.Bleed.War + Cost.War,
                            data=dta.Warf.C1)

MM.Warf.CV0 <- lm(cbind(TP.c, FN.c, TN.c, FP.c, TP.q, FN.q, TN.q, FP.q) ~
##########MAIN INPUTS#########################
  Hr.Stroke.LAABN + Hr.Stroke.C0 +Hr.Stroke.C1 + Hr.Stroke.C2 + Hr.Stroke.C3 + Hr.Stroke.C4 + Prob.Death.stroke +
  Prob.DepState.stroke + Prob.IndState.stroke + Prob.Gos2.Ich +Prob.Gos3.Ich + Prob.Gos4.Ich +  Prob.Gos5.Ich + Util.DepState.Stroke +
  Util.IndState.Stroke + Util.Gos3 + Util.Gos4 + Util.Gos5 +  Util.Nich + Cost.Death.Stroke + Cost.DepState.Inst + Cost.DepState.Cont +
  Cost.IndState.Inst + Cost.IndState.Cont + Cost.Gos2.Inst  + Cost.Gos2.Cont  + Cost.Gos3.Inst  + Cost.Gos3.Cont  + Cost.Gos4.Inst  +
  Cost.Gos5.Inst  + Cost.Death.Bleed + Cost.Nich  + Sens  + Spec  + 
###############################################                  
                            hrProp.CV0 +        
### Warf Inputs##################################
 Hr.Bleed.lt75.War + Hr.Bleed.75over.War + Hr.Stroke.War + Prob.Death.Bleed.War + Prob.Ich.Bleed.War  + Prob.Nich.Bleed.War + Cost.War,
                            data=dta.Warf.C0)

MM.Warf.CV1 <- lm(cbind(TP.c, FN.c, TN.c, FP.c, TP.q, FN.q, TN.q, FP.q) ~
##########MAIN INPUTS#########################
  Hr.Stroke.LAABN + Hr.Stroke.C0 +Hr.Stroke.C1 + Hr.Stroke.C2 + Hr.Stroke.C3 + Hr.Stroke.C4 + Prob.Death.stroke +
  Prob.DepState.stroke + Prob.IndState.stroke + Prob.Gos2.Ich +Prob.Gos3.Ich + Prob.Gos4.Ich +  Prob.Gos5.Ich + Util.DepState.Stroke +
  Util.IndState.Stroke + Util.Gos3 + Util.Gos4 + Util.Gos5 +  Util.Nich + Cost.Death.Stroke + Cost.DepState.Inst + Cost.DepState.Cont +
  Cost.IndState.Inst + Cost.IndState.Cont + Cost.Gos2.Inst  + Cost.Gos2.Cont  + Cost.Gos3.Inst  + Cost.Gos3.Cont  + Cost.Gos4.Inst  +
  Cost.Gos5.Inst  + Cost.Death.Bleed + Cost.Nich  + Sens  + Spec  +
########################################  
                                      hrProp.CV1 +          
### Warf Inputs##################################
 Hr.Bleed.lt75.War + Hr.Bleed.75over.War + Hr.Stroke.War + Prob.Death.Bleed.War + Prob.Ich.Bleed.War  + Prob.Nich.Bleed.War + Cost.War,
                            data=dta.Warf.C1)

EVPPI.Dab.C0 <- Calc.EVPPI(InData=PSA.inputs, mod=MM.Dab.C0, hrProp.option="C0", 
                           params=c(
                               "Hr.Stroke.LAABN", "Hr.Stroke.C0", "Hr.Stroke.C1", "Hr.Stroke.C2",
                               "Hr.Stroke.C3", "Hr.Stroke.C4", "Prob.Death.stroke",
                                "Prob.DepState.stroke","Prob.IndState.stroke", "Prob.Gos2.Ich", "Prob.Gos3.Ich","Prob.Gos4.Ich", "Prob.Gos5.Ich", 
                               "Util.DepState.Stroke", "Util.IndState.Stroke", "Util.Gos3", "Util.Gos4", "Util.Gos5", "Util.Nich", 
                               "Cost.Death.Stroke", "Cost.DepState.Inst",
                               "Cost.DepState.Cont", "Cost.IndState.Inst", "Cost.IndState.Cont", "Cost.Gos2.Inst", "Cost.Gos2.Cont",
                               "Cost.Gos3.Inst", "Cost.Gos3.Cont", "Cost.Gos4.Inst", "Cost.Gos5.Inst", "Cost.Death.Bleed",
                               "Cost.Nich", "hrProp.C0"
                               )
                           )

EVPPI.Dab.CV0 <- Calc.EVPPI(InData=PSA.inputs, mod=MM.Dab.CV0, hrProp.option="CV0", 
                           params=c(
                               "Hr.Stroke.LAABN", "Hr.Stroke.C0", "Hr.Stroke.C1", "Hr.Stroke.C2",
                               "Hr.Stroke.C3", "Hr.Stroke.C4", "Prob.Death.stroke",
                                "Prob.DepState.stroke","Prob.IndState.stroke", "Prob.Gos2.Ich", "Prob.Gos3.Ich","Prob.Gos4.Ich", "Prob.Gos5.Ich", 
                               "Util.DepState.Stroke", "Util.IndState.Stroke", "Util.Gos3", "Util.Gos4", "Util.Gos5", "Util.Nich", 
                               "Cost.Death.Stroke", "Cost.DepState.Inst",
                               "Cost.DepState.Cont", "Cost.IndState.Inst", "Cost.IndState.Cont", "Cost.Gos2.Inst", "Cost.Gos2.Cont",
                               "Cost.Gos3.Inst", "Cost.Gos3.Cont", "Cost.Gos4.Inst", "Cost.Gos5.Inst", "Cost.Death.Bleed",
                               "Cost.Nich",              "hrProp.CV0"
                               )
                            )

EVPPI.Dab.C0.hrProp <- Calc.EVPPI(InData=PSA.inputs, mod=MM.Dab.C0, hrProp.option="C0", 
                           params=c("hrProp.C0")
                           )

EVPPI.Dab.CV0.hrProp <- Calc.EVPPI(InData=PSA.inputs, mod=MM.Dab.CV0, hrProp.option="CV0", 
                           params=c("hrProp.CV0")
                            )

EVPPI.Dab.C0.JointSensSpec <- Calc.EVPPI(InData=PSA.inputs, mod=MM.Dab.C0, hrProp.option="C0", 
                           params=NULL, joint.SensSpec=T)

EVPPI.Dab.CV0.JointSensSpec <- Calc.EVPPI(InData=PSA.inputs, mod=MM.Dab.CV0, hrProp.option="CV0", 
                           params=NULL, joint.SensSpec=T)


EVPPI.Warf.C0 <- Calc.EVPPI(InData=PSA.inputs, mod=MM.Warf.C0, hrProp.option="C0", 
                           params=c(
                               "Hr.Stroke.LAABN", "Hr.Stroke.C0", "Hr.Stroke.C1", "Hr.Stroke.C2",
                               "Hr.Stroke.C3", "Hr.Stroke.C4", "Prob.Death.stroke",
                                "Prob.DepState.stroke","Prob.IndState.stroke", "Prob.Gos2.Ich", "Prob.Gos3.Ich","Prob.Gos4.Ich", "Prob.Gos5.Ich", 
                               "Util.DepState.Stroke", "Util.IndState.Stroke", "Util.Gos3", "Util.Gos4", "Util.Gos5", "Util.Nich", 
                               "Cost.Death.Stroke", "Cost.DepState.Inst",
                               "Cost.DepState.Cont", "Cost.IndState.Inst", "Cost.IndState.Cont", "Cost.Gos2.Inst", "Cost.Gos2.Cont",
                               "Cost.Gos3.Inst", "Cost.Gos3.Cont", "Cost.Gos4.Inst", "Cost.Gos5.Inst", "Cost.Death.Bleed",
                               "Cost.Nich", "hrProp.C0"
                               )
                           )

EVPPI.Warf.C1 <- Calc.EVPPI(InData=PSA.inputs, mod=MM.Warf.C1, hrProp.option="C1", 
                           params=c(
                               "Hr.Stroke.LAABN", "Hr.Stroke.C0", "Hr.Stroke.C1", "Hr.Stroke.C2",
                               "Hr.Stroke.C3", "Hr.Stroke.C4", "Prob.Death.stroke",
                                "Prob.DepState.stroke","Prob.IndState.stroke", "Prob.Gos2.Ich", "Prob.Gos3.Ich","Prob.Gos4.Ich", "Prob.Gos5.Ich", 
                               "Util.DepState.Stroke", "Util.IndState.Stroke", "Util.Gos3", "Util.Gos4", "Util.Gos5", "Util.Nich", 
                               "Cost.Death.Stroke", "Cost.DepState.Inst",
                               "Cost.DepState.Cont", "Cost.IndState.Inst", "Cost.IndState.Cont", "Cost.Gos2.Inst", "Cost.Gos2.Cont",
                               "Cost.Gos3.Inst", "Cost.Gos3.Cont", "Cost.Gos4.Inst", "Cost.Gos5.Inst", "Cost.Death.Bleed",
                               "Cost.Nich", "hrProp.C1"
                               )
                           )

EVPPI.Warf.CV0 <- Calc.EVPPI(InData=PSA.inputs, mod=MM.Warf.CV0, hrProp.option="CV0", 
                           params=c(
                               "Hr.Stroke.LAABN", "Hr.Stroke.C0", "Hr.Stroke.C1", "Hr.Stroke.C2",
                               "Hr.Stroke.C3", "Hr.Stroke.C4", "Prob.Death.stroke",
                                "Prob.DepState.stroke","Prob.IndState.stroke", "Prob.Gos2.Ich", "Prob.Gos3.Ich","Prob.Gos4.Ich", "Prob.Gos5.Ich", 
                               "Util.DepState.Stroke", "Util.IndState.Stroke", "Util.Gos3", "Util.Gos4", "Util.Gos5", "Util.Nich", 
                               "Cost.Death.Stroke", "Cost.DepState.Inst",
                               "Cost.DepState.Cont", "Cost.IndState.Inst", "Cost.IndState.Cont", "Cost.Gos2.Inst", "Cost.Gos2.Cont",
                               "Cost.Gos3.Inst", "Cost.Gos3.Cont", "Cost.Gos4.Inst", "Cost.Gos5.Inst", "Cost.Death.Bleed",
                               "Cost.Nich", "hrProp.CV0"
                               )
                            )

EVPPI.Warf.CV1 <- Calc.EVPPI(InData=PSA.inputs, mod=MM.Warf.CV1, hrProp.option="CV1", 
                           params=c(
                               "Hr.Stroke.LAABN",       "Hr.Stroke.C0",         "Hr.Stroke.C1",        "Hr.Stroke.C2",
                               "Hr.Stroke.C3",          "Hr.Stroke.C4",         "Prob.Death.stroke",
                               "Prob.DepState.stroke",  "Prob.IndState.stroke", "Prob.Gos2.Ich",       "Prob.Gos3.Ich",  "Prob.Gos4.Ich", "Prob.Gos5.Ich", 
                               "Util.DepState.Stroke",  "Util.IndState.Stroke", "Util.Gos3",           "Util.Gos4",      "Util.Gos5",     "Util.Nich", 
                               "Cost.Death.Stroke",     "Cost.DepState.Inst",
                               "Cost.DepState.Cont",    "Cost.IndState.Inst",    "Cost.IndState.Cont", "Cost.Gos2.Inst", "Cost.Gos2.Cont",
                               "Cost.Gos3.Inst",        "Cost.Gos3.Cont",        "Cost.Gos4.Inst",     "Cost.Gos5.Inst", "Cost.Death.Bleed",
                               "Cost.Nich",             "hrProp.CV1"
                               )
                            )


EVPPI.Dab.C0.hrProp <- Calc.EVPPI(InData=PSA.inputs, mod=MM.Dab.C0, hrProp.option="C0", 
                           params=c("hrProp.C0")
                           )

EVPPI.Dab.CV0.hrProp <- Calc.EVPPI(InData=PSA.inputs, mod=MM.Dab.CV0, hrProp.option="CV0", 
                           params=c("hrProp.CV0")
                            )

EVPPI.Dab.C0.JointSensSpec <- Calc.EVPPI(InData=PSA.inputs, mod=MM.Dab.C0, hrProp.option="C0", 
                           params=NULL, joint.SensSpec=T)

EVPPI.Dab.CV0.JointSensSpec <- Calc.EVPPI(InData=PSA.inputs, mod=MM.Dab.CV0, hrProp.option="CV0", 
                           params=NULL, joint.SensSpec=T)

##########################################
##########################################

EVPPI.Warf.C0.hrProp <- Calc.EVPPI(InData=PSA.inputs, mod=MM.Warf.C0, hrProp.option="C0", 
                           params=c("hrProp.C0")
                           )

EVPPI.Warf.C1.hrProp <- Calc.EVPPI(InData=PSA.inputs, mod=MM.Warf.C1, hrProp.option="C1", 
                           params=c("hrProp.C1")
                           )

EVPPI.Warf.CV0.hrProp <- Calc.EVPPI(InData=PSA.inputs, mod=MM.Warf.CV0, hrProp.option="CV0", 
                           params=c("hrProp.CV0")
                            )

EVPPI.Warf.CV1.hrProp <- Calc.EVPPI(InData=PSA.inputs, mod=MM.Warf.CV1, hrProp.option="CV1", 
                           params=c("hrProp.CV1")
                            )

EVPPI.Warf.C0.JointSensSpec <- Calc.EVPPI(InData=PSA.inputs, mod=MM.Warf.C0, hrProp.option="C0", 
                           params=NULL, joint.SensSpec=T)
EVPPI.Warf.C1.JointSensSpec <- Calc.EVPPI(InData=PSA.inputs, mod=MM.Warf.C0, hrProp.option="C1", 
                           params=NULL, joint.SensSpec=T)


EVPPI.Warf.CV0.JointSensSpec <- Calc.EVPPI(InData=PSA.inputs, mod=MM.Warf.CV0, hrProp.option="CV0", 
                           params=NULL, joint.SensSpec=T)

EVPPI.Warf.CV1.JointSensSpec <- Calc.EVPPI(InData=PSA.inputs, mod=MM.Warf.CV1, hrProp.option="CV1", 
                           params=NULL, joint.SensSpec=T)


save(EVPPI.Dab.C0, 
     EVPPI.Dab.CV0, 
     EVPPI.Dab.C0.JointSensSpec,
     EVPPI.Dab.CV0.JointSensSpec,
     file="EVPPI_Dabigatran.rData")

tmp <- Draw.EVPI(EVPPI.Dab.C0.JointSensSpec[["Sens_and_Spec"]])
tmp2 <- Draw.EVPI(EVPPI.Dab.CV0.JointSensSpec[["Sens_and_Spec"]])

tmp3 <- Draw.EVPI(EVPPI.Dab.C0.hrProp[["hrProp.C0"]])
tmp4 <- Draw.EVPI(EVPPI.Dab.CV0.hrProp[["hrProp.CV0"]])



#######################
Dab.hrIcer <-  HrProp.Icer.graph(
               D.C=  Data.DabC0,
               D.CV= Data.DabCV0,
               sens.PSA=sens.PSA,
               spec.PSA=spec.PSA,
               main=""
               )

Warf.0.hrIcer <-  HrProp.Icer.graph(
               D.C=  Data.WarfC0,
               D.CV= Data.WarfCV0,
               sens.PSA=sens.PSA,
               spec.PSA=spec.PSA,
               main=""
               )


Warf.1.hrIcer <-  HrProp.Icer.graph(
               D.C=  Data.WarfC1,
               D.CV= Data.WarfCV1,
               sens.PSA=sens.PSA,
               spec.PSA=spec.PSA,
               main=""
               )

Riv.0.hrIcer <-  HrProp.Icer.graph(
               D.C=  Data.RivC0,
               D.CV= Data.RivCV0,
               sens.PSA=sens.PSA,
               spec.PSA=spec.PSA,
               main=""
               )


Riv.1.hrIcer <-  HrProp.Icer.graph(
               D.C=  Data.RivC1,
               D.CV= Data.RivCV1,
               sens.PSA=sens.PSA,
               spec.PSA=spec.PSA,
               main=""
               )

with(Dab.hrIcer, which)


#####################################################################################################################################
#####################################################################################################################################

#plot(predict(MetaModel.Dab.C0.qaly) ~ Dab.C0.Analyses[["PSA"]]$qaly)
#plot(predict(MetaModel.Dab.C0.cost) ~ Dab.C0.Analyses[["PSA"]]$cost)

#cor(predict(MetaModel.Dab.C0.cost), Dab.C0.Analyses[["PSA"]]$cost)
# correlation of 0.94

#cor(predict(MetaModel.Dab.C0.qaly), Dab.C0.Analyses[["PSA"]]$qaly)
# correlation of 0.99


#Calc.EVPPI <- function(D, mod, params,  Nruns=1000, maicers=seq(0,50000, by=1000)){
#test <- D=dta, mod=MetaModel.Dab.C0.qaly  
  
############################################################################################
MetaModel.War.C0.qaly <- lm(qaly ~ 
##########MAIN INPUTS#########################
  Hr.Stroke.LAABN + Hr.Stroke.C0 +Hr.Stroke.C1 + Hr.Stroke.C2 + Hr.Stroke.C3 + Hr.Stroke.C4 + Prob.Death.stroke +
  Prob.DepState.stroke + Prob.IndState.stroke + Prob.Gos2.Ich +Prob.Gos3.Ich + Prob.Gos4.Ich +  Prob.Gos5.Ich + Util.DepState.Stroke +
  Util.IndState.Stroke + Util.Gos3 + Util.Gos4 + Util.Gos5 +  Util.Nich + Cost.Death.Stroke + Cost.DepState.Inst + Cost.DepState.Cont +
  Cost.IndState.Inst + Cost.IndState.Cont + Cost.Gos2.Inst  + Cost.Gos2.Cont  + Cost.Gos3.Inst  + Cost.Gos3.Cont  + Cost.Gos4.Inst  +
  Cost.Gos5.Inst  + Cost.Death.Bleed + Cost.Nich  + Sens  + Spec  + hrProp.C0 + hrProp.C1 + hrProp.CV0 + hrProp.CV1 +          
### War Inputs##################################
 Hr.Bleed.lt75.War + Hr.Bleed.75over.War + Hr.Stroke.War + Prob.Death.Bleed.War + Prob.Ich.Bleed.War  + Prob.Nich.Bleed.War + Cost.War,
                            data=dta)

dta <- data.frame(cost=WarfC1.PSA$cost, qaly=WarfC1.PSA.qaly, PSA.inputs)

MetaModel.War.C1.cost <- lm(cost ~ 
  ##########MAIN INPUTS#########################
  Hr.Stroke.LAABN + Hr.Stroke.C0 +Hr.Stroke.C1 + Hr.Stroke.C2 + Hr.Stroke.C3 + Hr.Stroke.C4 + Prob.Death.stroke +
  Prob.DepState.stroke + Prob.IndState.stroke + Prob.Gos2.Ich +Prob.Gos3.Ich + Prob.Gos4.Ich +  Prob.Gos5.Ich + Util.DepState.Stroke +
  Util.IndState.Stroke + Util.Gos3 + Util.Gos4 + Util.Gos5 +  Util.Nich + Cost.Death.Stroke + Cost.DepState.Inst + Cost.DepState.Cont +
  Cost.IndState.Inst + Cost.IndState.Cont + Cost.Gos2.Inst  + Cost.Gos2.Cont  + Cost.Gos3.Inst  + Cost.Gos3.Cont  + Cost.Gos4.Inst  +
  Cost.Gos5.Inst  + Cost.Death.Bleed + Cost.Nich  + Sens  + Spec  + hrProp.C0 + hrProp.C1 + hrProp.CV0 + hrProp.CV1 +          
### War Inputs##################################
 Hr.Bleed.lt75.War + Hr.Bleed.75over.War + Hr.Stroke.War + Prob.Death.Bleed.War + Prob.Ich.Bleed.War  + Prob.Nich.Bleed.War + Cost.War,
                          data=dta)
MetaModel.War.C1.qaly <- lm(qaly ~ 
##########MAIN INPUTS#########################
  Hr.Stroke.LAABN + Hr.Stroke.C0 +Hr.Stroke.C1 + Hr.Stroke.C2 + Hr.Stroke.C3 + Hr.Stroke.C4 + Prob.Death.stroke +
  Prob.DepState.stroke + Prob.IndState.stroke + Prob.Gos2.Ich +Prob.Gos3.Ich + Prob.Gos4.Ich +  Prob.Gos5.Ich + Util.DepState.Stroke +
  Util.IndState.Stroke + Util.Gos3 + Util.Gos4 + Util.Gos5 +  Util.Nich + Cost.Death.Stroke + Cost.DepState.Inst + Cost.DepState.Cont +
  Cost.IndState.Inst + Cost.IndState.Cont + Cost.Gos2.Inst  + Cost.Gos2.Cont  + Cost.Gos3.Inst  + Cost.Gos3.Cont  + Cost.Gos4.Inst  +
  Cost.Gos5.Inst  + Cost.Death.Bleed + Cost.Nich  + Sens  + Spec  + hrProp.C0 + hrProp.C1 + hrProp.CV0 + hrProp.CV1 +          
### War Inputs##################################
 Hr.Bleed.lt75.War + Hr.Bleed.75over.War + Hr.Stroke.War + Prob.Death.Bleed.War + Prob.Ich.Bleed.War  + Prob.Nich.Bleed.War + Cost.War,
                            data=dta)


dta <- data.frame(cost=WarfCV0.PSA$cost, qaly=WarfCV0.PSA.qaly, PSA.inputs)
MetaModel.War.CV0.cost <- lm(cost ~ 
##########MAIN INPUTS#########################
  Hr.Stroke.LAABN + Hr.Stroke.C0 +Hr.Stroke.C1 + Hr.Stroke.C2 + Hr.Stroke.C3 + Hr.Stroke.C4 + Prob.Death.stroke +
  Prob.DepState.stroke + Prob.IndState.stroke + Prob.Gos2.Ich +Prob.Gos3.Ich + Prob.Gos4.Ich +  Prob.Gos5.Ich + Util.DepState.Stroke +
  Util.IndState.Stroke + Util.Gos3 + Util.Gos4 + Util.Gos5 +  Util.Nich + Cost.Death.Stroke + Cost.DepState.Inst + Cost.DepState.Cont +
  Cost.IndState.Inst + Cost.IndState.Cont + Cost.Gos2.Inst  + Cost.Gos2.Cont  + Cost.Gos3.Inst  + Cost.Gos3.Cont  + Cost.Gos4.Inst  +
  Cost.Gos5.Inst  + Cost.Death.Bleed + Cost.Nich  + Sens  + Spec  + hrProp.C0 + hrProp.C1 + hrProp.CV0 + hrProp.CV1 +          
### War Inputs##################################
 Hr.Bleed.lt75.War + Hr.Bleed.75over.War + Hr.Stroke.War + Prob.Death.Bleed.War + Prob.Ich.Bleed.War  + Prob.Nich.Bleed.War + Cost.War,
                             data=dta)
MetaModel.War.CV0.qaly <- lm(qaly ~ 
##########MAIN INPUTS#########################
  Hr.Stroke.LAABN + Hr.Stroke.C0 +Hr.Stroke.C1 + Hr.Stroke.C2 + Hr.Stroke.C3 + Hr.Stroke.C4 + Prob.Death.stroke +
  Prob.DepState.stroke + Prob.IndState.stroke + Prob.Gos2.Ich +Prob.Gos3.Ich + Prob.Gos4.Ich +  Prob.Gos5.Ich + Util.DepState.Stroke +
  Util.IndState.Stroke + Util.Gos3 + Util.Gos4 + Util.Gos5 +  Util.Nich + Cost.Death.Stroke + Cost.DepState.Inst + Cost.DepState.Cont +
  Cost.IndState.Inst + Cost.IndState.Cont + Cost.Gos2.Inst  + Cost.Gos2.Cont  + Cost.Gos3.Inst  + Cost.Gos3.Cont  + Cost.Gos4.Inst  +
  Cost.Gos5.Inst  + Cost.Death.Bleed + Cost.Nich  + Sens  + Spec  + hrProp.C0 + hrProp.C1 + hrProp.CV0 + hrProp.CV1 +          
### War Inputs##################################
 Hr.Bleed.lt75.War + Hr.Bleed.75over.War + Hr.Stroke.War + Prob.Death.Bleed.War + Prob.Ich.Bleed.War  + Prob.Nich.Bleed.War + Cost.War,
                           data=dta)

dta <- data.frame(cost=WarfCV1.PSA$cost, qaly=WarfCV1.PSA.qaly, PSA.inputs)
MetaModel.War.CV1.cost <- lm(cost ~ 
  ##########MAIN INPUTS#########################
  Hr.Stroke.LAABN + Hr.Stroke.C0 +Hr.Stroke.C1 + Hr.Stroke.C2 + Hr.Stroke.C3 + Hr.Stroke.C4 + Prob.Death.stroke +
  Prob.DepState.stroke + Prob.IndState.stroke + Prob.Gos2.Ich +Prob.Gos3.Ich + Prob.Gos4.Ich +  Prob.Gos5.Ich + Util.DepState.Stroke +
  Util.IndState.Stroke + Util.Gos3 + Util.Gos4 + Util.Gos5 +  Util.Nich + Cost.Death.Stroke + Cost.DepState.Inst + Cost.DepState.Cont +
  Cost.IndState.Inst + Cost.IndState.Cont + Cost.Gos2.Inst  + Cost.Gos2.Cont  + Cost.Gos3.Inst  + Cost.Gos3.Cont  + Cost.Gos4.Inst  +
  Cost.Gos5.Inst  + Cost.Death.Bleed + Cost.Nich  + Sens  + Spec  + hrProp.C0 + hrProp.C1 + hrProp.CV0 + hrProp.CV1 +          
### War Inputs##################################
 Hr.Bleed.lt75.War + Hr.Bleed.75over.War + Hr.Stroke.War + Prob.Death.Bleed.War + Prob.Ich.Bleed.War  + Prob.Nich.Bleed.War + Cost.War,
                             data=dta)
MetaModel.War.CV1.qaly <- lm(qaly ~
##########MAIN INPUTS#########################
  Hr.Stroke.LAABN + Hr.Stroke.C0 +Hr.Stroke.C1 + Hr.Stroke.C2 + Hr.Stroke.C3 + Hr.Stroke.C4 + Prob.Death.stroke +
  Prob.DepState.stroke + Prob.IndState.stroke + Prob.Gos2.Ich +Prob.Gos3.Ich + Prob.Gos4.Ich +  Prob.Gos5.Ich + Util.DepState.Stroke +
  Util.IndState.Stroke + Util.Gos3 + Util.Gos4 + Util.Gos5 +  Util.Nich + Cost.Death.Stroke + Cost.DepState.Inst + Cost.DepState.Cont +
  Cost.IndState.Inst + Cost.IndState.Cont + Cost.Gos2.Inst  + Cost.Gos2.Cont  + Cost.Gos3.Inst  + Cost.Gos3.Cont  + Cost.Gos4.Inst  +
  Cost.Gos5.Inst  + Cost.Death.Bleed + Cost.Nich  + Sens  + Spec  + hrProp.C0 + hrProp.C1 + hrProp.CV0 + hrProp.CV1 +          
### War Inputs##################################
 Hr.Bleed.lt75.War + Hr.Bleed.75over.War + Hr.Stroke.War + Prob.Death.Bleed.War + Prob.Ich.Bleed.War  + Prob.Nich.Bleed.War + Cost.War,
                             data=dta)

#####################################################
#####################################################

main.inputs <- c(
  "Hr.Stroke.LAABN" +
  "Hr.Stroke.C0" +
  "Hr.Stroke.C1" +
  "Hr.Stroke.C2"  +
  "Hr.Stroke.C3"  +
  "Hr.Stroke.C4"  +
#  "Hr.Bleed.lt75.Dab" +
#  "Hr.Bleed.75over.Dab" +
# "Hr.Bleed.lt75.War" +
# "Hr.Bleed.75over.War"  +
# "Hr.Bleed.lt75.Riv" +
# "Hr.Bleed.75over.Riv" + 
#  "Hr.Stroke.Dab" +
#  "Hr.Stroke.War" +
#  "Hr.Stroke.Riv" +  
  "Prob.Death.stroke" +
  "Prob.DepState.stroke" +
  "Prob.IndState.stroke" +
# "Prob.Death.Bleed.Dab" +
# "Prob.Death.Bleed.War" +
# "Prob.Death.Bleed.Riv" + 
# "Prob.Ich.Bleed.Dab" +
# "Prob.Ich.Bleed.War"  
# "Prob.Ich.Bleed.Riv"  
# "Prob.Nich.Bleed.Dab" +
# "Prob.Nich.Bleed.War"  
# "Prob.Nich.Bleed.Riv"  
  "Prob.Gos2.Ich" + 
  "Prob.Gos3.Ich" +
  "Prob.Gos4.Ich" +       
  "Prob.Gos5.Ich" +       
  "Util.DepState.Stroke" +
  "Util.IndState.Stroke" +
  "Util.Gos3" +
  "Util.Gos4" + 
  "Util.Gos5" +          
  "Util.Nich" +
  "Cost.Death.Stroke" +
  "Cost.DepState.Inst" +
  "Cost.DepState.Cont"+
  "Cost.IndState.Inst" +
  "Cost.IndState.Cont" +
  "Cost.Gos2.Inst"  +
  "Cost.Gos2.Cont"  +
  "Cost.Gos3.Inst"  +
  "Cost.Gos3.Cont"  +
  "Cost.Gos4.Inst"  +
  "Cost.Gos5.Inst"  +    
 "Cost.Death.Bleed" +
 "Cost.Nich"  +
# "Cost.War" +
 "Sens"   +
 "Spec"   +
 "hrProp.C0" +
 "hrProp.C1" +
 "hrProp.CV0" +
 "hrProp.CV1"          

dab.inputs <- c(
  "Hr.Bleed.lt75.Dab" +
  "Hr.Bleed.75over.Dab" +
  "Hr.Stroke.Dab" +
 "Prob.Death.Bleed.Dab" +
 "Prob.Ich.Bleed.Dab" +
 "Prob.Nich.Bleed.Dab" +

war.inputs <- c(
  "Hr.Bleed.lt75.War" +
  "Hr.Bleed.75over.War" +
  "Hr.Stroke.War" +
 "Prob.Death.Bleed.War" +
 "Prob.Ich.Bleed.War" +
 "Prob.Nich.Bleed.War" +
  "Cost.War" +

  
riv.inputs <- c(
  "Hr.Bleed.lt75.Riv" +
  "Hr.Bleed.75over.Riv" +
  "Hr.Stroke.Riv" +
 "Prob.Death.Bleed.Riv" +
 "Prob.Ich.Bleed.Riv" +
 "Prob.Nich.Bleed.Riv" +

# uncertainty in proportions

 
 ###########################################
 ###########################################
 ##########MAIN INPUTS######################
  Hr.Stroke.LAABN + Hr.Stroke.C0 +Hr.Stroke.C1 + Hr.Stroke.C2 + Hr.Stroke.C3 + Hr.Stroke.C4 + Prob.Death.stroke +
  Prob.DepState.stroke + Prob.IndState.stroke + Prob.Gos2.Ich +Prob.Gos3.Ich + Prob.Gos4.Ich +  Prob.Gos5.Ich + Util.DepState.Stroke +
  Util.IndState.Stroke + Util.Gos3 + Util.Gos4 + Util.Gos5 +  Util.Nich + Cost.Death.Stroke + Cost.DepState.Inst + Cost.DepState.Cont +
  Cost.IndState.Inst + Cost.IndState.Cont + Cost.Gos2.Inst  + Cost.Gos2.Cont  + Cost.Gos3.Inst  + Cost.Gos3.Cont  + Cost.Gos4.Inst  +
  Cost.Gos5.Inst  + Cost.Death.Bleed + Cost.Nich  + Sens  + Spec  + hrProp.C0 + hrProp.C1 + hrProp.CV0 + hrProp.CV1          

### Dab Inputs##############################
 Hr.Bleed.lt75.Dab + Hr.Bleed.75over.Dab + Hr.Stroke.Dab + Prob.Death.Bleed.Dab + Prob.Ich.Bleed.Dab + Prob.Nich.Bleed.Dab
  
### War Inputs##############################
 Hr.Bleed.lt75.War + Hr.Bleed.75over.War + Hr.Stroke.War + Prob.Death.Bleed.War + Prob.Ich.Bleed.War  + Prob.Nich.Bleed.War + Cost.War
  
### Riv Inputs #############################
 Hr.Bleed.lt75.Riv + Hr.Bleed.75over.Riv + Hr.Stroke.Riv + Prob.Death.Bleed.Riv + Prob.Ich.Bleed.Riv + Prob.Nich.Bleed.Riv
  
  
###################################################################
###################################################################

# Calculate relationship between proportion of population in patient groups and true HR proportion 
# assuming sens and spec fixed at expected values
  
trueHrProp <- seq(from=0.00, to= 0.20, by=0.001)
results.matrix <- matrix(nrow=4, ncol=length(trueHrProp))
rownames(results.matrix) <- c("TP", "TN", "FP", "FN")
colnames(results.matrix) <- trueHrProp
  
for (i in 1:length(trueHrProp)){
  tmp <- Calc.Prop(trueHrProp[i], sens=sens.Exp, spec=spec.Exp)
  results.matrix["TP", i] <- tmp$prop.comp["TP"]
  results.matrix["TN", i] <- tmp$prop.comp["TN"]
  results.matrix["FP", i] <- tmp$prop.comp["FP"]
  results.matrix["FN", i] <- tmp$prop.comp["FN"]   
}  

  
  
#DabCV.Exp.SS <- Test.SensSpec(0.5/12, cost=this.cost, qaly=this.qaly)

# Should do this differently: want to know what the effect of different proportions of people 
# with structurally high risk in this population group

trueHRprop <- seq(0,0.2, by = 0.001)
icerVec <- vector("numeric", length(trueHRprop))

for (i in 1:length(trueHRprop)){
  X <- Calc.MeanCostQaly(Data.DabCV)
  this.cost <- X$cost
  this.qaly <- X$qaly
  X <- Calc.Prop(trueHRprop[i], sens=sens.Exp, spec=spec.Exp)
  this.prop.base <- X$prop.base
  this.prop.comp <- X$prop.comp
  DabCV.Exp.icer <- Calc.CostQALY(prop.base=this.prop.base, prop.comp=this.prop.comp, cost=this.cost, qaly=this.qaly)
  icerVec[i] <- DabCV.Exp.icer$icer
}

plot(icerVec ~ trueHRprop, type="l", ylab="ICER (£/QALY)", xlab="True proportion high risk", lwd=2)

abline(v=0, col="grey")

segments(x0=0, x1=2/24, y0=icerVec[max(which(trueHRprop < 2/24))], y1=icerVec[max(which(trueHRprop < 2/24))], lty=1)
segments(2/24, x1=2/24, y0=0, y1=icerVec[max(which(trueHRprop < 2/24))], lty=1)
# implies ICER of £29,033 / QALY

segments(x0=0, x1=0.5/12, y0=icerVec[max(which(trueHRprop < 0.5/12))], y1=icerVec[max(which(trueHRprop < 0.5/12))], lty=2)
segments(x0=0.5/12, x1=0.5/12, y0=0, y1=icerVec[max(which(trueHRprop < 0.5/12))], lty=2)
# implies ICER of £51,227 / QALY

segments(x0=0, x1=trueHRprop[max(which(icerVec > 20000))], y0=20000, y1=20000, lwd=2, lty=1)
segments(x0=trueHRprop[max(which(icerVec > 20000))], x1=trueHRprop[max(which(icerVec > 20000))], y0=0, y1=20000,lwd=2, lty=1)
# implies HR prop of 0.117

segments(x0=0, x1=trueHRprop[max(which(icerVec > 30000))], y0=30000, y1=30000, lwd=2, lty=2)
segments(x0=trueHRprop[max(which(icerVec > 30000))], x1=trueHRprop[max(which(icerVec > 30000))], y0=0, y1=30000, lwd=2, lty=2)
# implies HR prop of 0.08

legend("topright", lwd=c(1,1,NA, 2,2), col=c("black"), lty=c(1,2,NA, 1,2), 
       legend=c("ICER if true proportion is 2/24", 
                "ICER if true proportion is 0.5/12 (noninformative prior)", 
                NA,
                "True proportion required for ICER of £20,000 / QALY", 
                "True proportion required for ICER of £30,000 / QALY"))


########################################################################################################################
#[5]
hrProp <- CV.Data.prior$events[CV.Data.prior$score==0] / CV.Data.prior$total[CV.Data.prior$score==0] # proportion high risk | score
DabCV.PSA <- Run.PSA(Data.DabCV0, hrProp, sens=sens.PSA, spec=spec.PSA)
Draw.PSA(DabCV0.PSA, main="")
#Draw.PSA(DabCV.PSA, main="", contour=T)

#plot(density(DabCV.PSA$icer), main="Density of ICERs from PSA", xlab="ICER (£/QALY)", ylab="Proportion of draws", xlim=c(-500000, 500000), lwd=2)
#abline(v=0)

#quantile(DabCV.PSA$icer, c(0.025, 0.5, 0.975))
#length(which(DabCV.PSA$icer > 20000)) / length(DabCV.PSA$icer)
#length(which(DabCV.PSA$icer > 30000)) / length(DabCV.PSA$icer)



####################################################################
####################################################################




# Should do this differently: want to know what the effect of different proportions of people 
# with structurally high risk in this population group

# Dabigatran
  
  
trueHRprop <- seq(0,0.2, by = 0.001)
icerVec <- vector("numeric", length(trueHRprop))

for (i in 1:length(trueHRprop)){
  X <- Calc.MeanCostQaly(Data.DabC0)
  this.cost <- X$cost
  this.qaly <- X$qaly
  X <- Calc.Prop(trueHRprop[i], sens=sens.Exp, spec=spec.Exp)
  this.prop.base <- X$prop.base
  this.prop.comp <- X$prop.comp
  DabCV.Exp.icer <- Calc.CostQALY(prop.base=this.prop.base, prop.comp=this.prop.comp, cost=this.cost, qaly=this.qaly)
  icerVec[i] <- DabCV.Exp.icer$icer
}

plot(icerVec ~ trueHRprop, type="l", main="Relationship between mean ICER and true high risk proportion\nDabigatran", ylab="ICER (£/QALY)", ylim=c(0, 80000), xlab="True proportion high risk", lwd=2)

abline(v=0, col="grey")
abline(h=0, col="grey")

trueHRprop <- seq(0,0.2, by = 0.001)
icerVec <- vector("numeric", length(trueHRprop))

for (i in 1:length(trueHRprop)){
  X <- Calc.MeanCostQaly(Data.DabCV0)
  this.cost <- X$cost
  this.qaly <- X$qaly
  X <- Calc.Prop(trueHRprop[i], sens=sens.Exp, spec=spec.Exp)
  this.prop.base <- X$prop.base
  this.prop.comp <- X$prop.comp
  DabCV.Exp.icer <- Calc.CostQALY(prop.base=this.prop.base, prop.comp=this.prop.comp, cost=this.cost, qaly=this.qaly)
  icerVec[i] <- DabCV.Exp.icer$icer
}

lines(icerVec ~ trueHRprop, type="l", lwd=2, lty="dashed")


  legend("topright", lty=c("solid", "dashed"), bty="n", inset=0.01, col=c("black", "black", "white", "white"), lwd=c(2,2), legend=c(expression(CHADS[2]),
                                                                expression(paste(CHA[2], DS[2],-VASc), sep="")) )

###################################################################################################################
####################################################################################################################
# Dabigatran
  
  
trueHRprop <- seq(0,0.2, by = 0.001)
icerVec <- vector("numeric", length(trueHRprop))

for (i in 1:length(trueHRprop)){
  X <- Calc.MeanCostQaly(Data.DabC0)
  this.cost <- X$cost
  this.qaly <- X$qaly
  X <- Calc.Prop(trueHRprop[i], sens=sens.Exp, spec=spec.Exp)
  this.prop.base <- X$prop.base
  this.prop.comp <- X$prop.comp
  DabCV.Exp.icer <- Calc.CostQALY(prop.base=this.prop.base, prop.comp=this.prop.comp, cost=this.cost, qaly=this.qaly)
  icerVec[i] <- DabCV.Exp.icer$icer
}

plot(icerVec ~ trueHRprop, type="l", main="Relationship between mean ICER and true high risk proportion\nDabigatran", ylab="ICER (£/QALY)", ylim=c(0, 80000), xlab="True proportion high risk", lwd=2)

abline(v=0, col="grey")
abline(h=0, col="grey")

trueHRprop <- seq(0,0.2, by = 0.001)
icerVec <- vector("numeric", length(trueHRprop))

for (i in 1:length(trueHRprop)){
  X <- Calc.MeanCostQaly(Data.DabCV0)
  this.cost <- X$cost
  this.qaly <- X$qaly
  X <- Calc.Prop(trueHRprop[i], sens=sens.Exp, spec=spec.Exp)
  this.prop.base <- X$prop.base
  this.prop.comp <- X$prop.comp
  DabCV.Exp.icer <- Calc.CostQALY(prop.base=this.prop.base, prop.comp=this.prop.comp, cost=this.cost, qaly=this.qaly)
  icerVec[i] <- DabCV.Exp.icer$icer
}

lines(icerVec ~ trueHRprop, type="l", lwd=2, lty="dashed")


  legend("topright", lty=c("solid", "dashed"), bty="n", inset=0.01, col=c("black", "black", "white", "white"), lwd=c(2,2), legend=c(expression(CHADS[2]),
                                                                expression(paste(CHA[2], DS[2],-VASc), sep="")) )

###################################################################################################################
####################################################################################################################
  
  
  
segments(x0=0, x1=2/24, y0=icerVec[max(which(trueHRprop < 2/24))], y1=icerVec[max(which(trueHRprop < 2/24))], lty=1)
segments(2/24, x1=2/24, y0=0, y1=icerVec[max(which(trueHRprop < 2/24))], lty=1)
# implies ICER of £29,033 / QALY

segments(x0=0, x1=0.5/12, y0=icerVec[max(which(trueHRprop < 0.5/12))], y1=icerVec[max(which(trueHRprop < 0.5/12))], lty=2)
segments(x0=0.5/12, x1=0.5/12, y0=0, y1=icerVec[max(which(trueHRprop < 0.5/12))], lty=2)
# implies ICER of £51,227 / QALY

segments(x0=0, x1=trueHRprop[max(which(icerVec > 20000))], y0=20000, y1=20000, lwd=2, lty=1)
segments(x0=trueHRprop[max(which(icerVec > 20000))], x1=trueHRprop[max(which(icerVec > 20000))], y0=0, y1=20000,lwd=2, lty=1)
# implies HR prop of 0.117

segments(x0=0, x1=trueHRprop[max(which(icerVec > 30000))], y0=30000, y1=30000, lwd=2, lty=2)
segments(x0=trueHRprop[max(which(icerVec > 30000))], x1=trueHRprop[max(which(icerVec > 30000))], y0=0, y1=30000, lwd=2, lty=2)
# implies HR prop of 0.08

legend("topright", lwd=c(1,1,NA, 2,2), col=c("black"), lty=c(1,2,NA, 1,2), 
       legend=c("ICER if true proportion is 2/24", 
                "ICER if true proportion is 0.5/12 (noninformative prior)", 
                NA,
                "True proportion required for ICER of £20,000 / QALY", 
                "True proportion required for ICER of £30,000 / QALY"))

########################################################################################################################

