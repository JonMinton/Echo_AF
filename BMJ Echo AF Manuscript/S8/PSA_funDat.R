######################################################################################################################################
################################ F U N C T I O N S ###################################################################################
######################################################################################################################################

# Calculate proportions of TN, TP, FN and FP in both base and comparator
Calc.Prop <- function(hrProp, sens=1, spec=1){
  prop.base.TP <- 0
  prop.base.TN <- 1 - hrProp
  prop.base.FP <- 0
  prop.base.FN <-     hrProp

  prop.comp.TP <-      hrProp  *      sens
  prop.comp.TN <- (1 - hrProp) *      spec
  prop.comp.FP <- (1 - hrProp) * (1 - spec)
  prop.comp.FN <-      hrProp  * (1 - sens)

  prop.base <- c(TP=prop.base.TP, TN=prop.base.TN, FP=prop.base.FP, FN=prop.base.FN)
  prop.comp <- c(TP=prop.comp.TP, TN=prop.comp.TN, FP=prop.comp.FP, FN=prop.comp.FN)
  
  return(list(prop.base=prop.base, prop.comp=prop.comp))
}

# Expected value of perfect information
Calc.EVPI <- function(D, hrProp, Nruns=1000, sens, spec, maicers=seq(0,50000, by=100)){
  
  Niter <- length(maicers)
  EVPI.out <- vector("list", Niter)
  base.NetMB <- rep(NA,      Niter)
  comp.NetMB <- rep(NA,      Niter)
  
  # for each element in maicers
    # for each element in PSA
       # calculate cost of baseline
       # calculate QALY of baseline
       # calculate cost of comparator
       # calculate QALY of comparator
  
       # calculate netMB of baseline at maicer[i]
       # calculate netMB of comparator at maicer[i]
  
  tmp <- Run.PSA.evpi(D=D, hrProp=hrProp, Nruns=Nruns, sens=sens, spec=spec)

  for (i in 1:Niter){
    base.NetMB <- Calc.NetMB(cost=tmp$base.cost, qaly=tmp$base.qaly, maicer=maicers[i])
    comp.NetMB <- Calc.NetMB(cost=tmp$comp.cost, qaly=tmp$comp.qaly, maicer=maicers[i])
    
    base.Exp.NetMB <- mean(base.NetMB)
    comp.Exp.NetMB <- mean(comp.NetMB)
    
    Exp.optimalChoice <- ifelse(comp.Exp.NetMB > base.Exp.NetMB, "Comparator", "Baseline")
    
    tmp2 <- cbind(base.NetMB, comp.NetMB)
    Max.NetMB <- apply(tmp2,1, function(x) max(x))
    Dif.NetMB <- apply(tmp2,1, function(x) abs(x[1] - x[2]))
    Exp.Max.NetMB <- mean(Max.NetMB)
    
    Ind.optimalChoice <- ifelse(comp.NetMB > base.NetMB, "Comparator", "Baseline")
    
    Ind.opLoss <- ifelse(Ind.optimalChoice==Exp.optimalChoice, 0, Dif.NetMB)
    Sum.opLoss <- sum(Ind.opLoss)
    Exp.opLoss <- Sum.opLoss/ Nruns
    
    EVPI.out[[i]] <- list(maicer=maicers[i], base.NetMB=base.NetMB, comp.NetMB=comp.NetMB, 
                          ind.opLoss=Ind.opLoss, sum.opLoss=Sum.opLoss, exp.opLoss=Exp.opLoss,
                          ind.optimalChoice=Ind.optimalChoice, exp.optimalChoice=Exp.optimalChoice)
  }
  return(EVPI.out)
}

#######
# Attempt at building EVPPI function:

# Expected value of perfect parameter information
Calc.EVPPI <- function(InData, mod, params, joint.SensSpec=F, hrProp.option,  Nruns=1000, maicers=seq(0,50000, by=1000)){
  
  Niter <- length(maicers)
  Nparam <- ifelse(joint.SensSpec==F, length(params), 1)
  EVPPI <- vector("list", Nparam)
  names(EVPPI) <- ifelse(joint.SensSpec==F, params, c("Sens_and_Spec"))
  base.NetMB <- rep(NA, Niter)
  comp.NetMB <- rep(NA, Niter)
  
  # for each element in maicers
    # for each element in PSA
       # calculate cost of baseline
       # calculate QALY of baseline
       # calculate cost of comparator
       # calculate QALY of comparator
  
       # calculate netMB of baseline at maicer[i]
       # calculate netMB of comparator at maicer[i]
  print(paste("Niter=", Niter))
  print(paste("Nparam=", Nparam))
#  browser()
  

  for (j in 1:Nparam){
  print(paste("Entered loop j=", j))
#    browser()
    EVPPI.inner <- vector("list", Niter)
    for (i in 1:Niter){
      if (i%%10==0) {print(paste("Entered inner loop i = ", i))}
      tmpData <- InData
      if (joint.SensSpec==F){tmpData[,params[j]] <- InData[i,params[j]]} else {tmpData[,c("Sens", "Spec")] <- InData[i,c("Sens", "Spec")]}
      PredOutData <- as.data.frame(predict(mod, tmpData))
      if (hrProp.option=="C0")  {hrProp <- tmpData$hrProp.C0}
      if (hrProp.option=="C1")  {hrProp <- tmpData$hrProp.C1}
      if (hrProp.option=="CV0") {hrProp <- tmpData$hrProp.CV0}
      if (hrProp.option=="CV0") {hrProp <- tmpData$hrProp.CV1}

      tmp <- Run.PSA.evpi(D=PredOutData, hrProp=hrProp, Nruns=Nruns, sens=tmpData$Sens, spec=tmpData$Spec)
      
      base.NetMB <- Calc.NetMB(cost=tmp$base.cost, qaly=tmp$base.qaly, maicer=maicers[i])
      comp.NetMB <- Calc.NetMB(cost=tmp$comp.cost, qaly=tmp$comp.qaly, maicer=maicers[i])
      
      base.Exp.NetMB <- mean(base.NetMB)
      comp.Exp.NetMB <- mean(comp.NetMB)
      
      Exp.optimalChoice <- ifelse(comp.Exp.NetMB > base.Exp.NetMB, "Comparator", "Baseline")
      
      tmp2 <- cbind(base.NetMB, comp.NetMB)
      Max.NetMB <- apply(tmp2,1, function(x) max(x)             )
      Dif.NetMB <- apply(tmp2,1, function(x) abs(x[1] - x[2])   )
      Exp.Max.NetMB <- mean(Max.NetMB)
      
      Ind.optimalChoice <- ifelse(comp.NetMB > base.NetMB, "Comparator", "Baseline")
      
      Ind.opLoss <- ifelse(Ind.optimalChoice==Exp.optimalChoice, 0, Dif.NetMB)
      Sum.opLoss <- sum(Ind.opLoss)
      Exp.opLoss <- Sum.opLoss/ Nruns
      
      EVPPI.inner[[i]] <- list(
                            maicer=maicers[i], 
                            base.NetMB=base.NetMB, 
                            comp.NetMB=comp.NetMB, 
                            ind.opLoss=Ind.opLoss, 
                            sum.opLoss=Sum.opLoss, 
                            exp.opLoss=Exp.opLoss,
                            ind.optimalChoice=Ind.optimalChoice, 
                            exp.optimalChoice=Exp.optimalChoice
                          )      
    }
    EVPPI[[j]] <- EVPPI.inner  
  print(paste("Exited loop j=", j))
  }
  return(EVPPI)
}

#######################################################
# Produce PSA: cost & qaly for both baseline and comparator
Run.PSA.evpi <- function(D, hrProp, Nruns=1000, sens=1, spec=1){

 
  PSA.results <- data.frame(
    base.cost=rep(NA, Nruns), 
    base.qaly=rep(NA, Nruns), 
    comp.cost=rep(NA, Nruns),
    comp.qaly=rep(NA, Nruns)
    )
  
  for (i in 1:Nruns){

    tmp <- Calc.Prop(hrProp[i], sens=sens[i], spec=spec[i])
    
    groups <- c("TP", "TN", "FP", "FN")
    prop.base <- tmp$prop.base
    prop.comp <- tmp$prop.comp
    
    cost <- c(D$TP.c[i], D$TN.c[i], D$FP.c[i], D$FN.c[i])
    qaly <- c(D$TP.q[i], D$TN.q[i], D$FP.q[i], D$FN.q[i])
    
    tmp2 <- Calc.CostQALY.evpi(prop.base, prop.comp, cost, qaly)
    
    PSA.results[i, "base.cost"] <- tmp2$base$cost
    PSA.results[i, "base.qaly"] <- tmp2$base$qaly
    
    PSA.results[i, "comp.cost"] <- tmp2$comp$cost
    PSA.results[i, "comp.qaly"] <- tmp2$comp$qaly
  }

  return(PSA.results)
}
########################################################
Calc.CostQALY.evpi <- function(prop.base, prop.comp, cost, qaly, cost.test=66){
  
  c.base <- sum(prop.base * cost)
  c.comp <- cost.test + sum(prop.comp * cost)
  q.base <- sum(prop.base * qaly)
  q.comp <- sum(prop.comp * qaly)
  
  
  return(list(base=list(cost=c.base, qaly=q.base), 
              comp=list(cost=c.comp, qaly=q.comp)))
}


#ICER sort
# Now no longer using this function following an argument with Matt

# Sort.ICER <- function(D, qnt=c(0.025, 0.5, 0.975)){
#   icers <- D$icer
#   quadrant <- D$quadrant
#   icers[quadrant=="NW"] <- Inf
#   icers[quadrant=="SE"] <- -Inf  
#   if (length(which(quadrant=="SW"))!=0){
#     # do some code for when this occurs
#     # The approach below probably isn't correct'
#     mean.sw.icer <- mean(icers[quadrant=="SW"])
#     mean.ne.icer <- mean(icers[quadrant=="NE"])
#     weight <- length(which(quadrant=="SW")) / length(which(quadrant=="NE"))
#     diff <- mean.ne.icer * (1/weight) - mean.sw.icer * weight
#     icers[quadrant=="NE"] <- icers[quadrant=="NE"] - diff
#   }
#   
#   qi <- quantile(icers, qnt, na.rm=T)
#   return(list(icers=icers, qi=qi))
# }

Calc.MeanCostQaly <- function(D){
  tmp.c <- c(tp=mean(D$TP.c), tn=mean(D$TN.c), fp=mean(D$FP.c), fn=mean(D$FN.c))
  tmp.q <- c(tp=mean(D$TP.q), tn=mean(D$TN.q), fp=mean(D$FP.q), fn=mean(D$FN.q))
  qi <- list(cost=tmp.c, qaly=tmp.q)
  return(qi)
}

# Calculate net benefit
Calc.NetMB <- function(cost, qaly, maicer=20000){
  netBenefit <- maicer * qaly - cost
  return(netBenefit)
}

Draw.EVPI <- function(D, main="", max.y=NA, export=NULL, ylab="Individual EVPI (£)", xlab="Willingness-to-pay threshold (£/QALY)"){
  
  maicers <- sapply(D, function(x) x$maicer)
  opLoss <- sapply(D, function(x) x$exp.opLoss)
  
  Data.out <- data.frame(maicer=maicers, opLoss=opLoss)

  if (!is.null(export)){jpeg(filename=export, height=800, width=800)}

  if (is.na(max.y)){
      plot(opLoss ~ maicer, data=Data.out, type="l", ylab=ylab, xlab=xlab, main=main)
  } else {
    plot(opLoss ~ maicer, data=Data.out, type="l", ylab=ylab, xlab=xlab, main=main, ylim=c(0, max.y))  
  }
  if(!is.null(export)){dev.off()}
  
  return(Data.out)
}




# CEAF
# NOTE : THIS FUNCTION HAS NOW BEEN SUPERSEDED BY Draw.CEAF.NB as it is more reliable. 
#   To do at some point in the future:
#      Extend to more than two options.
#Draw.CEAF <- function(icers, D, main="", export=NULL, xlab="Threshold (£/QALY)", ylab="Probability cost effective", ylim=c(0,1)){
#  maicers <-  sapply(D, function(x) x$maicer)
#  optChoices <- sapply(D, function(x) x$exp.optimalChoice)
#  
#  NoptChoices <- length(unique(optChoices))
#  pce1 <- rep(NA, length(maicers))
#  pce2 <- rep(NA, length(maicers))
#
#  pce <- matrix(nrow=length(maicers), ncol=NoptChoices)
#  colnames(pce) <- unique(optChoices)
#  
#  for (i in 1:length(maicers)){
#    ifelse(optChoices[i]=="Comparator", #
#
#           pce1[i] <- length(which(icers < maicers[i]))/length(icers), 
#           pce2[i] <- 1 - length(which(icers < maicers[i]))/length(icers)
#    )
#  }
#  
#  Data <- data.frame(maicer=maicers, pce1=pce1, pce2=pce2)
#  if (!is.null(export)){jpeg(filename=export, width=800, height=800)}
#  plot(pce1 ~ maicer, data=Data, type="l", lwd=2, main=main, xlab=xlab, ylab=ylab, ylim=ylim)
#  lines(pce2 ~ maicer, data=Data, lwd=2, lty="dashed")
#  if (!is.null(export)){ dev.off()}
#  return(Data)  
#}

###################################  
Draw.CEAF.NB <- function(PSA, D, main="", export=NULL, xlab="Willingness-to-pay threshold (£/QALY)", 
                         ylab="Probability cost effective", legend=NULL, legpos="topright", ylim=c(0,1)){
  maicers <-  sapply(D, function(x) x$maicer)
  optChoices <- sapply(D, function(x) x$exp.optimalChoice)
  
  #  NoptChoices <- length(unique(optChoices))
  pce1 <- rep(NA, length(maicers))
  pce2 <- rep(NA, length(maicers))
  
  #  pce <- matrix(nrow=length(maicers), ncol=NoptChoices)
  #  colnames(pce) <- unique(optChoices)
  qaly <- PSA$qaly
  cost <- PSA$cost
  
  for (i in 1:length(maicers)){
    ifelse(optChoices[i]=="Comparator", 
           
           pce1[i] <- length(which(qaly * maicers[i] - cost > 0)) / length(qaly * maicers[i] - cost), 
           pce2[i] <- 1 - length(which(qaly * maicers[i] - cost > 0)) / length(qaly * maicers[i] - cost)
           )
  }
  
  Data <- data.frame(maicer=maicers, pce1=pce1, pce2=pce2)
  if (!is.null(export)){jpeg(filename=export, width=800, height=800)}
  plot(pce1 ~ maicer, data=Data, type="l", lwd=2, main=main, xlab=xlab, ylab=ylab, ylim=ylim)
  lines(pce2 ~ maicer, data=Data, lwd=2, lty="dashed")
  if (!is.null(legend)){
    legend(legpos, legend=legend, lty=c("solid", "dashed"), lwd=2)
  }
  if (!is.null(export)){ dev.off()}
  return(Data)  
}  
 
# Calculate cost and qalys (assuming cost of test is £66)
Calc.CostQALY <- function(prop.base, prop.comp, cost, qaly, cost.test=66, maicer1=20000, maicer2=30000){

  c.base <- sum(prop.base * cost)
  c.comp <- cost.test + sum(prop.comp * cost)
  q.base <- sum(prop.base * qaly)
  q.comp <- sum(prop.comp * qaly)

  dif.cost <- c.comp - c.base
  dif.qaly <- q.comp - q.base
  
  icer <- dif.cost / dif.qaly  
  
  if (dif.cost < 0 & dif.qaly < 0) quadrant <- "SW"
  if (dif.cost >= 0 & dif.qaly < 0) quadrant <- "NW"
  if (dif.cost >= 0 & dif.qaly >= 0) quadrant <- "NE"
  if (dif.cost < 0 & dif.qaly >= 0) quadrant <- "SE"
  netMB1 <- Calc.NetMB(cost=dif.cost, qaly=dif.qaly, maicer=maicer1)
  netMB2 <- Calc.NetMB(cost=dif.cost, qaly=dif.qaly, maicer=maicer2)
  
  return(list(base=list(cost=c.base, qaly=q.base), 
              comp=list(cost=c.comp, qaly=q.comp), 
              dif.cost=dif.cost, 
              dif.qaly=dif.qaly, 
              icer=icer, 
              quadrant=quadrant,
              netMB1=netMB1, 
              netMB2=netMB2)
         )
}

Calc.OtherQI <- function(prop.base, prop.comp, D){

  lyg             <- c(D$TP.lyg,            D$TN.lyg,             D$FP.lyg,            D$FN.lyg)
  death.other     <- c(D$TP.death.other,    D$TN.death.other,     D$FP.death.other,    D$FN.death.other)
  death.stroke    <- c(D$TP.death.stroke,   D$TN.death.stroke,    D$FP.death.stroke,   D$FN.death.stroke)
  death.bleed     <- c(D$TP.death.bleed,    D$TN.death.bleed,     D$FP.death.bleed,    D$FN.death.bleed)
  stroke.dep      <- c(D$TP.stroke.dep,     D$TN.stroke.dep,      D$FP.stroke.dep,     D$FN.stroke.dep)
  stroke.indep    <- c(D$TP.stroke.indep,   D$TN.stroke.indep,    D$FP.stroke.indep,   D$FN.stroke.indep)
  bleed.ich       <- c(D$TP.bleed.ich,      D$TN.bleed.ich,       D$FP.bleed.ich,      D$FN.bleed.ich)
  bleed.nich      <- c(D$TP.bleed.nich,     D$TN.bleed.nich,      D$FP.bleed.nich,     D$FN.bleed.nich)
  count.p.dep     <- c(D$TP.count.p.dep,    D$TN.count.p.dep,     D$FP.count.p.dep,    D$FN.count.p.dep)
  count.p.indep   <- c(D$TP.count.p.indep,  D$TN.count.p.indep,   D$FP.count.p.indep,  D$FN.count.p.indep)
  count.p.ich     <- c(D$TP.count.p.ich,    D$TN.count.p.ich,     D$FP.count.p.ich,    D$FN.count.p.ich)
  count.p.nich    <- c(D$TP.count.p.nich,   D$TN.count.p.nich,    D$FP.count.p.nich,   D$FN.count.p.nich)
  
  base.lyg            <- sum(prop.base * lyg)
  base.death.other    <- sum(prop.base * death.other)
  base.death.stroke   <- sum(prop.base * death.stroke)
  base.death.bleed    <- sum(prop.base * death.bleed)
  base.stroke.dep     <- sum(prop.base * stroke.dep)
  base.stroke.indep   <- sum(prop.base * stroke.indep)
  base.bleed.ich      <- sum(prop.base * bleed.ich)
  base.bleed.nich     <- sum(prop.base * bleed.nich)
  base.count.p.dep    <- sum(prop.base * count.p.dep)
  base.count.p.indep  <- sum(prop.base * count.p.indep)
  base.count.p.ich    <- sum(prop.base * count.p.ich)
  base.count.p.nich   <- sum(prop.base * count.p.nich)

  comp.lyg            <- sum(prop.comp * lyg)
  comp.death.other    <- sum(prop.comp * death.other) 
  comp.death.stroke   <- sum(prop.comp * death.stroke)  
  comp.death.bleed    <- sum(prop.comp * death.bleed)   
  comp.stroke.dep     <- sum(prop.comp * stroke.dep)    
  comp.stroke.indep   <- sum(prop.comp * stroke.indep)  
  comp.bleed.ich      <- sum(prop.comp * bleed.ich)     
  comp.bleed.nich     <- sum(prop.comp * bleed.nich)    
  comp.count.p.dep    <- sum(prop.comp * count.p.dep)   
  comp.count.p.indep  <- sum(prop.comp * count.p.indep) 
  comp.count.p.ich    <- sum(prop.comp * count.p.ich)   
  comp.count.p.nich   <- sum(prop.comp * count.p.nich)  
  
  dif.lyg            <- comp.lyg           -  base.lyg
  dif.death.other    <- comp.death.other   -  base.death.other
  dif.death.stroke   <- comp.death.stroke  -  base.death.stroke
  dif.death.bleed    <- comp.death.bleed   -  base.death.bleed 
  dif.stroke.dep     <- comp.stroke.dep    -  base.stroke.dep
  dif.stroke.indep   <- comp.stroke.indep  -  base.stroke.indep   
  dif.bleed.ich      <- comp.bleed.ich     -  base.bleed.ich      
  dif.bleed.nich     <- comp.bleed.nich    -  base.bleed.nich   
  dif.count.p.dep    <- comp.count.p.dep   -  base.count.p.dep    
  dif.count.p.indep  <- comp.count.p.indep -  base.count.p.indep 
  dif.count.p.ich    <- comp.count.p.ich   -  base.count.p.ich  
  dif.count.p.nich   <- comp.count.p.nich  -  base.count.p.nich 
  
  Base <- list(
    lyg           = base.lyg,
    death.other   = base.death.other,
    death.stroke  = base.death.stroke, 
    death.bleed   = base.death.bleed,
    stroke.dep    = base.stroke.dep,
    stroke.indep  = base.stroke.indep,
    bleed.ich     = base.bleed.ich,
    bleed.nich    = base.bleed.nich,
    count.p.dep   = base.count.p.dep,
    count.p.indep = base.count.p.indep,
    count.p.ich   = base.count.p.ich,
    count.p.nich  = base.count.p.nich
    )

  Comp <- list(
    lyg           = comp.lyg,
    death.other   = comp.death.other,
    death.stroke  = comp.death.stroke, 
    death.bleed   = comp.death.bleed,
    stroke.dep    = comp.stroke.dep,
    stroke.indep  = comp.stroke.indep,
    bleed.ich     = comp.bleed.ich,
    bleed.nich    = comp.bleed.nich,
    count.p.dep   = comp.count.p.dep,
    count.p.indep = comp.count.p.indep,
    count.p.ich   = comp.count.p.ich,
    count.p.nich  = comp.count.p.nich
    )

  Dif <- list(
    lyg           =  dif.lyg,
    death.other   =  dif.death.other,
    death.stroke  =  dif.death.stroke, 
    death.bleed   =  dif.death.bleed,
    stroke.dep    =  dif.stroke.dep,
    stroke.indep  =  dif.stroke.indep,
    bleed.ich     =  dif.bleed.ich,
    bleed.nich    =  dif.bleed.nich,
    count.p.dep   =  dif.count.p.dep,
    count.p.indep =  dif.count.p.indep,
    count.p.ich   =  dif.count.p.ich,
    count.p.nich  =  dif.count.p.nich
    )
  
  return(list(base=Base, comp=Comp, dif=Dif))
}



# test the sensitivity of the results to assumptions about the proportion of LAABNs 
# in the population under consideration

Test.SensSpec <- function(hrProp, sens.range=seq(0,1, by=0.1), spec.range=seq(0,1, by=0.1), cost, qaly){
  
  icer.block <- matrix(NA, nrow=length(sens.range), ncol=length(spec.range))
  rownames(icer.block) <- sens.range
  colnames(icer.block) <- spec.range
  
  for (i in 1:length(sens.range)){
    for (j in 1:length(spec.range)){
      sens <- sens.range[i]
      spec <- spec.range[j]
      
      tmp <- Calc.Prop(hrProp, sens=sens, spec=spec)
      
      groups <- c("TP", "TN", "FP", "FN")
      prop.base <- tmp$prop.base
      prop.comp <- tmp$prop.comp
  
      
      tmp2 <- Calc.CostQALY(prop.base, prop.comp, cost, qaly)
  
      icer.block[i,j] <- tmp2$icer
    }
  }
  return(icer.block)
}

# Produce PSA
Run.PSA <- function(D, hrProp, Nruns=1000, sens=1, spec=1, maicer1=20000, maicer2=30000){

  PSA.results <- data.frame(
    base.cost=rep(NA, Nruns),
    base.qaly=rep(NA, Nruns),
    comp.cost=rep(NA, Nruns),
    comp.qaly=rep(NA, Nruns),
    cost=rep(NA, Nruns), 
    qaly=rep(NA, Nruns), 
    icer=rep(NA, Nruns), 
    quadrant=rep(NA,Nruns), 
    netMB1=rep(NA, Nruns),
    netMB2=rep(NA, Nruns)
    )
 
  for (i in 1:Nruns){
    
    tmp <- Calc.Prop(hrProp[i], sens=sens[i], spec=spec[i])
    
    groups <- c("TP", "TN", "FP", "FN")
    prop.base <- tmp$prop.base
    prop.comp <- tmp$prop.comp
    
    cost <- c(D$TP.c[i], D$TN.c[i], D$FP.c[i], D$FN.c[i])
    qaly <- c(D$TP.q[i], D$TN.q[i], D$FP.q[i], D$FN.q[i])
    
    tmp2 <- Calc.CostQALY(prop.base, prop.comp, cost, qaly, maicer1=maicer1, maicer2=maicer2)
  

    PSA.results[i, "base.cost"] <- tmp2$base$cost
    PSA.results[i, "base.qaly"] <- tmp2$base$qaly
    PSA.results[i, "comp.cost"] <- tmp2$comp$cost
    PSA.results[i, "comp.qaly"] <- tmp2$comp$qaly
    PSA.results[i,"cost"]       <- tmp2$dif.cost
    PSA.results[i, "qaly"]      <- tmp2$dif.qaly
    PSA.results[i, "icer"]      <- tmp2$icer
    PSA.results[i, "quadrant"]  <- tmp2$quadrant
    PSA.results[i, "netMB1"]    <- tmp2$netMB1
    PSA.results[i, "netMB2"]   <- tmp2$netMB2
    
   }
  return(PSA.results)
}

# Produce PSA
Run.PSA.otherQI <- function(D, hrProp, Nruns=1000, sens=1, spec=1){
  

  PSA.results <- data.frame(
    base.lyg           = rep(NA, Nruns), 
    base.death.other   = rep(NA, Nruns), 
    base.death.stroke  = rep(NA, Nruns), 
    base.death.bleed   = rep(NA, Nruns), 
    base.stroke.dep    = rep(NA, Nruns),
    base.stroke.indep  = rep(NA, Nruns), 
    base.bleed.ich     = rep(NA, Nruns),
    base.bleed.nich    = rep(NA, Nruns),
    base.count.p.dep   = rep(NA, Nruns),
    base.count.p.indep = rep(NA, Nruns),
    base.count.p.ich   = rep(NA, Nruns),
    base.count.p.nich  = rep(NA, Nruns),
    
    comp.lyg           = rep(NA, Nruns), 
    comp.death.other   = rep(NA, Nruns), 
    comp.death.stroke  = rep(NA, Nruns), 
    comp.death.bleed   = rep(NA, Nruns), 
    comp.stroke.dep    = rep(NA, Nruns),
    comp.stroke.indep  = rep(NA, Nruns), 
    comp.bleed.ich     = rep(NA, Nruns),
    comp.bleed.nich    = rep(NA, Nruns),
    comp.count.p.dep   = rep(NA, Nruns),
    comp.count.p.indep = rep(NA, Nruns),
    comp.count.p.ich   = rep(NA, Nruns),
    comp.count.p.nich  = rep(NA, Nruns),
    
     dif.lyg           = rep(NA, Nruns), 
     dif.death.other   = rep(NA, Nruns), 
     dif.death.stroke  = rep(NA, Nruns), 
     dif.death.bleed   = rep(NA, Nruns), 
     dif.stroke.dep    = rep(NA, Nruns),
     dif.stroke.indep  = rep(NA, Nruns), 
     dif.bleed.ich     = rep(NA, Nruns),
     dif.bleed.nich    = rep(NA, Nruns),
     dif.count.p.dep   = rep(NA, Nruns),
     dif.count.p.indep = rep(NA, Nruns),
     dif.count.p.ich   = rep(NA, Nruns),
     dif.count.p.nich  = rep(NA, Nruns)
    )
    
 
  for (i in 1:Nruns){
    
    tmp <- Calc.Prop(hrProp[i], sens=sens[i], spec=spec[i])
    
    groups <- c("TP", "TN", "FP", "FN")
    prop.base <- tmp$prop.base
    prop.comp <- tmp$prop.comp
    
    
    tmp <- Calc.OtherQI(prop.base, prop.comp, D=D[i,])
  
    PSA.results[i,  "base.lyg"           ] <- tmp$base$lyg 
    PSA.results[i,  "base.death.other"   ] <- tmp$base$death.other 
    PSA.results[i,  "base.death.stroke"  ] <- tmp$base$death.stroke 
    PSA.results[i,  "base.death.bleed"   ] <- tmp$base$death.bleed 
    PSA.results[i,  "base.stroke.dep"    ] <- tmp$base$stroke.dep
    PSA.results[i,  "base.stroke.indep"  ] <- tmp$base$stroke.indep 
    PSA.results[i,  "base.bleed.ich"     ] <- tmp$base$bleed.ich
    PSA.results[i,  "base.bleed.nich"    ] <- tmp$base$bleed.nich
    PSA.results[i,  "base.count.p.dep"   ] <- tmp$base$count.p.dep
    PSA.results[i,  "base.count.p.indep" ] <- tmp$base$count.p.indep
    PSA.results[i,  "base.count.p.ich"   ] <- tmp$base$count.p.ich
    PSA.results[i,  "base.count.p.nich"  ] <- tmp$base$count.p.nich
    PSA.results[i,  "comp.lyg"           ] <- tmp$comp$lyg 
    PSA.results[i,  "comp.death.other"   ] <- tmp$comp$death.other 
    PSA.results[i,  "comp.death.stroke"  ] <- tmp$comp$death.stroke 
    PSA.results[i,  "comp.death.bleed"   ] <- tmp$comp$death.bleed 
    PSA.results[i,  "comp.stroke.dep"    ] <- tmp$comp$stroke.dep
    PSA.results[i,  "comp.stroke.indep"  ] <- tmp$comp$stroke.indep 
    PSA.results[i,  "comp.bleed.ich"     ] <- tmp$comp$bleed.ich
    PSA.results[i,  "comp.bleed.nich"    ] <- tmp$comp$bleed.nich
    PSA.results[i,  "comp.count.p.dep"   ] <- tmp$comp$count.p.dep
    PSA.results[i,  "comp.count.p.indep" ] <- tmp$comp$count.p.indep
    PSA.results[i,  "comp.count.p.ich"   ] <- tmp$comp$count.p.ich
    PSA.results[i,  "comp.count.p.nich"  ] <- tmp$comp$count.p.nich
    PSA.results[i,  "dif.lyg"           ]  <- tmp$dif$lyg 
    PSA.results[i,  "dif.death.other"   ]  <- tmp$dif$death.other 
    PSA.results[i,  "dif.death.stroke"  ]  <- tmp$dif$death.stroke 
    PSA.results[i,  "dif.death.bleed"   ]  <- tmp$dif$death.bleed 
    PSA.results[i,  "dif.stroke.dep"    ]  <- tmp$dif$stroke.dep
    PSA.results[i,  "dif.stroke.indep"  ]  <- tmp$dif$stroke.indep 
    PSA.results[i,  "dif.bleed.ich"     ]  <- tmp$dif$bleed.ich
    PSA.results[i,  "dif.bleed.nich"    ]  <- tmp$dif$bleed.nich
    PSA.results[i,  "dif.count.p.dep"   ]  <- tmp$dif$count.p.dep
    PSA.results[i,  "dif.count.p.indep" ]  <- tmp$dif$count.p.indep
    PSA.results[i,  "dif.count.p.ich"   ]  <- tmp$dif$count.p.ich
    PSA.results[i,  "dif.count.p.nich"  ]  <- tmp$dif$count.p.nich
    
   }
  return(PSA.results)
}

# Draw either scatterplots or contour plots for PSA
Draw.PSA <- function(Data, export=NULL, cost.swing=NULL, qaly.swing=NULL, contour=F, main="PSA Results"){
  if(is.null(cost.swing)){ cost.swing <- max(abs(Data$cost)) * 1.2 }
  if(is.null(qaly.swing)){ qaly.swing <- max(abs(Data$qaly)) * 1.2 }

  if(!is.null(export)){jpeg(filename=export, height=800, width=800)}
  if(contour==F){
    plot(cost ~ qaly, data=Data, ylim=c(-cost.swing, cost.swing), xlim=c(-qaly.swing,qaly.swing), 
         ylab="Difference in cost (£)", xlab="Difference in QALYs", main=main)
    abline(v=0)
    abline(h=0)
    legend("topleft", legend="More costly and less effective", bty="n")
    legend("topright", legend="More costly and more effective", bty="n")
    legend("bottomleft", legend="Less costly and less effective", bty="n")
    legend("bottomright", legend="Less costly and more effective", bty="n")
  }
  
  if(contour==T){
    require(MASS)
    contour(kde2d(Data$qaly, Data$cost), ylim=c(-cost.swing, cost.swing), xlim=c(-qaly.swing, qaly.swing), 
            ylab="Difference in cost", xlab="Difference in QALYs", main=main)
    abline(v=0)
    abline(h=0)
  }
  
  if(!is.null(export)) {dev.off()}
  
}

#Draw CEAC

# Draw.CEAC <- function(icers, export=NULL, maicers= seq(0, 50000, by=100), main="", xlab="Threshold (£/QALY)", ylab="Probability cost effective", ylim=c(0,1))
# {
#   pce <- rep(NA, length(maicers))
#   for (i in 1:length(maicers)){
#     pce[i] <- length(which(icers < maicers[i]))/length(icers) 
#   }
#   Data <- data.frame(maicer=maicers, pce=pce)
#   
#   if(!is.null(export)){jpeg(filename=export, height=800, width=800)}
#   plot(pce ~ maicer, data=Data, type="l", lwd=2, main=main, xlab=xlab, ylab=ylab, ylim=ylim)
#   if(!is.null(export)){dev.off()}
#   
#   # Want a breakdown separately at £20,000 and £30,000 / QALY
#   pce20k <- length(which(icers < 20000)) / length(icers)
#   pce30k <- length(which(icers < 30000)) / length(icers)
#   
#   qi <- list(pce20k=pce20k, pce30k=pce30k)
#   Data <- list(main=Data, qi=qi)
#   return(Data)
# }

Draw.CEAC.NB <- function(D, export=NULL, maicers= seq(0, 50000, by=100), main="", xlab="Willingness-to-Pay Threshold (£/QALY)", ylab="Probability cost effective", ylim=c(0,1))
{
  pce <- rep(NA, length(maicers))
  qaly <- D[["PSA"]]$qaly
  cost <- D[["PSA"]]$cost
  for (i in 1:length(maicers)){
    pce[i] <- length(which(qaly * maicers[i] - cost > 0)) / length(qaly * maicers[i] - cost)
  }
  Data <- data.frame(maicer=maicers, pce=pce)
  
  if(!is.null(export)){jpeg(filename=export, height=800, width=800)}
  plot(pce ~ maicer, data=Data, type="l", lwd=2, main=main, xlab=xlab, ylab=ylab, ylim=ylim)
  if(!is.null(export)){dev.off()}
  
  # Want a breakdown separately at £20,000 and £30,000 / QALY
  pce20k <- length(which(qaly * 20000 - cost > 0)) / length(qaly * 20000 - cost > 0)
  pce30k <- length(which(qaly * 30000 - cost > 0)) / length(qaly * 30000 - cost > 0)
  
  qi <- list(pce20k=pce20k, pce30k=pce30k)
  Data <- list(main=Data, qi=qi)
  return(Data)
}




# Weighted average of the proportion of LAABNs in C=0 and C=1 combined. 
#Calc.WgtFnProp <- function(Data){
#  proportions <- Data$events / Data$total
#  weights <- Data$total / sum(Data$total)
#  return(sum(proportions * weights))
#}

# Jackknife ICER
Jackknife.Icer <- function(D, qnt=c(0.025, 0.5, 0.975)){
  qaly <- D$qaly
  cost <- D$cost
  
  jack <- rep(NA, length(qaly))
  for (i in 1:length(qaly)){
    jack[i] <- mean(cost[-i]) / mean(qaly[-i])
  }
  qi =quantile(jack, qnt, na.rm=T)
  return(list(icers=jack, qi=qi))
}



# Meta function

Meta.Function <- function(Data, hrProp.Exp, hrProp.PSA, sens.Exp, spec.Exp, sens.PSA, spec.PSA){

  Data.cq.Exp <- Calc.MeanCostQaly(Data)
  this.cost <- Data.cq.Exp$cost
  this.qaly <- Data.cq.Exp$qaly
  X <- Calc.Prop(hrProp.Exp, sens=sens.Exp, spec=spec.Exp)
  # proportion high risk | score
  this.prop.base <- X$prop.base
  this.prop.comp <- X$prop.comp
  Data.Exp.icer <- Calc.CostQALY(prop.base=this.prop.base, prop.comp=this.prop.comp, cost=this.cost, qaly=this.qaly)
  sensto.SensSpec <- Test.SensSpec(hrProp=hrProp.Exp, cost=this.cost, qaly=this.qaly)
  # PSA
  Data.PSA <- Run.PSA(Data, hrProp.PSA, sens=sens.PSA, spec=spec.PSA)
  # Other quantities of interest
  Data.PSA.otherQI <- Run.PSA.otherQI(Data, hrProp.PSA, Nruns=1000, sens=sens.PSA, spec=sens.PSA)
 # Find credible intervals for icers

# Changing approach for calculating credible intervals around ICERs following another argument with Matt
# Now using jackknifing process
#  Data.icers <- Sort.ICER(Data.PSA)
  Data.icers <- Jackknife.Icer(Data.PSA)
  Data.evpi <- Calc.EVPI(D=Data, hrProp=hrProp.PSA, sens=sens.PSA, spec=spec.PSA)
  
  return(list(Expectations=Data.Exp.icer, PSA=Data.PSA, PSA.otherQI=Data.PSA.otherQI, Icers=Data.icers, EVPI=Data.evpi, sensSpec=sensto.SensSpec))
}

Extract.mainQI <- function(D){
  tmp <- c(
      base.cost = mean(D[["PSA"]]$base.cost),
      base.qaly = mean(D[["PSA"]]$base.qaly),
      comp.cost = mean(D[["PSA"]]$comp.cost),
      comp.qaly = mean(D[["PSA"]]$comp.qaly),
      mean.icer = mean(D[["PSA"]]$cost)   / mean(D[["PSA"]]$qaly),  
      D[["Icers"]]$qi[1:3]    
    )
return(tmp)  
}


############################
HrProp.Icer.graph <- function(trueHrProp=seq(0,0.2, by=0.01), main="Relationship between mean ICER and true high risk proportion", D.C, D.CV, sens.PSA, spec.PSA){

    icerVec.C <- vector("numeric", length(trueHrProp))
    icerVec.CV <- vector("numeric", length(trueHrProp))
   

    for (i in 1:length(trueHrProp)){
      mean.icer = mean(Dab.C0.Analyses[["PSA"]]$cost)   / mean(Dab.C0.Analyses[["PSA"]]$qaly)
      this.hr <- rep(trueHrProp[i],1000)
      
      this.PSA.C  <- Run.PSA(D.C,  this.hr, sens=sens.PSA, spec=spec.PSA)
      this.PSA.CV <- Run.PSA(D.CV, this.hr, sens=sens.PSA, spec=spec.PSA)
      
      this.icer.C  <- mean(this.PSA.C$cost)  / mean(this.PSA.C$qaly)
      this.icer.CV <- mean(this.PSA.CV$cost) / mean(this.PSA.CV$qaly) 

      icerVec.C[i] <- this.icer.C
      icerVec.CV[i] <- this.icer.CV
    }
    icerVec.C[icerVec.C<0] <- NA
    icerVec.CV[icerVec.CV<0] <- NA
    
    # add some code doing some stuff 
plot(icerVec.C ~ trueHrProp, type="l", 
     main=main, 
     ylab="ICER (£/QALY)", 
     ylim=c(0, 80000), 
     xlab="True proportion high risk", 
     lwd=2)
    
abline(v=0, col="grey")
abline(h=0, col="grey")
    
lines(icerVec.CV ~ trueHrProp, 
      type="l", 
      lwd=2, 
      lty="dashed")
    
legend("topright", 
       lty=c("solid", "dashed"), 
       bty="n", 
       inset=0.01, 
       col=c("black", "black", "white", "white"), 
       lwd=c(2,2), 
       legend=c(expression(CHADS[2]),
                expression(paste(CHA[2], DS[2],-VASc), sep="")) 
       )
    Data.out <- data.frame(trueHrProp=trueHrProp, icer.C=icerVec.C, icer.CV=icerVec.CV)
    return(Data.out)
}


#############################################################################################################################
######################## D A T A    A N D     V A R I A B L E S #############################################################
#############################################################################################################################


# From Providencia
# now with selectedly added prior
CHADS.Data <- data.frame(
  score = c(0, 1),
  events = c(2.5, 17),
  total = c(25, 79)
  )

CV.Data <- data.frame(
  score=c(0, 1),
  events = c(0.5, 1.5),
  total = c(12, 21)
  )


# true high risk proportions for PSA 

#attach(CHADS.Data.prior)
#hrProp.C0 <- rbeta(1000, events[score==0], total[score==0] - events[score==0])
#hrProp.C1 <- rbeta(1000, events[score==1], total[score==1] - events[score==1])
#detach(CHADS.Data.prior)

#attach(CV.Data.prior)
#hrProp.CV0 <- rbeta(1000, events[score==0], total[score==0] - events[score==0])
#hrProp.CV1 <- rbeta(1000, events[score==1], total[score==1] - events[score==1])
#detach(CV.Data.prior)

#write.csv(cbind(hrProp.C0, hrProp.C1, hrProp.CV0, hrProp.CV1), file="C:/temp/tmp10.csv")

#Calc.WgtFnProp(CV.Data)

#Calc.WgtFnProp(CHADS.Data)
##
###### Alternative sensitivity and specificity estimates:
# weighted average of sensitivity values reported in Table 3 ('structural defect') of HTA 
# sensitivity : 0.540
# specificity : 0.900 (poorly reported)

# CHADS = 0
# 2/24 patients identified as having LA ABN
# Possible sensitivity depends on assumption about true prevalence
# if true prevalence is 2/24, sensitivity is 100% (2/2)
# if true prevalence is 24/24, sensitivity is 2/24 (about 8.3%)
# I have greater believe that true answer is near 100% than 8.3%
# Use Beta distribution with 0.5 added both to event counts to calculate distribution ranging from 
# 100% to 0%.
# Then transform this onto region {8.3% to 100%}

# mean value:
#sens.Exp.prior <- 2.5/25 * (1 - 2/24) + 2/24

# PSA
#sens.PSA.prior <- rbeta(1000, 2.5, 0.5) * (1- 2/ 24) + 2/24




################################################ P S A  D A T A ##################################################################
# PSA Input data

#PSA.inputs <- read.csv("PSA_data/PSA_Input_Data.csv")
#PSA.inputs <- PSA.inputs[-c(1:2), -1]

# append uncertainty in sens and spec onto PSA.inputs
#PSA.inputs <- data.frame(PSA.inputs, sens=sens.PSA, spec=spec.PSA)
# Note: Sens and spec calculated afresh each time, so will differ each time this 
# R script is used. 
# To fix these values, they should be run once, then saved and loaded. 
#write.csv(PSA.inputs, "PSA_data/PSA_Inputs_simplified.csv")

# This is now done: use...

##### Sensitivity and specificity

# Using bottom left of Table 2 of Providencia

freqs <- c(fn=5, tn=83, tp=87, fp=159)

sens.Exp <- freqs["tp"] / (freqs["tp"] + freqs["fn"])
spec.Exp <- freqs["tn"] / (freqs["tn"] + freqs["fp"])

names(sens.Exp) <- NULL
names(spec.Exp) <- NULL 

#require(MCMCpack)
#freq.dist <- rdirichlet(1000, freqs)

#colnames(freq.dist) <- c("fn", "tn", "tp", "fp")

#sens.PSA <- apply(freq.dist, 1, function(x) x["tp"] / (x["tp"] + x["fn"]))
#spec.PSA <- apply(freq.dist, 1, function(x) x["tn"] / (x["tn"] + x["fp"]))

#plot(spec.PSA ~ sens.PSA, xlim=c(0.75,1), ylim=c(0.25,0.5))
#image(kde2d(sens.PSA, spec.PSA), xlim=c(0.75, 1), ylim=c(0.25,0.5), xlab=c("Sensitivity"), ylab="Specificity", main="Contour map of estimated sensitivity and specificity")


##########
##########

PSA.inputs <- read.csv("PSA_Inputs_simplified.csv")

sens.PSA <- PSA.inputs$Sens
spec.PSA <- PSA.inputs$Spec

hrProp.C0 <- PSA.inputs$hrProp.C0
hrProp.C1 <- PSA.inputs$hrProp.C1

hrProp.CV0 <- PSA.inputs$hrProp.CV0
hrProp.CV1 <- PSA.inputs$hrProp.CV1

# Warfarin

# 50_F
Data.WarfC0_50_F <- read.csv("W_50_F.csv")
# 50_M
Data.WarfC0_50_M <- read.csv("W_50_M.csv")
# 65_F
Data.WarfC0_65_F <- read.csv("W_65_F.csv")
# 65_M
Data.WarfC0_65_M <- read.csv("W_65_M.csv")

# Rivaroxaban

# 50_F
Data.RivC0_50_F <- read.csv("R_50_F.csv")
# 50_M
Data.RivC0_50_M <- read.csv("R_50_M.csv")
# 65_F
Data.RivC0_65_F <- read.csv("R_65_F.csv")
# 65_M
Data.RivC0_65_M <- read.csv("R_65_M.csv")

# Dabigatran

# 65_F
Data.DabC0_65_F <- read.csv("D_65_F.csv")
# 65_M
Data.DabC0_65_M <- read.csv("D_65_M.csv")


# #Warfarin
# #Data.WarfC0.old <- read.csv("PSA_data/Warf_C0_PSA.csv")
# 
# Data.WarfC0 <- read.csv("PSA_data/NewPSA_Warf_C0.csv")
# Data.WarfC1 <- read.csv("PSA_data/NewPSA_Warf_C1.csv")
# 
# Data.WarfCV0 <- read.csv("PSA_data/NewPSA_Warf_CV0.csv")
# Data.WarfCV1 <- read.csv("PSA_data/NewPSA_Warf_CV1.csv")
# 
# 
# #Data.WarfCV0 <- read.csv("PSA_data/Warf_CV0_PSA.csv") 
# 
# #Data.WarfC1 <- read.csv("PSA_data/Warf_C1_PSA.csv")
# # some issues in one line...
# 
# 
# 
# #Data.WarfCV1 <- read.csv("PSA_data/Warf_CV1_PSA.csv") 
# 
# Rivaroxaban
#Data.RivC0 <- read.csv("PSA_data/Riv_C0_PSA.csv") 
# #Data.RivCV0 <- read.csv("PSA_data/Riv_CV0_PSA.csv") 
# 
# Data.RivC0 <- read.csv("PSA_data/NewPSA_Riv_C0.csv") 
# Data.RivC1 <- read.csv("PSA_data/NewPSA_Riv_C1.csv")
# 
# #Data.RivC1 <- read.csv("PSA_data/Riv_C1_PSA.csv") 
# #Data.RivCV1 <- read.csv("PSA_data/Riv_CV1_PSA.csv") 
# 
# Data.RivCV0 <- read.csv("PSA_data/NewPSA_Riv_CV0.csv")
# Data.RivCV1 <- read.csv("PSA_data/NewPSA_Riv_CV1.csv")
# 
