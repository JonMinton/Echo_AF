




#############################################################################################################################
######################## D A T A    A N D     V A R I A B L E S #############################################################
#############################################################################################################################


# From Providencia
# now with selectedly added prior
chads_data <- data.frame(
    score = c(0, 1),
    events = c(2.5, 17),
    total = c(25, 79)
)

chadsvasc_data <- data.frame(
    score=c(0, 1),
    events = c(0.5, 1.5),
    total = c(12, 21)
)


# true high risk proportions for PSA 

attach(CHADS.Data.prior)

fn <- function(x, n_psa=1000){
    out <- rbeta(
        n=n_psa,
        shape1=x$total,
        shape2=x$total - x$events
        )
    out <- t(out)
    return(out)
}

hrprop_c <- ddply(chads_data , .(score), fn)
hrprop_c <- melt(hrprop_c, id.var="score", variable.name="psa")
hrprop <- hrprop_c
hrprop$instrument <- "chads2"
hrprop <- arrange(hrprop, instrument, score, psa, value)

hrprop_cv <- ddply(chadsvasc_data, .(score), fn)
hrprop_cv <- melt(hrprop_cv, id.var="score", variable.name="psa")
hrprop_cv$instrument <- "chadsvasc"
hrprop_cv <- arrange(hrprop_cv, instrument, score, psa, value)
hrprop <- rbind(hrprop, hrprop_cv)

write.csv(hrprop, file="Data/csv/proportion_high_risk.csv")


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
