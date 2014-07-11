




#############################################################################################################################
######################## D A T A    A N D     V A R I A B L E S #############################################################
#############################################################################################################################


# From Providencia
# now with selectedly added prior
providencia_data <- data.frame(
    instrument=c(
        "chads2", "chads2",
        "chadsvasc", "chadsvasc"),
    score=c(0, 1,
            0, 1),
    events=c(
        2.5, 17,
        0.5, 1.5
        ),
    total=c(
        25, 79,
        12, 21
        )
)




fn <- function(x, n_psa=1000){
    out <- rbeta(
        n=n_psa,
        shape1=x$total,
        shape2=x$total - x$events
        )
    return(out)
}

hrprop <- ddply(providencia_data , .(instrument, score), fn)
hrprop <- melt(hrprop, id.vars=c("instrument", "score"), variable.name="psa")
hrprop$psa <- str_replace(hrprop$psa, "V", "")

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

psa_inputs <- read.csv("Data/csv/PSA_Inputs_simplified.csv")
names(psa_inputs) <- tolower(names(psa_inputs))
psa_inputs$psa <- 1:dim(psa_inputs)[1]

psa_inputs <- melt(psa_inputs, id.var="psa")

tmp <- psa_inputs$variable
psa_inputs$type <- word(tmp, start=1, sep=".", end=str_locate(tmp, "[.$]")[,1])

write.csv(psa_inputs, "data/csv/psa_inputs_long.csv")


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

fn <- function(x){
    # The first . will separate the diagnostic population (e.g. tp) from other pieces of 
    # variable information
    tmp <- str_split(x, "[.]")
    tmp <- tmp[[1]] # it will produce a list by default, but always of length 1
    tmp <- tmp[-1] # remove the diag bit
    out <- paste(tmp, collapse="_") # combine the rest
    return(out)
}
# Warfarin

# 50_F
raw_data <- read.csv("Data/csv/W_50_F.csv")
names(raw_data) <- tolower(names(raw_data))
data_long <- melt(raw_data, id.var=c("run"))
tmp <- data_long$variable
data_long$diag <- word(tmp, start=1, sep=".", end=str_locate(tmp, "[.$]")[,1])
data_long$instrument="chads2"
data_long$risk_score=0
data_long$sex="female"
data_long$drug="warfarin"
data_long$age=50
data_long$measure <- sapply(data_long$variable, fn)
data_long$variable <- NULL
data_tidy <- dcast(data_long, drug + sex + age + diag + run ~ measure)



# 50_M
raw_data <- read.csv("Data/csv/W_50_M.csv")
names(raw_data) <- tolower(names(raw_data))
data_long <- melt(raw_data, id=c("run"))
tmp <- data_long$variable
data_long$diag <- word(tmp, start=1, sep=".", end=str_locate(tmp, "[.$]")[,1])
data_long$instrument="chads2"
data_long$risk_score=0
data_long$sex="male"
data_long$drug="warfarin"
data_long$age=50
data_long$measure <- sapply(data_long$variable, fn)
data_long$variable <- NULL

data_tidy.this <- dcast(data_long, drug + sex + age + diag + run ~ measure)
data_tidy <- rbind(data_tidy, data_tidy.this)


# 65_F
raw_data <- read.csv("Data/csv/W_65_F.csv")
names(raw_data) <- tolower(names(raw_data))
data_long <- melt(raw_data, id=c("run"))
tmp <- data_long$variable
data_long$diag <- word(tmp, start=1, sep=".", end=str_locate(tmp, "[.$]")[,1])
data_long$instrument="chads2"
data_long$risk_score=0
data_long$sex="female"
data_long$drug="warfarin"
data_long$age=65
data_long$measure <- sapply(data_long$variable, fn)
data_long$variable <- NULL

data_tidy.this <- dcast(data_long, drug + sex + age + diag + run ~ measure)
data_tidy <- rbind(data_tidy, data_tidy.this)

# 65_M
raw_data <- read.csv("Data/csv/W_65_M.csv")
names(raw_data) <- tolower(names(raw_data))
data_long <- melt(raw_data, id=c("run"))
tmp <- data_long$variable
data_long$diag <- word(tmp, start=1, sep=".", end=str_locate(tmp, "[.$]")[,1])
data_long$instrument="chads2"
data_long$risk_score=0
data_long$sex="male"
data_long$drug="warfarin"
data_long$age=65
data_long$measure <- sapply(data_long$variable, fn)
data_long$variable <- NULL

data_tidy.this <- dcast(data_long, drug + sex + age + diag + run ~ measure)
data_tidy <- rbind(data_tidy, data_tidy.this)


# Rivaroxaban

# 50_F
raw_data <- read.csv("Data/csv/R_50_F.csv")
names(raw_data) <- tolower(names(raw_data))
data_long <- melt(raw_data, id=c("run"))
tmp <- data_long$variable
data_long$diag <- word(tmp, start=1, sep=".", end=str_locate(tmp, "[.$]")[,1])
data_long$instrument="chads2"
data_long$risk_score=0
data_long$sex="female"
data_long$drug="rivaroxaban"
data_long$age=50
data_long$measure <- sapply(data_long$variable, fn)
data_long$variable <- NULL

data_tidy.this <- dcast(data_long, drug + sex + age + diag + run ~ measure)
data_tidy <- rbind(data_tidy, data_tidy.this)



# 50_M
raw_data <- read.csv("Data/csv/R_50_M.csv")
names(raw_data) <- tolower(names(raw_data))
data_long <- melt(raw_data, id=c("run"))
tmp <- data_long$variable
data_long$diag <- word(tmp, start=1, sep=".", end=str_locate(tmp, "[.$]")[,1])
data_long$instrument="chads2"
data_long$risk_score=0
data_long$sex="male"
data_long$drug="rivaroxaban"
data_long$age=50
data_long$measure <- sapply(data_long$variable, fn)
data_long$variable <- NULL

data_tidy.this <- dcast(data_long, drug + sex + age + diag + run ~ measure)
data_tidy <- rbind(data_tidy, data_tidy.this)

# 65_F
raw_data <- read.csv("Data/csv/R_65_F.csv")
names(raw_data) <- tolower(names(raw_data))
data_long <- melt(raw_data, id=c("run"))
tmp <- data_long$variable
data_long$diag <- word(tmp, start=1, sep=".", end=str_locate(tmp, "[.$]")[,1])
data_long$instrument="chads2"
data_long$risk_score=0
data_long$sex="female"
data_long$drug="rivaroxaban"
data_long$age=65
data_long$measure <- sapply(data_long$variable, fn)
data_long$variable <- NULL

data_tidy.this <- dcast(data_long, drug + sex + age + diag + run ~ measure)
data_tidy <- rbind(data_tidy, data_tidy.this)

# 65_M
raw_data <- read.csv("Data/csv/R_65_M.csv")
names(raw_data) <- tolower(names(raw_data))
data_long <- melt(raw_data, id=c("run"))
tmp <- data_long$variable
data_long$diag <- word(tmp, start=1, sep=".", end=str_locate(tmp, "[.$]")[,1])
data_long$instrument="chads2"
data_long$risk_score=0
data_long$sex="male"
data_long$drug="rivaroxaban"
data_long$age=65
data_long$measure <- sapply(data_long$variable, fn)
data_long$variable <- NULL
data_tidy.this <- dcast(data_long, drug + sex + age + diag + run ~ measure)
data_tidy <- rbind(data_tidy, data_tidy.this)

# Dabigatran

# 65_F
raw_data <- read.csv("Data/csv/D_65_F.csv")
names(raw_data) <- tolower(names(raw_data))
data_long <- melt(raw_data, id=c("run"))
tmp <- data_long$variable
data_long$diag <- word(tmp, start=1, sep=".", end=str_locate(tmp, "[.$]")[,1])
data_long$instrument="chads2"
data_long$risk_score=0
data_long$sex="female"
data_long$drug="dabigatran"
data_long$age=65
data_long$measure <- sapply(data_long$variable, fn)
data_long$variable <- NULL
data_tidy.this <- dcast(data_long, drug + sex + age + diag + run ~ measure)
data_tidy <- rbind(data_tidy, data_tidy.this)

# 65_M
raw_data <- read.csv("Data/csv/D_65_M.csv")
names(raw_data) <- tolower(names(raw_data))
data_long <- melt(raw_data, id=c("run"))
tmp <- data_long$variable
data_long$diag <- word(tmp, start=1, sep=".", end=str_locate(tmp, "[.$]")[,1])
data_long$instrument="chads2"
data_long$risk_score=0
data_long$sex="male"
data_long$drug="dabigatran"
data_long$age=65
data_long$measure <- sapply(data_long$variable, fn)
data_long$variable <- NULL
data_tidy.this <- dcast(data_long, drug + sex + age + diag + run ~ measure)
data_tidy <- rbind(data_tidy, data_tidy.this)



write.csv(data_tidy, file="Data/tidy/input_data_tidied.csv")
#### Now to save the tidied data


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
