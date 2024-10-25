#Intensive Care Unit Code#

set.seed(100)
rm(list = ls())
#dev.off()
devtools::source_url("https://github.com/mbounthavong/Decision_Analysis/blob/master/tornado_diagram_code.R?raw=TRUE")
par(cex.main=1)
library(ggplot2)
library(dplyr)

#Input values#
#Selecting which model to run - untag the desired analysis#

scenario = "basecase"

#Scenario analyses#
#scenario = "nodiffabx"
#scenario = "noduabx"
#scenario = "nodiffLOS"
#scenario = "noduLOS"
#scenario = "nodiffmort"
#scenario = "UKsepsismort"
#scenario = "PCT50"

#Subgroup analyses#
#scenario = "renal"
#scenario = "resp"
#scenario = "shock"
#scenario = "sepsis3"
#scenario = "SOFA06"
#scenario = "SOFA1024"

#SECTION 1 - Defining parameters#

#Base case standard care sample size
SC_n <- 2230
if(scenario == "renal"){
  SC_n = 767
} else if(scenario == "resp"){
  SC_n = 1434
} else if(scenario == "shock"){
  SC_n = 1593
} else if(scenario == "sepsis3"){
  SC_n = 1630
} else if(scenario == "SOFA06"){
  SC_n = 763
} else if(scenario == "SOFA1024"){
  SC_n = 486
}

# Base case mortality risk
Mort <- 0.237
if(scenario == "UKsepsismort"){
  Mort = 0.317
} else if(scenario == "renal"){
  Mort = 0.297
} else if(scenario == "resp"){
  Mort = 0.280
} else if(scenario == "shock"){
  Mort = 0.264
} else if(scenario == "sepsis3"){
  Mort = 0.244
} else if(scenario == "SOFA06"){
  Mort = 0.138
} else if(scenario == "SOFA1024"){
  Mort = 0.391
}
#Odds
MortOdds <- Mort / (1 - Mort)
#Standard error
Mort_SE <- sqrt(Mort * (1 - Mort) / SC_n)
#95% confidence interval
Mort_95LL <- Mort - qnorm(0.975) * Mort_SE
Mort_95UL <- Mort + qnorm(0.975) * Mort_SE
print(paste("95% Confidence Interval Lower Limit:", Mort_95LL))
print(paste("95% Confidence Interval Upper Limit:", Mort_95UL))
Mort_var <- (Mort-Mort_95LL)/1.96
mu <- Mort
var <- Mort_var*Mort_var

# Function to calculate alpha and beta for a Beta distribution
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(list(alpha = alpha, beta = beta))
}
Mort_betaEstimates <- estBetaParams(mu, var)
print(Mort_betaEstimates)
Mort_alpha <- Mort_betaEstimates$alpha
Mort_beta <- Mort_betaEstimates$beta
#Checking the alpha and beta parameters correspond with 95% confidence intervals
qbeta(c(0.025, 0.975), Mort_betaEstimates$alpha, Mort_betaEstimates$beta)

#Base case mortality adjusted odds ratio
MortOR <- 0.89
#95% confidence interval
MortOR_95LL <- 0.80
MortOR_95UL <- 0.99
if(scenario == "nodiffmort"){
  MortOR = 1
  MortOR_95LL = 1
  MortOR_95UL = 1
} else if(scenario == "renal"){
  MortOR = 0.96
  MortOR_95LL = 0.83
  MortOR_95UL = 1.11
} else if(scenario == "resp"){
  MortOR = 0.89
  MortOR_95LL = 0.79
  MortOR_95UL = 1
} else if(scenario == "shock"){
  MortOR = 0.89
  MortOR_95LL = 0.79
  MortOR_95UL = 1
} else if(scenario == "sepsis3"){
  MortOR = 0.86
  MortOR_95LL = 0.76
  MortOR_95UL = 0.98
} else if(scenario == "SOFA06"){
  MortOR = 0.85
  MortOR_95LL = 0.66
  MortOR_95UL = 1.1
} else if(scenario == "SOFA1024"){
  MortOR = 0.86
  MortOR_95LL = 0.72
  MortOR_95UL = 1.01
}
#Standard error
MortOR_SE <- (MortOR_95UL-MortOR)/qnorm(0.975)
#Mortality adjusted OR, margin of error and estimated standard deviation
MortOR_MoE <- (MortOR_95UL- MortOR_95LL)/2
MortOR_SD <- MortOR_MoE/qnorm(0.975)
MortOR_VAR <- MortOR_SD^2

#PCT mortality odds
PCTMortOdds <- MortOdds * MortOR
#PCT mortality risk
PCTMort <- PCTMortOdds/(PCTMortOdds+1)
#Standard care survival risk
Surv <- 1-Mort
#PCT survival risk
PCTSurv <- 1-(PCTMort)

#Antibiotic duration (days)
dAbx <- 10.4
dAbx_SD <- 9.7
if(scenario == "renal"){
  dAbx = 12.3
  dAbx_SD = 10.2
} else if(scenario == "resp"){
  dAbx = 10.2
  dAbx_SD = 9.0
} else if(scenario == "shock"){
  dAbx = 10.4
  dAbx_SD = 9.8
} else if(scenario == "sepsis3"){
  dAbx = 10.5
  dAbx_SD = 9.2
} else if(scenario == "SOFA06"){
  dAbx = 10.9
  dAbx_SD = 10.3
} else if(scenario == "SOFA1024"){
  dAbx = 10.7
  dAbx_SD = 9.1
}

#Antibiotic duration (days), standard error
dAbx_SE <- dAbx_SD/sqrt(SC_n)
#Antibiotic duration (days), variance (SE)
dAbx_VAR <- dAbx_SE**2
#Antibiotic duration (days), Gamma distribution parameters
dAbx_SCALE <- dAbx_VAR/dAbx
dAbx_SHAPE <- dAbx/dAbx_SCALE
dAbx_RATE <- 1/dAbx_SCALE
#Antibiotic duration (days), 95% confidence interval
dAbx_95LL <- dAbx - qnorm(0.975)*dAbx_SE
dAbx_95UL <- dAbx + qnorm(0.975)*dAbx_SE

#Difference in antibiotic duration (days)
diffAbx <- -1.19
#Difference in antibiotic duration (days), 95% confidence interval
diffAbx_95LL <- -1.73
diffAbx_95UL <- -0.66
if(scenario == "nodiffabx"){
  diffAbx = 0
  diffAbx_95LL = 0
  diffAbx_95UL = 0
} else if(scenario == "renal"){
  diffAbx = -1.45
  diffAbx_95LL = -2.44
  diffAbx_95UL = -0.46
} else if(scenario == "resp"){
  diffAbx = -1.09
  diffAbx_95LL = -1.73
  diffAbx_95UL = -0.45
} else if(scenario == "shock"){
  diffAbx = -1.03
  diffAbx_95LL = -1.68
  diffAbx_95UL = -0.39
} else if(scenario == "sepsis3"){
  diffAbx = -1.22
  diffAbx_95LL = -1.82
  diffAbx_95UL = -0.62
} else if(scenario == "SOFA06"){
  diffAbx = -2.62
  diffAbx_95LL = -3.51
  diffAbx_95UL = -1.73
} else if(scenario == "SOFA1024"){
  diffAbx = -0.63
  diffAbx_95LL = -1.71
  diffAbx_95UL = 0.45
}
#Difference in antibiotic duration (days), margin of error and estimated standard deviation
diffAbx_MoE <- (diffAbx_95UL- diffAbx_95LL)/2
diffAbx_SD <- diffAbx_MoE/qnorm(0.975)
diffAbx_VAR <- diffAbx_SD^2

#PCT antibiotic duration (days)
PCTdAbx <- dAbx + diffAbx

#Total length of stay (days)
LOS <- 28.7
#Total length of stay (days), SD
LOS_SD <- 27.9
if(scenario == "renal"){
  LOS = 26.6
  LOS_SD = 10.2
} else if(scenario == "resp"){
  LOS = 30.9
  LOS_SD = 30.9
} else if(scenario == "shock"){
  LOS = 30.1
  LOS_SD = 27.0
} else if(scenario == "sepsis3"){
  LOS = 29.5
  LOS_SD = 27.9
} else if(scenario == "SOFA06"){
  LOS = 28.6
  LOS_SD = 28.6
} else if(scenario == "SOFA1024"){
  LOS = 29.4
  LOS_SD = 27.5
}
#Total length of stay (days), SE
LOS_SE <- LOS_SD/sqrt(SC_n)
#Total length of stay (days), variance
LOS_VAR <- LOS_SE**2
#Total length of stay (days), gamma distribution parameters
LOS_SCALE <- LOS_VAR/LOS
LOS_SHAPE <- LOS/LOS_SCALE
LOS_RATE <- 1/LOS_SCALE
#Total length of stay (days), 95% confidence interval
LOS_95LL <- LOS - qnorm(0.975)*LOS_SE
LOS_95UL <- LOS + qnorm(0.975)*LOS_SE

#Difference in total length of stay (days)
diffLOS <- 0.09
#Difference in total length of stay (days), 95% confidence interval
diffLOS_95LL <- -1.51
diffLOS_95UL <- 1.7
if(scenario == "nodiffLOS"){
  diffLOS = 0
  diffLOS_95LL = 0
  diffLOS_95UL = 0
} else if(scenario == "renal"){
  diffLOS = 2.97
  diffLOS_95LL = 0.31
  diffLOS_95UL = 5.63
} else if(scenario == "resp"){
  diffLOS = -0.39
  diffLOS_95LL = -2.42
  diffLOS_95UL = 1.65
} else if(scenario == "shock"){
  diffLOS = 0.6
  diffLOS_95LL = -1.23
  diffLOS_95UL = 2.43
} else if(scenario == "sepsis3"){
  diffLOS = 0.07
  diffLOS_95LL = -1.89
  diffLOS_95UL = 2.03
} else if(scenario == "SOFA06"){
  diffLOS = -1.96
  diffLOS_95LL = -4.65
  diffLOS_95UL = 0.72
} else if(scenario == "SOFA1024"){
  diffLOS = 3.21
  diffLOS_95LL = -0.76
  diffLOS_95UL = 7.18
}
#Difference in total length of stay (days), margin of error and estimated standard deviation
diffLOS_MoE <- (diffLOS_95UL- diffLOS_95LL)/2
diffLOS_SD <- diffLOS_MoE/qnorm(0.975)
diffLOS_VAR <- diffLOS_SD^2
#PCT total length of stay (days)
PCTLOS <- LOS + diffLOS

#ICU length of stay (days)
ICU <- 14.7
#ICU length of stay (days), SD
ICU_SD <- 16.3
if(scenario == "renal"){
  ICU = 14.9
  ICU_SD = 14.0
} else if(scenario == "resp"){
  ICU = 16.2
  ICU_SD = 17.5
} else if(scenario == "shock"){
  ICU = 15
  ICU_SD = 16.1
} else if(scenario == "sepsis3"){
  ICU = 14.1
  ICU_SD = 15.5
} else if(scenario == "SOFA06"){
  ICU = 12.9
  ICU_SD = 16.7
} else if(scenario == "SOFA1024"){
  ICU = 15.6
  ICU_SD = 15.4
}
#ICU length of stay (days), SE
ICU_SE <- ICU_SD/sqrt(SC_n)
#ICU length of stay (days), variance
ICU_VAR <- ICU_SE**2
#ICU length of stay (days), gamma distribution parameters
ICU_SCALE <- ICU_VAR/ICU
ICU_SHAPE <- ICU/ICU_SCALE
ICU_RATE <- 1/ICU_SCALE
#ICU length of stay (days), 95% confidence interval
ICU_95LL <- ICU - qnorm(0.975)*ICU_SE
ICU_95UL <- ICU + qnorm(0.975)*ICU_SE

#Difference in ICU length of stay (days)
diffICU <- 0.04
#Difference in ICU length of stay (days), 95% confidence interval
diffICU_95LL <- -0.9
diffICU_95UL <- 0.99
if(scenario == "nodiffLOS"){
  diffICU = 0
  diffICU_95LL = 0
  diffLOS_95UL = 0
} else if(scenario == "renal"){
  diffICU = 1.43
  diffICU_95LL = -0.16
  diffICU_95UL = 3.02
} else if(scenario == "resp"){
  diffICU = -0.3
  diffICU_95LL = -1.56
  diffICU_95UL = 0.95
} else if(scenario == "shock"){
  diffICU = 0.34
  diffICU_95LL = -0.8
  diffICU_95UL = 1.49
} else if(scenario == "sepsis3"){
  diffICU = 0.37
  diffICU_95LL = -0.74
  diffICU_95UL = 1.48
} else if(scenario == "SOFA06"){
  diffICU = -0.92
  diffICU_95LL = -2.52
  diffICU_95UL = 0.69
} else if(scenario == "SOFA1024"){
  diffICU = 2.08
  diffICU_95LL = -0.08
  diffICU_95UL = 4.24
}
#Difference in ICU length of stay (days), margin of error and estimated standard deviation
diffICU_MoE <- (diffICU_95UL- diffICU_95LL)/2
diffICU_SD <- diffICU_MoE/qnorm(0.975)
diffICU_VAR <- diffICU_SD^2
#PCT ICU length of stay (days)
PCTICU <- ICU + diffICU
#Hospital length of stay (days)
Hos <- LOS - ICU
# PCT hospital length of stay (days)
PCTHos <- PCTLOS - PCTICU

#Hospital cost (day)
cHos <- 767.31
#ICU cost (day)
cICU <- 2095.76
#Antibiotic cost (day)
cAbx <- 6.38
#PCT cost
cPCT <- 18.23
if(scenario == "PCT50"){
  cPCT <- 50
}
#Number of PCT tests
nPCT <- 7

#Proportion of female sepsis patients
SepPF <- 0.45
#Proportion of female sepsis patients, sample size
SepPF_n <- 197142
if(scenario == "shock" | scenario == "renal" | scenario == "resp" | scenario == "SOFA1024"){
  SepPF = 0.446
  SepPF_n = 39262
}
#Proportion of female sepsis patients, SE
SepPF_SE <- sqrt(SepPF*(1-SepPF)/SepPF_n)
#Proportion of female sepsis patients, 95% confidence interval
SepPF_95LL <- SepPF - qnorm(0.975)*SepPF_SE
SepPF_95UL <- SepPF + qnorm(0.975)*SepPF_SE
print(paste("95% Confidence Interval Lower Limit:", SepPF_95LL))
print(paste("95% Confidence Interval Upper Limit:", SepPF_95UL))
SepPF_var <- ((SepPF-SepPF_95LL)/1.96)
mu <- SepPF
var <- SepPF_var*SepPF_var
# Function to calculate alpha and beta for a Beta distribution
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(list(alpha = alpha, beta = beta))
}
SepPF_betaEstimates <- estBetaParams(mu, var)
print(SepPF_betaEstimates)
SepPF_alpha <- SepPF_betaEstimates$alpha
SepPF_beta <- SepPF_betaEstimates$beta
#Checking the alpha and beta parameters correspond with 95% confidence intervals
qbeta(c(0.025, 0.975), SepPF_betaEstimates$alpha, SepPF_betaEstimates$beta)
#Proportion of male sepsis patients
SepPM <- 1-SepPF

#Female QALY 
uF <- 0.776
#Female QALY, sample size
uF_n <- 608
#Female QALY 95% confidence interval
uF_95LL <- 0.769
uF_95UL <- 0.797
#Female QALY standard error
uF_SE <- sqrt(uF*(1-uF)/uF_n)
if(scenario == "shock" | scenario == "renal" | scenario == "resp" | scenario == "SOFA1024"){
  uF = 0.775
  uF_n = 619
  uF_95LL = 0.770
  uF_95UL = 0.795
}
uF_var <- (uF-uF_95LL)/1.96
mu <- uF
var <- uF_var*uF_var

# Function to calculate alpha and beta for a Beta distribution
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(list(alpha = alpha, beta = beta))
}
uF_betaEstimates <- estBetaParams(mu, var)
print(uF_betaEstimates)
uF_alpha <- uF_betaEstimates$alpha
uF_beta <- uF_betaEstimates$beta
#Checking the alpha and beta parameters correspond with 95% confidence intervals
qbeta(c(0.025, 0.975), uF_betaEstimates$alpha, uF_betaEstimates$beta)


#Male QALY 
uM <- 0.803
#Male QALY, sample size
uM_n <- 532
#Male QALY 60-64 years old, 95% confidence interval
uM_95LL <- 0.798
uM_95UL <- 0.822
if(scenario == "shock" | scenario == "renal" | scenario == "resp" | scenario == "SOFA1024"){
  uM = 0.797
  uM_n = 568
  uM_95LL = 0.792
  uM_95UL = 0.818
}
#Standard care mortality rate, standard error
uM_SE <- sqrt(uM*(1-uM)/uM_n)
uM_var <- (uM-uM_95LL)/1.96
mu <- uM
var <- uM_var*uM_var

# Function to calculate alpha and beta for a Beta distribution
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(list(alpha = alpha, beta = beta))
}
uM_betaEstimates <- estBetaParams(mu, var)
print(uM_betaEstimates)
uM_alpha <- uM_betaEstimates$alpha
uM_beta <- uM_betaEstimates$beta
#Checking the alpha and beta parameters correspond with 95% confidence intervals
qbeta(c(0.025, 0.975), uM_betaEstimates$alpha, uM_betaEstimates$beta)

#Weighted population QALY 
uW <- (uF*SepPF) + (uM*SepPM)

#Sepsis utility
uSep <- 0.53
#Sepsis utility, sample size
uSep_n <- 701
#Sepsis utility, SE
uSep_SE <- sqrt(uSep*(1-uSep)/uSep_n)
#Sepsis utility, 95% confidence interval
uSep_95LL <- uSep - qnorm(0.975)*uSep_SE
uSep_95UL <- uSep + qnorm(0.975)*uSep_SE
print(paste("95% Confidence Interval Lower Limit:", uSep_95LL))
print(paste("95% Confidence Interval Upper Limit:", uSep_95UL))
uSep_var <- (uSep-uSep_95LL)/1.96
mu <- uSep
var <- uSep_var*uSep_var
# Function to calculate alpha and beta for a Beta distribution
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(list(alpha = alpha, beta = beta))
}
uSep_betaEstimates <- estBetaParams(mu, var)
print(uSep_betaEstimates)
uSep_alpha <- uSep_betaEstimates$alpha
uSep_beta <- uSep_betaEstimates$beta
#Checking the alpha and beta parameters correspond with 95% confidence intervals
qbeta(c(0.025, 0.975), uSep_betaEstimates$alpha, uSep_betaEstimates$beta)
#Sepsis disutility
duSep <- (uW-uSep)/12

#Hospital disutility (QALD/day)
duHosD <- 0.1
if(scenario == "noduLOS"){
  duHosD = 0
}
duHosD_n <- 159
duHosD_SE <- sqrt(duHosD*(1-duHosD)/duHosD_n)
duHosD_95LL <- duHosD - qnorm(0.975)*duHosD_SE
duHosD_95UL <- duHosD + qnorm(0.975)*duHosD_SE
print(paste("95% Confidence Interval Lower Limit:", duHosD_95LL))
print(paste("95% Confidence Interval Upper Limit:", duHosD_95UL))
duHosD_var <- (duHosD-duHosD_95LL)/1.96
mu <- duHosD
var <- duHosD_var*duHosD_var

# Function to calculate alpha and beta for a Beta distribution
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(list(alpha = alpha, beta = beta))
}
duHosD_betaEstimates <- estBetaParams(mu, var)
print(duHosD_betaEstimates)
duHosD_alpha <- duHosD_betaEstimates$alpha
duHosD_beta <- duHosD_betaEstimates$beta
#Checking the alpha and beta parameters correspond with 95% confidence intervals
qbeta(c(0.025, 0.975), duHosD_betaEstimates$alpha, duHosD_betaEstimates$beta)

#Antibiotic disutility (QALD/day)
duAbxD <- 0.057
if(scenario == "noduabx"){
  duAbxD = 0
}
duAbxD_n <- 349
duAbxD_SE <- sqrt(duAbxD*(1-duAbxD)/duAbxD_n)
duAbxD_95LL <- duAbxD - qnorm(0.975)*duAbxD_SE
duAbxD_95UL <- duAbxD + qnorm(0.975)*duAbxD_SE
print(paste("95% Confidence Interval Lower Limit:", duAbxD))
print(paste("95% Confidence Interval Upper Limit:", duAbxD))

duAbxD_var <- (duAbxD-duAbxD_95LL)/1.96
mu <- duAbxD
var <- duAbxD_var*duAbxD_var

# Function to calculate alpha and beta for a Beta distribution
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(list(alpha = alpha, beta = beta))
}
duAbxD_betaEstimates <- estBetaParams(mu, var)
print(duAbxD_betaEstimates)
duAbxD_alpha <- duAbxD_betaEstimates$alpha
duAbxD_beta <- duAbxD_betaEstimates$beta
#Checking the alpha and beta parameters correspond with 95% confidence intervals
qbeta(c(0.025, 0.975), duHosD_betaEstimates$alpha, duHosD_betaEstimates$beta)

#Total hospital disutility (QALY)
duHos <- duHosD*LOS/365.25
#Total antibiotic disutility (QALY)
duAbx <- duAbxD*dAbx/365.25
#PCT total hospital disutility (QALY)
PCTduHos <- duHosD*PCTLOS/365.25
#PCT total antibiotic disutility (QALY)
PCTduAbx <- duAbxD*PCTdAbx/365.25

#Female QALE, 3.5% discount
QALEF <- 11.64
#Male QALE, 3.5% discount
QALEM <- 11.22
if(scenario == "shock" | scenario == "renal" | scenario == "resp" | scenario == "SOFA1024"){
  QALEF = 10.99
  QALEM = 10.56
}
#Weighted QALE, 3.5% discount
QALEW <- (QALEF*SepPF) + (QALEM*SepPM)

#Female LE 
LEF <- 23
#Male LE 
LEM <- 20.51
if(scenario == "shock" | scenario == "renal" | scenario == "resp" | scenario == "SOFA1024"){
  LEF = 21.3
  LEM = 18.91
}
#Weighted LE 
LEW <- (LEF*SepPF) + (LEM*SepPM)

#SECTION 2 - Data frame for model inputs#

input <- data.frame(
  #1
  Mort <- Mort,
  MortOdds <- Mort/(1-Mort),
  MortOR <- MortOR,
  PCTMortOdds <- MortOdds * MortOR,
  PCTMort <- PCTMortOdds/(PCTMortOdds+1),
  #6
  Surv <- 1-Mort,
  PCTSurv <- 1-(PCTMort),
  dAbx <- dAbx,
  diffAbx <- diffAbx,
  PCTdAbx <- dAbx + diffAbx,
  #11
  LOS <- LOS,
  diffLOS <- diffLOS,
  PCTLOS <- LOS + diffLOS,
  ICU <- ICU,
  diffICU <- diffICU,
  #16
  PCTICU <- ICU + diffICU,
  Hos <- LOS - ICU,
  PCTHos <- PCTLOS - PCTICU,
  cHos <- cHos,
  cICU <- cICU, 
  #21
  cAbx <- cAbx,
  cPCT <- cPCT,
  nPCT <- nPCT,
  SepPF <- SepPF,
  SepPM <- 1-SepPF,
  #26
  uF <- uF,
  uM <- uM,
  uW <- (uF*SepPF) + (uM*SepPM),
  uSep <- uSep,
  duSep <- (uW-uSep)/12,
  #31
  duHos <- duHosD*LOS/365.25,
  duAbx <- duAbxD*dAbx/365.25,
  PCTduHos <- duHosD*PCTLOS/365.25,
  PCTduAbx <- duAbxD*PCTdAbx/365.25,
  QALEW <- (QALEF*SepPF) + (QALEM*SepPM),
  LEW <- (LEF*SepPF) + (LEM*SepPM)
)

#Section 3 - Decision Tree 
dec_tree <- function(params){
  with(
    as.list(params), 
    {
      #Input parameters
      Mort <- Mort
      MortOdds <- Mort/(1-Mort)
      MortOR <- MortOR
      PCTMortOdds <- MortOdds * MortOR
      PCTMort <- PCTMortOdds/(PCTMortOdds+1)
      Surv <- 1-Mort
      PCTSurv <- 1-(PCTMort)
      dAbx <- dAbx
      diffAbx <- diffAbx
      PCTdAbx <- dAbx + diffAbx
      LOS <- LOS
      diffLOS <- diffLOS
      PCTLOS <- LOS + diffLOS
      ICU <- ICU
      diffICU <- diffICU
      PCTICU <- ICU + diffICU
      Hos <- LOS - ICU
      PCTHos <- PCTLOS - PCTICU
      cHos <- cHos
      cICU <- cICU
      cAbx <- cAbx
      cPCT <- cPCT
      nPCT <- nPCT
      SepPF <- SepPF
      SepPM <- 1-SepPF
      uF <- uF
      uM <- uM
      uW <- (uF*SepPF) + (uM*SepPM)
      uSep <- uSep
      duSep <- (uW-uSep)/12
      duHos <- duHosD*LOS/365.25
      duAbx <- duAbxD*dAbx/365.25
      PCTduHos <- duHosD*PCTLOS/365.25
      PCTduAbx <- duAbxD*PCTdAbx/365.25
      uW <- (uF*SepPF) + (uM*SepPM)
      LEW <- (LEF*SepPF) + (LEM*SepPM)
      
      
      #Standard care costs
      C_SC <- (cAbx*dAbx) + (cHos*Hos) + (cICU*ICU)
      #PCT costs
      C_PCT <- (cPCT*nPCT) + (cAbx*PCTdAbx) + (cHos*PCTHos) + (cICU*PCTICU)
      #Standard care QALYs gained
      Q_SC <- (QALEW-duSep-duAbx-duHos)*Surv
      #PCT QALYs gained
      Q_PCT <- (QALEW-duSep-PCTduAbx-PCTduHos)*PCTSurv
      #Incremental Costs
      IC <- C_PCT - C_SC
      #Incremental QALYs
      IQALYx <- Q_PCT - Q_SC
      IQALY <- round(IQALYx, 4)
      
      if (IQALY == 0) {
        ICER <- NA
        ICER_flag <- "Excluded: IQALY of 0"
      } else {
        # Incremental Costs
        IC <- C_PCT - C_SC
        # ICER
        ICER <- IC / IQALY
        # Normal flag for non-excluded cases
        ICER_flag <- "Normal"
      }
      #Incremental Antibiotics
      IAbx <- PCTdAbx - dAbx 
      #Cost per antibiotic day avoided
      AbxAv <- (C_PCT - C_SC)/(dAbx-PCTdAbx)
      #Antibiotic value
      AbxV <- (IC-(IQALY*20000))/IAbx
      Inputs <- c(Mort, PCTMort, LOS, PCTLOS, ICU, PCTICU, Hos, PCTHos)
      BasicOutcomes <- c(C_SC, C_PCT, Q_SC, Q_PCT, dAbx, PCTdAbx)
      ICERs <- c(IC, IQALY, IAbx, ICER, AbxAv, AbxV)
      
      names(Inputs) <- paste ("Inputs", c("Mortality", "Procalcitonin mortality", "Length of stay", "Procalcitonin length of stay", "Length of ICU stay", "Procalcitonin length of ICU stay", "Length of regular ward stay", "Procalcitonin legnth of regular ward stay"))
      names(BasicOutcomes) <- paste ("Basic Outcomes", c("Standard care costs", "Procalcitonin costs", "Standard care QALYs", "Procalcitonin QALYs", "Standard care antibiotic duration", "Procalcitonin antibiotic duration"), sep = "_")
      names(ICERs)  <- paste("ICERs", c("Incremental costs", "Incremental QALYs", "Incremental antibiotics", "ICER", "Cost per antibiotic day avoided", "Antibiotic value"), sep = "_")
      
      return(c(Inputs, BasicOutcomes, ICERs))
    }
  )
}
options(scipen=999)
dec_tree(input)

#Section 4 - Tornado Plot#
if (scenario == "basecase"){
  Mort_range <- c(BaseCase = Mort, low = Mort_95LL, high = Mort_95UL)
  MortOR_range <- c(BaseCase = MortOR, low = MortOR_95LL, high = MortOR_95UL)
  dAbx_range <- c(BaseCase = dAbx, low = dAbx_95LL, high = dAbx_95UL)
  diffAbx_range <- c(BaseCase = diffAbx, low = diffAbx_95LL, high = diffAbx_95UL)
  LOS_range <- c(BaseCase = LOS, low = LOS_95LL, high = LOS_95UL)
  diffLOS_range <- c(BaseCase = diffLOS, low = diffLOS_95LL, high = diffLOS_95UL)
  ICU_range <- c(BaseCase= ICU, low = ICU_95LL, high = ICU_95UL)
  diffICU_range <- c(BaseCase = diffICU, low = diffICU_95LL, high = diffICU_95UL)
  cHos_range <- c(BaseCase = cHos, low = (cHos*0.75), high = (cHos*1.25))
  cICU_range <- c(BaseCase = cICU, low = (cICU*0.75), high = (cICU*1.25))
  cAbx_range <- c(BaseCase = cAbx, low = (cAbx*0.75), high = (cAbx*1.25))
  cPCT_range <- c(BaseCase = cPCT, low = 14.5, high = 90)
  nPCT_range <- c(BaseCase = nPCT, low = 1, high = 10)
  SepPF_range <- c(BaseCase = SepPF, low = SepPF_95LL, high = SepPF_95UL)
  uF_range <- c(BaseCase = uF, low = uF_95LL, high = uF_95UL)
  uM_range <- c(BaseCase = uM, low = uM_95LL, high = uM_95UL)
  duSep_range <- c(BaseCase = duSep, low = (duSep/30), high = duSep*2)
  
  #Plotting all parameters to find most influential
  paramNames <- c(
    "30-day mortality probability",
    "30-day mortality odds ratio",
    "Duration of antibiotics",
    "Antibiotic duration difference",
    "Hospital length of stay",
    "Hospital length of stay difference",
    "ICU length of stay",
    "ICU length of stay difference",
    "Hospital cost",
    "ICU cost",
    "Antibiotic cost",
    "Procalcitonin cost per test",
    "Number of procalcitonin tests used",
    "Proportion of female patients",
    "Female utility",
    "Male utility",
    "Sepsis disutility"
  )
  
  l.tor.in <- vector("list", 17)
  names(l.tor.in) <- paramNames
  
  l.tor.in$'30-day mortality probability' <- cbind(Mort = Mort_range, input [-1])
  l.tor.in$'30-day mortality odds ratio' <- cbind(MortOR = MortOR_range, input [-3])
  l.tor.in$'Duration of antibiotics' <- cbind(dAbx = dAbx_range, input [-8])
  l.tor.in$'Antibiotic duration difference' <- cbind(diffAbx = diffAbx_range, input [-9])
  l.tor.in$'Hospital length of stay' <-cbind(LOS = LOS_range, input [-11])
  l.tor.in$'Hospital length of stay difference' <-cbind(diffLOS = diffLOS_range, input [-12])
  l.tor.in$'ICU length of stay' <-cbind(ICU = ICU_range, input [-14])
  l.tor.in$'ICU length of stay difference' <-cbind(diffICU = diffICU_range, input [-15])
  l.tor.in$'Hospital cost' <-cbind(cHos = cHos_range, input [-19])
  l.tor.in$'ICU cost' <-cbind(cICU = cICU_range, input [-20])
  l.tor.in$'Antibiotic cost' <-cbind(cAbx = cAbx_range, input [-21])
  l.tor.in$'Procalcitonin cost per test' <-cbind(cPCT = cPCT_range, input [-22])
  l.tor.in$'Number of procalcitonin tests used' <-cbind(nPCT = nPCT_range, input [-23])
  l.tor.in$'Proportion of female patients' <-cbind(SepPF = SepPF_range, input [-24])
  l.tor.in$'Female utility' <-cbind(uF = uF_range, input [-26])
  l.tor.in$'Male utility' <-cbind(uM = uM_range, input [-27])
  l.tor.in$'Sepsis disutility' <-cbind(duSep = duSep_range, input [-30])
  
  #List of outputs
  l.tor.out <- vector("list", 17)
  names(l.tor.out) <- paramNames
  
  for(i in 1:17){
    l.tor.out[[i]] <- t(apply(l.tor.in[[i]], 1, dec_tree))[ , 18]
  }
  
  m.tor <- matrix(unlist(l.tor.out), nrow = 17, ncol = 3, byrow = TRUE, 
                  dimnames = list(paramNames, c("basecase", "low", "high")))
  
  TornadoPlot(main_title = "", Parms = paramNames, Outcomes = m.tor, 
              outcomeName = "", 
              xlab = "", 
              ylab = "", 
              col1="#3182bd", col2="#6baed6"
  )
  
  #Plotting Top 8 parameters
  paramNames <- c(
    "30-day mortality probability",
    "30-day mortality probability odds ratio",
    "Difference in total length of stay",
    "Difference in ICU length of stay",
    "Cost per regular ward day",
    "Cost per ICU day",
    "Cost per procalcitonin test",
    "Number of procalcitonin tests used"
  )
  
  l.tor.in <- vector("list", 8)
  names(l.tor.in) <- paramNames
  
  
  l.tor.in$'30-day mortality probability' <- cbind(Mort = Mort_range, input [-1])
  l.tor.in$'30-day mortality probability odds ratio' <- cbind(MortOR = MortOR_range, input [-3])
  l.tor.in$'Difference in total length of stay' <-cbind(diffLOS = diffLOS_range, input [-12])
  l.tor.in$'Difference in ICU length of stay' <-cbind(diffICU = diffICU_range, input [-15])
  l.tor.in$'Cost per regular ward day' <-cbind(cHos = cHos_range, input [-19])
  l.tor.in$'Cost per ICU day' <-cbind(cICU = cICU_range, input [-20])
  l.tor.in$'Cost per procalcitonin test' <-cbind(cPCT = cPCT_range, input [-22])
  l.tor.in$'Number of procalcitonin tests used' <-cbind(nPCT = nPCT_range, input [-23])
  
  #List of outputs
  l.tor.out <- vector("list", 8)
  names(l.tor.out) <- paramNames
  
  for(i in 1:8){
    l.tor.out[[i]] <- t(apply(l.tor.in[[i]], 1, dec_tree))[ , 18]
  }
  
  m.tor <- matrix(unlist(l.tor.out), nrow = 8, ncol = 3, byrow = TRUE, 
                  dimnames = list(paramNames, c("basecase", "low", "high")))
  
  TornadoPlot(main_title = "", Parms = paramNames, Outcomes = m.tor, 
              outcomeName = "", 
              xlab = "", 
              ylab = "", 
              col1="#3182bd", col2="#6baed6"
  )}

#Section 5 - Monte Carlo Simulation

#Defining parallel parameters to prevent simulation iterating on itself
#Base case mortality adjusted odds ratio
mcMortOR <- 0.89
#95% confidence interval
mcMortOR_95LL <- 0.80
mcMortOR_95UL <- 0.99
if(scenario == "nodiffmort"){
  mcMortOR = 1
  mcMortOR_95LL = 1
  mcMortOR_95UL = 1
} else if(scenario == "renal"){
  mcMortOR = 0.96
  mcMortOR_95LL = 0.83
  mcMortOR_95UL = 1.11
} else if(scenario == "resp"){
  mcMortOR = 0.89
  mcMortOR_95LL = 0.79
  mcMortOR_95UL = 1
} else if(scenario == "shock"){
  mcMortOR = 0.89
  mcMortOR_95LL = 0.79
  mcMortOR_95UL = 1
} else if(scenario == "sepsis3"){
  mcMortOR = 0.86
  mcMortOR_95LL = 0.76
  mcMortOR_95UL = 0.98
} else if(scenario == "SOFA06"){
  mcMortOR = 0.85
  mcMortOR_95LL = 0.66
  mcMortOR_95UL = 1.1
} else if(scenario == "SOFA1024"){
  mcMortOR = 0.86
  mcMortOR_95LL = 0.72
  mcMortOR_95UL = 1.01
}
#Standard error
mcMortOR_SE <- (mcMortOR_95UL-mcMortOR)/qnorm(0.975)
#Mortality adjusted OR, margin of error and estimated standard deviation
mcMortOR_MoE <- (mcMortOR_95UL- mcMortOR_95LL)/2
mcMortOR_SD <- mcMortOR_MoE/qnorm(0.975)
mcMortOR_VAR <- mcMortOR_SD^2


#Difference in antibiotic duration (days)
mcdiffAbx <- -1.19
#Difference in antibiotic duration (days), 95% confidence interval
mcdiffAbx_95LL <- -1.73
mcdiffAbx_95UL <- -0.66
if(scenario == "nodiffabx"){
  mcdiffAbx = 0
  mcdiffAbx_95LL = 0
  mcdiffAbx_95UL = 0
} else if(scenario == "renal"){
  mcdiffAbx = -1.45
  mcdiffAbx_95LL = -2.44
  mcdiffAbx_95UL = -0.46
} else if(scenario == "resp"){
  mcdiffAbx = -1.09
  mcdiffAbx_95LL = -1.73
  mcdiffAbx_95UL = -0.45
} else if(scenario == "shock"){
  mcdiffAbx = -1.03
  mcdiffAbx_95LL = -1.68
  mcdiffAbx_95UL = -0.39
} else if(scenario == "sepsis3"){
  mcdiffAbx = -1.22
  mcdiffAbx_95LL = -1.82
  mcdiffAbx_95UL = -0.62
} else if(scenario == "SOFA06"){
  mcdiffAbx = -2.62
  mcdiffAbx_95LL = -3.51
  mcdiffAbx_95UL = -1.73
} else if(scenario == "SOFA1024"){
  mcdiffAbx = -0.63
  mcdiffAbx_95LL = -1.71
  mcdiffAbx_95UL = 0.45
}
#Difference in antibiotic duration (days), margin of error and estimated standard deviation
mcdiffAbx_MoE <- (mcdiffAbx_95UL- mcdiffAbx_95LL)/2
mcdiffAbx_SD <- mcdiffAbx_MoE/qnorm(0.975)
mcdiffAbx_VAR <- mcdiffAbx_SD^2


#Difference in total length of stay (days)
mcdiffLOS <- 0.09
#Difference in total length of stay (days), 95% confidence interval
mcdiffLOS_95LL <- -1.51
mcdiffLOS_95UL <- 1.7
if(scenario == "nodiffLOS"){
  mcdiffLOS = 0
  mcdiffLOS_95LL = 0
  mcdiffLOS_95UL = 0
} else if(scenario == "renal"){
  mcdiffLOS = 2.97
  mcdiffLOS_95LL = 0.31
  mcdiffLOS_95UL = 5.63
} else if(scenario == "resp"){
  mcdiffLOS = -0.39
  mcdiffLOS_95LL = -2.42
  mcdiffLOS_95UL = 1.65
} else if(scenario == "shock"){
  mcdiffLOS = 0.6
  mcdiffLOS_95LL = -1.23
  mcdiffLOS_95UL = 2.43
} else if(scenario == "sepsis3"){
  mcdiffLOS = 0.07
  mcdiffLOS_95LL = -1.89
  mcdiffLOS_95UL = 2.03
} else if(scenario == "SOFA06"){
  mcdiffLOS = -1.96
  mcdiffLOS_95LL = -4.65
  mcdiffLOS_95UL = 0.72
} else if(scenario == "SOFA1024"){
  mcdiffLOS = 3.21
  mcdiffLOS_95LL = -0.76
  mcdiffLOS_95UL = 7.18
}
#Difference in total length of stay (days), margin of error and estimated standard deviation
mcdiffLOS_MoE <- (mcdiffLOS_95UL- mcdiffLOS_95LL)/2
mcdiffLOS_SD <- mcdiffLOS_MoE/qnorm(0.975)
mcdiffLOS_VAR <- mcdiffLOS_SD^2


#Difference in ICU length of stay (days)
mcdiffICU <- 0.04
#Difference in ICU length of stay (days), 95% confidence interval
mcdiffICU_95LL <- -0.9
mcdiffICU_95UL <- 0.99
if(scenario == "nodiffLOS"){
  mcdiffICU = 0
  mcdiffICU_95LL = 0
  mcdiffLOS_95UL = 0
} else if(scenario == "renal"){
  mcdiffICU = 1.43
  mcdiffICU_95LL = -0.16
  mcdiffICU_95UL = 3.02
} else if(scenario == "resp"){
  mcdiffICU = -0.3
  mcdiffICU_95LL = -1.56
  mcdiffICU_95UL = 0.95
} else if(scenario == "shock"){
  mcdiffICU = 0.34
  mcdiffICU_95LL = -0.8
  mcdiffICU_95UL = 1.49
} else if(scenario == "sepsis3"){
  mcdiffICU = 0.37
  mcdiffICU_95LL = -0.74
  mcdiffICU_95UL = 1.48
} else if(scenario == "SOFA06"){
  mcdiffICU = -0.92
  mcdiffICU_95LL = -2.52
  mcdiffICU_95UL = 0.69
} else if(scenario == "SOFA1024"){
  mcdiffICU = 2.08
  mcdiffICU_95LL = -0.08
  mcdiffICU_95UL = 4.24
}
#Difference in ICU length of stay (days), margin of error and estimated standard deviation
mcdiffICU_MoE <- (mcdiffICU_95UL- mcdiffICU_95LL)/2
mcdiffICU_SD <- mcdiffICU_MoE/qnorm(0.975)
mcdiffICU_VAR <- mcdiffICU_SD^2

# Number of iterations for the simulation
n <- 10000
results <- data.frame(matrix(ncol = 20, nrow = n))
names(results) <- c("Mortality", "Procalcitonin mortality", "Length of stay", "Procalcitonin length of stay", "Length of ICU stay", "Procalcitonin length of ICU stay", "Length of regular ward stay", "Procalcitonin length of regular ward stay", "Standard care costs", "Procalcitonin costs", "Standard care QALYs", "Procalcitonin QALYs", "Standard care antibiotic duration", "Procalcitonin antibiotic duration", "Incremental costs", "Incremental QALYs", "Incremental antibiotics", "ICER", "Cost per antibiotic day avoided", "Antibiotic value")

if(scenario == "basecase" | scenario == "nodiffabx" | scenario == "nodiffLOS" | scenario == "nodiffmort" | scenario == "UKsepsismort" | scenario == "PCT50" | scenario == "shock" | scenario == "renal" | scenario == "resp" | scenario == "sepsis3" | scenario == "SOFA06" | scenario == "SOFA1024"){
  for (i in 1:n) {
    Mort <- rbeta(1, Mort_alpha, Mort_beta)
    
    MortOR <- rlnorm(1, log(mcMortOR), (log(mcMortOR_95UL)-log(mcMortOR_95LL))/(2*qnorm(0.975)))
    
    dAbx <- rgamma(1, dAbx_SHAPE, dAbx_RATE)
    
    diffAbx <- rnorm(1, mcdiffAbx, mcdiffAbx_SD)
    
    LOS <- rgamma(1, LOS_SHAPE, LOS_RATE)
    
    diffLOS <- rnorm(1, mcdiffLOS, mcdiffLOS_SD)
    
    ICU <- rgamma(1, ICU_SHAPE, ICU_RATE)
    
    diffICU <- rnorm(1, mcdiffICU, mcdiffICU_SD)
    
    SepPF <- rbeta(1, SepPF_alpha, SepPF_beta)
    
    uF <- rbeta (1, uF_alpha, uF_beta)
    
    uM <- rbeta (1, uM_alpha, uM_beta)
    
    uSep <- rbeta (1, uSep_alpha, uSep_beta)
    
    duHosD <- rbeta(1, duHosD_alpha, duHosD_beta)
    
    duAbxD <- rbeta (1, duAbxD_alpha, duAbxD_beta)
    
    iteration_params <- c(Mort = Mort, MortOR = MortOR, dAbx = dAbx, diffAbx = diffAbx, LOS = LOS, diffLOS = diffLOS, ICU = ICU, diffICU= diffICU, SepPF = SepPF, uF = uF, uM = uM, uSep = uSep)
    
    iteration_results <- dec_tree(iteration_params)
    
    results[i,] <- iteration_results
  }} else if(scenario == "noduLOS"){
    for (i in 1:n) {
      Mort <- rbeta(1, Mort_alpha, Mort_beta)
      
      MortOR <- rlnorm(1, log(mcMortOR), (log(mcMortOR_95UL)-log(mcMortOR_95LL))/(2*qnorm(0.975)))
      
      dAbx <- rgamma(1, dAbx_SHAPE, dAbx_RATE)
      
      diffAbx <- rnorm(1, mcdiffAbx, mcdiffAbx_SD)
      
      LOS <- rgamma(1, LOS_SHAPE, LOS_RATE)
      
      diffLOS <- rnorm(1, mcdiffLOS, mcdiffLOS_SD)
      
      ICU <- rgamma(1, ICU_SHAPE, ICU_RATE)
      
      diffICU <- rnorm(1, mcdiffICU, mcdiffICU_SD)
      
      SepPF <- rbeta(1, SepPF_alpha, SepPF_beta)
      
      uF <- rbeta (1, uF_alpha, uF_beta)
      
      uM <- rbeta (1, uM_alpha, uM_beta)
      
      uSep <- rbeta (1, uSep_alpha, uSep_beta)
      
      duAbxD <- rbeta (1, duAbxD_alpha, duAbxD_beta)
      
      iteration_params <- c(Mort = Mort, MortOR = MortOR, dAbx = dAbx, diffAbx = diffAbx, LOS = LOS, diffLOS = diffLOS, ICU = ICU, diffICU= diffICU, SepPF = SepPF, uF = uF, uM = uM, uSep = uSep)
      
      iteration_results <- dec_tree(iteration_params)
      
      results[i,] <- iteration_results
    }} else if(scenario == "noduabx"){
      for (i in 1:n) {
        Mort <- rbeta(1, Mort_alpha, Mort_beta)
        
        MortOR <- rlnorm(1, log(mcMortOR), (log(mcMortOR_95UL)-log(mcMortOR_95LL))/(2*qnorm(0.975)))
        
        dAbx <- rgamma(1, dAbx_SHAPE, dAbx_RATE)
        
        diffAbx <- rnorm(1, mcdiffAbx, mcdiffAbx_SD)
        
        LOS <- rgamma(1, LOS_SHAPE, LOS_RATE)
        
        diffLOS <- rnorm(1, mcdiffLOS, mcdiffLOS_SD)
        
        ICU <- rgamma(1, ICU_SHAPE, ICU_RATE)
        
        diffICU <- rnorm(1, mcdiffICU, mcdiffICU_SD)
        
        SepPF <- rbeta(1, SepPF_alpha, SepPF_beta)
        
        uF <- rbeta (1, uF_alpha, uF_beta)
        
        uM <- rbeta (1, uM_alpha, uM_beta)
        
        uSep <- rbeta (1, uSep_alpha, uSep_beta)
        
        duHosD <- rbeta(1, duHosD_alpha, duHosD_beta)
        
        iteration_params <- c(Mort = Mort, MortOR = MortOR, dAbx = dAbx, diffAbx = diffAbx, LOS = LOS, diffLOS = diffLOS, ICU = ICU, diffICU= diffICU, SepPF = SepPF, uF = uF, uM = uM, uSep = uSep)
        
        iteration_results <- dec_tree(iteration_params)
        
        results[i,] <- iteration_results
      }}

# Calculate summary statistics
results_summary <- summary(results)

# View summary statistics
print(results_summary)

# For each column in the results dataframe
for (col_name in names(results)) {
  # Calculate the 95% confidence interval
  conf_interval <- quantile(results[[col_name]], probs = c(0.025, 0.975), na.rm = TRUE)
  
  # Print the confidence interval
  print(paste("95% confidence interval for", col_name, ":", conf_interval))
}

#Exporting results

results_with_ci <- data.frame(Column = character(),
                              Summary = character(),
                              stringsAsFactors = FALSE)
results_summary <- summary(results)
for (col_name in names(results)) {
  conf_interval <- quantile(results[[col_name]], probs = c(0.025, 0.975), na.rm = TRUE)
  result_value <- mean(results[[col_name]], na.rm = TRUE)
  if (grepl("cost", col_name, ignore.case = TRUE) || grepl("icer", col_name, ignore.case = TRUE)) {
    result_with_ci <- paste(round(result_value), " (", 
                            round(conf_interval[1]), " to ", 
                            round(conf_interval[2]), ")", sep = "")
  } else {
    result_with_ci <- paste(round(result_value, 3), " (", 
                            round(conf_interval[1], 3), " to ", 
                            round(conf_interval[2], 3), ")", sep = "")
  }
  results_with_ci <- rbind(results_with_ci, data.frame(Column = col_name, 
                                                       Summary = result_with_ci))
}
write.csv(results_with_ci, "results_with_confidence_intervals.csv", row.names = FALSE)


#Section 6 - Plotting results

#ICER Plot
mean_qalys <- mean(results$`Incremental QALYs`, na.rm = TRUE)
mean_costs <- mean(results$`Incremental costs`, na.rm = TRUE)

ggplot(data = results, aes(x = `Incremental QALYs`, y = `Incremental costs`)) +
  geom_point(color = 'grey', alpha = 0.5) +
  geom_point(aes(x = mean_qalys, y = mean_costs), color = 'red', size = 1.25) +
  geom_hline(yintercept = 0, color = "black") +
  geom_vline(xintercept = 0, color = "black") +
  geom_abline(slope = 20000, intercept = 0, linetype = "dashed", color = "blue") + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white")
  ) +
  scale_x_continuous(breaks = seq(-0.2, max(results$`Incremental QALYs`, na.rm = TRUE), 0.2)) +
  scale_y_continuous(breaks = seq(-3000, max(results$`Incremental costs`, na.rm = TRUE), 1000)) +
  ylab("Incremental costs (£GBP)") +
  xlab("Incremental QALYs")

ICUicer <- results$'ICER'
ICUicer <- ICUicer[!is.na(ICUicer)]

#Cost per antibiotic day avoided plot

ggplot(data = results, aes(x = `Incremental antibiotics`, y = `Incremental costs`)) +
  geom_point(color = 'grey', alpha = 0.5) 

#Emergency department code#

#Input values#
#Selecting which model to run - untag the desired analysis#

scenario = "basecase"

#Scenario analyses#
#scenario = "nodiffabx"
#scenario = "noduabx"
#scenario = "nodiffLOS"
#scenario = "noduLOS"
#scenario = "nodiffmort"
#scenario = "PCT50"

#Subgroup analyses#
#scenario = "bronchitis"
#scenario = "CAP"
#scenario = "COPD"


#SECTION 1 - Defining parameters#

#Base case standard care sample size
SC_n <- 1638
if(scenario == "bronchitis"){
  SC_n = 287
} else if(scenario == "CAP"){
  SC_n = 1468
} else if(scenario == "COPD"){
  SC_n = 631
}

# Base case mortality risk
Mort <- 0.038
if(scenario == "bronchitis"){
  Mort = 0
} else if(scenario == "CAP"){
  Mort = 0.14
} else if(scenario == "COPD"){
  Mort = 0.038
}
#Odds
MortOdds <- Mort / (1 - Mort)
#Standard error
Mort_SE <- sqrt(Mort * (1 - Mort) / SC_n)
#95% confidence interval
Mort_95LL <- Mort - qnorm(0.975) * Mort_SE
Mort_95UL <- Mort + qnorm(0.975) * Mort_SE
print(paste("95% Confidence Interval Lower Limit:", Mort_95LL))
print(paste("95% Confidence Interval Upper Limit:", Mort_95UL))
Mort_var <- (Mort-Mort_95LL)/1.96
mu <- Mort
var <- Mort_var*Mort_var

# Function to calculate alpha and beta for a Beta distribution
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(list(alpha = alpha, beta = beta))
}
Mort_betaEstimates <- estBetaParams(mu, var)
print(Mort_betaEstimates)
Mort_alpha <- Mort_betaEstimates$alpha
Mort_beta <- Mort_betaEstimates$beta
#Checking the alpha and beta parameters correspond with 95% confidence intervals
qbeta(c(0.025, 0.975), Mort_betaEstimates$alpha, Mort_betaEstimates$beta)

#Base case mortality adjusted odds ratio
MortOR <- 0.91
#95% confidence interval
MortOR_95LL <- 0.63
MortOR_95UL <- 1.33
if(scenario == "nodiffmort"| scenario == "bronchitis"){
  MortOR = 1
  MortOR_95LL = 1
  MortOR_95UL = 1
} else if(scenario == "CAP"){
  MortOR = 0.82
  MortOR_95LL = 0.66
  MortOR_95UL = 1.03
} else if(scenario == "COPD"){
  MortOR = 0.8
  MortOR_95LL = 0.43
  MortOR_95UL = 1.48
}
#Standard error
MortOR_SE <- (MortOR_95UL-MortOR)/qnorm(0.975)
#Mortality adjusted OR, margin of error and estimated standard deviation
MortOR_MoE <- (MortOR_95UL- MortOR_95LL)/2
MortOR_SD <- MortOR_MoE/qnorm(0.975)
MortOR_VAR <- MortOR_SD^2

#PCT mortality odds
PCTMortOdds <- MortOdds * MortOR
#PCT mortality risk
PCTMort <- PCTMortOdds/(PCTMortOdds+1)
#Standard care survival risk
Surv <- 1-Mort
#PCT survival risk
PCTSurv <- 1-(PCTMort)

#Initiation of antibiotic risk 
iAbx <- 0.832
if(scenario == "bronchitis"){
  iAbx = 0.66
} else if(scenario == "CAP"){
  iAbx = 0.99
} else if(scenario == "COPD"){
  iAbx = 0.72
}
#Odds
iAbxOdds <- iAbx / (1 - iAbx)
#Standard error
iAbx_SE <- sqrt(iAbx * (1 - iAbx) / SC_n)
#95% confidence interval
iAbx_95LL <- iAbx - qnorm(0.975) * iAbx_SE
iAbx_95UL <- iAbx + qnorm(0.975) * iAbx_SE
print(paste("95% Confidence Interval Lower Limit:", iAbx_95LL))
print(paste("95% Confidence Interval Upper Limit:", iAbx_95UL))
iAbx_var <- (iAbx-iAbx_95LL)/1.96
mu <- iAbx
var <- iAbx_var*iAbx_var

# Function to calculate alpha and beta for a Beta distribution
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(list(alpha = alpha, beta = beta))
}
iAbx_betaEstimates <- estBetaParams(mu, var)
print(iAbx_betaEstimates)
iAbx_alpha <- iAbx_betaEstimates$alpha
iAbx_beta <- iAbx_betaEstimates$beta
#Checking the alpha and beta parameters correspond with 95% confidence intervals
qbeta(c(0.025, 0.975), iAbx_betaEstimates$alpha, iAbx_betaEstimates$beta)

#Antibiotic initiation adjusted odds ratio
iAbxOR <- 0.49
iAbxOR_95LL <- 0.41
iAbxOR_95UL <- 0.58
if(scenario == "nodiffabx"){
  iAbxOR <- 1
  iAbxOR_95LL <- 1
  iAbxOR_95UL <- 1
} else if(scenario == "bronchitis"){
  iAbxOR <- 0.18
  iAbxOR_95LL <- 0.12
  iAbxOR_95UL <- 0.26
} else if(scenario == "CAP"){
  iAbxOR <- 0.49
  iAbxOR_95LL <- 0.41
  iAbxOR_95UL <- 0.58
} else if(scenario == "COPD"){
  iAbxOR <- 0.29
  iAbxOR_95LL <- 0.23
  iAbxOR_95UL <- 0.36
}
iAbxOR_SE <- (iAbxOR_95UL-iAbxOR)/qnorm(0.975)
iAbxOR_MoE <- (iAbxOR_95UL- iAbxOR_95LL)/2
iAbxOR_SD <- iAbxOR_MoE/qnorm(0.975)
iAbxOR_VAR <- iAbxOR_SD^2

PCTiAbxOdds <- iAbxOdds*iAbxOR
PCTiAbx <- PCTiAbxOdds/(PCTiAbxOdds+1)

#Antibiotic duration (days)
dAbx <- 9.8
dAbx_SD <- 5.4
if(scenario == "bronchitis"){
  dAbx <- 7.1
  dAbx_SD <- 3.0
} else if(scenario == "CAP"){
  dAbx <- 10.5
  dAbx_SD <- 6.2
} else if(scenario == "COPD"){
  dAbx <- 7.4
  dAbx_SD <- 5.3
}

#Antibiotic duration (days), standard error
dAbx_SE <- dAbx_SD/sqrt(SC_n)
#Antibiotic duration (days), variance (SE)
dAbx_VAR <- dAbx_SE**2
#Antibiotic duration (days), Gamma distribution parameters
dAbx_SCALE <- dAbx_VAR/dAbx
dAbx_SHAPE <- dAbx/dAbx_SCALE
dAbx_RATE <- 1/dAbx_SCALE
#Antibiotic duration (days), 95% confidence interval
dAbx_95LL <- dAbx - qnorm(0.975)*dAbx_SE
dAbx_95UL <- dAbx + qnorm(0.975)*dAbx_SE


#Difference in antibiotic duration (days)
diffAbx <- -2.45
#Difference in antibiotic duration (days), 95% confidence interval
diffAbx_95LL <- -2.86
diffAbx_95UL <- -2.08
if(scenario == "nodiffabx"){
  diffAbx = 0
  diffAbx_95LL = 0
  diffAbx_95UL = 0
} else if(scenario == "bronchitis"){
  diffAbx = -0.35
  diffAbx_95LL = -1.15
  diffAbx_95UL = 0.45
} else if(scenario == "CAP"){
  diffAbx = -2.45
  diffAbx_95LL = -2.86
  diffAbx_95UL = -2.05
} else if(scenario == "COPD"){
  diffAbx = -1.15
  diffAbx_95LL = -2
  diffAbx_95UL = -0.31
}

#Difference in antibiotic duration (days), margin of error and estimated standard deviation
diffAbx_MoE <- (diffAbx_95UL- diffAbx_95LL)/2
diffAbx_SD <- diffAbx_MoE/qnorm(0.975)
diffAbx_VAR <- diffAbx_SD^2

#PCT antibiotic duration (days)
PCTdAbx <- dAbx + diffAbx

#Total exposure of antibiotics
tAbx <- dAbx*iAbx
PCTtAbx <- PCTdAbx*PCTiAbx

#Total length of stay (days)
LOS <- 8.2
#Total length of stay (days), SD
LOS_SD <- 10.5
if(scenario == "bronchitis"){
  LOS <- 2.6
  LOS_SD <- 5.7
} else if(scenario == "CAP"){
  LOS <- 13.3
  LOS_SD <- 15.7
} else if(scenario == "COPD"){
  LOS <- 9.3
  LOS_SD <- 13.9
}
#Total length of stay (days), SE
LOS_SE <- LOS_SD/sqrt(SC_n)
#Total length of stay (days), variance
LOS_VAR <- LOS_SE**2
#Total length of stay (days), gamma distribution parameters
LOS_SCALE <- LOS_VAR/LOS
LOS_SHAPE <- LOS/LOS_SCALE
LOS_RATE <- 1/LOS_SCALE
#Total length of stay (days), 95% confidence interval
LOS_95LL <- LOS - qnorm(0.975)*LOS_SE
LOS_95UL <- LOS + qnorm(0.975)*LOS_SE

#Difference in total length of stay (days)
diffLOS <- -0.14
#Difference in total length of stay (days), 95% confidence interval
diffLOS_95LL <- -0.73
diffLOS_95UL <- 0.44
if(scenario == "nodiffLOS"){
  diffLOS = 0
  diffLOS_95LL = 0
  diffLOS_95UL = 0
} else if(scenario == "bronchitis"){
  diffLOS = -0.21
  diffLOS_95LL = -0.9
  diffLOS_95UL = 0.48
} else if(scenario == "CAP"){
  diffLOS = 0.74
  diffLOS_95LL = -0.25
  diffLOS_95UL = 1.73
} else if(scenario == "COPD"){
  diffLOS = -0.6
  diffLOS_95LL = -1.84
  diffLOS_95UL = 0.64
}
#Difference in total length of stay (days), margin of error and estimated standard deviation
diffLOS_MoE <- (diffLOS_95UL- diffLOS_95LL)/2
diffLOS_SD <- diffLOS_MoE/qnorm(0.975)
diffLOS_VAR <- diffLOS_SD^2
#PCT total length of stay (days)
PCTLOS <- LOS + diffLOS

#Hospital length of stay (days)
Hos <- LOS 
# PCT hospital length of stay (days)
PCTHos <- PCTLOS

#ED cost (admission)
cED <- 349.57
#Hospital cost (day)
cHos <- 767.31
#Antibiotic cost (day)
cAbx <- 0.29
#PCT cost
cPCT <- 18.23
if(scenario == "PCT50"){
  cPCT <- 50
}
#Number of PCT tests
nPCT <- 4

#Proportion of female LRTI patients
SepPF <- 0.49
SepPF_n <- 8010
if(scenario == "bronchitis"| scenario == "COPD"){
  SepPF <- 0.49
  SepPF_n <- 3005
} else if(scenario == "CAP"){
  SepPF <- 0.45
  SepPF_n <- 2402
} 
#Standard error
SepPF_SE <- sqrt(SepPF*(1-SepPF)/SepPF_n)
#95% confidence interval
SepPF_95LL <- SepPF - qnorm(0.975)*SepPF_SE
SepPF_95UL <- SepPF + qnorm(0.975)*SepPF_SE
print(paste("95% Confidence Interval Lower Limit:", SepPF_95LL))
print(paste("95% Confidence Interval Upper Limit:", SepPF_95UL))
SepPF_var <- ((SepPF-SepPF_95LL)/1.96)
mu <- SepPF
var <- SepPF_var*SepPF_var
# Function to calculate alpha and beta for a Beta distribution
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(list(alpha = alpha, beta = beta))
}
SepPF_betaEstimates <- estBetaParams(mu, var)
print(SepPF_betaEstimates)
SepPF_alpha <- SepPF_betaEstimates$alpha
SepPF_beta <- SepPF_betaEstimates$beta
#Checking the alpha and beta parameters correspond with 95% confidence intervals
qbeta(c(0.025, 0.975), SepPF_betaEstimates$alpha, SepPF_betaEstimates$beta)
#Proportion of male patients
SepPM <- 1-SepPF

#Female QALY 
uF <- 0.784
#Female QALY, sample size
uF_n <- 619
#Female QALY 95% confidence interval
uF_95LL <- 0.779
uF_95UL <- 0.801
#Female QALY standard error
uF_SE <- sqrt(uF*(1-uF)/uF_n)
if(scenario == "bronchitis" | scenario == "COPD"){
  uF = 0.776
  uF_n = 608
  uF_95LL = 0.769
  uF_95UL = 0.797
}
uF_var <- (uF-uF_95LL)/1.96
mu <- uF
var <- uF_var*uF_var

# Function to calculate alpha and beta for a Beta distribution
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(list(alpha = alpha, beta = beta))
}
uF_betaEstimates <- estBetaParams(mu, var)
print(uF_betaEstimates)
uF_alpha <- uF_betaEstimates$alpha
uF_beta <- uF_betaEstimates$beta
#Checking the alpha and beta parameters correspond with 95% confidence intervals
qbeta(c(0.025, 0.975), uF_betaEstimates$alpha, uF_betaEstimates$beta)


#Male QALY 
uM <- 0.801
#Male QALY, sample size
uM_n <- 505
#Male QALY 60-64 years old, 95% confidence interval
uM_95LL <- 0.794
uM_95UL <- 0.818
if(scenario == "bronchitis" | scenario == "COPD"){
  uM = 0.802
  uM_n = 532
  uM_95LL = 0.798
  uM_95UL = 0.822
}
#Standard care mortality rate, standard error
uM_SE <- sqrt(uM*(1-uM)/uM_n)
uM_var <- (uM-uM_95LL)/1.96
mu <- uM
var <- uM_var*uM_var

# Function to calculate alpha and beta for a Beta distribution
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(list(alpha = alpha, beta = beta))
}
uM_betaEstimates <- estBetaParams(mu, var)
print(uM_betaEstimates)
uM_alpha <- uM_betaEstimates$alpha
uM_beta <- uM_betaEstimates$beta
#Checking the alpha and beta parameters correspond with 95% confidence intervals
qbeta(c(0.025, 0.975), uM_betaEstimates$alpha, uM_betaEstimates$beta)

#Weighted population QALY 
uW <- (uF*SepPF) + (uM*SepPM)

#LRTI utility
uSep <- 0.705
#LRTI utility, sample size
uSep_n <- 349
#LRTI utility, SE
uSep_SE <- sqrt(uSep*(1-uSep)/uSep_n)
#LRTI utility, 95% confidence interval
uSep_95LL <- uSep - qnorm(0.975)*uSep_SE
uSep_95UL <- uSep + qnorm(0.975)*uSep_SE
print(paste("95% Confidence Interval Lower Limit:", uSep_95LL))
print(paste("95% Confidence Interval Upper Limit:", uSep_95UL))
uSep_var <- (uSep-uSep_95LL)/1.96
mu <- uSep
var <- uSep_var*uSep_var
# Function to calculate alpha and beta for a Beta distribution
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(list(alpha = alpha, beta = beta))
}
uSep_betaEstimates <- estBetaParams(mu, var)
print(uSep_betaEstimates)
uSep_alpha <- uSep_betaEstimates$alpha
uSep_beta <- uSep_betaEstimates$beta
#Checking the alpha and beta parameters correspond with 95% confidence intervals
qbeta(c(0.025, 0.975), uSep_betaEstimates$alpha, uSep_betaEstimates$beta)
#LRTI disutility
duSep <- (uW-uSep)/24

#Hospital disutility (QALD/day)
duHosD <- 0.1
if(scenario == "noduLOS"){
  duHosD = 0
}
duHosD_n <- 159
duHosD_SE <- sqrt(duHosD*(1-duHosD)/duHosD_n)
duHosD_95LL <- duHosD - qnorm(0.975)*duHosD_SE
duHosD_95UL <- duHosD + qnorm(0.975)*duHosD_SE
print(paste("95% Confidence Interval Lower Limit:", duHosD_95LL))
print(paste("95% Confidence Interval Upper Limit:", duHosD_95UL))
duHosD_var <- (duHosD-duHosD_95LL)/1.96
mu <- duHosD
var <- duHosD_var*duHosD_var

# Function to calculate alpha and beta for a Beta distribution
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(list(alpha = alpha, beta = beta))
}
duHosD_betaEstimates <- estBetaParams(mu, var)
print(duHosD_betaEstimates)
duHosD_alpha <- duHosD_betaEstimates$alpha
duHosD_beta <- duHosD_betaEstimates$beta
#Checking the alpha and beta parameters correspond with 95% confidence intervals
qbeta(c(0.025, 0.975), duHosD_betaEstimates$alpha, duHosD_betaEstimates$beta)

#Antibiotic disutility (QALD/day)
duAbxD <- 0.057
if(scenario == "noduabx"){
  duAbxD = 0
}
duAbxD_n <- 349
duAbxD_SE <- sqrt(duAbxD*(1-duAbxD)/duAbxD_n)
duAbxD_95LL <- duAbxD - qnorm(0.975)*duAbxD_SE
duAbxD_95UL <- duAbxD + qnorm(0.975)*duAbxD_SE
print(paste("95% Confidence Interval Lower Limit:", duAbxD))
print(paste("95% Confidence Interval Upper Limit:", duAbxD))

duAbxD_var <- (duAbxD-duAbxD_95LL)/1.96
mu <- duAbxD
var <- duAbxD_var*duAbxD_var

# Function to calculate alpha and beta for a Beta distribution
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(list(alpha = alpha, beta = beta))
}
duAbxD_betaEstimates <- estBetaParams(mu, var)
print(duAbxD_betaEstimates)
duAbxD_alpha <- duAbxD_betaEstimates$alpha
duAbxD_beta <- duAbxD_betaEstimates$beta
#Checking the alpha and beta parameters correspond with 95% confidence intervals
qbeta(c(0.025, 0.975), duAbxD_betaEstimates$alpha, duAbxD_betaEstimates$beta)


#Total hospital disutility (QALY)
duHos <- duHosD*LOS/365.25
#Total antibiotic disutility (QALY)
duAbx <- duAbxD*tAbx/365.25
#PCT total hospital disutility (QALY)
PCTduHos <- duHosD*PCTLOS/365.25
#PCT total antibiotic disutility (QALY)
PCTduAbx <- duAbxD*PCTtAbx/365.25

#Female QALE, 3.5% discount
QALEF <- 9.25
#Male QALE, 3.5% discount
QALEM <- 8.87
if(scenario == "bronchitis"| scenario == "COPD"){
  QALEF = 10.66
  QALEM = 10.23
} else if(scenario == "CAP"){
  QALEF = 8.88
  QALEM = 8.51}
#Weighted QALE, 3.5% discount
QALEW <- (QALEF*SepPF) + (QALEM*SepPM)

#Female LE 
LEF <- 17.17
#Male LE 
LEM <- 15.12
if(scenario == "bronchitis" | scenario == "COPD"){
  LEF = 20.46
  LEM = 18.13
} else if(scenario == "CAP"){
  LEF = 16.38
  LEM = 14.38
}
#Weighted LE 
LEW <- (LEF*SepPF) + (LEM*SepPM)

#SECTION 2 - Data frame for model inputs#

input <- data.frame(
  #1
  Mort <- Mort,
  MortOdds <- Mort/(1-Mort),
  MortOR <- MortOR,
  PCTMortOdds <- MortOdds * MortOR,
  PCTMort <- PCTMortOdds/(PCTMortOdds+1),
  #6
  Surv <- 1-Mort,
  PCTSurv <- 1-(PCTMort),
  iAbx <- iAbx,
  iAbxOdds <- iAbx/(1-iAbx),
  iAbxOR <- iAbxOR,
  #11
  PCTiAbxOdds <- iAbxOdds*iAbxOR,
  PCTiAbx <- PCTiAbxOdds/(PCTiAbxOdds+1),
  dAbx <- dAbx,
  diffAbx <- diffAbx,
  PCTdAbx <- dAbx + diffAbx,
  #16
  tAbx <- dAbx*iAbx,
  PCTtAbx <- PCTdAbx*PCTiAbx,
  LOS <- LOS,
  diffLOS <- diffLOS,
  PCTLOS <- LOS + diffLOS,
  #21
  cED <- cED,
  cHos <- cHos,
  cAbx <- cAbx,
  cPCT <- cPCT,
  nPCT <- nPCT,
  #26
  SepPF <- SepPF,
  SepPM <- 1-SepPF,
  uF <- uF,
  uM <- uM,
  uW <- (uF*SepPF) + (uM*SepPM),
  #31
  uSep <- uSep,
  duSep <- (uW-uSep)/24,
  duHos <- duHosD*LOS/365.25,
  duAbx <- duAbxD*tAbx/365.25,
  PCTduHos <- duHosD*PCTLOS/365.25,
  #36
  PCTduAbx <- duAbxD*PCTtAbx/365.25,
  QALEW <- (QALEF*SepPF) + (QALEM*SepPM),
  LEW <- (LEF*SepPF) + (LEM*SepPM)
)

#Section 3 - Decision Tree 
dec_tree <- function(params){
  with(
    as.list(params), 
    {
      #Input parameters
      Mort <- Mort
      MortOdds <- Mort/(1-Mort)
      MortOR <- MortOR
      PCTMortOdds <- MortOdds * MortOR
      PCTMort <- PCTMortOdds/(PCTMortOdds+1)
      Surv <- 1-Mort
      PCTSurv <- 1-(PCTMort)
      iAbx <- iAbx
      iAbxOdds <- iAbx/(1-iAbx)
      iAbxOR <- iAbxOR
      PCTiAbxOdds <- iAbxOdds*iAbxOR
      PCTiAbx <- PCTiAbxOdds/(PCTiAbxOdds+1)
      dAbx <- dAbx
      diffAbx <- diffAbx
      PCTdAbx <- dAbx + diffAbx
      tAbx <- dAbx*iAbx
      PCTtAbx <- PCTdAbx*PCTiAbx
      LOS <- LOS
      diffLOS <- diffLOS
      PCTLOS <- LOS + diffLOS
      cED <- cED
      cHos <- cHos
      cAbx <- cAbx
      cPCT <- cPCT
      nPCT <- nPCT
      SepPF <- SepPF
      SepPM <- 1-SepPF
      uF <- uF
      uM <- uM
      uW <- (uF*SepPF) + (uM*SepPM)
      uSep <- uSep
      duSep <- (uW-uSep)/24
      duHos <- duHosD*LOS/365.25
      duAbx <- duAbxD*tAbx/365.25
      PCTduHos <- duHosD*PCTLOS/365.25
      PCTduAbx <- duAbxD*PCTtAbx/365.25
      uW <- (uF*SepPF) + (uM*SepPM)
      LEW <- (LEF*SepPF) + (LEM*SepPM)
      
      #Standard care costs
      C_SC <- (cAbx*tAbx) + (cHos*LOS) + cED
      #PCT costs
      C_PCT <- (cPCT*nPCT) + (cAbx*PCTtAbx) + (cHos*PCTLOS) + cED
      #Standard care QALYs gained
      Q_SC <- (QALEW-duSep-duAbx-duHos)*Surv
      #PCT QALYs gained
      Q_PCT <- (QALEW-duSep-PCTduAbx-PCTduHos)*PCTSurv
      #Incremental Costs
      IC <- C_PCT - C_SC
      #Incremental QALYs
      IQALYx <- Q_PCT - Q_SC
      IQALY <- round(IQALYx, 4)
      
      if (IQALY == 0) {
        ICER <- NA
        ICER_flag <- "Excluded: IQALY of 0"
      } else {
        # Incremental Costs
        IC <- C_PCT - C_SC
        # ICER
        ICER <- IC / IQALY
        # Normal flag for non-excluded cases
        ICER_flag <- "Normal"
      }
      #Incremental Antibiotics
      IAbx <- PCTtAbx - tAbx
      #Cost per antibiotic day avoided
      AbxAv <- (C_PCT - C_SC)/(tAbx-PCTtAbx)
      #Antibiotic value
      AbxV <- (IC-(IQALY*20000))/IAbx
      Inputs <- c(Mort, PCTMort, iAbx, PCTiAbx, LOS, PCTLOS)
      BasicOutcomes <- c(C_SC, C_PCT, Q_SC, Q_PCT, tAbx, PCTtAbx)
      ICERs <- c(IC, IQALY, IAbx, ICER, AbxAv, AbxV)
      
      names(Inputs) <- paste ("Inputs", c("Mortality", "Procalcitonin mortality", "Standard care antibiotic initiation", "Procalcitonin antibiotic intiation", "Length of stay", "Procalcitonin length of stay"))
      names(BasicOutcomes) <- paste ("Basic Outcomes", c("Standard care costs", "Procalcitonin costs", "Standard care QALYs", "Procalcitonin QALYs", "Standard care antibiotic duration", "Procalcitonin antibiotic duration"), sep = "_")
      names(ICERs)  <- paste("ICERs", c("Incremental costs", "Incremental QALYs", "Incremental antibiotics", "ICER", "Cost per antibiotic day avoided", "Antibiotic value"), sep = "_")
      
      return(c(Inputs, BasicOutcomes, ICERs))
    }
  )
}
options(scipen=999)
dec_tree(input)

#Section 4 - Tornado Plot#
if (scenario == "basecase"){
  Mort_range <- c(BaseCase = Mort, low = Mort_95LL, high = Mort_95UL)
  MortOR_range <- c(BaseCase = MortOR, low = MortOR_95LL, high = MortOR_95UL)
  iAbx_range <- c(BaseCase = iAbx, low = iAbx_95LL, high = iAbx_95UL)
  iAbxOR_range <- c(BaseCase = iAbxOR, low = iAbxOR_95LL, high = iAbxOR_95UL)
  dAbx_range <- c(BaseCase = dAbx, low = dAbx_95LL, high = dAbx_95UL)
  diffAbx_range <- c(BaseCase = diffAbx, low = diffAbx_95LL, high = diffAbx_95UL)
  LOS_range <- c(BaseCase = LOS, low = LOS_95LL, high = LOS_95UL)
  diffLOS_range <- c(BaseCase = diffLOS, low = diffLOS_95LL, high = diffLOS_95UL)
  cED_range <- c(BaseCase = cED, low = (cED*0.75), high = (cED*1.25))
  cHos_range <- c(BaseCase = cHos, low = (cHos*0.75), high = (cHos*1.25))
  cAbx_range <- c(BaseCase = cAbx, low = (cAbx*0.75), high = (cAbx*1.25))
  cPCT_range <- c(BaseCase = cPCT, low = 14.5, high = 90)
  nPCT_range <- c(BaseCase = nPCT, low = 1, high = 8)
  SepPF_range <- c(BaseCase = SepPF, low = SepPF_95LL, high = SepPF_95UL)
  uF_range <- c(BaseCase = uF, low = uF_95LL, high = uF_95UL)
  uM_range <- c(BaseCase = uM, low = uM_95LL, high = uM_95UL)
  uSep_range <- c(BaseCase = uM, low = uSep_95LL, high = uSep_95UL)
  
  #Plotting all parameters to find most influential
  paramNames <- c(
    "30-day mortality probability",
    "30-day mortality odds ratio",
    "Antibiotic initiation probability",
    "Antibiotic initiation odds ratio",
    "Duration of antibiotics",
    "Antibiotic duration difference",
    "Hospital length of stay",
    "Hospital length of stay difference",
    "Admission cost",
    "Hospital cost",
    "Antibiotic cost",
    "Procalcitonin cost per test",
    "Number of procalcitonin tests used",
    "Proportion of female patients",
    "Female utility",
    "Male utility",
    "LRTI disutility"
  )
  
  l.tor.in <- vector("list", 17)
  names(l.tor.in) <- paramNames
  
  l.tor.in$'30-day mortality probability' <- cbind(Mort = Mort_range, input [-1])
  l.tor.in$'30-day mortality odds ratio' <- cbind(MortOR = MortOR_range, input [-3])
  l.tor.in$'Antibiotic initiation probability' <- cbind(iAbx = iAbx_range, input [-8])
  l.tor.in$'Antibiotic initiation odds ratio' <- cbind(iAbx = iAbx_range, input [-10])
  l.tor.in$'Duration of antibiotics' <- cbind(dAbx = dAbx_range, input [-13])
  l.tor.in$'Antibiotic duration difference' <- cbind(diffAbx = diffAbx_range, input [-14])
  l.tor.in$'Hospital length of stay' <-cbind(LOS = LOS_range, input [-18])
  l.tor.in$'Hospital length of stay difference' <-cbind(diffLOS = diffLOS_range, input [-19])
  l.tor.in$'Admission cost' <-cbind(cED = cED_range, input [-21])
  l.tor.in$'Hospital cost' <-cbind(cHos = cHos_range, input [-22])
  l.tor.in$'Antibiotic cost' <-cbind(cAbx = cAbx_range, input [-23])
  l.tor.in$'Procalcitonin cost per test' <-cbind(cPCT = cPCT_range, input [-24])
  l.tor.in$'Number of procalcitonin tests used' <-cbind(nPCT = nPCT_range, input [-25])
  l.tor.in$'Proportion of female patients' <-cbind(SepPF = SepPF_range, input [-26])
  l.tor.in$'Female utility' <-cbind(uF = uF_range, input [-28])
  l.tor.in$'Male utility' <-cbind(uM = uM_range, input [-29])
  l.tor.in$'LRTI disutility' <-cbind(uSep = uSep_range, input [-31])
  
  #List of outputs
  l.tor.out <- vector("list", 17)
  names(l.tor.out) <- paramNames
  
  for(i in 1:17){
    l.tor.out[[i]] <- t(apply(l.tor.in[[i]], 1, dec_tree))[ , 18]
  }
  
  m.tor <- matrix(unlist(l.tor.out), nrow = 17, ncol = 3, byrow = TRUE, 
                  dimnames = list(paramNames, c("basecase", "low", "high")))
  
  TornadoPlot(main_title = "", Parms = paramNames, Outcomes = m.tor, 
              outcomeName = "", 
              xlab = "", 
              ylab = "", 
              col1="#3182bd", col2="#6baed6"
  )
  
  #Plotting Top 8 parameters
  paramNames <- c(
    "Difference in total length of stay",
    "Cost per PCT test",
    "Number of PCT tests used",
    "Cost per regular ward day",
    "30-day mortality odds ratio",
    "30-day mortality probability odds ratio",
    "Cost per day of antibiotic therapy",
    "Difference in duration of antibiotic therapy"
  )
  
  
  l.tor.in <- vector("list", 8)
  names(l.tor.in) <- paramNames
  
  
  l.tor.in$'Difference in total length of stay' <- cbind(diffLOS = diffLOS_range, input [-19])
  l.tor.in$'Cost per PCT test' <- cbind(cPCT = cPCT_range, input [-24])
  l.tor.in$'Number of PCT tests used' <- cbind(nPCT = nPCT_range, input [-25])
  l.tor.in$'Cost per regular ward day' <- cbind(cHos = cHos_range, input [-22])
  l.tor.in$'30-day mortality odds ratio' <- cbind(MortOR = MortOR_range, input [-3])
  l.tor.in$'30-day mortality probability odds ratio' <- cbind(Mort = Mort_range, input [-1])
  l.tor.in$'Cost per day of antibiotic therapy' <- cbind(cAbx = cAbx_range, input [-23])
  l.tor.in$'Difference in duration of antibiotic therapy' <- cbind(diffAbx = diffAbx_range, input [-14])
  
  #List of outputs
  l.tor.out <- vector("list", 8)
  names(l.tor.out) <- paramNames
  
  for(i in 1:8){
    l.tor.out[[i]] <- t(apply(l.tor.in[[i]], 1, dec_tree))[ , 18]
  }
  
  m.tor <- matrix(unlist(l.tor.out), nrow = 8, ncol = 3, byrow = TRUE, 
                  dimnames = list(paramNames, c("basecase", "low", "high")))
  
  TornadoPlot(main_title = "", Parms = paramNames, Outcomes = m.tor, 
              outcomeName = "", 
              xlab = "", 
              ylab = "", 
              col1="#3182bd", col2="#6baed6"
  )}

#Section 5 - Monte Carlo Simulation

#Defining parallel parameters to prevent simulation iterating on itself

#Base case mortality adjusted odds ratio
mcMortOR <- 0.91
#95% confidence interval
mcMortOR_95LL <- 0.63
mcMortOR_95UL <- 1.33
if(scenario == "nodiffmort"| scenario == "bronchitis"){
  mcMortOR = 1
  mcMortOR_95LL = 1
  mcMortOR_95UL = 1
} else if(scenario == "CAP"){
  mcMortOR = 0.82
  mcMortOR_95LL = 0.66
  mcMortOR_95UL = 1.03
} else if(scenario == "COPD"){
  mcMortOR = 0.8
  mcMortOR_95LL = 0.43
  mcMortOR_95UL = 1.48
}
#Standard error
mcMortOR_SE <- (mcMortOR_95UL-mcMortOR)/qnorm(0.975)
#Mortality adjusted OR, margin of error and estimated standard deviation
mcMortOR_MoE <- (mcMortOR_95UL- mcMortOR_95LL)/2
mcMortOR_SD <- mcMortOR_MoE/qnorm(0.975)
mcMortOR_VAR <- mcMortOR_SD^2

#Antibiotic initiation adjusted odds ratio
mciAbxOR <- 0.49
mciAbxOR_95LL <- 0.41
mciAbxOR_95UL <- 0.58
if(scenario == "nodiffabx"){
  mciAbxOR <- 1
  mciAbxOR_95LL <- 1
  mciAbxOR_95UL <- 1
} else if(scenario == "bronchitis"){
  mciAbxOR <- 0.18
  mciAbxOR_95LL <- 0.12
  mciAbxOR_95UL <- 0.26
} else if(scenario == "CAP"){
  mciAbxOR <- 0.49
  mciAbxOR_95LL <- 0.41
  mciAbxOR_95UL <- 0.58
} else if(scenario == "COPD"){
  mciAbxOR <- 0.29
  mciAbxOR_95LL <- 0.23
  mciAbxOR_95UL <- 0.36
}

mciAbxOR_SE <- (mciAbxOR_95UL-mciAbxOR)/qnorm(0.975)
mciAbxOR_MoE <- (mciAbxOR_95UL- mciAbxOR_95LL)/2
mciAbxOR_SD <- mciAbxOR_MoE/qnorm(0.975)
mciAbxOR_VAR <- mciAbxOR_SD^2

#Difference in antibiotic duration (days)
mcdiffAbx <- -2.45
#Difference in antibiotic duration (days), 95% confidence interval
mcdiffAbx_95LL <- -2.86
mcdiffAbx_95UL <- -2.08
if(scenario == "nodiffabx"){
  mcdiffAbx = 0
  mcdiffAbx_95LL = 0
  mcdiffAbx_95UL = 0
} else if(scenario == "bronchitis"){
  mcdiffAbx = -0.35
  mcdiffAbx_95LL = -1.15
  mcdiffAbx_95UL = 0.45
} else if(scenario == "CAP"){
  mcdiffAbx = -2.45
  mcdiffAbx_95LL = -2.86
  mcdiffAbx_95UL = -2.05
} else if(scenario == "COPD"){
  mcdiffAbx = -1.15
  mcdiffAbx_95LL = -2
  mcdiffAbx_95UL = -0.31
}

#Difference in antibiotic duration (days), margin of error and estimated standard deviation
mcdiffAbx_MoE <- (mcdiffAbx_95UL- mcdiffAbx_95LL)/2
mcdiffAbx_SD <- mcdiffAbx_MoE/qnorm(0.975)
mcdiffAbx_VAR <- mcdiffAbx_SD^2

#Difference in total length of stay (days)
mcdiffLOS <- -0.14
#Difference in total length of stay (days), 95% confidence interval
mcdiffLOS_95LL <- -0.73
mcdiffLOS_95UL <- 0.44
if(scenario == "nodiffLOS"){
  mcdiffLOS = 0
  mcdiffLOS_95LL = 0
  mcdiffLOS_95UL = 0
} else if(scenario == "bronchitis"){
  mcdiffLOS = -0.21
  mcdiffLOS_95LL = -0.9
  mcdiffLOS_95UL = 0.48
} else if(scenario == "CAP"){
  mcdiffLOS = 0.74
  mcdiffLOS_95LL = -0.25
  mcdiffLOS_95UL = 1.73
} else if(scenario == "COPD"){
  mcdiffLOS = -0.6
  mcdiffLOS_95LL = -1.84
  mcdiffLOS_95UL = 0.64
}
#Difference in total length of stay (days), margin of error and estimated standard deviation
mcdiffLOS_MoE <- (mcdiffLOS_95UL- mcdiffLOS_95LL)/2
mcdiffLOS_SD <- mcdiffLOS_MoE/qnorm(0.975)
mcdiffLOS_VAR <- mcdiffLOS_SD^2

# Number of iterations for the simulation
n <- 10000
results <- data.frame(matrix(ncol = 18, nrow = n))
names(results) <- c("Mortality", "Procalcitonin mortality", "Standard care antibiotic initiation", "Procalcitonin antibiotic intiation", "Length of stay", "Procalcitonin length of stay", "Standard care costs", "Procalcitonin costs", "Standard care QALYs", "Procalcitonin QALYs", "Standard care antibiotic duration", "Procalcitonin antibiotic duration", "Incremental costs", "Incremental QALYs", "Incremental antibiotics", "ICER", "Cost per antibiotic day avoided", "Antibiotic value")

if(scenario == "basecase" | scenario == "nodiffabx" | scenario == "nodiffLOS" | scenario == "nodiffmort" | scenario == "PCT50" | scenario == "CAP" | scenario == "COPD"){
  for (i in 1:n) {
    Mort <- rbeta(1, Mort_alpha, Mort_beta)
    
    MortOR <- rlnorm(1, log(mcMortOR), (log(mcMortOR_95UL)-log(mcMortOR_95LL))/(2*qnorm(0.975)))
    
    iAbx <- rbeta(1, iAbx_alpha, iAbx_beta)
    
    iAbxOR <- rlnorm(1, log(mciAbxOR), (log(mciAbxOR_95UL)-log(mciAbxOR_95LL))/(2*qnorm(0.975)))
    
    dAbx <- rgamma(1, dAbx_SHAPE, dAbx_RATE)
    
    diffAbx <- rnorm(1, mcdiffAbx, mcdiffAbx_SD)
    
    LOS <- rgamma(1, LOS_SHAPE, LOS_RATE)
    
    diffLOS <- rnorm(1, mcdiffLOS, mcdiffLOS_SD)
    
    SepPF <- rbeta(1, SepPF_alpha, SepPF_beta)
    
    uF <- rbeta (1, uF_alpha, uF_beta)
    
    uM <- rbeta (1, uM_alpha, uM_beta)
    
    uSep <- rbeta (1, uSep_alpha, uSep_beta)
    
    duHosD <- rbeta(1, duHosD_alpha, duHosD_beta)
    
    duAbxD <- rbeta (1, duAbxD_alpha, duAbxD_beta)
    
    iteration_params <- c(Mort = Mort, MortOR = MortOR, iAbx = iAbx, iAbxOR = iAbxOR, dAbx = dAbx, diffAbx = diffAbx, LOS = LOS, diffLOS = diffLOS, SepPF = SepPF, uF = uF, uM = uM, uSep = uSep, duHosD = duHosD, duAbxD = duAbxD)
    
    iteration_results <- dec_tree(iteration_params)
    
    results[i,] <- iteration_results
  }} else if(scenario == "noduabx"){
    for (i in 1:n) {
      Mort <- rbeta(1, Mort_alpha, Mort_beta)
      
      MortOR <- rlnorm(1, log(mcMortOR), (log(mcMortOR_95UL)-log(mcMortOR_95LL))/(2*qnorm(0.975)))
      
      iAbx <- rbeta(1, iAbx_alpha, iAbx_beta)
      
      iAbxOR <- rlnorm(1, log(mciAbxOR), (log(mciAbxOR_95UL)-log(mciAbxOR_95LL))/(2*qnorm(0.975)))
      
      dAbx <- rgamma(1, dAbx_SHAPE, dAbx_RATE)
      
      diffAbx <- rnorm(1, mcdiffAbx, mcdiffAbx_SD)
      
      LOS <- rgamma(1, LOS_SHAPE, LOS_RATE)
      
      diffLOS <- rnorm(1, mcdiffLOS, mcdiffLOS_SD)
      
      SepPF <- rbeta(1, SepPF_alpha, SepPF_beta)
      
      uF <- rbeta (1, uF_alpha, uF_beta)
      
      uM <- rbeta (1, uM_alpha, uM_beta)
      
      uSep <- rbeta (1, uSep_alpha, uSep_beta)
      
      duHosD <- rbeta(1, duHosD_alpha, duHosD_beta)
      
      iteration_params <- c(Mort = Mort, MortOR = MortOR, iAbx = iAbx, iAbxOR = iAbxOR, dAbx = dAbx, diffAbx = diffAbx, LOS = LOS, diffLOS = diffLOS, SepPF = SepPF, uF = uF, uM = uM, uSep = uSep, duHosD = duHosD)
      
      iteration_results <- dec_tree(iteration_params)
      
      results[i,] <- iteration_results
    }} else if(scenario == "noduLOS"){
      for (i in 1:n) {
        Mort <- rbeta(1, Mort_alpha, Mort_beta)
        
        MortOR <- rlnorm(1, log(mcMortOR), (log(mcMortOR_95UL)-log(mcMortOR_95LL))/(2*qnorm(0.975)))
        
        iAbx <- rbeta(1, iAbx_alpha, iAbx_beta)
        
        iAbxOR <- rlnorm(1, log(mciAbxOR), (log(mciAbxOR_95UL)-log(mciAbxOR_95LL))/(2*qnorm(0.975)))
        
        dAbx <- rgamma(1, dAbx_SHAPE, dAbx_RATE)
        
        diffAbx <- rnorm(1, mcdiffAbx, mcdiffAbx_SD)
        
        LOS <- rgamma(1, LOS_SHAPE, LOS_RATE)
        
        diffLOS <- rnorm(1, mcdiffLOS, mcdiffLOS_SD)
        
        SepPF <- rbeta(1, SepPF_alpha, SepPF_beta)
        
        uF <- rbeta (1, uF_alpha, uF_beta)
        
        uM <- rbeta (1, uM_alpha, uM_beta)
        
        uSep <- rbeta (1, uSep_alpha, uSep_beta)
        
        duAbxD <- rbeta (1, duAbxD_alpha, duAbxD_beta)
        
        iteration_params <- c(Mort = Mort, MortOR = MortOR, iAbx = iAbx, iAbxOR = iAbxOR, dAbx = dAbx, diffAbx = diffAbx, LOS = LOS, diffLOS = diffLOS, SepPF = SepPF, uF = uF, uM = uM, uSep = uSep, duAbxD = duAbxD)
        
        iteration_results <- dec_tree(iteration_params)
        
        results[i,] <- iteration_results
      }} else if(scenario == "bronchitis"){
        for (i in 1:n) {
          iAbx <- rbeta(1, iAbx_alpha, iAbx_beta)
          
          iAbxOR <- rlnorm(1, log(mciAbxOR), (log(mciAbxOR_95UL)-log(mciAbxOR_95LL))/(2*qnorm(0.975)))
          
          dAbx <- rgamma(1, dAbx_SHAPE, dAbx_RATE)
          
          diffAbx <- rnorm(1, mcdiffAbx, mcdiffAbx_SD)
          
          LOS <- rgamma(1, LOS_SHAPE, LOS_RATE)
          
          diffLOS <- rnorm(1, mcdiffLOS, mcdiffLOS_SD)
          
          SepPF <- rbeta(1, SepPF_alpha, SepPF_beta)
          
          uF <- rbeta (1, uF_alpha, uF_beta)
          
          uM <- rbeta (1, uM_alpha, uM_beta)
          
          uSep <- rbeta (1, uSep_alpha, uSep_beta)
          
          duHosD <- rbeta(1, duHosD_alpha, duHosD_beta)
          
          duAbxD <- rbeta (1, duAbxD_alpha, duAbxD_beta)
          
          iteration_params <- c(iAbx = iAbx, iAbxOR = iAbxOR, dAbx = dAbx, diffAbx = diffAbx, LOS = LOS, diffLOS = diffLOS, SepPF = SepPF, uF = uF, uM = uM, uSep = uSep, duHosD = duHosD, duAbxD = duAbxD)
          
          iteration_results <- dec_tree(iteration_params)
          
          results[i,] <- iteration_results
        }}



# Calculate summary statistics
results_summary <- summary(results)

# View summary statistics
print(results_summary)

# For each column in the results dataframe
for (col_name in names(results)) {
  # Calculate the 95% confidence interval
  conf_interval <- quantile(results[[col_name]], probs = c(0.025, 0.975), na.rm = TRUE)
  
  # Print the confidence interval
  print(paste("95% confidence interval for", col_name, ":", conf_interval))
}

#Exporting results

results_with_ci <- data.frame(Column = character(),
                              Summary = character(),
                              stringsAsFactors = FALSE)
results_summary <- summary(results)
for (col_name in names(results)) {
  conf_interval <- quantile(results[[col_name]], probs = c(0.025, 0.975), na.rm = TRUE)
  result_value <- mean(results[[col_name]], na.rm = TRUE)
  if (grepl("cost", col_name, ignore.case = TRUE) || grepl("icer", col_name, ignore.case = TRUE)) {
    result_with_ci <- paste(round(result_value), " (", 
                            round(conf_interval[1]), " to ", 
                            round(conf_interval[2]), ")", sep = "")
  } else if (grepl("value", col_name, ignore.case = TRUE))
  {
    result_with_ci <- paste(round(result_value, 2), " (", 
                            round(conf_interval[1], 2), " to ", 
                            round(conf_interval[2], 2), ")", sep = "")
  }
  else {
    result_with_ci <- paste(round(result_value, 3), " (", 
                            round(conf_interval[1], 3), " to ", 
                            round(conf_interval[2], 3), ")", sep = "")
  }
  results_with_ci <- rbind(results_with_ci, data.frame(Column = col_name, 
                                                       Summary = result_with_ci))
}
write.csv(results_with_ci, "results_with_confidence_intervals.csv", row.names = FALSE)


#Section 6 - Plotting results

#ICER Plot
mean_qalys <- mean(results$`Incremental QALYs`, na.rm = TRUE)
mean_costs <- mean(results$`Incremental costs`, na.rm = TRUE)

ggplot(data = results, aes(x = `Incremental QALYs`, y = `Incremental costs`)) +
  geom_point(color = 'grey', alpha = 0.5) +
  geom_point(aes(x = mean_qalys, y = mean_costs), color = 'red', size = 1.25) +
  geom_hline(yintercept = 0, color = "black") +
  geom_vline(xintercept = 0, color = "black") +
  geom_abline(slope = 20000, intercept = 0, linetype = "dashed", color = "blue") + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white")
  ) +
  scale_x_continuous(breaks = seq(-0.8, max(results$`Incremental QALYs`, na.rm = TRUE), 0.2)) +
  scale_y_continuous(breaks = seq(-3000, max(results$`Incremental costs`, na.rm = TRUE), 1000)) +
  ylab("Incremental costs (£GBP)") +
  xlab("Incremental QALYs")

EDicer <- results$'ICER'
EDicer <- EDicer[!is.na(EDicer)]

#Cost per antibiotic day avoided plot

ggplot(data = results, aes(x = `Incremental antibiotics`, y = `Incremental costs`)) +
  geom_point(color = 'grey', alpha = 0.5) 

#General practice code#

#Input values#
#Selecting which model to run - untag the desired analysis#

scenario = "basecase"
#Scenario analyses#
#scenario = "nodiffabx"
#scenario = "noduabx"
#scenario = "nodiffadm"
#scenario = "nodifffail"
#scenario = "PCT50"

#SECTION 1 - Defining parameters#

#Standard care sample size
SC_n <- 501


#Initiation of antibiotic risk 
iAbx <- 0.631

#Odds
iAbxOdds <- iAbx / (1 - iAbx)
#Standard error
iAbx_SE <- sqrt(iAbx * (1 - iAbx) / SC_n)
#95% confidence interval
iAbx_95LL <- iAbx - qnorm(0.975) * iAbx_SE
iAbx_95UL <- iAbx + qnorm(0.975) * iAbx_SE
print(paste("95% Confidence Interval Lower Limit:", iAbx_95LL))
print(paste("95% Confidence Interval Upper Limit:", iAbx_95UL))
iAbx_var <- (iAbx-iAbx_95LL)/1.96
mu <- iAbx
var <- iAbx_var*iAbx_var

# Function to calculate alpha and beta for a Beta distribution
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(list(alpha = alpha, beta = beta))
}
iAbx_betaEstimates <- estBetaParams(mu, var)
print(iAbx_betaEstimates)
iAbx_alpha <- iAbx_betaEstimates$alpha
iAbx_beta <- iAbx_betaEstimates$beta
#Checking the alpha and beta parameters correspond with 95% confidence intervals
qbeta(c(0.025, 0.975), iAbx_betaEstimates$alpha, iAbx_betaEstimates$beta)

#Antibiotic initiation adjusted odds ratio
iAbxOR <- 0.13
iAbxOR_95LL <- 0.09
iAbxOR_95UL <- 0.18
if(scenario == "nodiffabx"){
  iAbxOR <- 1
  iAbxOR_95LL <- 1
  iAbxOR_95UL <- 1
} 
iAbxOR_SE <- (iAbxOR_95UL-iAbxOR)/qnorm(0.975)
iAbxOR_MoE <- (iAbxOR_95UL- iAbxOR_95LL)/2
iAbxOR_SD <- iAbxOR_MoE/qnorm(0.975)
iAbxOR_VAR <- iAbxOR_SD^2

PCTiAbxOdds <- iAbxOdds*iAbxOR
PCTiAbx <- PCTiAbxOdds/(PCTiAbxOdds+1)

#Antibiotic duration (days)
dAbx <- 7.3
dAbx_SD <- 2.5

#Antibiotic duration (days), standard error
dAbx_SE <- dAbx_SD/sqrt(SC_n)
#Antibiotic duration (days), variance (SE)
dAbx_VAR <- dAbx_SE**2
#Antibiotic duration (days), Gamma distribution parameters
dAbx_SCALE <- dAbx_VAR/dAbx
dAbx_SHAPE <- dAbx/dAbx_SCALE
dAbx_RATE <- 1/dAbx_SCALE
#Antibiotic duration (days), 95% confidence interval
dAbx_95LL <- dAbx - qnorm(0.975)*dAbx_SE
dAbx_95UL <- dAbx + qnorm(0.975)*dAbx_SE


#Difference in antibiotic duration (days)
diffAbx <- -0.52
#Difference in antibiotic duration (days), 95% confidence interval
diffAbx_95LL <- -1.07
diffAbx_95UL <- 0.04
if(scenario == "nodiffabx"){
  diffAbx = 0
  diffAbx_95LL = 0
  diffAbx_95UL = 0
} 

#Difference in antibiotic duration (days), margin of error and estimated standard deviation
diffAbx_MoE <- (diffAbx_95UL- diffAbx_95LL)/2
diffAbx_SD <- diffAbx_MoE/qnorm(0.975)
diffAbx_VAR <- diffAbx_SD^2

#PCT antibiotic duration (days)
PCTdAbx <- dAbx + diffAbx

#Total exposure of antibiotics
tAbx <- dAbx*iAbx
PCTtAbx <- PCTdAbx*PCTiAbx

#Days with restricted activities
LOS <- 8.9
#Days with restricted activities, SD
LOS_SD <- 4.2
#Days with restricted activities, SE
LOS_SE <- LOS_SD/sqrt(SC_n)
#Days with restricted activities, variance
LOS_VAR <- LOS_SE**2
#Days with restricted activities, Gamma distribution parameters
LOS_SCALE <- LOS_VAR/LOS
LOS_SHAPE <- LOS/LOS_SCALE
LOS_RATE <- 1/LOS_SCALE
#Days with restricted activities, 95% confidence interval
LOS_95LL <- LOS - qnorm(0.975)*LOS_SE
LOS_95UL <- LOS + qnorm(0.975)*LOS_SE

#Difference in days with restricted activities
diffLOS <- 0.07
#Difference in days with restricted activities, 95% confidence interval
diffLOS_95LL <- -0.44
diffLOS_95UL <- 0.59
if(scenario == "nodiffLOS"){
  diffLOS = 0
  diffLOS_95LL = 0
  diffLOS_95UL = 0
}
#Difference in days with restricted activities, margin of error and estimated standard deviation
diffLOS_MoE <- (diffLOS_95UL- diffLOS_95LL)/2
diffLOS_SD <- diffLOS_MoE/qnorm(0.975)
diffLOS_VAR <- diffLOS_SD^2

#PCT days with restricted activities
PCTLOS <- LOS + diffLOS

#Antibiotic cost (day)
cAbx <- 0.38
#PCT cost 16.26
cPCT <- 18.23
if(scenario == "PCT50"){
  cPCT <- 50
}
#Number of PCT tests
nPCT <- 2
#Cost of GP presentation
cGP <- 49
#Cost of LRTI admission
cAdmH <- (767.31*8.2)+(0.29*9.8)

#Treatment failure probability
pAdm <- 0.005306746
#Odds
pAdmOdds <- pAdm / (1 - pAdm)
#Standard error
pAdm_SD <- 0.0044185
#95% confidence interval
pAdm_95LL <- 0.0004084322
pAdm_95UL <- 0.01658112
print(paste("95% Confidence Interval Lower Limit:", pAdm_95LL))
print(paste("95% Confidence Interval Upper Limit:", pAdm_95UL))
pAdm_var <- ((pAdm-pAdm_95LL)/1.96)
mu <- pAdm
var <- pAdm_var*pAdm_var

# Function to calculate alpha and beta for a Beta distribution
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(list(alpha = alpha, beta = beta))
}
pAdm_betaEstimates <- estBetaParams(mu, var)
print(pAdm_betaEstimates)
pAdm_alpha <- pAdm_betaEstimates$alpha
pAdm_beta <- pAdm_betaEstimates$beta
#Checking the alpha and beta parameters correspond with 95% confidence intervals
qbeta(c(0.025, 0.975), pAdm_betaEstimates$alpha, pAdm_betaEstimates$beta)

pAdmOR <- 1.317062
pAdmOR_95LL <- 0.3921509
pAdmOR_95UL <- 3.32125
pAdmOR_SD <- 0.7659953 
if(scenario == "nodiffadm"){
  pAdmOR <- 1
  pAdmOR_95LL <- 1
  pAdmOR_95UL <- 1
  pAdmOR_SD <- 0
}
pAdmOR_SE <- sqrt(pAdmOR_SD/SC_n)
pAdmOR_VAR <- pAdmOR_SD^2
PCTpAdmOdds <- pAdmOdds*pAdmOR 
PCTpAdm <- PCTpAdmOdds/(PCTpAdmOdds+1)

pFail <- 0.327
#Odds
pFailOdds <- pFail / (1 - pFail)
#Standard error
pFail_SE <- sqrt(pFail * (1 - pFail) / SC_n)
#95% confidence interval
pFail_95LL <- pFail - qnorm(0.975) * pFail_SE
pFail_95UL <- pFail + qnorm(0.975) * pFail_SE
print(paste("95% Confidence Interval Lower Limit:", pFail_95LL))
print(paste("95% Confidence Interval Upper Limit:", pFail_95UL))
pFail_var <- (pFail-pFail_95LL)/1.96
mu <- pFail
var <- pFail_var*pFail_var

# Function to calculate alpha and beta for a Beta distribution
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(list(alpha = alpha, beta = beta))
}
pFail_betaEstimates <- estBetaParams(mu, var)
print(pFail_betaEstimates)
pFail_alpha <- pFail_betaEstimates$alpha
pFail_beta <- pFail_betaEstimates$beta
#Checking the alpha and beta parameters correspond with 95% confidence intervals
qbeta(c(0.025, 0.975), pFail_betaEstimates$alpha, pFail_betaEstimates$beta)

#Treatment failure odds ratio
pFailOR <- 0.96
pFailOR_95LL <- 0.73
pFailOR_95UL <- 1.25
if(scenario == "nodifffail"){
  pFailOR <- 1
  pFailOR_95LL <- 1
  pFailOR_95UL <- 1
}
#Mortality adjusted OR, standard error and distribution
pFailOR_SE <- (pFailOR_95UL-pFailOR)/qnorm(0.975)
#Mortality adjusted OR, margin of error and estimated standard deviation
pFailOR_MoE <- (pFailOR_95UL- pFailOR_95LL)/2
pFailOR_SD <- pFailOR_MoE/qnorm(0.975)
pFailOR_VAR <- pFailOR_SD^2

PCTpFailOdds <- pFailOdds*pFailOR
PCTpFail <- PCTpFailOdds/(PCTpFailOdds+1)

#Proportion of female LRTI patients
SepPF <- 0.49
SepPF_n <- 2905
#Standard error
SepPF_SE <- sqrt(SepPF*(1-SepPF)/SepPF_n)
#95% confidence interval
SepPF_95LL <- SepPF - qnorm(0.975)*SepPF_SE
SepPF_95UL <- SepPF + qnorm(0.975)*SepPF_SE
print(paste("95% Confidence Interval Lower Limit:", SepPF_95LL))
print(paste("95% Confidence Interval Upper Limit:", SepPF_95UL))
SepPF_var <- ((SepPF-SepPF_95LL)/1.96)
mu <- SepPF
var <- SepPF_var*SepPF_var
# Function to calculate alpha and beta for a Beta distribution
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(list(alpha = alpha, beta = beta))
}
SepPF_betaEstimates <- estBetaParams(mu, var)
print(SepPF_betaEstimates)
SepPF_alpha <- SepPF_betaEstimates$alpha
SepPF_beta <- SepPF_betaEstimates$beta
#Checking the alpha and beta parameters correspond with 95% confidence intervals
qbeta(c(0.025, 0.975), SepPF_betaEstimates$alpha, SepPF_betaEstimates$beta)
#Proportion of male patients
SepPM <- 1-SepPF

#Female QALY 
uF <- 0.775
#Female QALY, sample size
uF_n <- 619
#Female QALY 95% confidence interval
uF_95LL <- 0.770
uF_95UL <- 0.795
#Female QALY standard error
uF_SE <- sqrt(uF*(1-uF)/uF_n)
uF_var <- (uF-uF_95LL)/1.96
mu <- uF
var <- uF_var*uF_var

# Function to calculate alpha and beta for a Beta distribution
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(list(alpha = alpha, beta = beta))
}
uF_betaEstimates <- estBetaParams(mu, var)
print(uF_betaEstimates)
uF_alpha <- uF_betaEstimates$alpha
uF_beta <- uF_betaEstimates$beta
#Checking the alpha and beta parameters correspond with 95% confidence intervals
qbeta(c(0.025, 0.975), uF_betaEstimates$alpha, uF_betaEstimates$beta)

#Male QALY 
uM <- 0.797
#Male QALY, sample size
uM_n <- 568
#Male QALY 60-64 years old, 95% confidence interval
uM_95LL <- 0.792
uM_95UL <- 0.818
#Standard care mortality rate, standard error
uM_SE <- sqrt(uM*(1-uM)/uM_n)
uM_var <- (uM-uM_95LL)/1.96
mu <- uM
var <- uM_var*uM_var

# Function to calculate alpha and beta for a Beta distribution
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(list(alpha = alpha, beta = beta))
}
uM_betaEstimates <- estBetaParams(mu, var)
print(uM_betaEstimates)
uM_alpha <- uM_betaEstimates$alpha
uM_beta <- uM_betaEstimates$beta
#Checking the alpha and beta parameters correspond with 95% confidence intervals
qbeta(c(0.025, 0.975), uM_betaEstimates$alpha, uM_betaEstimates$beta)

#Weighted population QALY 
uW <- (uF*SepPF) + (uM*SepPM)

#LRTI utility
uSep <- 0.705
#LRTI utility, sample size
uSep_n <- 349
#LRTI utility, SE
uSep_SE <- sqrt(uSep*(1-uSep)/uSep_n)
#LRTI utility, 95% confidence interval
uSep_95LL <- uSep - qnorm(0.975)*uSep_SE
uSep_95UL <- uSep + qnorm(0.975)*uSep_SE
print(paste("95% Confidence Interval Lower Limit:", uSep_95LL))
print(paste("95% Confidence Interval Upper Limit:", uSep_95UL))
uSep_var <- (uSep-uSep_95LL)/1.96
mu <- uSep
var <- uSep_var*uSep_var
# Function to calculate alpha and beta for a Beta distribution
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(list(alpha = alpha, beta = beta))
}
uSep_betaEstimates <- estBetaParams(mu, var)
print(uSep_betaEstimates)
uSep_alpha <- uSep_betaEstimates$alpha
uSep_beta <- uSep_betaEstimates$beta
#Checking the alpha and beta parameters correspond with 95% confidence intervals
qbeta(c(0.025, 0.975), uSep_betaEstimates$alpha, uSep_betaEstimates$beta)
#LRTI disutility
duSep <- (uW-uSep)/24

#Antibiotic disutility (QALD/day)
duAbxD <- 0.057
if(scenario == "noduabx"){
  duAbxD = 0
}
duAbxD_n <- 349
duAbxD_SE <- sqrt(duAbxD*(1-duAbxD)/duAbxD_n)
duAbxD_95LL <- duAbxD - qnorm(0.975)*duAbxD_SE
duAbxD_95UL <- duAbxD + qnorm(0.975)*duAbxD_SE
print(paste("95% Confidence Interval Lower Limit:", duAbxD))
print(paste("95% Confidence Interval Upper Limit:", duAbxD))

duAbxD_var <- (duAbxD-duAbxD_95LL)/1.96
mu <- duAbxD
var <- duAbxD_var*duAbxD_var

# Function to calculate alpha and beta for a Beta distribution
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(list(alpha = alpha, beta = beta))
}
duAbxD_betaEstimates <- estBetaParams(mu, var)
print(duAbxD_betaEstimates)
duAbxD_alpha <- duAbxD_betaEstimates$alpha
duAbxD_beta <- duAbxD_betaEstimates$beta
#Checking the alpha and beta parameters correspond with 95% confidence intervals
qbeta(c(0.025, 0.975), duAbxD_betaEstimates$alpha, duAbxD_betaEstimates$beta)

#Total antibiotic disutility (QALY)
duAbx <- duAbxD*tAbx/365.25
#PCT total antibiotic disutility (QALY)
PCTduAbx <- duAbxD*PCTtAbx/365.25

#Female QALE, 3.5% discount
QALEF <- 10.66
#Male QALE, 3.5% discount 
QALEM <- 10.23
#Weighted QALE, 3.5% discount
QALEW <- (QALEF*SepPF) + (QALEM*SepPM)

#Female LE 
LEF <- 20.46
#Male LE 
LEM <- 18.13
#Weighted LE 
LEW <- (LEF*SepPF) + (LEM*SepPM)

#SECTION 2 - Data frame for model inputs#

input <- data.frame(
  #1
  iAbx <- iAbx,
  iAbxOdds <- iAbx/(1-iAbx),
  iAbxOR <- iAbxOR,
  PCTiAbxOdds <- iAbxOdds*iAbxOR,
  PCTiAbx <- PCTiAbxOdds/(PCTiAbxOdds+1),
  #6
  dAbx <- dAbx,
  diffAbx <- diffAbx,
  PCTdAbx <- dAbx + diffAbx,
  tAbx <- dAbx*iAbx,
  PCTtAbx <- PCTdAbx*PCTiAbx,
  #11
  LOS <- LOS,
  diffLOS <- diffLOS,
  PCTLOS <- LOS + diffLOS,
  cAbx <- cAbx,
  cPCT <- cPCT,
  #16
  nPCT <- nPCT,
  cGP <- cGP,
  cAdmH <- cAdmH,
  pAdm <- pAdm,
  pAdmOdds <- pAdm/(1-pAdm),
  #21
  pAdmOR <- pAdmOR,
  PCTpAdmOdds <- pAdmOdds*pAdmOR,
  PCTpAdm <- PCTpAdmOdds/(PCTpAdmOdds+1),
  pFail <- pFail,
  pFailOdds <- pFail/(1-pFail),
  #26
  pFailOR <- pFailOR,
  PCTpFailOdds <- pFailOdds*pFailOR,
  PCTpFail <- PCTpFailOdds/(PCTpFailOdds+1),
  SepPF <- SepPF,
  SepPM <- 1-SepPF,
  #31
  uF <- uF,
  uM <- uM,
  uW <- (uF*SepPF) + (uM*SepPM),
  uSep <- uSep,
  duSep <- (uW-uSep)/24,
  #36
  duAbx <- duAbxD*tAbx/365.25,
  PCTduAbx <- duAbxD*PCTtAbx/365.25,
  QALEW <- (QALEF*SepPF) + (QALEM*SepPM),
  LEW <- (LEF*SepPF) + (LEM*SepPM)
)

#Section 3 - Decision Tree 
dec_tree <- function(params){
  with(
    as.list(params), 
    {
      iAbx <- iAbx
      iAbxOdds <- iAbx/(1-iAbx)
      iAbxOR <- iAbxOR
      PCTiAbxOdds <- iAbxOdds*iAbxOR
      PCTiAbx <- PCTiAbxOdds/(PCTiAbxOdds+1)
      dAbx <- dAbx
      diffAbx <- diffAbx
      PCTdAbx <- dAbx + diffAbx
      tAbx <- dAbx*iAbx
      PCTtAbx <- PCTdAbx*PCTiAbx
      LOS <- LOS
      diffLOS <- diffLOS
      PCTLOS <- LOS + diffLOS
      cAbx <- cAbx
      cPCT <- cPCT
      nPCT <- nPCT
      cGP <- cGP
      cAdmH <- cAdmH
      pAdm <- pAdm
      pAdmOdds <- pAdm/(1-pAdm)
      pAdmOR <- pAdmOR
      PCTpAdmOdds <- pAdmOdds*pAdmOR
      PCTpAdm <- PCTpAdmOdds/(PCTpAdmOdds+1)
      pFail <- pFail
      pFailOdds <- pFail/(1-pFail)
      pFailOR <- pFailOR
      PCTpFailOdds <- pFailOdds*pFailOR
      PCTpFail <- PCTpFailOdds/(PCTpFailOdds+1)
      SepPF <- SepPF
      SepPM <- 1-SepPF
      uF <- uF
      uM <- uM
      uW <- (uF*SepPF) + (uM*SepPM)
      uSep <- uSep
      duSep <- (uW-uSep)/24
      duAbx <- duAbxD*tAbx/365.25
      PCTduAbx <- duAbxD*PCTtAbx/365.25
      QALEW <- (QALEF*SepPF) + (QALEM*SepPM)
      LEW <- (LEF*SepPF) + (LEM*SepPM)
      
      #Standard care costs
      C_SC <- (cAbx*tAbx) + (cGP+(cGP*pFail)) + (pAdm*cAdmH)
      #PCT costs
      C_PCT <- (cPCT*nPCT) + (cAbx*PCTtAbx) + (cGP+(cGP*PCTpFail)) + (PCTpAdm*cAdmH)
      #Standard care QALYs gained
      Q_SC <- (QALEW-duSep-duAbx-(duSep*pFail))
      #PCT QALYs gained
      Q_PCT <- (QALEW-duSep-PCTduAbx)
      #Incremental Costs
      IC <- C_PCT - C_SC
      #Incremental QALYs
      IQALYx <- Q_PCT - Q_SC
      IQALY <- IQALYx
      #IQALY <- round(IQALYx, 7)
      
      if (scenario == "basecase" | scenario == "nodiffadm" | scenario == "nodifffail" | scenario == "PCT50"){
        if (IQALY == 0) {
          ICER <- NA
          ICER_flag <- "Excluded: IQALY of 0"
        } else {
          # Incremental Costs
          IC <- C_PCT - C_SC
          # ICER
          ICER <- IC / IQALY
          # Normal flag for non-excluded cases
          ICER_flag <- "Normal"
        }} else if (scenario == "nodiffabx" | scenario == "noduabx"){
          # Incremental Costs
          IC <- C_PCT - C_SC
          # ICER
          ICER <- IC / IQALY
        }
      #Incremental Antibiotics
      IAbx <- PCTtAbx - tAbx
      #Cost per antibiotic day avoided
      AbxAv <- (C_PCT - C_SC)/(tAbx-PCTtAbx)
      #Antibiotic value
      AbxV <- (IC-(IQALY*20000))/(tAbx-PCTtAbx)
      Inputs <- c(iAbx, PCTiAbx, LOS, PCTLOS, pAdm, PCTpAdm, pFail, PCTpFail)
      BasicOutcomes <- c(C_SC, C_PCT, Q_SC, Q_PCT, tAbx, PCTtAbx)
      ICERs <- c(IC, IQALY, IAbx, ICER, AbxAv, AbxV)
      
      names(Inputs) <- paste ("Inputs", c("Standard care antibiotic initiation", "Procalcitonin antibiotic intiation", "Days with restricted activities", "Procalcitonin days with restricted activities", "Admission probability", "Procalcitonin admission probability", "Treatment failure probability", "Procalcitonin treatment failure probability"))
      names(BasicOutcomes) <- paste ("Basic Outcomes", c("Standard care costs", "Procalcitonin costs", "Standard care QALYs", "Procalcitonin QALYs", "Standard care antibiotic duration", "Procalcitonin antibiotic duration"), sep = "_")
      names(ICERs)  <- paste("ICERs", c("Incremental costs", "Incremental QALYs", "Incremental antibiotics", "ICER", "Cost per antibiotic day avoided", "Antibiotic value"), sep = "_")
      
      return(c(Inputs, BasicOutcomes, ICERs))
    }
  )
}
options(scipen=999)
dec_tree(input)

#Section 4 - Tornado Plot#
if (scenario == "basecase"){
  iAbx_range <- c(BaseCase = iAbx, low = iAbx_95LL, high = iAbx_95UL)
  iAbxOR_range <- c(BaseCase = iAbxOR, low = iAbxOR_95LL, high = iAbxOR_95UL)
  dAbx_range <- c(BaseCase = dAbx, low = dAbx_95LL, high = dAbx_95UL)
  diffAbx_range <- c(BaseCase = diffAbx, low = diffAbx_95LL, high = diffAbx_95UL)
  LOS_range <- c(BaseCase = LOS, low = LOS_95LL, high = LOS_95UL)
  diffLOS_range <- c(BaseCase = diffLOS, low = diffLOS_95LL, high = diffLOS_95UL)
  cAbx_range <- c(BaseCase = cAbx, low = cAbx*0.75, high = cAbx*1.25)
  cPCT_range <- c(BaseCase = cPCT, low = 14.5, high = 90)
  nPCT_range <- c(BaseCase = nPCT, low = 1, high = 3)
  cGP_range <- c(BaseCase = cGP, low = cGP*0.75, high = cGP*1.25)
  cAdmH_range <- c(BaseCase = cAdmH, low = cAdmH*0.75, high = cAdmH*1.25)
  pAdm_range <- c(BaseCase = pAdm, low = pAdm_95LL, high = pAdm_95UL)
  pAdmOR_range <- c(BaseCase = pAdmOR, low = pAdmOR_95LL, high = pAdmOR_95UL)
  pFail_range <- c(BaseCase = pFail, low = pFail_95LL, high = pFail_95UL)
  pFailOR_range <- c(BaseCase = pFailOR, low = pFailOR_95LL, high = pFailOR_95UL)
  SepPF_range <- c(BaseCase = SepPF, low = SepPF_95LL, high = SepPF_95UL)
  uF_range <- c(BaseCase = uF, low = uF_95LL, high = uF_95UL)
  uM_range <- c(BaseCase = uM, low = uM_95LL, high = uM_95UL)
  uSep_range <- c(BaseCase = uSep, low = uSep_95LL, high = uSep_95UL)
  
  #Plotting all parameters to find most influential
  paramNames <- c(
    "Antibiotic initiation probability",
    "Antibiotic initiation odds ratio",
    "Duration of antibiotics",
    "Antibiotic duration difference",
    "Days with restricted activity",
    "Days with restricted activity difference",
    "Antibiotic cost",
    "Procalcitonin cost per test",
    "Number of procalcitonin tests used",
    "GP visit cost",
    "Admission cost",
    "Admission probability",
    "Admission probability odds ratio",
    "Treatment failure probability",
    "Treatment failure probability odds ratio",
    "Proportion of female patients",
    "Female utility",
    "Male utility",
    "LRTI disutility"
  )
  
  l.tor.in <- vector("list", 19)
  names(l.tor.in) <- paramNames
  l.tor.in$'Antibiotic initiation probability' <- cbind(iAbx = iAbx_range, input [-1])
  l.tor.in$'Antibiotic initiation odds ratio' <- cbind(iAbxOR = iAbxOR_range, input [-3])
  l.tor.in$'Duration of antibiotics' <- cbind(dAbx = dAbx_range, input [-6])
  l.tor.in$'Antibiotic duration difference' <- cbind(diffAbx = diffAbx_range, input [-7])
  l.tor.in$'Days with restricted activity' <- cbind(LOS = LOS_range, input [-11])
  l.tor.in$'Days with restricted activity difference' <- cbind(diffLOS = diffLOS_range, input [-12])
  l.tor.in$'Antibiotic cost' <- cbind(cAbx = cAbx_range, input [-14])
  l.tor.in$'Procalcitonin cost per test' <- cbind(cPCT = cPCT_range, input [-15])
  l.tor.in$'Number of procalcitonin tests used' <- cbind(nPCT = nPCT_range, input [-16])
  l.tor.in$'GP visit cost' <- cbind(cGP = cGP_range, input [-17])
  l.tor.in$'Admission cost' <- cbind(cAdmH = cAdmH_range, input [-18])
  l.tor.in$'Admission probability' <- cbind(pAdm = pAdm_range, input [-19])
  l.tor.in$'Admission probability odds ratio' <- cbind(pAdmOR = pAdmOR_range, input [-21])
  l.tor.in$'Treatment failure probability' <- cbind(pFail = pFail_range, input [-24])
  l.tor.in$'Treatment failure probability odds ratio' <- cbind(pFailOR = pFailOR_range, input [-26])
  l.tor.in$'Proportion of female patients' <- cbind(SepPF = SepPF_range, input [-29])
  l.tor.in$'Female utility' <- cbind(uF = uF_range, input [-31])
  l.tor.in$'Male utility' <- cbind(uM = uM_range, input [-32])
  l.tor.in$'LRTI disutility' <- cbind(uSep = uSep_range, input [-34])
  
  #List of outputs
  l.tor.out <- vector("list", 19)
  names(l.tor.out) <- paramNames
  
  for(i in 1:19){
    l.tor.out[[i]] <- t(apply(l.tor.in[[i]], 1, dec_tree))[ , 18]
  }
  
  m.tor <- matrix(unlist(l.tor.out), nrow = 19, ncol = 3, byrow = TRUE, 
                  dimnames = list(paramNames, c("basecase", "low", "high")))
  
  TornadoPlot(main_title = "", Parms = paramNames, Outcomes = m.tor, 
              outcomeName = "", 
              xlab = "", 
              ylab = "", 
              col1="#3182bd", col2="#6baed6"
  )
  
  #Plotting top 8 parameters
  
  paramNames <- c(
    "Procalcitonin cost per test",
    "Number of procalcitonin tests used",
    "Admission cost",
    "Admission probability",
    "Admission probability odds ratio",
    "Treatment failure probability",
    "Treatment failure probability odds ratio",
    "LRTI disutility"
  )
  
  l.tor.in <- vector("list", 8)
  names(l.tor.in) <- paramNames
  l.tor.in$'Procalcitonin cost per test' <- cbind(cPCT = cPCT_range, input [-15])
  l.tor.in$'Number of procalcitonin tests used' <- cbind(nPCT = nPCT_range, input [-16])
  l.tor.in$'Admission cost' <- cbind(cAdmH = cAdmH_range, input [-18])
  l.tor.in$'Admission probability' <- cbind(pAdm = pAdm_range, input [-19])
  l.tor.in$'Admission probability odds ratio' <- cbind(pAdmOR = pAdmOR_range, input [-21])
  l.tor.in$'Treatment failure probability odds ratio' <- cbind(pFailOR = pFailOR_range, input [-26])
  l.tor.in$'Treatment failure probability' <- cbind(pFail = pFail_range, input [-24])
  l.tor.in$'LRTI disutility' <- cbind(uSep = uSep_range, input [-34])
  
  
  #List of outputs
  l.tor.out <- vector("list", 8)
  names(l.tor.out) <- paramNames
  
  for(i in 1:8){
    l.tor.out[[i]] <- t(apply(l.tor.in[[i]], 1, dec_tree))[ , 18]
  }
  
  m.tor <- matrix(unlist(l.tor.out), nrow = 8, ncol = 3, byrow = TRUE, 
                  dimnames = list(paramNames, c("basecase", "low", "high")))
  
  TornadoPlot(main_title = "", Parms = paramNames, Outcomes = m.tor, 
              outcomeName = "", 
              xlab = "", 
              ylab = "", 
              col1="#3182bd", col2="#6baed6"
  )}


#Section 5 - Monte Carlo Simulation

#Defining parallel parameters to prevent simulation iterating on itself

#Antibiotic initiation adjusted odds ratio
mciAbxOR <- 0.13
mciAbxOR_95LL <- 0.09
mciAbxOR_95UL <- 0.18
if(scenario == "nodiffabx"){
  mciAbxOR <- 1
  mciAbxOR_95LL <- 1
  mciAbxOR_95UL <- 1
} 
mciAbxOR_SE <- (mciAbxOR_95UL-mciAbxOR)/qnorm(0.975)
mciAbxOR_MoE <- (mciAbxOR_95UL- mciAbxOR_95LL)/2
mciAbxOR_SD <- mciAbxOR_MoE/qnorm(0.975)
mciAbxOR_VAR <- mciAbxOR_SD^2

#Difference in antibiotic duration (days)
mcdiffAbx <- -0.52
#Difference in antibiotic duration (days), 95% confidence interval
mcdiffAbx_95LL <- -1.07
mcdiffAbx_95UL <- 0.04
if(scenario == "nodiffabx"){
  mcdiffAbx = 0
  mcdiffAbx_95LL = 0
  mcdiffAbx_95UL = 0
} 

#Difference in antibiotic duration (days), margin of error and estimated standard deviation
mcdiffAbx_MoE <- (mcdiffAbx_95UL- mcdiffAbx_95LL)/2
mcdiffAbx_SD <- mcdiffAbx_MoE/qnorm(0.975)
mcdiffAbx_VAR <- mcdiffAbx_SD^2

#Difference in days with restricted activities
mcdiffLOS <- 0.07
#Difference in days with restricted activities, 95% confidence interval
mcdiffLOS_95LL <- -0.44
mcdiffLOS_95UL <- 0.59
if(scenario == "nodiffLOS"){
  mcdiffLOS = 0
  mcdiffLOS_95LL = 0
  mcdiffLOS_95UL = 0
}
#Difference in days with restricted activities, margin of error and estimated standard deviation
mcdiffLOS_MoE <- (mcdiffLOS_95UL- mcdiffLOS_95LL)/2
mcdiffLOS_SD <- mcdiffLOS_MoE/qnorm(0.975)
mcdiffLOS_VAR <- mcdiffLOS_SD^2

mcpAdmOR <- 1.317062
mcpAdmOR_95LL <- 0.3921509
mcpAdmOR_95UL <- 3.32125
mcpAdmOR_SD <- 0.7659953 
if(scenario == "nodiffadm"){
  mcpAdmOR <- 1
  mcpAdmOR_95LL <- 1
  mcpAdmOR_95UL <- 1
  mcpAdmOR_SD <- 0
}
mcpAdmOR_SE <- sqrt(mcpAdmOR_SD/SC_n)
mcpAdmOR_VAR <- mcpAdmOR_SD^2

#Treatment failure odds ratio
mcpFailOR <- 0.96
mcpFailOR_95LL <- 0.73
mcpFailOR_95UL <- 1.25
if(scenario == "nodifffail"){
  mcpFailOR <- 1
  mcpFailOR_95LL <- 1
  mcpFailOR_95UL <- 1
}
#Mortality adjusted OR, standard error and distribution
mcpFailOR_SE <- (mcpFailOR_95UL-mcpFailOR)/qnorm(0.975)
#Mortality adjusted OR, margin of error and estimated standard deviation
mcpFailOR_MoE <- (mcpFailOR_95UL- mcpFailOR_95LL)/2
mcpFailOR_SD <- mcpFailOR_MoE/qnorm(0.975)
mcpFailOR_VAR <- mcpFailOR_SD^2

# Number of iterations for the simulation
n <- 10000

# Create a data frame to store the results
results <- data.frame(matrix(ncol = 20, nrow = n))
names(results) <- c("Standard care antibiotic initiation", "Procalcitonin antibiotic initiation", "Days with restricted activities", "Procalcitonin days with restricted activities", "Admission probability", "Procalcitonin admission probability", "Treatment failure probability", "Procalcitonin treatment failure probability", "Standard care costs", "Procalcitonin costs", "Standard care QALYs", "Procalcitonin QALYs", "Standard care antibiotic duration", "Procalcitonin antibiotic duration", "Incremental costs", "Incremental QALYs", "Incremental antibiotics", "ICER", "Cost per antibiotic day avoided", "Antibiotic value")

# Perform the simulation

# Perform the simulation
if(scenario == "basecase"| scenario == "nodiffabx" | scenario == "nodiffadm" | scenario == "nodifffail" | scenario == "PCT50"){
  for (i in 1:n) {
    
    iAbx <- rbeta(1, iAbx_alpha, iAbx_beta)
    
    iAbxOR <- rlnorm(1, log(mciAbxOR), (log(mciAbxOR_95UL)-log(mciAbxOR_95LL))/(2*qnorm(0.975)))
    
    dAbx <- rgamma(1, dAbx_SHAPE, dAbx_RATE)
    
    diffAbx <- rnorm(1, mcdiffAbx, mcdiffAbx_SD)
    
    LOS <- rgamma(1, LOS_SHAPE, LOS_RATE)
    
    diffLOS <- rnorm(1, mcdiffLOS, mcdiffLOS_SD)
    
    pAdm <- rbeta(1, pAdm_alpha, pAdm_beta)
    
    pAdmOR <- rlnorm(1, log(mcpAdmOR), (log(mcpAdmOR_95UL)-log(mcpAdmOR_95LL))/(2*qnorm(0.975)))
    
    pFail <- rbeta(1, pFail_alpha, pFail_beta)
    
    pFailOR <- rlnorm(1, log(mcpFailOR), (log(mcpFailOR_95UL)-log(mcpFailOR_95LL))/(2*qnorm(0.975)))
    
    SepPF <- rbeta(1, SepPF_alpha, SepPF_beta)
    
    uF <- rbeta (1, uF_alpha, uF_beta)
    
    uM <- rbeta (1, uM_alpha, uM_beta)
    
    uSep <- rbeta (1, uSep_alpha, uSep_beta)
    
    duAbxD <- rbeta (1, duAbxD_alpha, duAbxD_beta)
    
    # Create the parameter vector for this iteration
    iteration_params <- c(iAbx = iAbx, iAbxOR = iAbxOR, dAbx = dAbx, diffAbx = diffAbx, LOS = LOS, diffLOS = diffLOS, pAdm = pAdm, pAdmOR = pAdmOR, pFail = pFail, pFailOR = pFailOR, SepPF = SepPF, uF = uF, uM = uM, uSep = uSep, duAbxD = duAbxD)
    
    # Calculate the decision tree for this iteration
    iteration_results <- dec_tree(iteration_params)
    
    # Store the results
    results[i,] <- iteration_results
  }} else if(scenario == "noduabx"){
    for (i in 1:n) {
      
      iAbx <- rbeta(1, iAbx_alpha, iAbx_beta)
      
      iAbxOR <- rlnorm(1, log(mciAbxOR), (log(mciAbxOR_95UL)-log(mciAbxOR_95LL))/(2*qnorm(0.975)))
      
      dAbx <- rgamma(1, dAbx_SHAPE, dAbx_RATE)
      
      diffAbx <- rnorm(1, mcdiffAbx, mcdiffAbx_SD)
      
      LOS <- rgamma(1, LOS_SHAPE, LOS_RATE)
      
      diffLOS <- rnorm(1, mcdiffLOS, mcdiffLOS_SD)
      
      pAdm <- rbeta(1, pAdm_alpha, pAdm_beta)
      
      pAdmOR <- rlnorm(1, log(mcpAdmOR), (log(mcpAdmOR_95UL)-log(mcpAdmOR_95LL))/(2*qnorm(0.975)))
      
      pFail <- rbeta(1, pFail_alpha, pFail_beta)
      
      pFailOR <- rlnorm(1, log(mcpFailOR), (log(mcpFailOR_95UL)-log(mcpFailOR_95LL))/(2*qnorm(0.975)))
      
      SepPF <- rbeta(1, SepPF_alpha, SepPF_beta)
      
      uF <- rbeta (1, uF_alpha, uF_beta)
      
      uM <- rbeta (1, uM_alpha, uM_beta)
      
      uSep <- rbeta (1, uSep_alpha, uSep_beta)
      
      # Create the parameter vector for this iteration
      iteration_params <- c(iAbx = iAbx, iAbxOR = iAbxOR, dAbx = dAbx, diffAbx = diffAbx, LOS = LOS, diffLOS = diffLOS, pAdm = pAdm, pAdmOR = pAdmOR, pFail = pFail, pFailOR = pFailOR, SepPF = SepPF, uF = uF, uM = uM, uSep = uSep, duAbxD = duAbxD)
      
      # Calculate the decision tree for this iteration
      iteration_results <- dec_tree(iteration_params)
      
      # Store the results
      results[i,] <- iteration_results
    }
  }

# Calculate summary statistics
results_summary <- summary(results)

# View summary statistics
print(results_summary)

# For each column in the results dataframe
for (col_name in names(results)) {
  # Calculate the 95% confidence interval
  conf_interval <- quantile(results[[col_name]], probs = c(0.025, 0.975), na.rm = TRUE)
  
  # Print the confidence interval
  print(paste("95% confidence interval for", col_name, ":", conf_interval))
}

#Exporting results
results_with_ci <- data.frame(Column = character(),
                              Summary = character(),
                              stringsAsFactors = FALSE)
results_summary <- summary(results)
for (col_name in names(results)) {
  conf_interval <- quantile(results[[col_name]], probs = c(0.025, 0.975), na.rm = TRUE)
  result_value <- mean(results[[col_name]], na.rm = TRUE)
  if (grepl("cost", col_name, ignore.case = TRUE) || grepl("icer", col_name, ignore.case = TRUE)) {
    result_with_ci <- paste(round(result_value), " (", 
                            round(conf_interval[1]), " to ", 
                            round(conf_interval[2]), ")", sep = "")
  } else if (grepl("value", col_name, ignore.case = TRUE))
  {
    result_with_ci <- paste(round(result_value, 2), " (", 
                            round(conf_interval[1], 2), " to ", 
                            round(conf_interval[2], 2), ")", sep = "")
  }
  else {
    result_with_ci <- paste(round(result_value, 3), " (", 
                            round(conf_interval[1], 3), " to ", 
                            round(conf_interval[2], 3), ")", sep = "")
  }
  results_with_ci <- rbind(results_with_ci, data.frame(Column = col_name, 
                                                       Summary = result_with_ci))
}
write.csv(results_with_ci, "results_with_confidence_intervals.csv", row.names = FALSE)


#Section 6 - Plotting results

#ICER Plot
mean_qalys <- mean(results$`Incremental QALYs`, na.rm = TRUE)
mean_costs <- mean(results$`Incremental costs`, na.rm = TRUE)

ggplot(data = results, aes(x = `Incremental QALYs`, y = `Incremental costs`)) +
  geom_point(color = 'grey', alpha = 0.5) +
  geom_point(aes(x = mean_qalys, y = mean_costs), color = 'red', size = 1.25) +
  geom_hline(yintercept = 0, color = "black") +
  geom_vline(xintercept = 0, color = "black") +
  geom_abline(slope = 20000, intercept = 0, linetype = "dashed", color = "blue") +
  scale_y_continuous(labels = c(0, 100, 200, 300, 400, 500, 600), breaks = c(0, 100, 200, 300, 400, 500, 600)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white")
  ) +
  ylab("Incremental costs (£GBP)") +
  xlab("Incremental QALYs")

GPicer <- results$'ICER'
GPicer <- GPicer[!is.na(GPicer)]


#Cost per antibiotic day avoided plot

ggplot(data = results, aes(x = `Incremental antibiotics`, y = `Incremental costs`)) +
  geom_point(color = 'grey', alpha = 0.5) 



####Combined CEAC plot####

max_wtp <- 50000

# Function to calculate CEAC data
calculate_ceac <- function(ICERs, max_wtp) {
  wtp_values <- seq(0, max_wtp, by = 1000)
  prob_values <- sapply(wtp_values, function(x) sum(ICERs <= x) / length(ICERs))
  ceac_data <- data.frame(WTP = wtp_values, Prob = prob_values)
  return(ceac_data)
}

# Calculate CEAC data for GPicers
ceac_data_gp <- calculate_ceac(GPicer, max_wtp)
ceac_data_gp$Source <- 'GP'

# Calculate CEAC data for EDicers
ceac_data_ed <- calculate_ceac(EDicer, max_wtp)
ceac_data_ed$Source <- 'ED'

# Calculate CEAC data for ICUicers
ceac_data_icu <- calculate_ceac(ICUicer, max_wtp)
ceac_data_icu$Source <- 'ICU'

# Combine all CEAC data into one data frame
combined_ceac_data <- bind_rows(ceac_data_gp, ceac_data_ed, ceac_data_icu)

# Plotting using ggplot2
ggplot() +
  geom_line(data = combined_ceac_data, aes(x = WTP, y = Prob, color = Source), linetype = "solid", size = 0.55) +
  labs(x = "Willingness to pay (£000/QALY)", 
       y = "Probability of cost-effectiveness (%)") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = function(x) paste0(x / 1000), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, max_wtp)) +
  theme_minimal() +
  theme(plot.title = element_blank(), # Remove the title
        axis.title.x = element_text(size = 10, color = "black"),  
        axis.title.y = element_text(size = 10, color = "black"),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        axis.line = element_line(color = "black", size = 0.5, linetype = "solid"),
        axis.ticks = element_line(color = "black"),
        legend.position = "right",
        legend.title = element_blank()) +
  scale_color_manual(values = c('GP' = "#d7191c", 'ED' = "#fdae61", 'ICU' = "#2c7bb6"), breaks = c('GP', 'ED', 'ICU')) +
  scale_linetype_manual(values = c('solid'))

# Print the data points for specific WTP values for GP, ED, and ICU
print(combined_ceac_data[combined_ceac_data$WTP %in% c(20000, 30000), ])

#END OF PRIMARY ANALYSIS#      

#Bayesian meta-analysis to generate values for GP#

rm(list = ls())

data <- data.frame(
  author = c("Burkhardt", "Briel"),
  event_intervention = c(1, 5),
  n_intervention = c(275, 232),
  event_control = c(1, 4),
  n_control = c(275, 226)
)

library(brms)
library(dplyr)


# Define the model formula
formula <- bf(event | trials(n) ~ group + (1 | author))

# Prepare the data for modeling
data_long <- tidyr::pivot_longer(
  data,
  cols = c(event_intervention, n_intervention, event_control, n_control),
  names_to = c(".value", "group"),
  names_pattern = "(.*)_(.*)"
) %>%
  mutate(group = factor(group, levels = c("control", "intervention")))

# Fit the model
model <- brm(
  formula = formula,
  data = data_long,
  family = binomial(),
  priors <- c(
    set_prior("normal(0, 1)", class = "b", coef = "groupintervention"),  
    set_prior("normal(-6.9, 1.2)", class = "Intercept"),                      
    set_prior("cauchy(0, 0.5)", class = "sd")                             
  ),
  chains = 4,
  iter = 2000,
  control = list(adapt_delta = 0.95),
  seed = 123  
)

# Extract posterior samples for the intercept (control group probability)
intercept_samples <- posterior_samples(model, pars = "b_Intercept")

# Calculate the control group probability (probability of success in the control group)
odds_control <- exp(intercept_samples$b_Intercept)
prob_control <- odds_control / (1 + odds_control)

# Calculate the mean, 95% credible interval, and standard deviation for the control group probability
mean_prob_control <- mean(prob_control)
ci_lower_prob_control <- quantile(prob_control, probs = 0.025)
ci_upper_prob_control <- quantile(prob_control, probs = 0.975)
sd_prob_control <- sd(prob_control)

# Print the results for control group probability
cat("Control Group Probability:\n")
cat("  Mean:", mean_prob_control, "\n")
cat("  95% Credible Interval: [", ci_lower_prob_control, ", ", ci_upper_prob_control, "]\n")
cat("  Standard Deviation:", sd_prob_control, "\n\n")

# Extract posterior samples for the coefficient representing the intervention effect
intervention_samples <- posterior_samples(model, pars = "b_groupintervention")

# Calculate the odds ratio from the posterior samples
odds_ratios <- exp(intervention_samples$b_groupintervention)

# Calculate the mean, 95% credible interval, and standard deviation for the odds ratio
mean_odds_ratio <- mean(odds_ratios)
ci_lower_odds_ratio <- quantile(odds_ratios, probs = 0.025)
ci_upper_odds_ratio <- quantile(odds_ratios, probs = 0.975)
sd_odds_ratio <- sd(odds_ratios)

# Print the results for odds ratio
cat("Odds Ratio:\n")
cat("  Mean:", mean_odds_ratio, "\n")
cat("  95% Credible Interval: [", ci_lower_odds_ratio, ", ", ci_upper_odds_ratio, "]\n")
cat("  Standard Deviation:", sd_odds_ratio, "\n")
