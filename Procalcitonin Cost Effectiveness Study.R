###Intensive care code###

rm(list = ls())
#dev.off()
devtools::source_url("https://github.com/mbounthavong/Decision_Analysis/blob/master/tornado_diagram_code.R?raw=TRUE")
par(cex.main=1)
library(ggplot2)
library(dplyr)

####SECTION 1 - Input values####
#Standard care sample size
SC_n <- 2230

#Standard care mortality risk
Mort <- 0.237
#Standard care mortality odds
MortOdds <- Mort/(1-Mort)
#Standard care mortality risk, SE
Mort_SE <- sqrt(Mort*(1-Mort)/SC_n)
#Standard care mortality risk, beta distribution parameters
Mort_alpha <- ((1 - Mort) / Mort_SE) - (1 / Mort) * (Mort ** 2)
Mort_beta <- Mort_alpha * (1 / Mort - 1)
#Standard care mortality risk, 95% confidence interval
Mort_95LL <- Mort - qnorm(0.975)*Mort_SE
Mort_95UL <- Mort + qnorm(0.975)*Mort_SE

#Mortality adjusted OR
MortOR <- 0.89
#Mortality adjusted OR, 95% confidence interval
MortOR_95LL <- 0.80
MortOR_95UL <- 0.99
#Mortality adjusted OR, standard error
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
#Antibiotic duration (days), SD
dAbx_SD <- 9.7
#Antibiotic duration (days), SE
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
#Difference in antibiotic duration (days), margin of error and estimated standard deviation
diffAbx_MoE <- (diffAbx_95UL- diffAbx_95LL)/2
diffAbx_SD <- diffAbx_MoE/qnorm(0.975)
diffAbx_VAR <- diffAbx_SD^2

#PCT antibiotic duration (days)
PCTdAbx <- dAbx + diffAbx

#Total length of stay (days)
LOS <- 28.6
#Total length of stay (days), SD
LOS_SD <- 27.9
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
cPCT <- 17.81
#Number of PCT tests
nPCT <- 7

#Proportion of female sepsis patients
SepPF <- 0.45
#Proportion of female sepsis patients, sample size
SepPF_n <- 197142
#Proportion of female sepsis patients, SE
SepPF_SE <- sqrt(SepPF*(1-SepPF)/SepPF_n)
#Proportion of female sepsis patients, Beta distribution parameters
SepPF_alpha <- ((1 - SepPF) / SepPF_SE) - (1 / SepPF) * (SepPF ** 2)
SepPF_beta <- SepPF_alpha * (1 / SepPF - 1)
#Proportion of female sepsis patients, 95% confidence interval
SepPF_95LL <- SepPF - qnorm(0.975)*SepPF_SE
SepPF_95UL <- SepPF + qnorm(0.975)*SepPF_SE
#Proportion of male sepsis patients
SepPM <- 1-SepPF

#Female QALY 60-64 years old
u60_64F <- 0.776
#Female QALY 60-64 years old, sample size
u60_64F_n <- 608
#Female QALY 60-64 years old, SE
u60_64F_SE <- sqrt(u60_64F*(1-u60_64F)/u60_64F_n)
#Female QALY 60-64 years old, Beta distribution parameters
u60_64F_alpha <- ((1 - u60_64F) / u60_64F_SE - 1 / u60_64F) * (u60_64F ** 2)
u60_64F_beta <- u60_64F_alpha * (1 / u60_64F - 1)
#Female QALY 60-64 years old, 95% confidence interval
u60_64F_95LL <- 0.769
u60_64F_95UL <- 0.797

#Male QALY 60-64 years old
u60_64M <- 0.803
#Male QALY 60-64 years old, sample size
u60_64M_n <- 532
#Standard care mortality rate, SE
u60_64M_SE <- sqrt(u60_64M*(1-u60_64M)/u60_64M_n)
#Standard care mortality rate, Beta distribution parameters
u60_64M_alpha <- ((1 - u60_64M) / u60_64M_SE - 1 / u60_64M) * (u60_64M ** 2)
u60_64M_beta <- u60_64F_alpha * (1 / u60_64F - 1)
#Male QALY 60-64 years old, 95% confidence interval
u60_64M_95LL <- 0.798
u60_64M_95UL <- 0.822

#Weighted QALY 60-64 years old
u60_64W <- (u60_64F*SepPF) + (u60_64M*SepPM)

#Sepsis utility
uSep <- 0.53
#Sepsis utility, sample size
uSep_n <- 701
#Sepsis utility, SE
uSep_SE <- sqrt(uSep*(1-uSep)/uSep_n)
#Sepsis utility, Beta distribution parameters
uSep_alpha <- ((1 - uSep) / uSep_SE - 1 / uSep) * (uSep ** 2)
uSep_beta <- uSep_alpha * (1 / uSep - 1)
#Sepsis utility, 95% confidence interval
uSep_95LL <- uSep - qnorm(0.975)*uSep_SE
uSep_95UL <- uSep + qnorm(0.975)*uSep_SE

#Sepsis 30 day disutility
duSep <- (u60_64W-uSep)/12
#Hospital disutility (QALD/day)
duHosD <- 0.1
duHosD_n <- 159
duHosD_SE <- sqrt(duHosD*(1-duHosD)/duHosD_n)
duHosD_95LL <- duHosD - qnorm(0.975)*duHosD_SE
duHosD_95UL <- duHosD + qnorm(0.975)*duHosD_SE
duHos_alpha <- ((1 - duHosD) / duHosD_SE) - (1 / duHosD) * (duHosD ** 2)
duHos_beta <- duHos_alpha * (1 / duHosD - 1)

#Antibiotic disutility (QALD/day)
duAbxD <- 0.057
duAbxD_n <- 349
duAbxD_SE <- sqrt(duAbxD*(1-duAbxD)/duAbxD_n)
duAbxD_95LL <- duAbxD - qnorm(0.975)*duAbxD_SE
duAbxD_95UL <- duAbxD + qnorm(0.975)*duAbxD_SE

#Total hospital disutility (QALY)
duHos <- duHosD*LOS/365.25
#Total antibiotic disutility (QALY)
duAbx <- duAbxD*dAbx/365.25
#PCT total hospital disutility (QALY)
PCTduHos <- duHosD*PCTLOS/365.25
#PCT total antibiotic disutility (QALY)
PCTduAbx <- duAbxD*PCTdAbx/365.25

#Female QALE 63 years old, 3.5% discount
u63F <- 11.64
#Male QALE 63 years old, 3.5% discount
u63M <- 11.24
#Weighted QALE 63 years old, 3.5% discount
u63W <- (u63F*SepPF) + (u63M*SepPM)

#Female LE 63 years old
le63F <- 23
#Male LE 63 years old
le63M <- 20.51
#Weighted LE 63 years old
le63W <- (le63F*SepPF) + (le63M*SepPM)

##SECTION 2 - input to allow values to be used in the model

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
  u60_64F <- u60_64F,
  u60_64M <- u60_64M,
  u60_64W <- (u60_64F*SepPF) + (u60_64M*SepPM),
  uSep <- uSep,
  duSep <- (u60_64W-uSep)/12,
  #31
  duHos <- duHosD*LOS/365.25,
  duAbx <- duAbxD*dAbx/365.25,
  PCTduHos <- duHosD*PCTLOS/365.25,
  PCTduAbx <- duAbxD*PCTdAbx/365.25,
  u63W <- (u63F*SepPF) + (u63M*SepPM),
  #36
  le63W <- (le63F*SepPF) + (le63M*SepPM)
)

#### SECTION 3 - Decision tree function 
dec_tree <- function(params){
  with(
    as.list(params), 
    {
      
      Mort = Mort
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
      u60_64F <- u60_64F
      u60_64M <- u60_64M
      u60_64W <- (u60_64F*SepPF) + (u60_64M*SepPM)
      uSep <- uSep
      duSep <- (u60_64W-uSep)/12
      duHos <- duHosD*LOS/365.25
      duAbx <- duAbxD*dAbx/365.25
      PCTduHos <- duHosD*PCTLOS/365.25
      PCTduAbx <- duAbxD*PCTdAbx/365.25
      u63W <- (u63F*SepPF) + (u63M*SepPM)
      le63W <- (le63F*SepPF) + (le63M*SepPM)
      
      #Standard care costs, hospital perspective
      hC_SC <- (cAbx*dAbx) + (cHos*Hos) + (cICU*ICU)
      #PCT costs, hospital perspective
      hC_PCT <- (cPCT*nPCT) + (cAbx*PCTdAbx) + (cHos*PCTHos) + (cICU*PCTICU)
      
      #Standard care QALYs gained
      Q_SC <- (u63W-duSep-duAbx-duHos)*Surv
      #PCT QALYs gained
      Q_PCT <- (u63W-duSep-PCTduAbx-PCTduHos)*PCTSurv
      
      #Standard care LYs gained
      LY_SC <- le63W*Surv
      #PCT LYs gained
      LY_PCT <- le63W*PCTSurv
      
      #Antibiotics avoided
      IAbx <- dAbx-PCTdAbx
      
      baseepi <- c(Mort, PCTMort, LOS, PCTLOS, ICU, PCTICU, Hos, PCTHos)
      C <- c(hC_SC, hC_PCT)
      QALY <- c(Q_SC, Q_PCT, LY_SC, LY_PCT)
      AbxU <- c(dAbx, PCTdAbx)
      hIC <- hC_PCT - hC_SC
      IQALY <- Q_PCT - Q_SC
      ILY <- LY_PCT - LY_SC
      IAbx <- dAbx-PCTdAbx
      hICER <- (hC_PCT - hC_SC)/(Q_PCT - Q_SC)
      hleICER <- (hC_PCT - hC_SC)/(LY_PCT - LY_SC)
      haICER <- (hC_PCT - hC_SC)/(dAbx-PCTdAbx)
      vAbxH <- (hIC-(IQALY*20000))/IAbx
      
      names(baseepi) <- paste("baseepi", c("Mort", "PCTMort", "LOS", "PCTLOS", "ICU", "PCTICU", "Hos", "PCTHos"), sep = "_")
      names(C) <- paste("C", c("hC_SC","hC_PCT"), sep = "_")
      names(QALY) <- paste("QALY", c("Q_SC","Q_PCT", "LY_SC", "LY_PCT"), sep = "_")
      names (AbxU) <- paste("Antibiotic use", c("dAbx", "PCTdAbx"), sep = "_")
      names(hIC)   <- paste("Healthcare Incremental Costs")
      names(IQALY)   <- paste("Incremental QALYs")
      names(ILY) <- paste("Incremental LYs")
      names(IAbx) <- paste("Incremental Antibiotics")
      names(hICER) <- paste("Healthcare ICER, QALY")
      names(hleICER) <- paste("Healthcare ICER, LY")
      names(haICER) <- paste("Cost per antibiotic day avoided, healthcare")
      names(vAbxH) <- paste("Antibiotic value, NHS")
      
      return(c(baseepi, C, QALY, AbxU, hIC, IQALY, ILY, IAbx, hICER, hleICER, haICER, vAbxH))
    }
  )
}
options(scipen=999)
dec_tree(input)


#### SECTION 5 - Tornado Plot 
########################
#Define ranges
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
cAbx_range <- c(BaseCase = cAbx, low = (cAbx*0.50), high = (cAbx*1.50))
cPCT_range <- c(BaseCase = cPCT, low = 10, high = 90)
nPCT_range <- c(BaseCase = nPCT, low = 1, high = 10)
SepPF_range <- c(BaseCase = SepPF, low = SepPF_95LL, high = SepPF_95UL)
u60_64F_range <- c(BaseCase = u60_64F, low = u60_64F_95LL, high = u60_64F_95UL)
u60_64M_range <- c(BaseCase = u60_64M, low = u60_64M_95LL, high = u60_64M_95UL)
uSep_range <- c(BaseCase = uSep, low = uSep_95LL, high = uSep_95UL)


####ICER Health system####
#Parameter names
paramNames <- c(
  "30-day mortality probability",
  "30-day mortality OR",
  
  "Difference in total length of stay",
  "Difference in ICU length of stay",
  "Cost per day on regular ward",
  "Cost per day in ICU",
  
  "Cost per PCT test",
  "Number of PCT tests used"
  
)

## List of inputs 
l.tor.in <- vector("list", 8)
names(l.tor.in) <- paramNames

l.tor.in$`30-day mortality probability`   <- cbind(Mort  = Mort_range,    input[-1])
l.tor.in$`30-day mortality OR`   <- cbind(MortOR  = MortOR_range,    input[-3])

l.tor.in$`Difference in total length of stay`  <- cbind(diffLOS = diffLOS_range,   input[-12])
l.tor.in$`Difference in ICU length of stay`  <- cbind(diffICU = diffICU_range,   input[-15])
l.tor.in$`Cost per day on regular ward`  <- cbind(cHos = cHos_range,   input[-19])
l.tor.in$`Cost per day in ICU`  <- cbind(cICU = cICU_range,   input[-20])

l.tor.in$`Cost per PCT test`  <- cbind(cPCT = cPCT_range,   input[-22])
l.tor.in$`Number of PCT tests used`  <- cbind(nPCT = nPCT_range,   input[-23])

## List of outputs
l.tor.out <- vector("list", 8)
names(l.tor.out) <- paramNames

for(i in 1:8){
  l.tor.out[[i]] <- t(apply(l.tor.in[[i]], 1, dec_tree))[ , 16] 
}

## Data structure: ymean, ymin, ymax
m.tor <- matrix(unlist(l.tor.out), nrow = 8, ncol = 3, byrow = TRUE, 
                dimnames = list(paramNames, c("basecase", "low", "high")))

TornadoPlot(main_title = "", Parms = paramNames, Outcomes = m.tor, 
            outcomeName = "", 
            xlab = "", 
            ylab = "", 
            col1="#3182bd", col2="#6baed6"
)

#Incremental costs, Healthcare
#Parameter names
paramNames <- c(
  
  
  "Difference in duration of antibiotic therapy",
  
  "Difference in total length of stay",
  "Difference in ICU length of stay",
  "Cost per day on regular ward",
  "Cost per ICU day",
  "Cost of antibiotic therapy",
  "Cost of PCT test",
  "Number of PCT tests used"
  
  
  
)

## List of inputs 
l.tor.in <- vector("list", 8)
names(l.tor.in) <- paramNames



l.tor.in$`Difference in duration of antibiotic therapy`  <- cbind(diffAbx = diffAbx_range,   input[-9])

l.tor.in$`Difference in total length of stay`  <- cbind(diffLOS = diffLOS_range,   input[-12])

l.tor.in$`Difference in ICU length of stay`  <- cbind(diffICU = diffICU_range,   input[-15])
l.tor.in$`Cost per day on regular ward`  <- cbind(cHos = cHos_range,   input[-19])
l.tor.in$`Cost per ICU day`  <- cbind(cICU = cICU_range,   input[-20])
l.tor.in$`Cost of antibiotic therapy`  <- cbind(cAbx = cAbx_range,   input[-21])
l.tor.in$`Cost of PCT test`  <- cbind(cPCT = cPCT_range,   input[-22])
l.tor.in$`Number of PCT tests used`  <- cbind(nPCT = nPCT_range,   input[-23])



## List of outputs
l.tor.out <- vector("list", 8)
names(l.tor.out) <- paramNames

for(i in 1:8){
  l.tor.out[[i]] <- t(apply(l.tor.in[[i]], 1, dec_tree))[ , 11] 
}

## Data structure: ymean, ymin, ymax
m.tor <- matrix(unlist(l.tor.out), nrow = 8, ncol = 3, byrow = TRUE, 
                dimnames = list(paramNames, c("basecase", "low", "high")))


TornadoPlot(main_title = "", Parms = paramNames, Outcomes = m.tor, 
            outcomeName = "", 
            xlab = "", 
            ylab = "", 
            col1="#3182bd", col2="#6baed6")
#######################################


#SECTION 6 - Monte Carlo simulation
###################################

# Number of iterations for the simulation
n <- 10000

# Create a data frame to store the results
results <- data.frame(matrix(ncol = 24, nrow = n))
names(results) <- c("Mort", "PCTMort", "LOS", "PCTLOS", "ICU", "PCTICU", "Hos", "PCTHos", "C_hC_SC", "C_hC_PCT", "QALY_Q_SC", "QALY_Q_PCT", "QALY_LY_SC", "QALY_LY_PCT", "Antibiotic use_dAbx", "Antibiotic use_PCTdAbx", "Healthcare Incr Costs", "Incr QALYs", "Incr LYs", "Incr Antibiotics", "Healthcare ICER, QALY", "Healthcare ICER, LY", "Cost per antibiotic day avoided, Healthcare", "Antibiotic value, NHS")

# Perform the simulation
for (i in 1:n) {
  ##### Probabilistic re-parameterisation
  
  Mort <- rbeta(1, Mort_alpha, Mort_beta)
  
  MortOR <- rlnorm(1, log(0.89), (log(0.99)-log(0.8))/(2*qnorm(0.975)))
  
  dAbx <- rgamma(1, dAbx_SHAPE, dAbx_RATE)
  
  diffAbx <- rnorm(1, -1.19, diffAbx_SD)
  
  LOS <- rgamma(1, LOS_SHAPE, LOS_RATE)
  
  diffLOS <- rnorm(1, 0.09, diffLOS_SD)
  
  ICU <- rgamma(1, ICU_SHAPE, ICU_RATE)
  
  diffICU <- rnorm(1, 0.04, diffICU_SD)
  
  SepPF <- rbeta(1, SepPF_alpha, SepPF_beta)
  
  u60_64F <- rbeta (1, u60_64F_alpha, u60_64F_beta)
  
  u60_64M <- rbeta (1, u60_64M_alpha, u60_64M_beta)
  
  uSep <- rbeta (1, uSep_alpha, uSep_beta)
  
  
  
  # Create the parameter vector for this iteration
  iteration_params <- c(Mort = Mort, MortOR = MortOR, dAbx = dAbx, diffAbx = diffAbx, LOS = LOS, diffLOS = diffLOS, ICU = ICU, diffICU= diffICU, SepPF = SepPF, u60_64F = u60_64F, u60_64M = u60_64M, uSep = uSep)
  
  # Calculate the decision tree for this iteration
  iteration_results <- dec_tree(iteration_params)
  
  # Store the results
  results[i,] <- iteration_results
}

# Calculate summary statistics
results_summary <- summary(results)

# View summary statistics
print(results_summary)

#install.packages("openxlsx")
library(openxlsx)
# Write the dataframe to an Excel file
#write.xlsx(results, "ICU results.xlsx")

# For each column in the results dataframe
for (col_name in names(results)) {
  # Calculate the 95% confidence interval
  conf_interval <- quantile(results[[col_name]], probs = c(0.025, 0.975))
  
  # Print the confidence interval
  print(paste("95% confidence interval for", col_name, ":", conf_interval))
}

#ICER plot#
# Calculate mean values for the red dot
mean_costs <- mean(results$`Healthcare Incr Costs`)
mean_qalys <- mean(results$`Incr QALYs`)

# Load necessary libraries
library(ggplot2)

# Generate the ICER plot

ggplot(data = results, aes(x = `Incr QALYs`, y = `Healthcare Incr Costs`)) +
  geom_point(color = 'grey', alpha = 0.5) +  # Grey dots for the Monte Carlo simulation results
  geom_point(aes(x = mean_qalys, y = mean_costs), color = 'red', size = 1.25) +  # Red dot for the mean value
  geom_hline(yintercept = 0, color = "black") +  # Horizontal line at y = 0
  geom_vline(xintercept = 0, color = "black") +  # Vertical line at x = 0
  geom_abline(slope = 20000, intercept = 0, linetype = "dashed", color = "blue") +  # WTP threshold line at £20,000/QALY
  scale_x_continuous(breaks = seq(from = floor(min(results$`Incr QALYs`)), to = ceiling(max(results$`Incr QALYs`)), by = 0.2),
                     labels = function(x) ifelse(x == 0, "", x)) +  # Omitting 0 from X axis labels, increments of 0.2
  scale_y_continuous(breaks = seq(from = floor(min(results$`Healthcare Incr Costs`)/1000) * 1000, 
                                  to = ceiling(max(results$`Healthcare Incr Costs`)/1000) * 1000, 
                                  by = 1000),
                     labels = function(y) ifelse(y == 0, "", scales::number(y))) +  # Omitting 0 from Y axis labels, increments of 1000
  theme_minimal() +  # Minimal theme for a clean look
  theme(
    axis.text.x = element_text(size = 10, margin = margin(t = -10)),  # Moving X axis text closer
    axis.text.y = element_text(size = 10, margin = margin(r = -10)),  # Moving Y axis text closer
    panel.grid = element_blank(),  # Removing grid lines
    axis.line.x = element_blank(),  # Removing default X axis line
    axis.line.y = element_blank(),  # Removing default Y axis line
    axis.ticks = element_blank()  # Removing ticks
  )


###Cost per antibiotic day plot###
# Calculate mean values for the red dot
mean_costs <- mean(results$`Incr Antibiotics`)
mean_qalys <- mean(results$`Incr QALYs`)

# Load necessary libraries
library(ggplot2)

# Generate the ICER plot
ggplot(data = results, aes(x = `Incr QALYs`, y = `Incr Antibiotics`)) +
  geom_point(color = 'grey', alpha = 0.5) +  # Grey dots for the Monte Carlo simulation results
  geom_point(aes(x = mean_qalys, y = mean_costs), color = 'red', size = 1.25) +  # Red dot for the mean value
  geom_hline(yintercept = 0, color = "black") +  # Horizontal line at y = 0
  geom_vline(xintercept = 0, color = "black") +  # Vertical line at x = 0
  scale_x_continuous(labels = function(x) ifelse(x == 0, "", x)) +  # Omitting 0 from X axis labels
  scale_y_continuous(labels = scales::number_format(accuracy = 1)) +  # Natural scaling for Y axis
  theme_minimal() +  # Minimal theme for a clean look
  theme(
    axis.text.x = element_text(size = 10, margin = margin(t = -10)),  # Moving X axis text closer
    axis.text.y = element_text(size = 10, margin = margin(r = -15)),  # Moving Y axis text much closer
    panel.grid = element_blank(),  # Removing grid lines
    axis.line.x = element_blank(),  # Removing default X axis line
    axis.line.y = element_blank(),  # Removing default Y axis line
    axis.ticks = element_blank()  # Removing ticks
  )



# Load necessary libraries
library(ggplot2)

# Extract the incremental costs, incremental QALYs, and incremental antibiotics
incremental_costs <- results$`Healthcare Incr Costs`
incremental_QALYs <- results$`Incr QALYs`
incremental_antibiotics <- results$`Incr Antibiotics`

# Define a sequence of willingness to pay thresholds
willingness_to_pay <- seq(0, 50000, by = 1000)
wtp <- 20000
# Calculate the value per antibiotic day avoided
value_per_antibiotic_day <- sapply(willingness_to_pay, function(wtp) {
  value_for_QALYs <- mean(incremental_QALYs) * wtp
  remaining_value <- mean(incremental_costs) - value_for_QALYs
  value_per_day <- remaining_value / mean(incremental_antibiotics)
  return(value_per_day)
})

value_for_QALYs <- mean(incremental_QALYs) * wtp
remaining_value <- mean(incremental_costs) - value_for_QALYs
value_per_antibiotic_day_at_20000 <- remaining_value / mean(incremental_antibiotics)
print(paste("The value needed per antibiotic day at a willingness to pay threshold of £20,000 is:", value_per_antibiotic_day_at_20000))
# Create a data frame for plotting
data_to_plot <- data.frame(WillingnessToPay = willingness_to_pay, 
                           ValuePerAntibioticDay = value_per_antibiotic_day)

# Plot the data
ggplot(data_to_plot, aes(x = WillingnessToPay, y = ValuePerAntibioticDay)) +
  geom_line() +
  ggtitle("Value per Antibiotic Day Avoided vs Willingness to Pay") +
  xlab("Willingness to Pay per QALY") +
  ylab("Value per Antibiotic Day Avoided")


ICUicers <- results$`Healthcare ICER, QALY`

###Emergency department code###


####SECTION 1 - Input values####
#Standard care sample size
SC_n <- 1638


#Standard care mortality risk
Mort <- 0.038
#Standard care mortality odds
MortOdds <- Mort/(1-Mort)
#Standard care mortality risk, SE
Mort_SE <- sqrt(Mort*(1-Mort)/SC_n)
#Standard care mortality risk, Beta distribution parameters
Mort_alpha <- ((1 - Mort) / Mort_SE) - (1 / Mort) * (Mort ** 2)
Mort_beta <- Mort_alpha * (1 / Mort - 1)
#Standard care mortality risk, 95% confidence interval
Mort_95LL <- Mort - qnorm(0.975)*Mort_SE
Mort_95UL <- Mort + qnorm(0.975)*Mort_SE

#Mortality adjusted OR
MortOR <- 0.91
#Mortality adjusted OR, 95% confidence interval
MortOR_95LL <- 0.63
MortOR_95UL <- 1.33
#Mortality adjusted OR, standard error and distribution
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
#Initiation of antibiotic odds
iAbxOdds <- iAbx/(1-iAbx)
#Initiation of antibiotic risk, SE
iAbx_SE <- sqrt(iAbx*(1-iAbx)/SC_n)
#Standard care mortality risk, Beta distribution parameters
iAbx_alpha <- ((1 - iAbx) / iAbx_SE) - (1 / iAbx) * (iAbx ** 2)
iAbx_beta <- iAbx_alpha * (1 / iAbx - 1)
#Standard care mortality risk, 95% confidence interval
iAbx_95LL <- iAbx - qnorm(0.975)*iAbx_SE
iAbx_95UL <- iAbx + qnorm(0.975)*iAbx_SE

iAbxOR <- 0.49
iAbxOR_95LL <- 0.41
iAbxOR_95UL <- 0.58
PCTiAbxOdds <- iAbxOdds*iAbxOR
PCTiAbx <- PCTiAbxOdds/(PCTiAbxOdds+1)
#Mortality adjusted OR, standard error and distribution
iAbxOR_SE <- (iAbxOR_95UL-iAbxOR)/qnorm(0.975)
#Mortality adjusted OR, margin of error and estimated standard deviation
iAbxOR_MoE <- (iAbxOR_95UL- iAbxOR_95LL)/2
iAbxOR_SD <- iAbxOR_MoE/qnorm(0.975)
iAbxOR_VAR <- iAbxOR_SD^2

#Antibiotic duration (days)
dAbx <- 9.8
#Antibiotic duration (days), SD
dAbx_SD <- 5.4
#Antibiotic duration (days), SE
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
diffAbx_95UL <- -2.05
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
#Total length of stay (days), SE
LOS_SE <- LOS_SD/sqrt(SC_n)
#Total length of stay (days), variance
LOS_VAR <- LOS_SE**2
#Total length of stay (days), Gamma distribution parameters
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
#Difference in total length of stay (days), margin of error and estimated standard deviation
diffLOS_MoE <- (diffLOS_95UL- diffLOS_95LL)/2
diffLOS_SD <- diffLOS_MoE/qnorm(0.975)
diffLOS_VAR <- diffLOS_SD^2

#PCT total length of stay (days)
PCTLOS <- LOS + diffLOS

#ICU length of stay (days)
ICU <- 0
#ICU length of stay (days), SD
ICU_SD <- 0
#ICU length of stay (days), SE
ICU_SE <- ICU_SD/sqrt(SC_n)
#ICU length of stay (days), variance
ICU_VAR <- ICU_SE**2
#ICU length of stay (days), Gamma distribution parameters
ICU_SCALE <- ICU_VAR/ICU
ICU_SHAPE <- ICU/ICU_SCALE
ICU_RATE <- 1/ICU_SCALE
#ICU length of stay (days), 95% confidence interval
ICU_95LL <- ICU - qnorm(0.975)*ICU_SE
ICU_95UL <- ICU + qnorm(0.975)*ICU_SE

#Difference in ICU length of stay (days)
diffICU <- 0
#Difference in ICU length of stay (days), 95% confidence interval
diffICU_95LL <- 0
diffICU_95UL <- 0
#Difference in ICU length of stay (days), margin of error and estimated standard deviation
diffICU_MoE <- (diffICU_95UL- diffICU_95LL)/2
diffICU_SD <- diffICU_MoE/qnorm(0.975)

#PCT ICU length of stay (days)
PCTICU <- ICU + diffICU

#Hospital length of stay (days)
Hos <- LOS - ICU

# PCT hospital length of stay (days)
PCTHos <- PCTLOS - PCTICU

#ED cost (admission)
cED <- 349.57
#Hospital cost (day)
cHos <- 767.31
#ICU cost (day)
cICU <- 2095.76
#Antibiotic cost (day)
cAbx <- 0.29
#PCT cost 16.26
cPCT <- 17.81
#Number of PCT tests
nPCT <- 4

#Proportion of female LRTI patients
SepPF <- 0.49
#Proportion of female LRTI patients, sample size
SepPF_n <- 8010
#Proportion of female sepsis patients, SE
SepPF_SE <- sqrt(SepPF*(1-SepPF)/SepPF_n)
#Proportion of female sepsis patients, Beta distribution parameters
SepPF_alpha <- ((1 - SepPF) / SepPF_SE) - (1 / SepPF) * (SepPF ** 2)
SepPF_beta <- SepPF_alpha * (1 / SepPF - 1)
#Proportion of female sepsis patients, 95% confidence interval
SepPF_95LL <- SepPF - qnorm(0.975)*SepPF_SE
SepPF_95UL <- SepPF + qnorm(0.975)*SepPF_SE
#Proportion of male sepsis patients
SepPM <- 1-SepPF

#Female QALY 70-74 years old
u60_64F <- 0.784
#Female QALY 70-74 years old, sample size
u60_64F_n <- 619
#Female QALY 70-74 years old, SE
u60_64F_SE <- sqrt(u60_64F*(1-u60_64F)/u60_64F_n)
#Female QALY 70-74 years old, Beta distribution parameters
u60_64F_alpha <- ((1 - u60_64F) / u60_64F_SE - 1 / u60_64F) * (u60_64F ** 2)
u60_64F_beta <- u60_64F_alpha * (1 / u60_64F - 1)
#Female QALY 70-74 years old, 95% confidence interval
u60_64F_95LL <- 0.779
u60_64F_95UL <- 0.801

#Male QALY 70-74 years old
u60_64M <- 0.801
#Male QALY 70-74 years old, sample size
u60_64M_n <- 505
#Standard care mortality rate, SE
u60_64M_SE <- sqrt(u60_64M*(1-u60_64M)/u60_64M_n)
#Standard care mortality rate, Beta distribution parameters
u60_64M_alpha <- ((1 - u60_64M) / u60_64M_SE - 1 / u60_64M) * (u60_64M ** 2)
u60_64M_beta <- u60_64F_alpha * (1 / u60_64F - 1)
#Male QALY 70-74 years old, 95% confidence interval
u60_64M_95LL <- 0.794
u60_64M_95UL <- 0.818

#Weighted QALY 60-64 years old
u60_64W <- (u60_64F*SepPF) + (u60_64M*SepPM)

#LRTI utility
uSep <- 0.705
#LRTI utility, sample size
uSep_n <- 349
#LRTI utility, SE
uSep_SE <- sqrt(uSep*(1-uSep)/uSep_n)
#LRTI utility, Beta distribution parameters
uSep_alpha <- ((1 - uSep) / uSep_SE - 1 / uSep) * (uSep ** 2)
uSep_beta <- uSep_alpha * (1 / uSep - 1)
#Sepsis utility, 95% confidence interval
uSep_95LL <- uSep - qnorm(0.975)*uSep_SE
uSep_95UL <- uSep + qnorm(0.975)*uSep_SE

#LRTI 30 day disutility
duSep <- (u60_64W-uSep)/12
#Hospital disutility (QALD/day)
duHosD <- 0.1
#Antibiotic disutility (QALD/day)
duAbx_n <- 349
duAbxD <- 0.057
duAbxD_SE <- sqrt(iAbx*(1-iAbx)/duAbx_n)
duAbxD_alpha <- ((1 - duAbxD) / duAbxD_SE) - (1 / duAbxD) * (duAbxD ** 2)
duAbxD_beta <- duAbxD_alpha * (1 / duAbxD - 1)
#Total hospital disutility (QALY)
duHos <- duHosD*LOS/365.25
#Total antibiotic disutility (QALY)
duAbx <- duAbxD*tAbx/365.25
#PCT total hospital disutility (QALY)
PCTduHos <- duHosD*PCTLOS/365.25
#PCT total antibiotic disutility (QALY)
PCTduAbx <- duAbxD*PCTtAbx/365.25

#Female QALE 70 years old, 3.5% discount
u63F <- 9.25
#Male QALE 63 years old, 3.5% discount
u63M <- 8.87
#Weighted QALE 63 years old, 3.5% discount
u63W <- (u63F*SepPF) + (u63M*SepPM)

#Female LE 70 years old
le63F <- 17.17
#Male LE 70 years old
le63M <- 15.12
#Weighted LE 70 years old
le63W <- (le63F*SepPF) + (le63M*SepPM)


##SECTION 2 - input to allow values to be used in the model

input <- data.frame(
  Mort <- Mort,
  MortOdds <- Mort/(1-Mort),
  MortOR <- MortOR,
  PCTMortOdds <- MortOdds * MortOR,
  PCTMort <- PCTMortOdds/(PCTMortOdds+1),
  Surv <- 1-Mort,
  PCTSurv <- 1-(PCTMort),
  iAbx <- iAbx,
  iAbxOdds <- iAbx/(1-iAbx),
  iAbxOR <- iAbxOR,
  PCTiAbxOdds <- iAbxOdds*iAbxOR,
  PCTiAbx <- PCTiAbxOdds/(PCTiAbxOdds+1),
  dAbx <- dAbx,
  diffAbx <- diffAbx,
  PCTdAbx <- dAbx + diffAbx,
  tAbx <- dAbx*iAbx,
  PCTtAbx <- PCTdAbx*PCTiAbx,
  LOS <- LOS,
  diffLOS <- diffLOS,
  PCTLOS <- LOS + diffLOS,
  ICU <- ICU,
  diffICU <- diffICU,
  PCTICU <- ICU + diffICU,
  Hos <- LOS - ICU,
  PCTHos <- PCTLOS - PCTICU,
  cED <- cED,
  cHos <- cHos,
  cICU <- cICU,
  cAbx <- cAbx,
  cPCT <- cPCT,
  nPCT <- nPCT,
  SepPF <- SepPF,
  SepPM <- 1-SepPF,
  u60_64F <- u60_64F,
  u60_64M <- u60_64M,
  u60_64W <- (u60_64F*SepPF) + (u60_64M*SepPM),
  uSep <- uSep,
  duSep <- (u60_64W-uSep)/12,
  duHos <- duHosD*LOS/365.25,
  duAbx <- duAbxD*tAbx/365.25,
  PCTduHos <- duHosD*PCTLOS/365.25,
  PCTduAbx <- duAbxD*PCTtAbx/365.25,
  u63W <- (u63F*SepPF) + (u63M*SepPM),
  le63W <- (le63F*SepPF) + (le63M*SepPM)
)

#### SECTION 3 - Decision tree function 
dec_tree <- function(params){
  with(
    as.list(params), 
    {
      
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
      ICU <- ICU
      diffICU <- diffICU
      PCTICU <- ICU + diffICU
      Hos <- LOS - ICU
      PCTHos <- PCTLOS - PCTICU
      cED <- cED
      cHos <- cHos
      cICU <- cICU
      cAbx <- cAbx
      cPCT <- cPCT
      nPCT <- nPCT
      SepPF <- SepPF
      SepPM <- 1-SepPF
      u60_64F <- u60_64F
      u60_64M <- u60_64M
      u60_64W <- (u60_64F*SepPF) + (u60_64M*SepPM)
      uSep <- uSep
      duSep <- (u60_64W-uSep)/12
      duHos <- duHosD*LOS/365.25
      duAbx <- duAbxD*tAbx/365.25
      PCTduHos <- duHosD*PCTLOS/365.25
      PCTduAbx <- duAbxD*PCTtAbx/365.25
      u63W <- (u63F*SepPF) + (u63M*SepPM)
      le63W <- (le63F*SepPF) + (le63M*SepPM)
      
      
      #Standard care costs, hospital perspective
      hC_SC <- (cAbx*tAbx) + (cHos*Hos) + (cICU*ICU) +cED
      #PCT costs, hospital perspective
      hC_PCT <- (cPCT*nPCT) + (cAbx*PCTtAbx) + (cHos*PCTHos) + (cICU*PCTICU) + cED
      
      #Standard care QALYs gained
      Q_SC <- (u63W-duSep-duAbx-duHos)*Surv
      #PCT QALYs gained
      Q_PCT <- (u63W-duSep-PCTduAbx-PCTduHos)*PCTSurv
      
      #Standard care LYs gained
      LY_SC <- le63W*Surv
      #PCT LYs gained
      LY_PCT <- le63W*PCTSurv
      
      #Antibiotics avoided
      IAbx <- tAbx-PCTtAbx
      
      baseepi <- c(Mort, PCTMort, iAbx, PCTiAbx, LOS, PCTLOS, ICU, PCTICU, Hos, PCTHos)
      C <- c(hC_SC, hC_PCT)
      QALY <- c(Q_SC, Q_PCT, LY_SC, LY_PCT)
      AbxU <- c(dAbx, PCTdAbx)
      hIC <- hC_PCT - hC_SC
      IQALY <- Q_PCT - Q_SC
      ILY <- LY_PCT - LY_SC
      IAbx <- dAbx-PCTdAbx
      hICER <- (hC_PCT - hC_SC)/(Q_PCT - Q_SC)
      hleICER <- (hC_PCT - hC_SC)/(LY_PCT - LY_SC)
      haICER <- (hC_PCT - hC_SC)/(dAbx-PCTdAbx)
      
      names(baseepi) <- paste("baseepi", c("Mort", "PCTMort", "iAbx", "PCTiAbx", "LOS", "PCTLOS", "ICU", "PCTICU", "Hos", "PCTHos"), sep = "_")
      names(C) <- paste("C", c("hC_SC","hC_PCT"), sep = "_")
      names(QALY) <- paste("QALY", c("Q_SC","Q_PCT", "LY_SC", "LY_PCT"), sep = "_")
      names (AbxU) <- paste("Antibiotic use", c("dAbx", "PCTdAbx"), sep = "_")
      names(hIC)   <- paste("Healthcare Incremental Costs")
      names(IQALY)   <- paste("Incremental QALYs")
      names(ILY) <- paste("Incremental LYs")
      names(IAbx) <- paste("Incremental Antibiotics")
      names(hICER) <- paste("Healthcare ICER, QALY")
      names(hleICER) <- paste("Healthcare ICER, LY")
      names(haICER) <- paste("Cost per antibiotic day avoided, healthcare")
      
      return(c(baseepi, C, QALY, AbxU, hIC, IQALY, ILY, IAbx, hICER, hleICER, haICER))
    }
  )
}
options(scipen=999)
dec_tree(input)

##Sensitivity analysis##
nPCT_range  <- seq(0, 9, length.out=8)    # Perform 11 calculations between 0.00 and 0.100.
nPCT_range
m.owsa.input <- cbind(nPCT = nPCT_range, input[-31])
m.owsa.input
outcomes_TC <- t(apply(m.owsa.input, 1, dec_tree))[ , 11:12]
outcomes_TC

#### SECTION 4 - Tornado Plot 
########################
#Define ranges
Mort_range <- c(BaseCase = Mort, low = Mort_95LL, high = Mort_95UL)
MortOR_range <- c(BaseCase = MortOR, low = MortOR_95LL, high = MortOR_95UL)
iAbx_range <- c(BaseCase = iAbx, low = iAbx_95LL, high = iAbx_95UL)
iAbxOR_range <- c(BaseCase = iAbxOR, low = iAbxOR_95LL, high = iAbxOR_95UL)
dAbx_range <- c(BaseCase = dAbx, low = dAbx_95LL, high = dAbx_95UL)
diffAbx_range <- c(BaseCase = diffAbx, low = diffAbx_95LL, high = diffAbx_95UL)
LOS_range <- c(BaseCase = LOS, low = LOS_95LL, high = LOS_95UL)
diffLOS_range <- c(BaseCase = diffLOS, low = diffLOS_95LL, high = diffLOS_95UL)
cED_range <- c(BaseCase = cED, low = cED*0.75, high = cED*0.75)
cHos_range <- c(BaseCase = cHos, low = (cHos*0.75), high = (cHos*1.25))
cAbx_range <- c(BaseCase = cAbx, low = (cAbx*0.75), high = (cAbx*1.25))
cPCT_range <- c(BaseCase = cPCT, low = 10, high = 90)
nPCT_range <- c(BaseCase = nPCT, low = 1, high = 8)
SepPF_range <- c(BaseCase = SepPF, low = SepPF_95LL, high = SepPF_95UL)

#Parameter names

paramNames <- c(
  "30-day mortality probability",
  "30-day mortality OR",
  "Difference in antibiotic therapy duration",
  "Difference in total length of stay",
  "Cost per day on regular ward",
  "Cost per day of antibiotic therapy",
  "Cost per PCT test",
  "Number of PCT tests used"
  
)

## List of inputs 
l.tor.in <- vector("list", 8)
names(l.tor.in) <- paramNames

l.tor.in$"30-day mortality probability"   <- cbind(Mort  = Mort_range,    input[-1])
l.tor.in$"30-day mortality OR"   <- cbind(MortOR  = MortOR_range,    input[-3])
l.tor.in$'Difference in antibiotic therapy duration'  <- cbind(diffAbx = diffAbx_range,   input[-14])
l.tor.in$'Difference in total length of stay'  <- cbind(diffLOS = diffLOS_range,   input[-19])
l.tor.in$`Cost per day on regular ward`  <- cbind(cHos = cHos_range,   input[-27])
l.tor.in$`Cost per day of antibiotic therapy`  <- cbind(cAbx = cAbx_range,   input[-29])
l.tor.in$`Cost per PCT test`  <- cbind(cPCT = cPCT_range,   input[-30])
l.tor.in$`Number of PCT tests used`  <- cbind(nPCT = nPCT_range,   input[-31])

## List of outputs
l.tor.out <- vector("list", 8)
names(l.tor.out) <- paramNames

for(i in 1:8){
  l.tor.out[[i]] <- t(apply(l.tor.in[[i]], 1, dec_tree))[ , 14] 
}

## Data structure: ymean, ymin, ymax
m.tor <- matrix(unlist(l.tor.out), nrow = 8, ncol = 3, byrow = TRUE, 
                dimnames = list(paramNames, c("basecase", "low", "high")))

TornadoPlot(main_title = "", Parms = paramNames, Outcomes = m.tor, 
            outcomeName = "", 
            xlab = "", 
            ylab = "", 
            col1="#3182bd", col2="#6baed6")


#SECTION 6 - Monte Carlo simulation
###################################

# Number of iterations for the simulation
n <- 10000

# Create a data frame to store the results
results <- data.frame(matrix(ncol = 25, nrow = n))
names(results) <- c("Mort", "PCTMort", "iAbx", "PCTiAbx", "LOS", "PCTLOS", "ICU", "PCTICU", "Hos", "PCTHos", "C_hC_SC", "C_hC_PCT", "QALY_Q_SC", "QALY_Q_PCT", "QALY_LY_SC", "QALY_LY_PCT", "AbxU_dAbx", "AbxU_PCTdAbx", "Healthcare Incremental Costs", "Incremental QALYs", "Incremental LYs", "Incremental Antibiotics", "Healthcare ICER, QALY", "Healthcare ICER, LY", "Cost per antibiotic day avoided, Healthcare")

# Perform the simulation
for (i in 1:n) {
  ##### Probabilistic re-parameterisation
  
  Mort <- rbeta(1, Mort_alpha, Mort_beta)
  
  MortOR <- rlnorm(1, log(0.91), (log(1.33)-log(0.63))/(2*qnorm(0.975)))
  
  iAbx <- rbeta(1, iAbx_alpha, iAbx_beta)
  
  iAbxOR <- rlnorm(1, log(0.49), (log(0.58)-log(0.41))/(2*qnorm(0.975)))
  
  dAbx <- rgamma(1, dAbx_SHAPE, dAbx_RATE)
  
  diffAbx <- rnorm(1, -0.52, diffAbx_SD)
  
  LOS <- rgamma(1, LOS_SHAPE, LOS_RATE)
  
  diffLOS <- rnorm(1, -0.14, diffLOS_SD)
  
  SepPF <- rbeta(1, SepPF_alpha, SepPF_beta)
  
  u60_64F <- rbeta (1, u60_64F_alpha, u60_64F_beta)
  
  u60_64M <- rbeta (1, u60_64M_alpha, u60_64M_beta)
  
  uSep <- rbeta (1, uSep_alpha, uSep_beta)
  
  # Create the parameter vector for this iteration
  iteration_params <- c(Mort = Mort, MortOR = MortOR, iAbx = iAbx, iAbxOR = iAbxOR, dAbx = dAbx, diffAbx = diffAbx, LOS = LOS, diffLOS = diffLOS, SepPF = SepPF, u60_64F = u60_64F, u60_64M = u60_64M, uSep = uSep)
  
  # Calculate the decision tree for this iteration
  iteration_results <- dec_tree(iteration_params)
  
  # Store the results
  results[i,] <- iteration_results
}

# Calculate summary statistics
results_summary <- summary(results)

# View summary statistics
print(results_summary)
range(results$'Healthcare ICER, QALY')

library(openxlsx)
# Write the dataframe to an Excel file
#write.xlsx(results, "ED results.xlsx")

# Create a histogram of the ICER values
hist(results$'Healthcare ICER, QALY', breaks = 100, main = "Distribution of ICER values", 
     xlab = "Healthcare ICER, QALY", col = "blue")

q1 <- quantile(results$'Healthcare ICER, QALY', 0.01)
q99 <- quantile(results$'Healthcare ICER, QALY', 0.99)

subset <- results$'Healthcare ICER, QALY'[results$'Healthcare ICER, QALY' > q1 & results$'Healthcare ICER, QALY' < q99]

hist(subset, breaks = 100, main = "PCT testing for LRTI", xlab = "ICER, Health system perspective", col = "blue")

# For each column in the results dataframe
for (col_name in names(results)) {
  # Calculate the 95% confidence interval
  conf_interval <- quantile(results[[col_name]], probs = c(0.025, 0.975))
  
  # Print the confidence interval
  print(paste("95% confidence interval for", col_name, ":", conf_interval))
}

# Save the results to a CSV file
write.csv(results, file = "ED_results.csv", row.names = FALSE)

# Save the summary statistics to another CSV file
write.csv(results_summary, file = "ED_summary.csv", row.names = FALSE)


####
#ICER plot#
# Calculate mean values for the red dot
mean_costs <- mean(results$`Healthcare Incremental Costs`)
mean_qalys <- mean(results$`Incremental QALYs`)

# Load necessary libraries
library(ggplot2)

# Generate the ICER plot
ggplot(data = results, aes(x = `Incremental QALYs`, y = `Healthcare Incremental Costs`)) +
  geom_point(color = 'grey', alpha = 0.5) +  # Grey dots for the Monte Carlo simulation results
  geom_point(aes(x = mean_qalys, y = mean_costs), color = 'red', size = 1.25) +  # Red dot for the mean value
  geom_hline(yintercept = 0, color = "black") +  # Horizontal line at y = 0
  geom_vline(xintercept = 0, color = "black") +  # Vertical line at x = 0
  geom_abline(slope = 20000, intercept = 0, linetype = "dashed", color = "blue") +  # WTP threshold line at £20,000/QALY
  scale_x_continuous(breaks = seq(from = floor(min(results$`Incremental QALYs`)), to = ceiling(max(results$`Incremental QALYs`)), by = 0.1),
                     labels = function(x) ifelse(x == 0, "", x)) +  # Omitting 0 from X axis labels, increments of 0.1
  scale_y_continuous(breaks = seq(from = floor(min(results$`Healthcare Incremental Costs`)/500) * 500, 
                                  to = ceiling(max(results$`Healthcare Incremental Costs`)/500) * 500, 
                                  by = 500),
                     labels = function(y) ifelse(y == 0, "", scales::number(y))) +  # Omitting 0 from Y axis labels, increments of 500
  theme_minimal() +  # Minimal theme for a clean look
  theme(
    axis.text.x = element_text(size = 10, margin = margin(t = -10)),  # Moving X axis text closer
    axis.text.y = element_text(size = 10, margin = margin(r = -15)),  # Moving Y axis text much closer
    panel.grid = element_blank(),  # Removing grid lines
    axis.line.x = element_blank(),  # Removing default X axis line
    axis.line.y = element_blank(),  # Removing default Y axis line
    axis.ticks = element_blank()  # Removing ticks
  )


###Cost per antibiotic day###
# Calculate mean values for the red dot
mean_costs <- mean(results$`Incremental Antibiotics`)
mean_qalys <- mean(results$`Incremental QALYs`)

# Load necessary libraries
library(ggplot2)

# Generate the ICER plot
ggplot(data = results, aes(x = `Incremental QALYs`, y = `Incremental Antibiotics`)) +
  geom_point(color = 'grey', alpha = 0.5) +  # Grey dots for the Monte Carlo simulation results
  geom_point(aes(x = mean_qalys, y = mean_costs), color = 'red', size = 1.25) +  # Red dot for the mean value
  geom_hline(yintercept = 0, color = "black") +  # Horizontal line at y = 0
  geom_vline(xintercept = 0, color = "black") +  # Vertical line at x = 0
  scale_x_continuous(labels = function(x) ifelse(x == 0, "", x)) +  # Omitting 0 from X axis labels
  scale_y_continuous(labels = scales::number_format(accuracy = 1)) +  # Natural scaling for Y axis
  theme_minimal() +  # Minimal theme for a clean look
  theme(
    axis.text.x = element_text(size = 10, margin = margin(t = -10)),  # Moving X axis text closer
    axis.text.y = element_text(size = 10, margin = margin(r = -15)),  # Moving Y axis text much closer
    panel.grid = element_blank(),  # Removing grid lines
    axis.line.x = element_blank(),  # Removing default X axis line
    axis.line.y = element_blank(),  # Removing default Y axis line
    axis.ticks = element_blank()  # Removing ticks
  )

EDicers <- results$`Healthcare ICER, QALY`

###General practice code###


####SECTION 1 - Input values####
#Standard care sample size
SC_n <- 501

#Initiation of antibiotic risk 
iAbx <- 0.631
#Initiation of antibiotic odds
iAbxOdds <- iAbx/(1-iAbx)
#Initiation of antibiotic risk, SE
iAbx_SE <- sqrt(iAbx*(1-iAbx)/SC_n)
#Standard care mortality risk, Beta distribution parameters
iAbx_alpha <- ((1 - iAbx) / iAbx_SE) - (1 / iAbx) * (iAbx ** 2)
iAbx_beta <- iAbx_alpha * (1 / iAbx - 1)
#Standard care mortality risk, 95% confidence interval
iAbx_95LL <- iAbx - qnorm(0.975)*iAbx_SE
iAbx_95UL <- iAbx + qnorm(0.975)*iAbx_SE

iAbxOR <- 0.13
iAbxOR_95LL <- 0.09
iAbxOR_95UL <- 0.18
#Mortality adjusted OR, standard error and distribution
iAbxOR_SE <- (iAbxOR_95UL-iAbxOR)/qnorm(0.975)
#Mortality adjusted OR, margin of error and estimated standard deviation
iAbxOR_MoE <- (iAbxOR_95UL- iAbxOR_95LL)/2
iAbxOR_SD <- iAbxOR_MoE/qnorm(0.975)
iAbxOR_VAR <- iAbxOR_SD^2

PCTiAbxOdds <- iAbxOdds*iAbxOR
PCTiAbx <- PCTiAbxOdds/(PCTiAbxOdds+1)

#Antibiotic duration (days)
dAbx <- 7.3
#Antibiotic duration (days), SD
dAbx_SD <- 2.5
#Antibiotic duration (days), SE
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
#Difference in days with restricted activities, margin of error and estimated standard deviation
diffLOS_MoE <- (diffLOS_95UL- diffLOS_95LL)/2
diffLOS_SD <- diffLOS_MoE/qnorm(0.975)
diffLOS_VAR <- diffLOS_SD^2

#PCT days with restricted activities
PCTLOS <- LOS + diffLOS

#Antibiotic cost (day)
cAbx <- 0.38
#PCT cost 16.26
cPCT <- 17.81
#Number of PCT tests
nPCT <- 2
#Cost of GP presentation
cGP <- 42
#Cost of LRTI admission
cAdmH <- (767.31*8.2)+(0.29*9.8)
#cAdmH_95LL <- 6260.68
#cAdmH_95UL <- 7037.02
#cAdm_n <- 1638
#cAdmH_MoE <- (cAdmH_95UL- cAdmH_95LL)/2
#cAdmH_SD <- cAdmH_MoE/qnorm(0.975)

#Risk of LRTI admission
pAdm <- 0.005306746
pAdmOdds <- pAdm/(1-pAdm)
pAdm_SD <- 0.0044185
pAdm_SE <- sqrt(pAdm_SD/SC_n)
pAdm_alpha <- ((1 - pAdm) / pAdm_SE) - (1 / pAdm) * (pAdm ** 2)
pAdm_beta <- pAdm_alpha * (1 / pAdm - 1)
pAdm_95LL <- 0.0004084322
pAdm_95UL <- 0.01658112

pAdmOR <- 1.317062
pAdmOR_95LL <- 0.3921509
pAdmOR_95UL <- 3.32125
pAdmOR_SD <- 0.7659953 
pAdmOR_SE <- sqrt(pAdmOR_SD/SC_n)
PCTpAdmOdds <- pAdmOdds*pAdmOR 
PCTpAdm <- PCTpAdmOdds/(PCTpAdmOdds+1)
pAdmOR_VAR <- pAdmOR_SD^2

#Probability of treatment failure
pFail <- 0.327
#Treatment failure odds
pFailOdds <- pFail/(1-pFail)
#Treatment failure, SE
pFail_SE <- sqrt(pFail*(1-pFail)/SC_n)
#Treatment failure, Beta distribution parameters
pFail_alpha <- ((1 - pFail) / pFail_SE) - (1 / pFail) * (pFail ** 2)
pFail_beta <- pFail_alpha * (1 / pFail - 1)
#Treatment failure, 95% confidence interval
pFail_95LL <- pFail - qnorm(0.975)*pFail_SE
pFail_95UL <- pFail + qnorm(0.975)*pFail_SE

#Treatment failure odds ratio
pFailOR <- 0.96
pFailOR_95LL <- 0.73
pFailOR_95UL <- 1.25
PCTpFailOdds <- pFailOdds*pFailOR
PCTpFail <- PCTpFailOdds/(PCTpFailOdds+1)
#Mortality adjusted OR, standard error and distribution
pFailOR_SE <- (pFailOR_95UL-pFailOR)/qnorm(0.975)
#Mortality adjusted OR, margin of error and estimated standard deviation
pFailOR_MoE <- (pFailOR_95UL- pFailOR_95LL)/2
pFailOR_SD <- pFailOR_MoE/qnorm(0.975)
pFailOR_VAR <- pFailOR_SD^2

#Proportion of female sepsis patients
SepPF <- 0.49
#Proportion of female sepsis patients, sample size
SepPF_n <- 2905
#Proportion of female sepsis patients, SE
SepPF_SE <- sqrt(SepPF*(1-SepPF)/SepPF_n)
#Proportion of female sepsis patients, Beta distribution parameters
SepPF_alpha <- ((1 - SepPF) / SepPF_SE) - (1 / SepPF) * (SepPF ** 2)
SepPF_beta <- SepPF_alpha * (1 / SepPF - 1)
#Proportion of female sepsis patients, 95% confidence interval
SepPF_95LL <- SepPF - qnorm(0.975)*SepPF_SE
SepPF_95UL <- SepPF + qnorm(0.975)*SepPF_SE
#Proportion of male sepsis patients
SepPM <- 1-SepPF

#Female QALY 65-69 years old
u60_64F <- 0.775
#Female QALY 60-64 years old, sample size
u60_64F_n <- 619
#Female QALY 60-64 years old, SE
u60_64F_SE <- sqrt(u60_64F*(1-u60_64F)/u60_64F_n)
#Female QALY 60-64 years old, Beta distribution parameters
u60_64F_alpha <- ((1 - u60_64F) / u60_64F_SE - 1 / u60_64F) * (u60_64F ** 2)
u60_64F_beta <- u60_64F_alpha * (1 / u60_64F - 1)
#Female QALY 60-64 years old, 95% confidence interval
u60_64F_95LL <- 0.770
u60_64F_95UL <- 0.795

#Male QALY 60-64 years old
u60_64M <- 0.797
#Male QALY 60-64 years old, sample size
u60_64M_n <- 568
#Standard care mortality rate, SE
u60_64M_SE <- sqrt(u60_64M*(1-u60_64M)/u60_64M_n)
#Standard care mortality rate, Beta distribution parameters
u60_64M_alpha <- ((1 - u60_64M) / u60_64M_SE - 1 / u60_64M) * (u60_64M ** 2)
u60_64M_beta <- u60_64F_alpha * (1 / u60_64F - 1)
#Male QALY 60-64 years old, 95% confidence interval
u60_64M_95LL <- 0.792
u60_64M_95UL <- 0.818

#Weighted QALY 60-64 years old
u60_64W <- (u60_64F*SepPF) + (u60_64M*SepPM)

#LRTI utility
uSep <- 0.705
#LRTI utility, sample size
uSep_n <- 349
#LRTI utility, SE
uSep_SE <- sqrt(uSep*(1-uSep)/uSep_n)
#LRTI utility, Beta distribution parameters
uSep_alpha <- ((1 - uSep) / uSep_SE - 1 / uSep) * (uSep ** 2)
uSep_beta <- uSep_alpha * (1 / uSep - 1)
#LRTI utility, 95% confidence interval
uSep_95LL <- uSep - qnorm(0.975)*uSep_SE
uSep_95UL <- uSep + qnorm(0.975)*uSep_SE

#Sepsis 30 day disutility
duSep <- (u60_64W-uSep)/12
#Antibiotic disutility (QALD/day)
duAbxD <- 0.057
#Total antibiotic disutility (QALY)
duAbx <- duAbxD*tAbx/365.25
#PCT total antibiotic disutility (QALY)
PCTduAbx <- duAbxD*PCTtAbx/365.25

#Female QALE 66 years old, 3.5% discount
u63F <- 10.66
#Male QALE 66 years old, 3.5% discount
u63M <- 10.23
#Weighted QALE 63 years old, 3.5% discount
u63W <- (u63F*SepPF) + (u63M*SepPM)

#Female LE 66 years old
le63F <- 20.46
#Male LE 66 years old
le63M <- 18.13
#Weighted LE 63 years old
le63W <- (le63F*SepPF) + (le63M*SepPM)

WTP <- 20000

##SECTION 2 - input to allow values to be used in the model

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
  pAdmOdds <- pAdmOdds,
  #21
  pAdmOR <- pAdmOR,
  PCTpAdmOdds <- PCTpAdmOdds,
  PCTpAdm <- PCTpAdm,
  pFail <- pFail,
  pFailOdds <- pFail/(1-pFail),
  #26
  pFailOR <- pFailOR,
  PCTpFailOdds <- pFailOdds*pFailOR,
  PCTpFail <- PCTpFailOdds/(PCTpFailOdds+1),
  SepPF <- SepPF,
  SepPM <- 1-SepPF,
  #31
  u60_64F <- u60_64F,
  u60_64M <- u60_64M,
  u60_64W <- (u60_64F*SepPF) + (u60_64M*SepPM),
  uSep <- uSep,
  duSep <- (u60_64W-uSep)/12,
  #36
  duAbx <- duAbxD*tAbx/365.25,
  PCTduAbx <- duAbxD*PCTtAbx/365.25,
  u63W <- (u63F*SepPF) + (u63M*SepPM),
  le63W <- (le63F*SepPF) + (le63M*SepPM),
  WTP <- WTP
)

#### SECTION 3 - Decision tree function 
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
      pAdmOdds <- pAdmOdds
      pAdmOR <- pAdmOR
      PCTpAdmOdds <- PCTpAdmOdds
      PCTpAdm <- PCTpAdm
      pFail <- pFail
      pFailOdds <- pFail/(1-pFail)
      pFailOR <- pFailOR
      PCTpFailOdds <- pFailOdds*pFailOR
      PCTpFail <- PCTpFailOdds/(PCTpFailOdds+1)
      SepPF <- SepPF
      SepPM <- 1-SepPF
      u60_64F <- u60_64F
      u60_64M <- u60_64M
      u60_64W <- (u60_64F*SepPF) + (u60_64M*SepPM)
      uSep <- uSep
      duSep <- (u60_64W-uSep)/12
      duAbx <- duAbxD*tAbx/365.25
      PCTduAbx <- duAbxD*PCTtAbx/365.25
      u63W <- (u63F*SepPF) + (u63M*SepPM)
      le63W <- (le63F*SepPF) + (le63M*SepPM)
      
      #Standard care costs, health system perspective
      hC_SC <- (cAbx*tAbx) + (cGP+(cGP*pFail)) + (pAdm*cAdmH)
      #PCT costs, hospital perspective
      hC_PCT <- (cPCT*nPCT) + (cAbx*PCTtAbx) + (cGP+(cGP*PCTpFail)) + (PCTpAdm*cAdmH)
      
      #Standard care QALYs gained
      Q_SC <- (u63W-duSep-duAbx)
      #PCT QALYs gained
      Q_PCT <- (u63W-duSep-PCTduAbx)
      
      #Antibiotics avoided
      IAbx <- tAbx-PCTtAbx
      
      C <- c(hC_SC, hC_PCT)
      QALY <- c(Q_SC, Q_PCT)
      AbxU <- c(dAbx, PCTdAbx)
      hIC <- hC_PCT - hC_SC
      IQALY <- Q_PCT - Q_SC
      IAbx <- dAbx-PCTdAbx
      hICER <- (hC_PCT - hC_SC)/(Q_PCT - Q_SC)
      haICER <- (hC_PCT - hC_SC)/(dAbx-PCTdAbx)
      vAbxH <- (hIC-(IQALY*WTP))/IAbx
      vAbxT <- (hIC-(IQALY*WTP))
      
      names(LOS) <- paste("LOS")
      names(PCTLOS) <- paste("PCTLOS")
      names(iAbx) <- paste("iAbx")
      names(PCTiAbx) <- paste("PCTiAbx")
      names(pAdm) <- paste("pAdm")
      names(PCTpAdm) <- paste("PCTpAdm")
      names(pFail) <- paste("pFail")
      names(PCTpFail) <- paste("PCTpFail")
      names(C) <- paste("C", c("hC_SC","hC_PCT"), sep = "_")
      names(QALY) <- paste("QALY", c("Q_SC","Q_PCT"), sep = "_")
      names(AbxU) <- paste("Antibiotic use", c("dAbx", "PCTdAbx"), sep = "_")
      names(hIC)   <- paste("Healthcare Incremental Costs")
      names(IQALY)   <- paste("Incremental QALYs")
      names(IAbx) <- paste("Incremental Antibiotics")
      names(hICER) <- paste("Healthcare ICER, QALY")
      names(haICER) <- paste("Cost per antibiotic day avoided, healthcare")
      names(vAbxH) <- paste("Antibiotic value, NHS")
      names(vAbxT) <- paste("Antibiotic numerator")
      
      return(c(LOS, PCTLOS, iAbx, PCTiAbx, pAdm, PCTpAdm, pFail, PCTpFail, C, QALY, AbxU, hIC, IQALY, IAbx, hICER, haICER, vAbxH, vAbxT))
    }
  )
}
options(scipen=999)
dec_tree(input)

#### SECTION 4 - Tornado Plot 
########################
#Define ranges

iAbx_range <- c(BaseCase = iAbx, low = iAbx_95LL, high = iAbx_95UL)
iAbxOR_range <- c(BaseCase = iAbxOR, low = iAbxOR_95LL, high = iAbxOR_95UL)
dAbx_range <- c(BaseCase = dAbx, low = dAbx_95LL, high = dAbx_95UL)
diffAbx_range <- c(BaseCase = diffAbx, low = diffAbx_95LL, high = diffAbx_95UL)
LOS_range <- c(BaseCase = LOS, low = LOS_95LL, high = LOS_95UL)
diffLOS_range <- c(BaseCase = diffLOS, low = diffLOS_95LL, high = diffLOS_95UL)
cAbx_range <- c(BaseCase = cAbx, low = cAbx*0.75, high = cAbx*1.25)
cPCT_range <- c(BaseCase = cPCT, low = cPCT*0.75, high = cPCT*1.25)
nPCT_range <- c(BaseCase = nPCT, low = 1, high = 3)
cGP_range <- c(BaseCase = cGP, low = cGP*0.75, high = cGP*1.25)
cAdmH_range <- c(BaseCase = cAdmH, low = cAdmH*0.75, high = cAdmH*1.25)
pAdm_range <- c(BaseCase = pAdm, low = pAdm_95LL, high = pAdm_95UL)
pAdmOR_range <- c(BaseCase = pAdmOR, low = pAdmOR_95LL, high = pAdmOR_95UL)
pFail_range <- c(BaseCase = pFail, low = pFail_95LL, high = pFail_95UL)
pFailOR_range <- c(BaseCase = pFailOR, low = pFailOR_95LL, high = pFailOR_95UL)
SepPF_range <- c(BaseCase = SepPF, low = SepPF_95LL, high = SepPF_95UL)

#Parameter names
paramNames <- c(
  
  "Antibiotic initiation probability",
  "Antibiotic initiation OR",
  "Antibiotic therapy duration",
  "Difference in antibiotic therapy duration",
  "Cost per day of antibiotic therapy",
  "Cost per PCT test",
  "Hospital admission probability",
  "Treatment failure OR"
  
)

## List of inputs 
l.tor.in <- vector("list", 8)
names(l.tor.in) <- paramNames


l.tor.in$"Antibiotic initiation probability"   <- cbind(iAbx  = iAbx_range,    input[-1])
l.tor.in$"Antibiotic initiation OR"   <- cbind(iAbxOR  = iAbxOR_range,    input[-3])
l.tor.in$"Antibiotic therapy duration"  <- cbind(dAbx = dAbx_range,   input[-6])
l.tor.in$"Difference in antibiotic therapy duration"  <- cbind(diffAbx = diffAbx_range,   input[-7])
l.tor.in$"Cost per day of antibiotic therapy"  <- cbind(cAbx = cAbx_range,   input[-14])
l.tor.in$"Cost per PCT test"  <- cbind(cPCT = cPCT_range,   input[-15])
l.tor.in$"Hospital admission probability"  <- cbind(pAdm = pAdm_range,   input[-20])
l.tor.in$"Treatment failure OR"  <- cbind(pFailOR = pFailOR_range,   input[-27])



## List of outputs
l.tor.out <- vector("list", 8)
names(l.tor.out) <- paramNames

for(i in 1:8){
  l.tor.out[[i]] <- t(apply(l.tor.in[[i]], 1, dec_tree))[ , 10] 
}

## Data structure: ymean, ymin, ymax
m.tor <- matrix(unlist(l.tor.out), nrow = 8, ncol = 3, byrow = TRUE, 
                dimnames = list(paramNames, c("basecase", "low", "high")))

TornadoPlot(main_title = "", Parms = paramNames, Outcomes = m.tor, 
            outcomeName = "", 
            xlab = "", 
            ylab = "", 
            col1="#3182bd", col2="#6baed6")


#SECTION 5 - Monte Carlo simulation
###################################

# Number of iterations for the simulation
n <- 10000

# Create a data frame to store the results
results <- data.frame(matrix(ncol = 21, nrow = n))
names(results) <- c("LOS", "PCTLOS", "iAbx", "PCTiAbx", "pAdm", "PCTpAdm", "pFail", "PCTpFail", "C_hC_SC", "C_hC_PCT", "QALY_Q_SC", "QALY_Q_PCT", "dAbx", "PCTdAbx", "Healthcare Incremental Costs", "Incremental QALYs", "Incremental Antibiotics", "Healthcare ICER, QALY", "Cost per antibiotic day avoided, Healthcare", "Antibiotic value, NHS", "Antibiotic numerator")

# Perform the simulation
for (i in 1:n) {
  ##### Probabilistic re-parameterisation
  
  
  iAbx <- rbeta(1, iAbx_alpha, iAbx_beta)
  
  iAbxOR <- rlnorm(1, log(0.49), (log(0.58)-log(0.41))/(2*qnorm(0.975)))
  
  dAbx <- rgamma(1, dAbx_SHAPE, dAbx_RATE)
  
  diffAbx <- rnorm(1, -0.52, diffAbx_SD)
  
  LOS <- rgamma(1, LOS_SHAPE, LOS_RATE)
  
  diffLOS <- rnorm(1, -0.14, diffLOS_SD)
  
  pAdm <- rbeta(1, pAdm_alpha, pAdm_beta)
  
  pAdmOR <- rlnorm(1, log(0.988), (log(3.434)-log(0.285))/(2*qnorm(0.975)))
  
  pFail <- rbeta(1, pFail_alpha, pFail_beta)
  
  pFailOR <- rlnorm(1, log(0.96), (log(1.25)-log(0.73))/(2*qnorm(0.975)))
  
  SepPF <- rbeta(1, SepPF_alpha, SepPF_beta)
  
  u60_64F <- rbeta (1, u60_64F_alpha, u60_64F_beta)
  
  u60_64M <- rbeta (1, u60_64M_alpha, u60_64M_beta)
  
  uSep <- rbeta (1, uSep_alpha, uSep_beta)
  
  # Create the parameter vector for this iteration
  iteration_params <- c(iAbx = iAbx, iAbxOR = iAbxOR, dAbx = dAbx, diffAbx = diffAbx, LOS = LOS, diffLOS = diffLOS, pAdm = pAdm, pAdmOR = pAdmOR, pFail = pFail, pFailOR = pFailOR, SepPF = SepPF, u60_64F = u60_64F, u60_64M = u60_64M, uSep = uSep)
  
  # Calculate the decision tree for this iteration
  iteration_results <- dec_tree(iteration_params)
  
  # Store the results
  results[i,] <- iteration_results
}

# Calculate summary statistics
results_summary <- summary(results)

# View summary statistics
print(results_summary)

standard_deviations <- apply(results, 2, sd)
print(standard_deviations)

# For each column in the results dataframe
for (col_name in names(results)) {
  # Calculate the 95% confidence interval
  conf_interval <- quantile(results[[col_name]], probs = c(0.025, 0.975))
  
  # Print the confidence interval
  print(paste("95% confidence interval for", col_name, ":", conf_interval))
}
###########################################################
#ICER plot 

# Calculate mean values for the red dot
mean_costs <- mean(results$`Healthcare Incremental Costs`)
mean_qalys <- mean(results$`Incremental QALYs`)

# Load necessary libraries
library(ggplot2)

# Generate the ICER plot
ggplot(data = results, aes(x = `Incremental QALYs`, y = `Healthcare Incremental Costs`)) +
  geom_point(color = 'grey', alpha = 0.5) +  # Grey dots for the Monte Carlo simulation results
  geom_point(aes(x = mean_qalys, y = mean_costs), color = 'red', size = 1.25) +  # Red dot for the mean value
  geom_abline(intercept = 0, slope = 20000, linetype = "dashed", color = "blue") +  # WTP threshold line
  geom_hline(yintercept = 0, color = "black") +  # Horizontal line at y = 0
  geom_vline(xintercept = 0, color = "black") +  # Vertical line at x = 0
  scale_x_continuous(
    breaks = seq(0.0001, 0.00035, by = 0.0001), 
    limits = c(0, 0.00035), 
    labels = function(x) ifelse(x == 0, "", x)
  ) +  # Labeling X axis in intervals of 0.0001
  scale_y_continuous(breaks = seq(0, 55, by = 5), limits = c(0, 55)) +  # Extending Y axis to 55 and labeling in increments of 5
  theme_minimal() +  # Minimal theme for a clean look
  theme(
    axis.text.x = element_text(size = 10, margin = margin(t = -10)),  # Moving X axis text closer
    axis.text.y = element_text(size = 10, margin = margin(r = -10)),  # Moving Y axis text closer
    panel.grid = element_blank(),  # Removing grid lines
    axis.line.x = element_blank(),  # Removing default X axis line
    axis.line.y = element_blank(),  # Removing default Y axis line
    axis.ticks = element_blank()  # Removing ticks
  )

########
# Calculate mean values for the red dot
mean_costs <- mean(results$`Antibiotic numerator`)
mean_qalys <- mean(results$`Incremental Antibiotics`)

# Load necessary libraries
library(ggplot2)

# Generate the ICER plot
ggplot(data = results, aes(x = `Incremental Antibiotics`, y = `Antibiotic numerator`)) +
  geom_point(color = 'grey', alpha = 0.5) +  # Grey dots for the Monte Carlo simulation results
  geom_point(aes(x = mean_qalys, y = mean_costs), color = 'red', size = 1.25) +  # Red dot for the mean value
  geom_hline(yintercept = 0, color = "black") +  # Horizontal line at y = 0
  geom_vline(xintercept = 0, color = "black") +  # Vertical line at x = 0
  scale_x_continuous(labels = function(x) ifelse(x == 0, "", x)) +  # Omitting 0 from X axis labels
  scale_y_continuous() +  # Natural scaling for Y axis
  theme_minimal() +  # Minimal theme for a clean look
  theme(
    axis.text.x = element_text(size = 10, margin = margin(t = -10)),  # Moving X axis text closer
    axis.text.y = element_text(size = 10, margin = margin(r = -15)),  # Moving Y axis text much closer
    panel.grid = element_blank(),  # Removing grid lines
    axis.line.x = element_blank(),  # Removing default X axis line
    axis.line.y = element_blank(),  # Removing default Y axis line
    axis.ticks = element_blank()  # Removing ticks
  )

###Cost per antibiotic day##
mean_costs <- mean(results$`Healthcare Incremental Costs`)
mean_qalys <- mean(results$`Incremental Antibiotics`)

# Load necessary libraries
library(ggplot2)

# Generate the ICER plot
ggplot(data = results, aes(x = `Incremental Antibiotics`, y = `Healthcare Incremental Costs`)) +
  geom_point(color = 'grey', alpha = 0.5) +  # Grey dots for the Monte Carlo simulation results
  geom_point(aes(x = mean_qalys, y = mean_costs), color = 'red', size = 1.25) +  # Red dot for the mean value
  geom_hline(yintercept = 0, color = "black") +  # Horizontal line at y = 0
  geom_vline(xintercept = 0, color = "black") +  # Vertical line at x = 0
  scale_x_continuous(labels = function(x) ifelse(x == 0, "", x)) +  # Omitting 0 from X axis labels
  scale_y_continuous() +  # Natural scaling for Y axis
  theme_minimal() +  # Minimal theme for a clean look
  theme(
    axis.text.x = element_text(size = 10, margin = margin(t = -10)),  # Moving X axis text closer
    axis.text.y = element_text(size = 10, margin = margin(r = -15)),  # Moving Y axis text much closer
    panel.grid = element_blank(),  # Removing grid lines
    axis.line.x = element_blank(),  # Removing default X axis line
    axis.line.y = element_blank(),  # Removing default Y axis line
    axis.ticks = element_blank()  # Removing ticks
  )

###Value per antibiotic day plot###


# vector of incremental costs (sub in real values from PSA)
vec_c_i = rnorm(1000, 45.22, 2.21149019491)

# vector of incremental qalys (sub in real values from PSA)
vec_q_i = rnorm(1000, 0.00022939, 0.00003229733)

# vector of incremental abx days (sub in real values from PSA)
vec_a_i = rnorm(1000, 0.5233, 0.28473361935)

# function to calculate abx value
f_abx_value = function(cost, # incremental costs
                       qaly, # incremental qalys
                       abx, # incremental abx days
                       cet # cost-effectiveness threshold
){
  value_per_abx_day = (cost - (qaly*cet))/abx
  return(value_per_abx_day)
}

# empty vectors for outputs
vec_value_abx_day = c()
vec_cet = c()

# loop through each run of psa 
for(psa_i in 1:1000){
  print(paste0("on psa simulation ", psa_i, " of 1000"))
  
  # loop through considered CETs
  for(cet_i in seq(1000, 50000, by = 1000)){
    
    # for each run of psa, extract corresponding outcomes from vectors
    c_i = vec_c_i[psa_i]
    q_i = vec_q_i[psa_i]
    a_i = vec_a_i[psa_i]
    
    value_per_abx_day_i = f_abx_value(c_i, q_i, a_i, cet = cet_i)
    
    # fill vectors with results
    vec_cet = append(vec_cet, cet_i)
    vec_value_abx_day = append(vec_value_abx_day, value_per_abx_day_i)
    
  }
}

df_value_abx_day = data.frame(cet = vec_cet,
                              value_abx_day = vec_value_abx_day)

df_value_abx_day_summarized = df_value_abx_day%>%
  group_by(cet)%>%
  summarise(mean = mean(value_abx_day), lower = quantile(value_abx_day, 0.025), upper = quantile(value_abx_day, 0.975))

df_value_abx_day_summarized%>%
  ggplot(aes(x = cet, y = mean, ymin = lower, ymax = upper))+
  geom_line()+
  geom_ribbon(alpha = 0.3)

GPicers <- results$`Healthcare ICER, QALY`

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
ceac_data_gp <- calculate_ceac(GPicers, max_wtp)
ceac_data_gp$Source <- 'GP'

# Calculate CEAC data for EDicers
ceac_data_ed <- calculate_ceac(EDicers, max_wtp)
ceac_data_ed$Source <- 'ED'

# Calculate CEAC data for ICUicers
ceac_data_icu <- calculate_ceac(ICUicers, max_wtp)
ceac_data_icu$Source <- 'ICU'

# Combine all CEAC data into one data frame
combined_ceac_data <- bind_rows(ceac_data_gp, ceac_data_ed, ceac_data_icu)

# Plotting using ggplot2
ggplot() +
  geom_line(data = combined_ceac_data, aes(x = WTP, y = Prob, color = Source, linetype = "dashed"), size = 0.55) +
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
  scale_color_manual(values = c('GP' = "green", 'ED' = "blue", 'ICU' = "red"), breaks = c('GP', 'ED', 'ICU')) +
  scale_linetype_manual(values = c('dashed'))

# Print the data points for specific WTP values for GP, ED, and ICU
print(combined_ceac_data[combined_ceac_data$WTP %in% c(20000, 30000), ])

