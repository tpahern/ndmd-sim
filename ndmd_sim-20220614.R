# Simulation: Robustness of NDMD w/ perfect Sp to:
# - small departures from Sp=1
# - stochastic departures from ND
# - study size of 10,000 fixed for all simulations

# Programmed by Ya Tuo and Tom Ahern (02tahern@med.uvm.edu)
# June 2022

# set working directory
wd <- ("/Users/02tahern/OneDrive - UVM Larner College of Medicine/projects/ndmd-sim/results/simulation")
setwd(wd)

# load required packages
library(tidyverse)
library(writexl)

# seed the random number generator
seed <- 2718
set.seed(seed)

# number of iterations desired per unique parameter set
niter <- 1e4

# maximum number of iterations to attempt per parameter set
niter.max <- 2e5

# set the cohort size
total <- 1e4

# parameter vectors
Pe.vals   <- c(0.10, 0.50) # exposure prevalence
IPRt.vals <- c(1.0, 1.5, 2.0, 5.0) # true risk ratio
IP0.vals  <- c(0.01, 0.10) # risk in unexposed group
Se.vals   <- c(1.0, 0.8, 0.5) # sensitivity for outcome
#Sp.vals   <- c(0.990, 0.992, 0.994, 0.996, 0.998, 1.000) # specificity for outcome
Sp.vals <- c(1.000, 0.998, 0.996, 0.994, 0.992, 0.990) # specificity for outcome

# calculate the number of unique parameter sets
nParmSets <- length(IP0.vals)*length(Pe.vals)*length(IPRt.vals)*length(Sp.vals)*length(Se.vals)

# populate the simulation parameter matrix
sim.parms <- matrix(NA,nrow=nParmSets,ncol=6)
colnames(sim.parms) <- c('scenario', 'Pe', 'IPRt', 'IP0', 'Se', 'Sp')
parmset <- 0
for (i1 in 1:length(Pe.vals)) {
  for (i2 in 1:length(IPRt.vals)) {
    for (i3 in 1:length(IP0.vals)) {
      for (i4 in 1:length(Se.vals)) {
        for (i5 in 1:length(Sp.vals)) {
            parmset <- parmset+1
            sim.parms[parmset,] <- c(parmset, Pe.vals[i1], IPRt.vals[i2],
                                     IP0.vals[i3], Se.vals[i4], Sp.vals[i5])
           }
        }
     }
  }
}

# convert sim.parms to data frame and save as Rdata
sim.parms <- data.frame(sim.parms)
save(sim.parms, file="sim_parms-20220614.Rdata")

# initialize sim.data
n.row <- nParmSets*niter
sim.data <- matrix(NA,nrow=n.row,ncol=26)

# add column names to sim.data
colnames(sim.data) <- c('scenario',
                        'se_set', 'sp_set', 'pe_set', 'ip_set', 'r_set',
                        'flag',
                        'se1_s', 'se0_s', 'sp1_s', 'sp0_s',
                        'se1_d', 'se0_d', 'sp1_d', 'sp0_d', 
                        'delta_sp_s', 'delta_se_s', 'delta_sp_d', 'delta_se_d',
                        'iprc', 'iprm_s', 'iprm_d',
                        'v_s', 'v_d',
                        'obf_s', 'obf_d'
                        )

# loop over the parameter sets and write output to text file
sink("sim-console-20220614.txt")
i.row <- 0
for (parmset in 1:nParmSets) {

   # retrieve parameter values
   ip     <- sim.parms$IP0[parmset]
   pe     <- sim.parms$Pe[parmset]
   r      <- sim.parms$IPRt[parmset]
   se     <- sim.parms$Se[parmset]
   sp     <- sim.parms$Sp[parmset]
   
   # calculate risk in the exposed
   prob   <- min(ip*r,1)

   # compute number of exposed (N1) and unexposed (N0) cases
   N1 <- round(total*pe, digits=0)
   N0 <- total-N1

   # record number of successes and counter per scenario
   nsuccess <- 0
   ctr <- 0

   # want up to 'niter' iterations, but stop after trying 'niter.max' times
   while((nsuccess < niter)&(ctr < niter.max)) {
     
     # uptick the counter
      ctr <- ctr+1
      
     # CORRECTLY CLASSIFIED DATA
     # generate table cells for correctly-classified scenario
     A <- rbinom(1,N1,prob)
     B <- rbinom(1,N0,ip)
     C <- N1-A
     D <- N0-B
     # compute risk ratio for correctly-classified data
     ipr_c <- (A/N1)/(B/N0)

     # STOCHASTIC MISCLASSIFICATION
     # false-negative cases
     FN1s <- A-rbinom(1,A,se)
     FN0s <- B-rbinom(1,B,se)
     # false-positive cases
     FP1s <- rbinom(1,N1-A,1-sp)
     FP0s <- rbinom(1,N0-B,1-sp)
     # compute cell frequencies
     a_s <- A - FN1s + FP1s
     b_s <- B - FN0s + FP0s
     c_s <- N1 - a_s
     d_s <- N0 - b_s
     # compute risk ratio for stochastically misclassified data
     iprm_sto <- (a_s / N1) / (b_s / N0)
     # calculate var[ln(iprm_sto)]
     v_sto <- 1 / a_s - 1 / N1 + 1 / b_s - 1 / N0
     # calculate actual values of Se and Sp
     se_exp_s <- (A - FN1s) / A
     se_unexp_s <- (B - FN0s) / B
     sp_exp_s <- (C - FP1s) / C
     sp_unexp_s <- (D - FP0s) / D
     # calculate degree of non-differentiality for Se and Sp
     delta_sp_sto <- sp_exp_s - sp_unexp_s
     delta_se_sto <- se_exp_s - se_unexp_s
     # calculate overall bias factor
     obf_sto <- iprm_sto / ipr_c

     
     # DETERMINISTIC MISCLASSIFICATION
     # false-negative cases
     FN1d <- A-round(se*A, digits=0)
     FN0d <- B-round(se*B,digits=0)
     # false-positive cases
     FP1d = round((N1 - A) * (1 - sp),digits=0)
     FP0d = round((N0 - B) * (1 - sp),digits=0)
     # compute cell frequencies
     a_d <- A - FN1d + FP1d
     b_d <- B - FN0d + FP0d
     c_d <- N1 - a_d
     d_d <- N0 - b_d
     # compute risk ratio for deterministically misclassified data
     iprm_det <- (a_d / N1) / (b_d / N0)
     # calculate var[ln(iprm_det)]
     v_det <- 1 / a_d - 1 / N1 + 1 / b_d - 1 / N0
     # calculate actual values of Se and Sp
     se_exp_d <- (A - FN1d) / A
     se_unexp_d <- (B - FN0d) / B
     sp_exp_d <- (C - FP1d) / C
     sp_unexp_d <- (D - FP0d) / D
     # calculate degree of non-differentiality for Se and Sp
     delta_sp_det <- sp_exp_d - sp_unexp_d
     delta_se_det <- se_exp_d - se_unexp_d
     # calculate overall bias factor
     obf_det <- iprm_det / ipr_c
     
     
     # check for NaN or Inf results (invalid iterations)
     flag <- min(is.finite(c(ipr_c, iprm_sto, iprm_det, 
                             se_exp_s, se_unexp_s, sp_exp_s, sp_unexp_s, 
                             se_exp_d, se_unexp_d, sp_exp_d, sp_unexp_d, 
                             v_sto, v_det, 
                             obf_sto, obf_det)))
     
     # flag==0 means there is Nan or Inf
     if(flag == 1){
       nsuccess <- nsuccess+1
       i.row <- i.row+1
        
        # add valid results to sim.data
        result <- matrix(c(parmset, 
                           se, sp, pe, ip, r, 
                           flag,
                           se_exp_s, se_unexp_s, sp_exp_s, sp_unexp_s, 
                           se_exp_d, se_unexp_d, sp_exp_d, sp_unexp_d, 
                           delta_sp_sto, delta_se_sto, delta_sp_det, delta_se_det, 
                           ipr_c, iprm_sto, iprm_det, 
                           v_sto, v_det, 
                           obf_sto, obf_det
                           ), nrow=1,ncol=26
                         )
        
        sim.data[i.row,] <- result
     }
   }

   # print message if niter.max reached
   cat("Scenario: ",parmset," - Iterations: ",ctr,"\n")
   flush.console()
   if((nsuccess < niter)&(ctr == niter.max)) {
     cat("Maximum iterations reached with ",nsuccess," successes.\n")
     flush.console()
   }
}

# remove unfilled rows from sim.data
i.row <- i.row+1
if(i.row <= n.row) {
   sim.data <- sim.data[-(i.row:n.row),]
}

# convert sim.data to a data frame
sim.data.df <- data.frame(sim.data)

# generate derived variables
sim.data.df <- sim.data.df |> 
  mutate(
         # confidence limits for misclassified risk ratio
         ul_s = exp(log(iprm_s) + 1.96 * sqrt(v_s)),
         ll_s = exp(log(iprm_s) - 1.96 * sqrt(v_s)),
         ul_d = exp(log(iprm_d) + 1.96 * sqrt(v_d)),
         ll_d = exp(log(iprm_d) - 1.96 * sqrt(v_d)),
         
         # log bias factors (base 10)
         logobf_s = log10(obf_s),
         logobf_d = log10(obf_d),
         
         # indicator: observed true risk ratio within interval of misclassified risk ratio
         iprc_in_ci_s = ifelse(ll_s <= iprc & iprc <= ul_s,1,0),
         iprc_in_ci_d = ifelse(ll_d <= iprc & iprc <= ul_d,1,0),
         
         # indicator: set true risk ratio within interval of misclassified risk ratio
         iprt_in_ci_s = ifelse(ll_s <= r_set & r_set <= ul_s,1,0),  
         iprt_in_ci_d = ifelse(ll_d <= r_set & r_set <= ul_d,1,0),
         
         # indicator: misclassified risk ratio greater than set true risk ratio
         iprm_gt_iprt_s = ifelse(iprm_s > r_set,1,0),
         iprm_gt_iprt_d = ifelse(iprm_d > r_set,1,0),
         
         # indicator: misclassified risk ratio greater than observed true risk ratio
         iprm_gt_iprc_s = ifelse(iprm_s > iprc,1,0),
         iprm_gt_iprc_d = ifelse(iprm_d > iprc,1,0),
         
         # indicator: misclassified risk ratio between null and set true risk ratio
         iprm_btw1_iprt_s = ifelse(1 <= iprm_s & iprm_s <= r_set,1,0),
         iprm_btw1_iprt_d = ifelse(1 <= iprm_d & iprm_d <= r_set,1,0),
         
         # indicator: misclassified risk ratio between null and observed true risk ratio
         iprm_btw1_iprc_s = ifelse(1 <= iprm_s & iprm_s <= iprc,1,0),
         iprm_btw1_iprc_d = ifelse(1 <= iprm_d & iprm_d <= iprc,1,0),
         
         # indicator: misclassified risk ratio less than null
         iprm_lt_1_s = ifelse(iprm_s < 1,1,0),
         iprm_lt_1_d = ifelse(iprm_d < 1,1,0),
         
         # define set specificity as difference from 1
         sp_setm1=round(1-sp_set,3), 
         
         # re-scale deltas to unit=thousandths for regression
         delta_se_s_1k=delta_se_s*1000, 
         delta_sp_s_1k=delta_sp_s*1000,
         
         # checksums for mutually exclusive indicator sets
         checksum_s=iprm_lt_1_s+iprm_btw1_iprt_s+iprm_gt_iprt_s,
         checksum_d=iprm_lt_1_d+iprm_btw1_iprt_d+iprm_gt_iprt_d
         )

# write complete simulation output to Rdata file
save(sim.data.df, file="ndmd_sim-20220614.Rdata")
rm(list=ls())
sink()