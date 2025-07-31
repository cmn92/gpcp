# full_sim_and_plot.R
# ------------------------------
# 1.  Load libraries & C++ code
# ------------------------------
library(AlphaSimR)
library(sommer)
library(AGHmatrix)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Rcpp)
source("~/gpcp_simulation_functions.R")
# point this to wherever you put your CalcCrossMeans.cpp
Rcpp::sourceCpp("~/gpcp/src/CalcCrossMeans.cpp")  

# ------------------------------
# 2. functions
#    - makeCrossPlan()
#    - introduceGenotypingErrors()
#    - computeGeneticDiversity()
#    - runSimulationCycle()
#    - runSimulation()

# ------------------------------
# 5. One‐run pipeline across pop sizes, meanDD, methods & p_sel
# ------------------------------

# Parameters
pop_sizes   <- c(250, 500, 750, 1000)
meanDD_vals <- c(0.01, 1, 2, 4)            # your meanDD
p_list      <- c(0.25, 0.5, 0.75, 1)   # proportion selected
nBurnIn     <- 20
nCycles     <- 50
nReps       <- 2

nChr = 18
segSites = 5400
nQTL=56
# Storage
res_list <- list()

for (N in pop_sizes) {
  # 5.2.1. Generate founders
  founders <- runMacs(nInd = N, nChr = nChr, segSites = segSites)
  
  for (m in meanDD_vals) {
    cat("Population =", N, " meanDD =", m, "\n")
    
    # 5.2.2. Set up SimParam
    SP <- SimParam$new(founders)
    SP$addSnpChip(1000)
    SP$addTraitAD(nQtlPerChr = nQTL, meanDD = m, varDD = 0.001)
    
    # 5.2.3. Single burn-in to get true d/a ratio
    pop0 <- newPop(founders, simParam = SP)
    for (i in seq_len(nBurnIn)) {
      pop0 <- setPheno(pop0, simParam = SP)
      plan <- makeCrossPlan(ids = pop0@id, nProgeny = N)
      pop0 <- makeCross(pop0, plan, simParam = SP)
    }
    varD0   <- varD(pop0)
    varG0   <- varG(pop0)
    da_prop <- if (m == 0) 0 else varD0 / varG0
    
    for (method in c("GEBV", "GPCP")) {
      addDom <- (method == "GPCP")
      cat("Population =", N, " meanDD =", m, "method ==", method, "\n")
      
      for (p in p_list) {
        cat("Population =", N, " meanDD =", m, "method ==", method, "prop==",p, "\n")
        # 5.2.4. Convert p → truncation intensity i
        z         <- qnorm(1 - p)
        sel_int   <- dnorm(z) / p
        sel_label <- paste0("i=", format(round(sel_int, 3), nsmall = 3))
        
        # totalselc <- round(N * p)
        totalselc <- N
        nProgeny <- round(p * 1000)
        
        # 5.2.5. Matrices to hold full (burn-in + selection) replicates
        Ttot    <- nBurnIn + nCycles
        gains   <- matrix(NA, nrow = Ttot, ncol = nReps)
        divs    <- matrix(NA, nrow = Ttot, ncol = nReps)
        
        for (rep in seq_len(nReps)) {
          cat("rep =",rep, "\n")
          # --- fresh burn-in per replicate ---
          pop_rep   <- newPop(founders, simParam = SP)
          burn_gv   <- numeric(nBurnIn)
          burn_div  <- numeric(nBurnIn)
          for (i in seq_len(nBurnIn)) {
            pop_rep   <- setPheno(pop_rep, simParam = SP)
            plan      <- makeCrossPlan(ids = pop_rep@id, nProgeny = N)
            pop_rep   <- makeCross(pop_rep, plan, simParam = SP)
            burn_gv[i]  <- mean(AlphaSimR::gv(pop_rep))
            snpdata     <- pullSnpGeno(pop_rep, simParam = SP)
            burn_div[i] <- mean(computeGeneticDiversity(snpdata)$heterozygosity)
          }
          
          # --- selection cycles via existing runSimulation() ---
          sim <- runSimulation(
            initial_pop         = pop_rep,
            SP                  = SP,
            nCycles             = nCycles,
            selection_intensity = sel_int,
            p=p,
            totalselc           = totalselc,
            nProgeny            = nProgeny,
            addDominance        = addDom
          )
          
          # 5.2.6. Combine burn-in + selection
          gains[, rep] <- c(burn_gv,   sim$genetic_gain)
          divs[,  rep] <- c(burn_div,  sim$diversity$Heterozygosity)
        }
        
        # 5.2.7. Build the result data.frame
        df <- data.frame(
          Population    = factor(N, levels = pop_sizes),
          Method        = factor(method, levels = c("GEBV", "GPCP")),
          d_a           = da_prop,
          p_sel         = p,
          sel_intensity = sel_int,
          sel_label     = sel_label,
          nProgeny = round(p * 1000),
          Cycle         = seq(-nBurnIn + 1, nCycles),
          gain_mean     = rowMeans(gains, na.rm = TRUE),
          gain_sd       = apply(gains, 1, function(x) sd(x, na.rm = TRUE)),
          div_mean      = rowMeans(divs, na.rm = TRUE)
        )
        
        res_list[[paste(N, m, method, p)]] <- df
      }
    }
  }
}