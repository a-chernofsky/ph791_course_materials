#################################################################
#
# PH791 R Code: for STANDARDIZATION using parametric models (g-formula)
#           for a saturated model (no model assumption)
# author: Ariel Chernofsky
# created: September 2, 2020
#
# Saturated model 
# Goals: to estimate the average causal effect of 
# in the aspirin data using a simulation method that does not require 
# specifying the distribution of P[L=l]
# We will see that the estimate will be identical to one we computed in
# Week 1 when we defined the standardized risk as a weighted average
#
#################################################################

set.seed(30459584)

# load libraries ----------------------------------------------------------

library(readxl)

bp <- read_excel("data/dbp_1.xlsx")

names(bp) <- tolower(names(bp))

# create a blood pressure change variable ---------------------------------

bp$change <- bp$dbp4 - bp$dbp0

bp$asp <- ifelse(bp$trt == "A", 1, 0)

bp$sex <- relevel(factor(bp$sex), ref = "M")


# mean outcome model ------------------------------------------------------

out <- lm(change ~ asp + bad_progn + asp * bad_progn,
          data = bp)


# step 1: Expansion of the dataset ----------------------------------------

onesample <- bp[rep(bp$subject, each = 3),]
onesample$interv <- rep(c(-1, 0, 1), nrow(bp))
onesample$change <- ifelse(onesample$interv == -1, onesample$change, NA)
onesample$asp <- ifelse(onesample$interv == 0, 0,
                        ifelse(onesample$interv == 1, 1, onesample$asp))


# step 2: outcome model and prediction ------------------------------------

onesample$predicted_mean <- predict(out, newdata = onesample)

results <- tapply(onesample$predicted_mean, onesample$interv, mean)

ace <- c(EYA_0 = results[["0"]], 
         EYA_1 = results[["1"]], 
         ACE = results[["1"]] - results[["0"]])

ace

# Bootstrap confidence intervals ------------------------------------------

boot_gformula <- function(){
  
  #bootstrap sample
  bindex <- sample(1:nrow(bp), replace = T)
  bsample <- bp[bindex, ] 
  
  # step 1
  onesample <- bsample[rep(bsample$subject, each = 3),]
  onesample$interv <- rep(c(-1, 0, 1), nrow(bsample))
  onesample$change <- ifelse(onesample$interv == -1, 
                             onesample$change, NA)
  onesample$asp <- ifelse(onesample$interv == 0, 0,
                          ifelse(onesample$interv == 1, 1, onesample$asp))
  
  # step 2
  
  bsout <- lm(change ~ asp + bad_progn + asp * bad_progn,
            data = bsample)
  
  onesample$predicted_mean <- predict(bsout, newdata = onesample)
  
  results <- tapply(onesample$predicted_mean, onesample$interv, mean)
  
  ace <- c(EYA_0 = results[["0"]], 
           EYA_1 = results[["1"]], 
           ACE = results[["1"]] - results[["0"]])
  ace 
}

#Number of bootstrap samples
B <- 1000

#initialize dataframe to store results from bootstrap
boot <- data.frame(sample = 1:B, EYA_0 = rep(NA, B), 
                   EYA_1 = rep(NA, B), 
                   ACE = rep(NA, B))

#run repeated bootstrap samples and estimates
for(i in 1:B){
  boot[i,2:4] <- boot_gformula()
}

c(EYA_0 = mean(boot$EYA_0), 
  quantile(boot$EYA_0, probs = 0.025),
  quantile(boot$EYA_0, probs = 0.975))

c(EYA_1 = mean(boot$EYA_1), 
  quantile(boot$EYA_1, probs = 0.025),
  quantile(boot$EYA_1, probs = 0.975))

c(ACE = mean(boot$ACE), 
  quantile(boot$ACE, probs = 0.025),
  quantile(boot$ACE, probs = 0.975))


