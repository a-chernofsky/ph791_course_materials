#################################################################
#
# PH791 R Code: Assignment 1
# author: Ari Chernofsky
# date created: August 25, 2020
# version: 001
#
#################################################################



# load libraries ----------------------------------------------------------

library(readxl)

acup <- read_excel("data/acup.xlsx")

names(acup) <- tolower(names(acup))

######################## Question 1 ######################################

# step 1: Expansion of the dataset ----------------------------------------

onesample <- acup[rep(1:nrow(acup), each = 3),]
onesample$interv <- rep(c(-1, 0, 1), nrow(acup))
onesample$response <- ifelse(onesample$interv == -1, onesample$response, NA)
onesample$group <- ifelse(onesample$interv == 0, 0,
                        ifelse(onesample$interv == 1, 1, onesample$group))


# step 2: outcome model and prediction ------------------------------------

out_model <- glm(response ~ group + sex + age + I(age^2) + migraine + chronicity,
                 data = onesample, family = binomial())

onesample$p_response <- predict(out_model, newdata = onesample,
                                type = "response")

results <- tapply(onesample$p_response, onesample$interv, mean)

ace <- c(EYA_0 = results[["0"]], 
         EYA_1 = results[["1"]], 
         ACE = results[["1"]] / results[["0"]])

ace

# Bootstrap confidence intervals ------------------------------------------

set.seed(234)

boot_gformula <- function(){
  
  #bootstrap sample
  bindex <- sample(1:nrow(acup), replace = T)
  bsample <- acup[bindex, ] 
  
  onesample <- bsample[rep(1:nrow(bsample), each = 3),]
  onesample$interv <- rep(c(-1, 0, 1), nrow(bsample))
  onesample$response <- ifelse(onesample$interv == -1, onesample$response, NA)
  onesample$group <- ifelse(onesample$interv == 0, 0,
                            ifelse(onesample$interv == 1, 1, onesample$group))
  
  
  # step 2: outcome model and prediction ------------------------------------
  
  out_model <- glm(response ~ group + sex + age + I(age^2) + migraine + chronicity,
                   data = onesample, family = binomial())
  
  onesample$p_response <- predict(out_model, newdata = onesample,
                                  type = "response")
  
  results <- tapply(onesample$p_response, onesample$interv, mean)
  
  ace <- c(EYA_0 = results[["0"]], 
           EYA_1 = results[["1"]], 
           ACE = results[["1"]] / results[["0"]])
  
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

############################# Question 2 ##########################################

# estimation of denominator probs for IP weights
acup$pd_group <- glm(group ~ sex + age + I(age^2) + migraine + chronicity , 
              family = binomial(), data = acup)$fitted.values

# estimation of numerator probs for IP weights
acup$pn_group <- glm(group ~ 1, 
              family = binomial(), data = acup)$fitted.values

# estimation of stabilized weights
acup$swt <- ifelse(acup$group == 1, acup$pn_group/acup$pd_group,
                (1 - acup$pn_group)/(1 - acup$pd_group))


results_sipw <- lapply(split(acup, acup$group), 
                  function(z) weighted.mean(z$response, z$swt))

ace_sipw <- c(EYA_0 = results_sipw$`0`, 
             EYA_1 = results_sipw$`1`, 
             ACE = results_sipw$`1` / results_sipw$`0`)

ace_sipw


# IPW non-parametric bootstrap confidence intervals -----------------------

boot_ipw <- function(){
  
  #bootstrap sample
  bindex <- sample(1:nrow(acup), replace = T)
  bsample <- acup[bindex, ] 
  
  # estimation of denominator probs for IP weights
  bsample$pd_group <- glm(group ~ sex + age + I(age^2) + migraine + chronicity , 
                       family = binomial(), data = bsample)$fitted.values
  
  # estimation of numerator probs for IP weights
  bsample$pn_group <- glm(group ~ 1, 
                       family = binomial(), data = bsample)$fitted.values
  
  # estimation of stabilized weights
  bsample$swt <- ifelse(bsample$group == 1, bsample$pn_group/bsample$pd_group,
                     (1 - bsample$pn_group)/(1 - bsample$pd_group))
  
  
  results_sipw <- lapply(split(bsample, bsample$group), 
                         function(z) weighted.mean(z$response, z$swt))
  
  ace_sipw <- c(EYA_0 = results_sipw$`0`, 
                EYA_1 = results_sipw$`1`, 
                ACE = results_sipw$`1` / results_sipw$`0`)
  
  ace_sipw
}

#Number of bootstrap samples
B <- 1000

#initialize dataframe to store results from bootstrap
boot <- data.frame(sample = 1:B, EYA_0 = rep(NA, B), 
                   EYA_1 = rep(NA, B), 
                   ACE = rep(NA, B))

#run repeated bootstrap samples and estimates
for(i in 1:B){
  boot[i,2:4] <- boot_ipw()
}

sipw_ey0 <- c(mean = mean(boot$EYA_0), 
              q = quantile(boot$EYA_0, probs = 0.025),
              q = quantile(boot$EYA_0, probs = 0.975))
sipw_ey1 <- c(mean = mean(boot$EYA_1), 
              q = quantile(boot$EYA_1, probs = 0.025),
              q = quantile(boot$EYA_1, probs = 0.975))

sipw_ace <- c(mean = mean(boot$ACE), 
              q = quantile(boot$ACE, probs = 0.025),
              q = quantile(boot$ACE, probs = 0.975))

rbind(sipw_ey0, sipw_ey1, sipw_ace)

