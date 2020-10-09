#################################################################
#
# PH791 R Code: Inverse Probability Weighting
# author: Ariel Chernofsky
# created: August 19, 2020
#
#
#################################################################

# load libraries ----------------------------------------------------------

library(readxl)

bp <- read_excel("data/dbp_1.xlsx")


# create a blood pressure change variable ---------------------------------

bp$change <- bp$DBP4 - bp$DBP0

bp$asp <- ifelse(bp$TRT == "A", 1, 0)


# 2 by 2 table ------------------------------------------------------------

addmargins(table(bp$asp, bp$bad_progn),c(1,2))

prop.table(table(bp$asp, bp$bad_progn), margin = 2)


# stratified means by bad prognosis ---------------------------------------

mean_clm <- function(X){
  EX <- mean(X)
  upper <- EX + qt(p = .975, df = length(X)-1)*(sd(X)/sqrt(length(X)))
  lower <- EX - qt(p = .975, df = length(X)-1)*(sd(X)/sqrt(length(X)))
  c(mean = EX, lower = lower, upper = upper)
}

mean_clm(bp$change[which(bp$asp == 0 & bp$bad_progn == 0)])
mean_clm(bp$change[which(bp$asp == 0 & bp$bad_progn == 1)])
mean_clm(bp$change[which(bp$asp == 1 & bp$bad_progn == 0)])
mean_clm(bp$change[which(bp$asp == 1 & bp$bad_progn == 1)])


# Inverse Probabilty Weight -----------------------------------------------

#estimate probability of treatments
p <- glm(asp ~ bad_progn, data = bp, family = binomial())$fitted.values
#for A = 0 compute 1-prob of treatment
prob <- ifelse(bp$asp == 0,1 - p, p)
#IP weights
wt <- 1/prob

EYA_1 <- weighted.mean(bp$change[which(bp$asp == 1)], 
                       w = wt[which(bp$asp == 1)])

EYA_0 <- weighted.mean(bp$change[which(bp$asp == 0)], 
                       w = wt[which(bp$asp == 0)])

ACE <- EYA_1 - EYA_0


# Bootstrap confidence intervals ------------------------------------------

set.seed(79112)

boot_ipw <- function(data){
  bindex <- sample(1:nrow(data), replace = T)
  bsample <- data[bindex, ] 
  p <- glm(asp ~ bad_progn, family = binomial(),
           data = bsample)$fitted.values
  prob <- ifelse(bsample$asp == 0,1 - p, p)
  wt <- 1/prob
  
  eya_1 <- weighted.mean(bsample$change[which(bsample$asp == 1)], 
                         w = wt[which(bsample$asp == 1)])
  
  eya_0 <- weighted.mean(bsample$change[which(bsample$asp == 0)], 
                         w = wt[which(bsample$asp == 0)])
  c(EYA_0 = eya_0, EYA_1 = eya_1)#, ACE = eya_1 - eya_0)
}

#Number of bootstrap samples
B <- 1000

#initialize dataframe to store results from bootstrap
boot <- data.frame(sample = 1:B, EYA_0 = rep(NA, B), 
                   EYA_1 = rep(NA, B), 
                   ACE = rep(NA, B))

#run repeated bootstrap samples and estimates
for(i in 1:B){
  boot[i, 2:3] <- boot_ipw(bp)
}

boot$ACE <- boot$EYA_1 - boot$EYA_0

c(EYA_0 = mean(boot$EYA_0), 
  quantile(boot$EYA_0, probs = 0.025),
  quantile(boot$EYA_0, probs = 0.975))

c(EYA_1 = mean(boot$EYA_1), 
  quantile(boot$EYA_1, probs = 0.025),
  quantile(boot$EYA_1, probs = 0.975))

c(ACE = mean(boot$ACE), 
  quantile(boot$ACE, probs = 0.025),
  quantile(boot$ACE, probs = 0.975))
