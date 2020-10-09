#################################################################
#
# PH791 R Code: Standardization
#     Assume these data were collected in randomized trial
#     with conditional randomization
#     Randomization is conditional on baseline prognosis 
#     bad_progn 
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


# Standardization ---------------------------------------------------------

means <- data.frame(asp = c(0, 0, 1, 1),
                    bad_progn = c(0, 1, 0, 1),
                    mean = c(mean(bp$change[which(bp$asp == 0 & bp$bad_progn == 0)]),
                             mean(bp$change[which(bp$asp == 0 & bp$bad_progn == 1)]),
                             mean(bp$change[which(bp$asp == 1 & bp$bad_progn == 0)]),
                             mean(bp$change[which(bp$asp == 1 & bp$bad_progn == 1)])),
                    count = c(length(which(bp$asp == 0 & bp$bad_progn == 0)),
                              length(which(bp$asp == 0 & bp$bad_progn == 1)),
                              length(which(bp$asp == 1 & bp$bad_progn == 0)),
                              length(which(bp$asp == 1 & bp$bad_progn == 1))))

means$pr_l <- ifelse(means$bad_progn == 1, length(which(bp$bad_progn == 1))/nrow(bp),
                     length(which(bp$bad_progn == 0))/nrow(bp))
means$total <- means$mean * means$pr_l

EYA_1 <- sum(means$total[which(means$asp == 1)])
EYA_0 <- sum(means$total[which(means$asp == 0)])
ACE <- EYA_1 - EYA_0


# Bootstrap Confidence Intervals ------------------------------------------

set.seed(818)

boot_standard <- function(data){
  bindex <- sample(1:nrow(data), replace = T)
  bsample <- data[bindex, ]
  
  means <- data.frame(asp = c(0, 0, 1, 1),
                      bad_progn = c(0, 1, 0, 1),
                      mean = c(mean(bsample$change[which(bsample$asp == 0 & bsample$bad_progn == 0)]),
                               mean(bsample$change[which(bsample$asp == 0 & bsample$bad_progn == 1)]),
                               mean(bsample$change[which(bsample$asp == 1 & bsample$bad_progn == 0)]),
                               mean(bsample$change[which(bsample$asp == 1 & bsample$bad_progn == 1)])),
                      count = c(length(which(bsample$asp == 0 & bsample$bad_progn == 0)),
                                length(which(bsample$asp == 0 & bsample$bad_progn == 1)),
                                length(which(bsample$asp == 1 & bsample$bad_progn == 0)),
                                length(which(bsample$asp == 1 & bsample$bad_progn == 1))))
  
  means$pr_l <- ifelse(means$bad_progn == 1, length(which(bsample$bad_progn == 1))/nrow(bsample),
                       length(which(bsample$bad_progn == 0))/nrow(bsample))
  means$total <- means$mean * means$pr_l
  
  eya_1 <- sum(means$total[which(means$asp == 1)])
  eya_0 <- sum(means$total[which(means$asp == 0)])
  #ACE <- EYA_1 - EYA_0
  list(E = c(EYA_0 = eya_0, EYA_1 = eya_1), data = bsample)
}

#Number of bootstrap samples
B <- 1000

#initialize dataframe to store results from bootstrap
boot <- data.frame(sample = 1:B, EYA_0 = rep(NA, B), 
                   EYA_1 = rep(NA, B), 
                   ACE = rep(NA, B))

datalist <- vector(mode = "list", length = B)

#run repeated bootstrap samples and estimates
for(i in 1:B){
  bs <- boot_standard(bp)
  boot[i, 2:3] <- bs$E
  datalist[[i]] <- bs$data
}

boot$ACE <- boot$EYA_1 - boot$EYA_0

c(EYA_0 = mean(boot$EYA_0, na.rm = T), 
  quantile(boot$EYA_0, probs = 0.025),
  quantile(boot$EYA_0, probs = 0.975))

c(EYA_1 = mean(boot$EYA_1, na.rm = T), 
  quantile(boot$EYA_1, probs = 0.025, na.rm = T),
  quantile(boot$EYA_1, probs = 0.975, na.rm = T))

c(ACE = mean(boot$ACE, na.rm = T), 
  quantile(boot$ACE, probs = 0.025, na.rm = T),
  quantile(boot$ACE, probs = 0.975, na.rm = T))
