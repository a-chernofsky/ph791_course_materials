#################################################################
#
# PH791 R Code: for STANDARDIZATION with time to event outcome and 
#     time fixed covariates
# author: Ariel Chernofsky
# created: October 6, 2020
#
# Saturated model 
# Goals: to estimate the causal effect of fish consumption at 
#     baseline on time to death
#
#################################################################

# load libraries ----------------------------------------------------------

library(haven)
library(tidyverse)

# read in data ------------------------------------------------------------

fish <- read_sas("data/fishdata_1.sas7bdat")

summary(fish$time)
summary(fish$fish_bl)

# POTENTIAL CONFOUNDERS: baseline meet consumption 
# (continous, number of serving per week,
# high blood pressure (binary, yes versus no), 
# race (categorical 1,2, and 3), number of cigarettes per day

##################################################################################
######## STEP 1 POOLED LOGISTIC REGRESSION AND PREDICTION
######## STEP 2 PREDICTION OF DISCRETE HAZARD AT EACH TIME FOR EACH INDIVIDUAL
##################################################################################

inmodel <- glm(death ~ fish_bl + time + fish_bl:time + 
                 race_1 + race_2 + 
                 mt_bl + I(mt_bl^2) +
                 cig_bl + I(cig_bl^2) + 
                 hbp_bl,
               family = binomial(),
               data = fish)

##################################################################################
####### STEP 3 COMPUTE SURVIVAL ESTIMATES 
####### OPTIONAL: COMPUTE THE HAZARD RATIOS 
####### Make copies for each individual under each intervention
####### Each individuals will have a follow-up as long as 8 visits;
##################################################################################

#create copies 
fish_copies <- fish[fish$time == 0,]
fish_copies <- fish_copies[rep(1:10000, each = 18),]
fish_copies$time <- rep(0:8, each =2, times = 10000)
fish_copies$fish_bl <- rep(c(2,6), times = 9*10000)
fish_copies$timesq  <- fish_copies$time^2

#predict probabilities of event = 0 on the copied data
p_0 = 1 - predict(inmodel, newdata = fish_copies, type = "response")
fish_tf_pred <- data.frame(id = fish_copies$id, 
                           time = fish_copies$time,
                           fish_bl = fish_copies$fish_bl,
                           p_noevent = p_0)
#reorder rows 
fish_tf_pred <- fish_tf_pred[order(fish_tf_pred$id, 
                                   fish_tf_pred$fish_bl,
                                   fish_tf_pred$time), ]

#create new dataframe to store survival estimates
surv <- fish_tf_pred
surv$p_event = 1 - surv$p_noevent

#estimate survival function for each id for each fish
s_t <- tapply(surv$p_noevent, 
              list(surv$fish_bl, surv$id), 
              cumprod)
surv$surv <- unlist(s_t)


#calculate mean per time point per type of fish
mean_surv <- tapply(surv$surv, list(surv$fish_bl, surv$time), mean)

mysurv <- data.frame(expand.grid(time = 0:8, fish_bl = c(2,6)), 
                     surv = c(mean_surv[1,], mean_surv[2,]))

effects <- reshape(mysurv, v.names = "surv", 
                   timevar = "fish_bl",
                   idvar = "time",
                   direction = "wide")
effects$ci_6 <- 1 - effects$surv.6
effects$ci_2 <- 1 - effects$surv.2
effects$risk_diff <- effects$ci_6 - effects$ci_2
effects$risk_ratio <- effects$ci_6/effects$ci_2

#plot survival functions 
plot(mysurv$time[mysurv$fish_bl == 2], mysurv$surv[mysurv$fish_bl == 2], 
     type = "l",
     col = "blue", 
     ylim = c(0.8,1),
     xlab = "time",
     ylab = "S(t)")
lines(mysurv$time[mysurv$fish_bl == 6], 
      mysurv$surv[mysurv$fish_bl == 6], 
      type = "l",
      col = "red", ylim = c(0.8,1))
legend("bottomleft",
       c("A = 2 serving/week","A = 6 serving/week"),
       fill=c("blue","red"))

#############################################################################
#### OPTIONAL COMPUTE THE HAZARDS AT EACH TIME AND UNDER EACH STRATEGY 
#### THE HR is average of the ratio of the ratio of these hazards over time
#############################################################################

mean_p_event <- tapply(surv$p_event, list(surv$fish_bl, surv$time), mean)
myhazard <- data.frame(expand.grid(time = 0:8, fish_bl = c(2,6)), 
                       hazard = c(mean_p_event[1,], mean_p_event[2,]))

hazard_ratio <- reshape(myhazard, v.names = "hazard", 
                        timevar = "fish_bl",
                        idvar = "time",
                        direction = "wide")
hazard_ratio$hazard_ratio <- hazard_ratio$hazard.6 / hazard_ratio$hazard.2

mean(hazard_ratio$hazard_ratio)

##################################################################################
#### Step 4: Estimate Variability using bootstrap
##################################################################################

boot <- function(){
  bid <- sample(unique(fish$id), replace = TRUE)
  rows <-  sapply(bid, function(x) rownames(fish)[fish$id == x])
  bsample <- fish[as.numeric(do.call(c, rows)), ]
  
  inmodel <- glm(death ~ fish_bl + time + fish_bl:time + 
                   race_1 + race_2 + 
                   mt_bl + I(mt_bl^2) +
                   cig_bl + I(cig_bl^2) + 
                   hbp_bl,
                 family = binomial(),
                 data = bsample)
  
  
  
  bsample_copies <- bsample[bsample$time == 0,]
  bsample_copies <- bsample_copies[rep(1:10000, each = 18),]
  bsample_copies$time <- rep(0:8, each =2, times = 10000)
  bsample_copies$fish_bl <- rep(c(2,6), times = 9*10000)
  bsample_copies$timesq  <- bsample_copies$time^2
  bsample_copies$newid <- rep(1:10000, each = 18)
  
  p_0 = 1 - predict(inmodel, newdata = bsample_copies, type = "response")
  bsample_pred <- data.frame(id = bsample_copies$id, 
                             newid = bsample_copies$newid,
                             time = bsample_copies$time,
                             fish_bl = bsample_copies$fish_bl,
                             p_noevent = p_0)
  bsample_pred <- bsample_pred[order(bsample_pred$newid, 
                                     bsample_pred$fish_bl,
                                     bsample_pred$time), ]
  
  
  surv <- bsample_pred
  surv$p_event = 1 - surv$p_noevent
  
  s_t <- tapply(surv$p_noevent,
                list(surv$fish_bl, surv$newid),
                cumprod)
  surv$surv <- unlist(s_t)
  
  mean_surv <- tapply(surv$surv, list(surv$fish_bl, surv$time), mean)
  
  mysurv <- data.frame(expand.grid(time = 0:8, fish_bl = c(2,6)), 
                       surv = c(mean_surv[1,], mean_surv[2,]))
  
  effects <- reshape(mysurv, v.names = "surv", 
                     timevar = "fish_bl",
                     idvar = "time",
                     direction = "wide")
  effects$ci_6 <- 1 - effects$surv.6
  effects$ci_2 <- 1 - effects$surv.2
  effects$risk_diff <- effects$ci_6 - effects$ci_2
  effects$risk_ratio <- effects$ci_6/effects$ci_2
  
  effects_8 <- effects[effects$time == 8,]
  
  list(ci_2 = effects_8$ci_2, 
       ci_6 = effects_8$ci_6,
       risk_diff = effects_8$risk_diff,
       risk_ratio = effects_8$risk_ratio)
}


#Number of bootstrap samples
B <- 10

#initialize dataframe to store results from bootstrap
boots <- data.frame(sample = 1:B, ci_2 = rep(NA, B), 
                   ci_6 = rep(NA, B), 
                   risk_diff = rep(NA, B),
                   risk_ratio = rep(NA, B))

#run repeated bootstrap samples and estimates
s <- Sys.time()
for(i in 1:B){
  bs <- boot()
  boots$ci_2[i] <- bs$ci_2
  boots$ci_6[i] <- bs$ci_6
  boots$risk_diff[i] <- bs$risk_diff
  boots$risk_ratio[i] <- bs$risk_ratio
}
e <- Sys.time()
e-s

c(ci_2 = mean(boots$ci_2), 
  quantile(boots$ci_2, probs = 0.025),
  quantile(boots$ci_2, probs = 0.975))

c(ci_6 = mean(boots$ci_6), 
  quantile(boots$ci_6, probs = 0.025),
  quantile(boots$ci_6, probs = 0.975))

c(risk_diff = mean(boots$risk_diff), 
  quantile(boots$risk_diff, probs = 0.025),
  quantile(boots$risk_diff, probs = 0.975))

c(risk_ratio = mean(boots$risk_ratio), 
  quantile(boots$risk_ratio, probs = 0.025),
  quantile(boots$risk_ratio, probs = 0.975))
