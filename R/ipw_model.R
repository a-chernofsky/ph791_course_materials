#################################################################
#
# PH791 R Code for IPW using parametric models
#           using non-saturated model (making model assumptions)
# author: Ariel Chernofsky
# created: September 15, 2020
#
# Non-saturated model 
# Goal: to estimate the average causal effect of 
# in the aspirin data using a simulation method that does not require 
# specifying the distribution of P[L=l] while adjusted for
# several confounders incluing continous ones
#
#################################################################



# load libraries ----------------------------------------------------------

library(readxl)
library(geepack)
library(doBy)
library(weights)

# read in data ------------------------------------------------------------

bp <- read_excel("data/dbp_1.xlsx")

names(bp) <- tolower(names(bp))

# create a blood pressure change variable ---------------------------------

bp$change <- bp$dbp4 - bp$dbp0

bp$asp <- ifelse(bp$trt == "A", 1, 0)

bp$sex <- relevel(factor(bp$sex), ref = "M")


# estimation of ip weights via a logistic regression model ------------------------------------------------------

#treatment model
trt_model <- glm(asp ~ bad_progn + sex + 
                   age + I(age^2) + 
                   dbp0 + I(dbp0^2), 
            family = binomial(), data = bp)

#estimated probability of treatment
p_asp <- trt_model$fitted.values

#calculate weights
bp$wt <- ifelse(bp$asp == 1, 1/p_asp, 1/(1-p_asp))

#check weight values
summary(bp$wt)

#  saturated model to estimate the average treatment effect --------

# Now, use a saturated model to estimate the average 
# treatment effect in the two populations;
# To compute the variance I am using robust standard error;

#outcome model
w_out_model <- geeglm(change ~ asp, data = bp, 
                  weights = wt, 
                  corstr = "independence",
                  id = subject)

#contrast matrix
l <- matrix(c(1, 1, 1, 0), nrow = 2)

# estimate IPW means with confidence intervals with 
# robust standard errors
(ipw_means <- data.frame(esticon(w_out_model, l)))

# Show that the IPW make A and L independent in 
# the pseudopopulation

wtd.chi.sq(bp$asp, bp$bad_progn,weight=bp$wt)
chisq.test(bp$asp, bp$bad_progn, correct = F)

# stabilized weights ------------------------------------------------------

# estimation of denominator probs for IP weights
pd_asp <- trt_model$fitted.values
# estimation of numerator probs for IP weights
pn_asp <- glm(asp ~ 1, 
                 family = binomial(), data = bp)$fitted.values
# estimation of stabilized weights
bp$sw <- ifelse(bp$asp == 1, pn_asp/pd_asp,
                (1 - pn_asp)/(1 - pd_asp))

summary(bp$sw)

# Now, use a saturated model to estimate the 
# average treatment effect in the two populations
# To compute the variance I am using robust standard error
sw_out_model <- geeglm(change ~ asp, data = bp, 
                       weights = sw, 
                       corstr = "independence",
                       id = subject)

# contrast matrix
l <- matrix(c(1, 1, 1, 0), nrow = 2)

# estimate IPW means with stabilized weights with confidence intervals with 
# robust standard errors
(ipsw_means <- data.frame(esticon(sw_out_model, l)))


# check positivity --------------------------------------------------------

#index for females without a bad prognosis
index <- which(bp$sex == 'F', bp$bad_progn == 0)
#check positivity
table(bp$age[index], bp$asp[index])
