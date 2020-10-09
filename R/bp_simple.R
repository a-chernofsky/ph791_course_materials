#################################################################
#
# PH791 R Code: Simple Randomization
# Study participants are randomized 1:1 to treatment (say aspirin) or placebo 
# Primary outcome: change in blood pressure between end of the study 
# and baseline (DBP4-DBP0);
# Secondary outcome: coronary heart disease (Y) 
# author: Ariel Chernofsky
# created: August 19, 2020
#
#
#
#################################################################


# load libraries ----------------------------------------------------------

library(readxl)

bp <- read_excel("data/dbp_1.xlsx")


# create a blood pressure change variable ---------------------------------

bp$change <- bp$DBP4 - bp$DBP0

bp$asp <- ifelse(bp$TRT == "A", 1, 0)


# means and 95% CI --------------------------------------------------------

Y1 <- bp$change[which(bp$asp == 1)]
EY1 <- mean(Y1)
upper1 <- EY1 + qt(p = .975, df = length(Y1)-1)*(sd(Y1)/sqrt(length(Y1)))
lower1 <- EY1 - qt(p = .975, df = length(Y1)-1)*(sd(Y1)/sqrt(length(Y1)))

Y0 <- bp$change[which(bp$asp == 0)]
EY0 <- mean(Y0)
upper0 <- EY0 + qt(p = .975, df = length(Y0)-1)*(sd(Y0)/sqrt(length(Y0)))
lower0 <- EY0 - qt(p = .975, df = length(Y0)-1)*(sd(Y0)/sqrt(length(Y0)))

c(mean1 = EY1, lower1 = lower1, upper1 = upper1)
c(mean0 = EY0, lower0 = lower0, upper0 = upper0)


# similar output using T-Test ---------------------------------------------

t.test(Y1)$conf.int
t.test(Y0)$conf.int

t.test(change ~ asp, data = bp)


# similar output using lm -------------------------------------------------

summary(lm(change ~ asp, data = bp))

t.test(Y ~ asp, data = bp)
