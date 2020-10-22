#################################################################
#
# PH791 R Code: Fish IPW
# author: Ariel Chernofsky
# created: October 19, 2020
#
#
#################################################################

# load libraries ----------------------------------------------------------

library(haven)
library(broom)
library(geepack)

# read in data ------------------------------------------------------------

fish <- read_sas("data/fish.sas7bdat")

fish <- fish[order(fish$id, fish$time), ]

table(fish$time)
prop.table(table(fish$time))

############################### Baseline Analysis ####################################################

# visit count variable ----------------------------------------------------

for(i in seq_along(unique(fish$id))){
  fish[fish$id == i, "vis_count"] <- 1:sum(fish$id == i)
}


fish$fish_bl <- ifelse(fish$fishprebl <= 5, 0,
                       ifelse(fish$fishprebl >= 6, 1, NA))


# maximum number of visits variable ---------------------------------------
maxvisits <- data.frame(id = seq_along(unique(fish$id)), maxvisit = NA)
maxvisits$maxvisit <- tapply(fish$vis_count, fish$id, max)

fish <- merge(fish, maxvisits)

# maximum number of time variable ---------------------------------------
fish$maxtime <- fish$maxvisit - 1


# death overall variable --------------------------------------------------
deaths <- data.frame(id = seq_along(unique(fish$id)), deathoverall = NA)
deaths$deathoverall <- tapply(fish$death, fish$id, function(x) max(x, na.rm = T) )

fish <- merge(fish, deaths)

# fit the weights ---------------------------------------------------------

#create a baseline dataset
fish0 <- fish[fish$time == 0,]

#fit denominator model for weights
dfit_0 <- glm(fish_bl ~ race_1 + race_2 + mtprebl + I(mtprebl^2) + 
                cigprebl + I(cigprebl^2) + hbpprebl, 
              data = fish0, 
              family = binomial())
fish0$pdenom <- dfit_0$fitted.values

#fit numerator model for weights
nfit_0 <- glm(fish_bl ~ 1, family = binomial(), data = fish0)
fish0$pnum <- nfit_0$fitted.values

#calculate weights
fish0$numcont <- fish0$fish_bl*fish0$pnum + (1-fish0$fish_bl)*(1-fish0$pnum)
fish0$dencont <- fish0$fish_bl*fish0$pdenom + (1-fish0$fish_bl)*(1-fish0$pdenom) 

fish0$unstabw <- 1/fish0$dencont
fish0$stabw <- fish0$numcont/fish0$dencont

#summary statistics for weights
c(mean = mean(fish0$unstabw), 
  min = min(fish0$unstabw), 
  max = max(fish0$unstabw), 
  q95 = quantile(fish0$unstabw, probs = 0.95), 
  q99 = quantile(fish0$unstabw, probs = 0.99))
c(mean = mean(fish0$stabw), 
  min = min(fish0$stabw), 
  max = max(fish0$stabw), 
  q95 = quantile(fish0$stabw, probs = 0.95), 
  q99 = quantile(fish0$stabw, probs = 0.99))

#baseline unstabalized weight model
unstabw_base_model <- glm(deathoverall ~ fish_bl, data = fish0, family = quasibinomial(), weights = unstabw)
#coefficient estimates
tidy(unstabw_base_model)
#odds ratios
tidy(unstabw_base_model, exp = T)

#baseline stabilized weight model
stabw_base_model <- glm(deathoverall ~ fish_bl, data = fish0, family = quasibinomial(), weights = stabw)
#coefficient estimates
tidy(stabw_base_model)
#odds ratios
tidy(stabw_base_model, exp = T)

############################### 2 vs. 6 over time ####################################################

# 2 vs. 6 over time  ------------------------------------------------------

#create binary versions of fish and fish_l1 variable
fish$fish_bin <- ifelse(fish$fish <= 5, 0,
                        ifelse(fish$fish >= 6, 1, NA))

fish$fish_bin_l1 <- ifelse(fish$fish_l1 <= 5, 0,
                        ifelse(fish$fish_l1 >= 6, 1, NA))

#time squared variable
fish$time2 <- (fish$time)^2

fish_t <- fish

#create a censoring variable for fish
fish_t$cens_fish <- NA

for(i in seq_along(unique(fish_t$id))){
  temp <- fish_t[fish_t$id == i, ]
  for(j in 1:nrow(temp)){
    if(j == 1) temp$cens_fish[j] <- 0
    else if (j > 1 && (temp$cens_fish[j-1] == 1) || is.na(temp$cens_fish[j-1])) temp$cens_fish[j] <- NA
    else if(temp$fish_bin[j] == temp$fish_bin_l1[j]) temp$cens_fish[j] <- 0
    else if(temp$fish_bin[j] != temp$fish_bin_l1[j]) temp$cens_fish[j] <- 1
  }
  fish_t[fish_t$id == i, ]$cens_fish <- temp$cens_fish
}

#subset data to time after baseline
fish_t_tgt0 <- fish_t[fish_t$time > 0, ]

#denominator model for weights
dfit_t <- glm(fish_bin ~ time + time2 + fishprebl + race_1 + race_2 + mt + 
                I(mt^2) + mtprebl + I(mtprebl^2) + cig + I(cig^2) + 
                cigprebl + I(cigprebl^2) + hbp + hbpprebl, 
              family = binomial(),
              data = fish_t_tgt0)

#numerator model for weights
nfit_t <- glm(fish_bin ~ time + time2 + fishprebl + race_1 + race_2  + 
                mtprebl + I(mtprebl^2) +  
                cigprebl + I(cigprebl^2) + hbpprebl, 
              family = binomial(),
              data = fish_t_tgt0)

fish_t_tgt0$pdenom_t <- dfit_t$fitted.values
fish_t_tgt0$pnum_t <- nfit_t$fitted.values

fish_t <- merge(fish_t, fish_t_tgt0, all.x = T)

fish_t$k1_0 <- 1
fish_t$k1_w <- 1

fish_t$pnum_1 <- ifelse(is.na(fish_t$pnum_t), 1, NA)
fish_t$pdenom_1 <- ifelse(is.na(fish_t$pdenom_t), 1, NA)

fish_t$numcont <- fish_t$fish_bin*fish_t$pnum_t +(1-fish_t$fish_bin)*(1-fish_t$pnum_t)
fish_t$dencont <- fish_t$fish_bin*fish_t$pdenom_t +(1-fish_t$fish_bin)*(1-fish_t$pdenom_t)

for(i in seq_along(unique(fish_t$id))){
  temp <- fish_t[fish_t$id == i, ]
  for(j in 1:nrow(temp)){
    if(j == 1) {
      temp$k1_0[j] <- 1
      temp$k1_w[j] <- 1
      }
    else if (j > 1) {
      temp$k1_0[j] <- temp$k1_0[j - 1] * temp$numcont[j]
      temp$k1_w[j] <- temp$k1_w[j - 1] * temp$dencont[j]
    } 
  }
  fish_t[fish_t$id == i, ]$k1_0 <- temp$k1_0
  fish_t[fish_t$id == i, ]$k1_w <- temp$k1_w
}

#calculate weights
fish_t$unstabw <- 1/fish_t$k1_w
fish_t$stabw <- fish_t$k1_0/fish_t$k1_w

stabw_p99 <- quantile(fish_t$stabw, probs = 0.99)

#truncate weights above 99 percentile
fish_t$stabw_t <- ifelse(fish_t$stabw > stabw_p99, stabw_p99, fish_t$stabw)
summary(fish_t$stabw_t)
summary(fish_t$stabw)


# Marginal structural models ----------------------------------------------

### Unstabalized weights ###
unstabw_out_model <- geeglm(death ~ time + time2 + fish_bl, 
                       data = fish_t,
                       family = binomial(),
                       weights = unstabw,
                       corstr = "independence",
                       id = id,
                       subset = which(cens_fish == 0))

tidy(unstabw_out_model)
tidy(unstabw_out_model, conf.int = T, exp = T)[4,]

### Stabalized weights ###
stabw_out_model <- geeglm(death ~ time + time2 + fish_bl +
                            race_1 + race_2 + 
                            mtprebl + I(mtprebl^2) +
                            cigprebl + I(cigprebl^2) + 
                            hbpprebl, 
                            data = fish_t,
                            family = binomial(),
                            weights = stabw,
                            corstr = "independence",
                            id = id,
                            subset = which(cens_fish == 0))
tidy(stabw_out_model)
tidy(stabw_out_model, conf.int = T, exp = T)[4,]

### Stabalized truncated weights ###
stabwt_out_model <- geeglm(death ~ time + time2 + fish_bl +
                            race_1 + race_2 + 
                            mtprebl + I(mtprebl^2) +
                            cigprebl + I(cigprebl^2) + 
                            hbpprebl, 
                          data = fish_t,
                          family = binomial(),
                          weights = stabw_t,
                          corstr = "independence",
                          id = id,
                          subset = which(cens_fish == 0))
tidy(stabwt_out_model)
tidy(stabwt_out_model, conf.int = T, exp = T)[4,] 

### Unweighted ###
unw_out_model <- geeglm(death ~ time + time2 + fish_bl +
                             race_1 + race_2 + 
                             mtprebl + I(mtprebl^2) +
                             cigprebl + I(cigprebl^2) + 
                             hbpprebl, 
                           data = fish_t,
                           family = binomial(),
                           corstr = "independence",
                           id = id,
                           subset = which(cens_fish == 0))
tidy(unw_out_model)
tidy(unw_out_model, conf.int = T, exp = T)[4,]
