##

require(tidyverse)
library(dplyr)

library(survival)
library(survminer)

require(rstan)
require(brms)
require(splines2)
require(rstanarm)

rstan_options(auto_write = TRUE)


# Import the ovarian cancer dataset and have a look at it

data(cancer, package="survival")

# ovarian <- cancer %>%
#               filter()
            
glimpse(ovarian)

# Dichotomize age and change data labels
ovarian$rx <- factor(ovarian$rx, 
                     levels = c("1", "2"), 
                     labels = c("A", "B"))
ovarian$resid.ds <- factor(ovarian$resid.ds, 
                           levels = c("1", "2"), 
                           labels = c("no", "yes"))
ovarian$ecog.ps <- factor(ovarian$ecog.ps, 
                          levels = c("1", "2"), 
                          labels = c("good", "bad"))

# Data seems to be bimodal
hist(ovarian$age) 


ovarian <- ovarian %>% mutate(age_group = ifelse(age >=50, "old", "young"))
ovarian$age_group <- factor(ovarian$age_group)

surv_object <- Surv(time = ovarian$futime, event = ovarian$fustat)
surv_object 

fit1 <- survfit(surv_object ~ rx, data = ovarian)
summary(fit1)



ggsurvplot(fit1, data = ovarian, pval = TRUE)



# Fit a Cox proportional hazards model

fit.coxph1 <- coxph(surv_object ~ rx + resid.ds + age_group + ecog.ps, 
                   data = ovarian)

summary(fit.coxph1)

fit.coxph2 <- coxph(surv_object ~ rx , 
                    data = ovarian)

summary(fit.coxph2)



ggforest(fit.coxph1, data = ovarian)

### Bayesian cox model

fit1.cox <- brm(data = ovarian,
                  family = brmsfamily(family = "cox"),
                  formula = futime | cens(1 - fustat) ~ 1 + rx,
                  iter = 2000, warmup = 1000, chains = 4, cores = 4,
                  seed = 14,
                  file = "fits/fit1.cox",
                  save_model = "stan_cox1.stan")

print(fit1.cox)


fit14.1 <- brm(data = ovarian,
               family = brmsfamily(family = "exponential"),
               formula = futime | cens(1 - fustat) ~ 1 + rx,
               iter = 2000, warmup = 1000, chains = 4, cores = 4,
               seed = 14,
               file = "fits/fit14.01",
               save_model = "stan_expon1.stan")



print(fit14.1)

fit14.1$model


