##############################################################################################################################
#
# Piecewise Exponential Model in Stan
# Shamil Sadikhov
#
# Source: https://rpubs.com/kaz_yos/surv_stan_piecewise1
#
#
## Plotting predicted results: https://rstudio-pubs-static.s3.amazonaws.com/435225_07b4ab5afa824342a4680c9fb2de6098.html
#
##############################################################################################################################

require(dplyr)
require(rstan)
library(tidyverse)
library(survival)
library(rstanarm)
library(broom)
library(directlabels)
library(bayesplot)
require(tidybayes)
require(survminer)



data(cancer, package="survival")
glimpse(leukemia)

leukemia <- tibble(leukemia) %>%
                mutate(id = seq_len(n())) %>%
                select(id, everything())
leukemia


# plot KM

require("survival")
fit <- survfit(Surv(time, status) ~ x, data = leukemia)

ggsurvplot(fit, data = leukemia)


## One piece
cut_one <- max(leukemia$time) + 1
cut_one


## Two pieces
cut_two <- c(median(leukemia$time), cut_one) %>% round()
cut_two


## Three pieces
cut_three <- c(quantile(leukemia$time, probs = c(1/3, 2/3)), cut_one) %>% round()
cut_three


## At all event times
cut_events <- c(leukemia %>%
                  filter(status == 1) %>%
                  select(time) %>%
                  arrange(time) %>%
                  magrittr::extract2("time") %>%
                  unique,
                cut_one) %>% round()
cut_events

## At all event or censoring
cut_times <- c(sort(unique(leukemia$time)),
               cut_one) %>% round()
cut_times

##############################################################################################################################
## No cut for all same constant hazard
##############################################################################################################################
leukemia_one <- survival::survSplit(formula = Surv(time, status) ~ ., data = leukemia, cut = cut_one) %>%
  mutate(interval = factor(tstart),
         interval_length = time - tstart) %>%
  as_data_frame
leukemia_one

## Split into two observations
leukemia_two <- survival::survSplit(formula = Surv(time, status) ~ ., data = leukemia, cut = cut_two) %>%
  mutate(interval = factor(tstart),
         interval_length = time - tstart) %>%
  as_data_frame
leukemia_two

##############################################################################################################################
## Split into three observations
leukemia_three <- survival::survSplit(formula = Surv(time, status) ~ ., data = leukemia, cut = cut_three) %>%
  mutate(interval = factor(tstart),
         interval_length = time - tstart) %>%
  as_data_frame
leukemia_three

## Split at event times
leukemia_events <- survival::survSplit(formula = Surv(time, status) ~ ., data = leukemia, cut = cut_events) %>%
                   mutate(interval = factor(tstart),
                   interval_length = time - tstart) %>%
                    as_data_frame
leukemia_events

# add timepoint id
# leukemia_events %<>% group_by(id) %>% mutate(tid = seq_along(id))

leukemia_events$tid <-  with(leukemia_events, match(time, sort(unique(time))))

##############################################################################################################################
## Split at event and censoring times
leukemia_times <- survival::survSplit(formula = Surv(time, status) ~ ., data = leukemia, cut = cut_times) %>%
  mutate(interval = factor(tstart),
         interval_length = time - tstart) %>%
  as_data_frame
leukemia_times


##############################################################################################################################

##############################################################################################################################

coxph(formula = Surv(time, status) ~ x,
      data    = leukemia,
      ties    = c("efron","breslow","exact")[1]) %>% summary

glm(formula = status ~ x + offset(log(interval_length)),
    family  = poisson(link = "log"),
    data    = leukemia_one) %>% summary

## Drop intercept to show all interval-specific log rates
glm(formula = status ~ -1 + interval + x + offset(log(interval_length)),
    data = leukemia_two,
    family = poisson(link = "log")) %>% summary

glm(formula = status ~ -1 + interval + x + offset(log(interval_length)),
    data = leukemia_three,
    family = poisson(link = "log")) %>% summary

glm(formula = status ~ -1 + interval + x + offset(log(interval_length)),
    data = leukemia_events,
    family = poisson(link = "log")) %>% summary

glm(formula = status ~ -1 + interval + x + offset(log(interval_length)),
    data = leukemia_times,
    family = poisson(link = "log")) %>% summary



## One Cut
fit_leukemia_one <- rstanarm::stan_glm(formula = status ~ x + offset(log(interval_length)),
                                       data = leukemia_one,
                                       family = poisson(link = "log"),
                                       prior_intercept = normal(location = 0, scale = 4),
                                       prior = normal(location = 0, scale = 4))

fit_leukemia_one <- rstanarm::stan_glm(formula = status ~ x + offset(log(interval_length)),
                                       data = leukemia_one,
                                       family = poisson(link = "log"),
                                       prior_intercept = norm(0, 4),
                                       prior = norm(0, 4))

fit_leukemia_one <- rstanarm::stan_glm(formula = status ~ x + offset(log(interval_length)),
                                       data = leukemia_one,
                                       family = poisson(link = "log")
                                      )


rstanarm::prior_summary(fit_leukemia_one)

summary(fit_leukemia_one)

stan_code <- get_stancode(fit_leukemia_one$stanfit)

# show readable stan code
cat(stan_code )





##############################################################################################################################
## AT events
##############################################################################################################################

# Piecewise exponential https://rpubs.com/kaz_yos/surv_stan_piecewise1


fit_leukemia_events <- rstanarm::stan_glm(formula = status ~ -1 + interval + x + offset(log(interval_length)),
                                          data = leukemia_events,
                                          family = poisson(link = "log"))
rstanarm::prior_summary(fit_leukemia_events)
fit_leukemia_events


stan_leuk_events <- stan('model/pem_leuk_stan_SS.stan',
                         data=list(N=NROW(unique(leukemia_events$id)),
                                   J=NROW(cut_events),
                                   NR = NROW(leukemia_events),
                                   x = as.integer(leukemia_events$x)-1,
                                   status=as.integer(leukemia_events$status),
                                   interval=as.integer(leukemia_events$interval),
                                   interval_fu_time = as.integer(leukemia_events$interval_length)),
                         chains = 4,
                         warmup = 2000,
                         iter = 4000,
                         cores = 4,
                         refresh = 0,
                         control = list(adapt_delta = 0.99))

stan_leuk_events 

# with hierarchical prior on baseline hazards

stan_leuk_events2 <- stan('model/pem_leuk_stan_SS2.stan',
                         data=list(N=NROW(unique(leukemia_events$id)),
                                   J=NROW(cut_events),
                                   NR = NROW(leukemia_events),
                                   x = as.integer(leukemia_events$x)-1,
                                   status = as.integer(leukemia_events$status),
                                   interval = as.integer(leukemia_events$interval),
                                   interval_fu_time = as.integer(leukemia_events$interval_length)),
                         chains = 4,
                         warmup = 4000,
                         iter = 8000,
                         cores = 4,
                         refresh = 0,
                         control = list(adapt_delta = 0.95))

summary(stan_leuk_events2)


##############################################################################################################################
##
##       PEM from Survivalstan python package
##      https://github.com/hammerlab/survivalstan/blob/master/survivalstan/stan/pem_survival_model_unstructured.stan
##
##############################################################################################################################


## Unstructured model pem_survival_model_unstructured.stan
t_obs1 <-  leukemia_events %>% group_by(tid) %>% distinct(time)
t_obs <-  sort(as.numeric(as.character( pull(t_obs1, time))))


t_dur <-  c(t_obs[1], diff(t_obs))   # duration at first timepoint = t_obs[1] ( implicit t0 = 0 )



stan_leuk_events_pem_un <- stan('model/survivalstan/pem_survival_model_unstructured.stan',
                          data=list(N = dim(leukemia_events)[1],          # total number of observations
                                    S = NROW(unique(leukemia_events$id)), # number of ids
                                    T = max(leukemia_events$tid),         # max timepoint (number of timepoint ids) # max(leukemia_events$time)
                                    M = 1,                                # number of covatiates
                                    s = leukemia_events$id,               # sample id for each obs
                                    t = leukemia_events$tid,              # timepoint id for each obs
                                    event = as.integer(leukemia_events$status),
                                    x =   as.matrix(as.numeric(leukemia_events$x) - 1),  # covariates matrix
                                    t_obs =  t_obs,       #observed time since origin for each timepoint id (end of period)
                                    t_dur = t_dur  # duration of each timepoint period (first diff of t_obs)
                                ),
                          chains = 4,
                          warmup = 4000,
                          iter = 8000,
                          cores = 4,
                          refresh = 0,
                          control = list(adapt_delta = 0.95))

summary(stan_leuk_events_pem_un)



##############################################################################################################################
##  pem_survival_model.stan
##############################################################################################################################

stan_leuk_events_pem <- stan('model/survivalstan/pem_survival_model.stan',
                                data=list(N = dim(leukemia_events)[1],          # total number of observations
                                          S = NROW(unique(leukemia_events$id)), # number of ids
                                          T = max(leukemia_events$tid),         # max timepoint (number of timepoint ids) # max(leukemia_events$time)
                                          M = 1,                                # number of covatiates
                                          s = leukemia_events$id,               # sample id for each obs
                                          t = leukemia_events$tid,              # timepoint id for each obs
                                          event = as.integer(leukemia_events$status),
                                          x =   as.matrix(as.numeric(leukemia_events$x) - 1),  # covariates matrix
                                          obs_t =  as.numeric(leukemia_events$time)       # observed end time for each obs
                                ),
                                chains = 1,
                                warmup = 4000,
                                iter = 8000,
                                cores = 4,
                                refresh = 0,
                                control = list(adapt_delta = 0.99))

summary(stan_leuk_events_pem)

##############################################################################################################################
##  pem_survival_model_randomwalk.stan
##############################################################################################################################

t_obs1 <-  leukemia_events %>% group_by(tid) %>% distinct(time)
t_obs <-  sort(as.numeric(as.character( pull(t_obs1, time))))


t_dur <-  c(t_obs[1], diff(t_obs))   # duration at first timepoint = t_obs[1] ( implicit t0 = 0 )


stan_leuk_events_pem_rw <- stan('model/survivalstan/pem_survival_model_randomwalk.stan',
                             data=list(N = dim(leukemia_events)[1],          # total number of observations
                                       S = NROW(unique(leukemia_events$id)), # number of ids
                                       T = max(leukemia_events$tid),         # max timepoint (number of timepoint ids) # max(leukemia_events$time)
                                       M = 1,                                # number of covatiates
                                       s = leukemia_events$id,               # sample id for each obs
                                       t = leukemia_events$tid,              # timepoint id for each obs
                                       event = as.integer(leukemia_events$status),
                                       x =   as.matrix(as.numeric(leukemia_events$x) - 1),  # covariates matrix
                                       t_obs =  t_obs,     #observed time since origin for each timepoint id (end of period)
                                       t_dur = t_dur  # duration of each timepoint period (first diff of t_obs)
                             ),
                             chains = 4,
                             warmup = 4000,
                             iter = 8000,
                             cores = 4,
                             refresh = 0,
                             control = list(adapt_delta = 0.99))

summary(stan_leuk_events_pem_rw)

mcmc_trace(stan_leuk_events_pem_rw, pars=c("beta[1]", "baseline_sigma"))

bayesplot::mcmc_areas(as.matrix(stan_leuk_events_pem_rw), pars = c("beta[1]", "baseline_sigma"), prob = 0.95)

stan_weibull_survival_model_draws <- tidybayes::tidy_draws(stan_leuk_events_pem_rw)

# posterior_predict(stan_leuk_events_pem_rw, draws = 500)

yhat_time <- rstan::extract(stan_leuk_events_pem_rw)[["y_hat_time"]]
yhat_event <- rstan::extract(stan_leuk_events_pem_rw)[["y_hat_event"]] 

plot(apply(yhat_time, 2, mean), apply(yhat_event, 2, mean))

plot(yhat_time, yhat_event)
