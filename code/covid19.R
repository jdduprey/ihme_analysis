# =================================================
# Author: Joe Duprey
# assessment for IHME researcher position
# the following script: 
# > visualizes covid cases, deaths, and hospitalizations
# > fits a logistic growth model to cumulative deaths by minimizing RMSE
# > fits two linear models with death as the independent variable and
# cases/hospitalizations as the dependent variable
# > forecasts deaths out 14 days based on the best-fit logistic growth model, then
# predicts cases and hospitalizations based on the linear models. 

library(tidyverse)
library(ggplot2)
library(dplyr)

# load data and create more readable colums
dat = read_csv("../data/covid_data_cases_deaths_hosp.csv")
colnames(dat)[2] = "state"
colnames(dat)[5] = "confirmed_cases"

# get daily deaths, cases and hospitalizations
# get state level data 
dat <- dat %>%
  filter(location_id != 102)

# get daily indicators from cumulative indicators
dat <- dat %>%
  group_by(location_id) %>%
  mutate(daily_deaths = c(Deaths[1], diff(Deaths))) %>%
  mutate(daily_cases = c(confirmed_cases[1], diff(confirmed_cases))) %>%
  mutate(daily_hosp = c(Hospitalizations[1], diff(Hospitalizations)))

# plot every state on one page 
pdf("../viz/US_cases_hosp_plot.pdf", width=18, height=9)  
ggplot(data=dat, aes(x=daily_cases, y=daily_hosp, color=state)) +
  geom_point() 
dev.off()

pdf("../viz/US_deaths_hosp_plot.pdf", width=12, height=9)  
ggplot(data=dat, aes(x=daily_deaths, y=daily_hosp, color=state)) +
  geom_point() +
  xlim(0, 300)
dev.off()

pdf("../viz/US_deaths_cases_plot.pdf", width=12, height=9)  
ggplot(data=dat, aes(x=daily_deaths, y=daily_cases, color=state)) +
  geom_point() +
  xlim(0, 300)
dev.off()

# Data exploration by state, plot relationships for each state individually 
pdf("../viz/cases_hosp_plots.pdf", onefile = TRUE)
for (loc in unique(dat$location_id)){
  print(loc)
  
  plot_dat <- dat %>%
    filter(location_id %in% loc)
  
  my_plot <- ggplot(data=plot_dat, aes(x=daily_cases, y=daily_hosp)) +
    geom_point() +
    labs(title = plot_dat$state) +
    theme_minimal()
  print(my_plot)
  
}
dev.off()

# each state deaths vs cases
pdf("../viz/deaths_cases_plots.pdf", onefile = TRUE)
for (loc in unique(dat$location_id)){
  print(loc)
  
  plot_dat <- dat %>%
    filter(location_id %in% loc)
  
  my_plot <- ggplot(data=plot_dat, aes(x=daily_deaths, y=daily_cases)) +
    geom_point() +
    labs(title = plot_dat$state) +
    theme_minimal()
  print(my_plot)
  
}
dev.off()

# daily death curve fitting 
# =======================================
# first get cumulative deaths at the national level

nat <- dat %>%
  ungroup() %>%
  select(state, Date, Deaths) %>%
  pivot_wider(id_cols=Date, names_from=state, values_from=Deaths) %>%
  arrange(Date) %>%
  replace(is.na(.), 0) # replace NAs with zero since they only appear before 1st death

nat_us <- nat %>%
  select(-Date) %>% 
  mutate(cumulative_US_deaths=rowSums(.)) # sum all state columns to get US col 

national_df <- cbind(nat$Date, nat_us)
colnames(national_df)[1] = "Date"

# TODO is this an appropriate approach for fitting logistic model
# might want to cut it off even later? 
# start when cases are larger to avoid low case oddities
national_df_cut <- national_df %>%
  filter(Date > 43896)   

# fit curve to cumulative deaths 
# calculate RMSE for logistic growth model +++++++++++++++++++++++++++
Ntobs <- national_df_cut$cumulative_US_deaths
ndata <- length(Ntobs)


logistic_model_rmse<- function(pars, Ntobs) {
  #accepts logistic growth params and observered data, returns RMSE,
  #to be used by optim to minimize RMSE
  r <- pars[1]
  K <- pars[2]
  init_cases <- pars[3]
  #sigma <- pars[4] 
  
  ndata <- length(Ntobs)
  Nt <- rep(NA, ndata)
  Nt[1] <- init_cases
  
  for (i in 1:(ndata-1)) Nt[i+1] <- Nt[i] + Nt[i]* r * ( 1 - Nt[i] / K)
  
  RMSE <- sqrt((sum((Nt - Ntobs)^2)) / ndata)
  #epst <- Ntobs - Nt
  #nll <- -dnorm(x = epst, mean = 0, sd = sigma, log = T)
  
  return(RMSE)
  #return(sum(nll, na.rm = TRUE))
}

logistic_preds<- function(pars, Ntobs) {
  r <- pars[1]
  K <- pars[2]
  init_cases <- pars[3]
  
  ndata <- length(Ntobs)
  Nt <- rep(NA, ndata)
  Nt[1] <- init_cases
  
  for (i in 1:(ndata-1)) Nt[i+1] <- Nt[i] + Nt[i]* r * ( 1 - Nt[i] / K)
  
  return(Nt)
}

logistic_forecast<- function(pars, Ntobs, forecast_days) {
  r <- pars[1]
  K <- pars[2]
  init_cases <- pars[3]

  ndata <- length(Ntobs)
  Nt <- rep(NA, ndata)
  Nt[1] <- init_cases
  
  for (i in 1:(ndata+forecast_days)) Nt[i+1] <- Nt[i] + Nt[i]* r * ( 1 - Nt[i] / K)
  
  #epst <- Ntobs - Nt
  #nll <- -dnorm(x = epst, mean = 0, sd = sigma, log = T)
  
  return(Nt)
}

exp_RMSE <- function(pars, Ntobs) {
  r <- pars[1]
  init_cases <- pars[2]
  
  ndata <- length(Ntobs)
  Nt <- rep(NA, ndata)
  Nt[1] <- init_cases
  
  for (i in 1:(ndata-1)) Nt[i+1] <- Nt[i] + Nt[i] * r 
  
  RMSE <- sqrt((sum((Nt - Ntobs)^2)) / ndata)
  
  return(RMSE)
}

exp_preds <- function(pars, Ntobs) {
  r <- pars[1]
  init_cases <- pars[2]
  
  ndata <- length(Ntobs)
  Nt <- rep(NA, ndata)
  Nt[1] <- init_cases
  
  for (i in 1:(ndata-1)) Nt[i+1] <- Nt[i] + Nt[i] * r
  
  return(Nt)
}


# optimization to find params LOGISTIC +++++++++++++++++
start.pars <- c(0.2, 200000, 1) 
sol <- optim(start.pars, 
            fn = logistic_model_rmse,
            method = "Nelder-Mead",
            Ntobs = Ntobs)

preds <- logistic_preds(pars=sol$par, Ntobs=Ntobs)
print(preds)

# optimization to find params EXPONENTIAL +++++++++++++++++
exp.pars <- c(0.02, 1)
exp.sol <- optim(exp.pars, 
             fn = exp_RMSE,
             method = "Nelder-Mead",
             Ntobs = Ntobs)

exp_preds <- exp_preds(pars=exp.sol$par, Ntobs=Ntobs)
print(exp_preds)

# add predictions for logistic and exp to dataframe 
national_df_cut$my_preds <- preds
national_df_cut$exp_preds <- exp_preds

# sanity check visualization LOGISTIC
pdf("../viz/US_cumulative_model_fit.pdf", onefile = TRUE)
ggplot(data=national_df_cut, aes(x=Date, y=cumulative_US_deaths)) +
  geom_point() +
  geom_line(data=national_df_cut, aes(x=Date, y=preds), color="red") +
  theme_minimal()
dev.off() 

# sanity check visualization EXP
pdf("../viz/US_cumulative_model_fit_exp.pdf", onefile = TRUE)
ggplot(data=national_df_cut, aes(x=Date, y=cumulative_US_deaths)) +
  geom_point() +
  geom_line(data=national_df_cut, aes(x=Date, y=exp_preds), color="blue") +
  theme_minimal()
dev.off() 

# convert back to daily numbers
national_df_cut <- national_df_cut %>% 
  mutate(daily_preds = c(245, diff(preds))) %>%
  mutate(daily_obs = c(cumulative_US_deaths[1], diff(cumulative_US_deaths)))

# visualization of daily model fit
pdf("../viz/US_daily_model_fit.pdf", onefile = TRUE)
ggplot(data=national_df_cut, aes(x=Date, y=daily_obs)) +
  geom_point() +
  geom_line(data=national_df_cut, aes(x=Date, y=daily_preds), color="red") +
  theme_minimal()
dev.off() 

# forecast 14 days into future
my_forecast <- logistic_forecast(pars=sol$par, forecast_days=14)
my_forecast <- as.data.frame(my_forecast)
daily_forecast <- my_forecast %>% 
  mutate(daily_forecast = c(245, diff(my_forecast)))

# get dataframe to plot forecasts and to fit/predict with linear models
forecast_date <- c(national_df_cut$Date, seq(44050, 44064, by=1))
daily_forecast_df <- cbind(daily_forecast, forecast_date)
colnames(daily_forecast_df)[2] <- "daily_deaths"

# Model hospitalizations and deaths as a function of deaths
hosp_lm <- lm(formula = daily_hosp ~ daily_deaths, data=dat)
summary(hosp_lm)

case_lm <- lm(formula = daily_cases ~ daily_deaths, data=dat)
summary(case_lm)

# use these linear models to model past hosp and cases, and to
# predict 14 days into the future 
daily_forecast_df$forecast_hosp <- predict(hosp_lm, newdata=daily_forecast_df)
daily_forecast_df$forecast_cases <- predict(case_lm, newdata=daily_forecast_df)

pdf("../viz/US_daily_forecast.pdf", onefile = TRUE)
ggplot(data=national_df_cut, aes(x=Date, y=daily_obs)) +
  geom_point() +
  geom_line(data=daily_forecast_df, aes(x=forecast_date, y=daily_deaths), color="dark green") +
  geom_line(data=daily_forecast_df, aes(x=forecast_date, y=forecast_hosp), color="purple") +
  geom_line(data=daily_forecast_df, aes(x=forecast_date, y=forecast_cases), color="dark orange") +
  theme_minimal() +
  ylim(0, 2000) + 
  geom_vline(xintercept=44050, linetype="dotted") +
  ylab("daily observed and predicted")
dev.off() 

pdf("../viz/US_daily_forecast_all_indicators.pdf", onefile = TRUE)
ggplot(data=national_df_cut, aes(x=Date, y=daily_obs)) +
  geom_point() +
  geom_line(data=daily_forecast_df, aes(x=forecast_date, y=daily_deaths), color="dark green") +
  geom_line(data=daily_forecast_df, aes(x=forecast_date, y=forecast_hosp), color="purple") +
  geom_line(data=daily_forecast_df, aes(x=forecast_date, y=forecast_cases), color="dark orange") +
  theme_minimal() +
  geom_vline(xintercept=44050, linetype="dotted") +
  ylab("daily observed and predicted")
dev.off() 

#TODO 
#fit logistic model for each state individually  ================
# cal_dat <- dat %>%
#   filter(state %in% c("California")) %>%
#   filter(Date > 43915)
# 
# cal_obs <- cal_dat$Deaths
# 
# cal_start.pars <- c(0.02, 20000, 80) 
# cal_sol <- optim(cal_start.pars, 
#              fn = logistic_model_rmse,
#              method = "Nelder-Mead",
#              Ntobs = cal_obs)
# 
# cal_sol$par
# 
# cal_preds <- logistic_preds(pars=cal_sol$par, Ntobs=cal_obs)
# print(cal_preds)
# 
# cal_dat$my_preds <- cal_preds
# 
# ggplot(data=cal_dat, aes(x=Date, y=Deaths)) +
#   geom_point() +
#   geom_line(data=cal_dat, aes(x=Date, y=my_preds, color="red"))
