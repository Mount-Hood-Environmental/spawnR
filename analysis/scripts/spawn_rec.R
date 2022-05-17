# Author: Mike Ackerman
#  - using code from Validation_SpawnRec.Rmd in KevinSee/QRFcapacity
# Purpose: 
# Created: May 17, 2022
# Last Modified: 
# Notes: 

# clear environment
rm(list = ls())

#-----------------------
# load packages
library(FSA)
library(here)
library(tidyverse)
library(tidyr)
library(readr)
library(magrittr)

#-----------------------
# load spawner-recruit data
load(here('analysis/data/spawn_recr_data.rda'))
#load(here("analysis/data/spawn_rec_params.rda"))

#-----------------------
# estimate capacities
sr_p = spawn_recr_data %>%
  filter(!is.na(Spawners), !is.na(Parr)) %>%
  ggplot(aes(x = Spawners,
             y = Parr,
             shape = Spawner_type)) +
  facet_wrap(~ Population, 
             scales = "free") +
  stat_smooth(method = 'nls',
              se = F,
              formula = y ~ x / (a + b * x),
              method.args = list(start = list(a = 2e-3, b = 4e-5),
                                 lower = c(0, 0),
                                 algorithm = 'port'),
              fullrange = T,
              lwd = 1.5,
              aes(color = 'Beverton Holt')) +
  stat_smooth(method = 'nls',
              se = F,
              formula = y ~ (a*x) * exp(-b * x),
              method.args = list(start = list(a = 200, b = 1/2e5),
                                 lower = c(0, 0),
                                 algorithm = 'port'),
              fullrange = T,
              lwd = 1.5,
              aes(color = 'Ricker')) +
  stat_smooth(method = 'nls',
              se = F,
              formula = y ~ (a*x) / (1 + (x/b)^c),
              method.args = list(start = list(a = 600, b = 200, c = 1.5),
                                 lower = c(0, 0, 0),
                                 control = list(maxiter = 100, tol = 2e-05),
                                 algorithm = 'port'),
              fullrange = T,
              lwd = 1.5,
              aes(color = 'Shepard')) +
  stat_smooth(method = 'nls',
              se = F,
              formula = y ~ b * (1 - exp(-a / b * x)),
              method.args = list(start = list(a = 1.5e3, b = 1e5),
                                 lower = c(0, 0),
                                 algorithm = 'port'),
              fullrange = T,
              lwd = 1.5,
              aes(color = 'Hockey Stick')) +
  stat_smooth(method = lm,
              formula = y ~ -1 + x,
              fullrange = T,
              lty = 2,
              se = F,
              aes(color = 'Linear')) +
  geom_point(size = 3) +
  scale_color_brewer(palette = 'Set1',
                     breaks = c('Beverton Holt',
                                'Ricker',
                                'Hockey Stick',
                                'Shepard',
                                'Linear'),
                     name = 'Model') +
  labs(shape = 'Spawner\nType')
sr_p

#------------------------------------
# fit various spawner-recruit curves
#------------------------------------
library(msm)
# the Beverton-Holt function from FSA
bh1 = srFuns("BevertonHolt", param = 3)
# fit Beverton-Holt for all populations
bevHolt_params = spawn_recr_data %>%
  group_by(Population) %>%
  nest() %>%
  mutate(form = 'BevertonHolt',
         inits = map(data,
                     .f = function(x) {
                       inits = srStarts(Parr ~ Spawners,
                                        data = x,
                                        type = 'BevertonHolt',
                                        param = 3)
                       
                       inits$a = if_else(inits$a < 0, 
                                         2e-3, 
                                         inits$a)
                       inits$b = if_else(inits$b < 0, 
                                         1 / max(x$Parr), 
                                         inits$b)
                       return(inits)
                     }),
         mod_fit = map2(.x = data,
                        .y = inits,
                        .f = function(x, y) {
                          fit = try(nls(log(Parr) ~ log(bh1(Spawners, a, b)),
                                        data = x,
                                        start = y,
                                        algorithm = 'port',
                                        lower = c(0, 0)))
                        }),
         coefs = map(mod_fit,
                     .f = function(x) {
                       if(class(x) == 'try-error') {
                         return(as.numeric(NA))
                       } else coef(x)
                     }),
         est = map_dbl(coefs,
                       .f = function(x) {
                         1 / (x)['b']
                       }),
         se = map_dbl(mod_fit,
                      .f = function(x) {
                        if(class(x) == 'try-error') {
                          return(NA)
                        } else deltamethod(~ 1 / x1,
                                           coef(x)['b'],
                                           vcov(x)[2,2])
                      }),
         preds = map2(.x = data,
                      .y = mod_fit,
                      .f = function(x, y) {
                        if(class(y) != 'try-error') {
                          tibble(Spawners = 1:max(x$Spawners)) %>%
                            mutate(Parr = predict(y,
                                                  newdata = .),
                                   Parr = exp(Parr))
                        } else return(tibble(Spawners = NA,
                                             Parr = NA))
                      })) %>%
  mutate(cv = se / est)

# fit Ricker for all populations
# the Ricker function from FSA
ric1 = srFuns("Ricker", param = 1)
ricker_params = spawn_recr_data %>%
  group_by(Population) %>%
  nest() %>%
  mutate(form = 'Ricker',
         inits = map(data,
                     .f = function(x) {
                       inits = srStarts(Parr ~ Spawners,
                                        data = x,
                                        type = 'Ricker',
                                        param = 1)
                       inits$b = if_else(inits$b < 0, 
                                         1 / max(x$Parr), 
                                         inits$b)
                       return(inits)
                     }),
         mod_fit = map2(.x = data,
                        .y = inits,
                        .f = function(x, y) {
                          fit = try(nls(log(Parr) ~ log(ric1(Spawners, a, b)),
                                        data = x,
                                        start = y,
                                        algorithm = 'port',
                                        lower = c(0, 0)))
                        }),
         coefs = map(mod_fit,
                     .f = function(x) {
                       if(class(x) == 'try-error') {
                         return(as.numeric(NA))
                       } else coef(x)
                     }),
         est = map_dbl(coefs,
                       .f = function(x) {
                         x[1] / x[2] * exp(-1)
                       }),
         se = map_dbl(mod_fit,
                      .f = function(x) {
                        if(class(x) == 'try-error') {
                          return(NA)
                        } else deltamethod(~ x1 / x2 * exp(-1),
                                           coef(x),
                                           vcov(x))
                      }),
         preds = map2(.x = data,
                      .y = mod_fit,
                      .f = function(x, y) {
                        if(class(y) != 'try-error') {
                          tibble(Spawners = 1:max(x$Spawners)) %>%
                            mutate(Parr = predict(y,
                                                  newdata = .),
                                   Parr = exp(Parr))
                        } else return(tibble(Spawners = NA,
                                             Parr = NA))
                      })) %>%
  mutate(cv = se / est)

# fit hockey stick for all populations
hoc1 = function(S, a, b) {
  b * (1 - exp(-a / b * S))
}
hockey_params = spawn_recr_data %>%
  group_by(Population) %>%
  nest() %>%
  mutate(form = 'Hockey',
         inits = map(data,
                     .f = function(x) {
                       inits = list(a = coef(lm(Parr ~ -1 + Spawners, x)),
                                    b = max(x$Parr, na.rm = T))
                       return(inits)
                     }),
         mod_fit = map2(.x = data,
                        .y = inits,
                        .f = function(x, y) {
                          fit = try(nls(log(Parr) ~ log(hoc1(Spawners, a, b)),
                                        data = x,
                                        start = y,
                                        algorithm = 'port'))
                        }),
         coefs = map(mod_fit,
                     .f = function(x) {
                       if(class(x) == 'try-error') {
                         return(as.numeric(NA))
                       } else coef(x)
                     }),
         est = map_dbl(coefs,
                       .f = function(x) {
                         if(is.na(x[1])) {
                           return(as.numeric(NA))
                         } else x[2]
                       }),
         se = map_dbl(mod_fit,
                      .f = function(x) {
                        if(class(x) == 'try-error') {
                          return(as.numeric(NA))
                        } else summary(x)$coefficients['b', 'Std. Error']
                      }),
         preds = map2(.x = data,
                      .y = mod_fit,
                      .f = function(x, y) {
                        if(class(y) != 'try-error') {
                          tibble(Spawners = 1:max(x$Spawners)) %>%
                            mutate(Parr = predict(y,
                                                  newdata = .),
                                   Parr = exp(Parr))
                        } else return(tibble(Spawners = NA,
                                             Parr = NA))
                      })) %>%
  mutate(cv = se / est)

# compile all fitted parameters
spawn_rec_params = bevHolt_params %>%
  bind_rows(ricker_params) %>%
  bind_rows(hockey_params) %>%
  ungroup() %>%
  # drop this because it's causing problems with package build
  select(-mod_fit)

# save to be used later
save(spawn_rec_params,
     file = here("analysis/data/spawn_rec_params.rda"))

# pull out estimates of capacity
cap_tbl = spawn_rec_params %>%
  select(Population, form, est, se, cv) %>%
  arrange(Population, form)

# pull out fitted spawner recruit curves
curve_fits = spawn_rec_params %>%
  select(Population, form, preds) %>%
  unnest() %>%
  arrange(Population, form)

#------------------------------------
# try making a plot
#------------------------------------
# define polygons for uncertainty shading
plot_max = Inf
plot_max = 5e5
cap_poly = cap_tbl %>%
  mutate(lwrCI = est + se * qnorm(0.025),
         uprCI = est + se * qnorm(0.975)) %>%
  mutate(lwrCI = if_else(lwrCI < 0, 0, lwrCI),
         uprCI = if_else(uprCI > plot_max, plot_max, uprCI)) %>%
  left_join(spawn_rec_params %>%
              select(Population, data) %>%
              unnest() %>%
              group_by(Population) %>%
              summarise_at(vars(Spawners, Parr),
                           list(min = min,
                                max = max),
                           na.rm = T)) %>%
  # select(Population, form, ends_with('min'), ends_with('max'), ends_with('CI'))
  group_by(Population, form) %>%
  summarise(coords = list(tibble(x = c(0, 0, Spawners_max, Spawners_max),
                                 y = c(lwrCI, uprCI, uprCI, lwrCI)))) %>%
  ungroup() %>%
  unnest()

spawn_rec_params %>%
  # filter out a couple populations with super poor fits
  filter(!Population %in% c('Methow R.', 'Tucannon River')) %>%
  select(Population, data) %>%
  unnest() %>%
  distinct() %>%
  ggplot(aes(x = Spawners,
             y = Parr)) +
  geom_polygon(data = cap_poly %>%
                 # filter out a couple populations with super poor fits
                 filter(!Population %in% c('Methow R.', 'Tucannon River')),
               aes(x = x,
                   y = y,
                   fill = form),
               alpha = 0.2) +
  geom_line(data = curve_fits %>%
              # filter out a couple populations with super poor fits
              filter(!Population %in% c('Methow R.', 'Tucannon River')),
            aes(color = form),
            lwd = 1.5) +
  geom_hline(data = cap_tbl %>%
               # filter out a couple populations with super poor fits
               filter(!Population %in% c('Methow R.', 'Tucannon River')),
             aes(yintercept = est,
                 color = form),
             linetype = 2) +
  geom_point(aes(shape = Spawner_type),
             size = 3) +
  facet_wrap(~ Population,
             scales = 'free') +
  scale_color_brewer(palette = 'Set1',
                     breaks = c('BevertonHolt',
                                'Ricker',
                                'Hockey'),
                     labels = c('Beverton Holt',
                                'Ricker',
                                'Hockey Stick'),
                     name = 'Model') +
  labs(shape = 'Spawner\nType')
