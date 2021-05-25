setwd("~/GitHub/EMSI-2021/epidemiology-lab")

renv::init()

## Using R package: EpiModel
##   http://www.epimodel.org/
install.packages("EpiModel")

library(EpiModel)

#######################################################
####    Model 1: “Childhood Disease” SIR Model
#######################################################

param <- param.dcm(inf.prob = 0.5, act.rate = 0.25, rec.rate = 0.01)
init <- init.dcm(s.num = 600, i.num = 1, r.num = 400)
control <- control.dcm(type = "SIR", nsteps = 600)

mod1 <- dcm(param, init, control)

mod1
# EpiModel Simulation
# =======================
#   Model class: dcm
# 
# Simulation Summary
# -----------------------
#   Model type: SIR
# No. runs: 1
# No. time steps: 600
# No. groups: 1
# 
# Model Parameters
# -----------------------
#   inf.prob = 0.5
# act.rate = 0.25
# rec.rate = 0.01
# 
# Model Output
# -----------------------
#   Variables: s.num i.num r.num si.flow ir.flow num

plot(mod1)  ## F1

## changing the model
param <- param.dcm(inf.prob = 0.5, act.rate = 0.25, rec.rate = seq(0.01, 0.03, 0.005))
init <- init.dcm(s.num = 500, i.num = 1, r.num = 400)
control <- control.dcm(type = "SIR", nsteps = 600)

mod2 <- dcm(param, init, control)

mod2
# EpiModel Simulation
# =======================
#   Model class: dcm
# 
# Simulation Summary
# -----------------------
#   Model type: SIR
# No. runs: 5
# No. time steps: 600
# No. groups: 1
# 
# Model Parameters
# -----------------------
#   inf.prob = 0.5
# act.rate = 0.25
# rec.rate = 0.01 0.015 0.02 0.025 0.03
# 
# Model Output
# -----------------------
#   Variables: s.num i.num r.num si.flow ir.flow num

plot(mod2) ## F2


## setting up multi-panel plots
## F3
par(mfrow = c(1,2)) # this sets up the plotting window, giving two panels
plot(mod2, y = "i.num", popfrac = TRUE, alpha = 1, col = "Blues", legend = "full", main = "Disease Prevalence")
plot(mod2, y = "si.flow", col = "Greens", alpha = 0.8, main = "Disease Incidence", legend = "full")


## creating a figure of the model diagram at step 33
## F4
par(mfrow = c(1, 1))
comp_plot(mod2, at = 33, digits = 2)

mod2$epi$si.flow[33,1]
# [1] 0.4661749
mod2$epi$si.flow[66,1]
# [1] 2.720477

#######################################################
####   Demographic Models: Births and Deaths
#######################################################

param <- param.dcm(inf.prob = 1/2, act.rate = 1/4, rec.rate = 1/100, 
                   a.rate = 1/80, ds.rate = 1/100, di.rate = 1/100, dr.rate = 1/100)
init <- init.dcm(s.num = 600, i.num = 1, r.num = 400)
control <- control.dcm(type = "SIR", nsteps = 600)

mod3 <- dcm(param, init, control)

plot(mod3) ## F5

## changing virulence
param <- param.dcm(inf.prob = 1/2, act.rate = 1/4, rec.rate = 1/100, 
                   a.rate = 1/80, ds.rate = 1/100, di.rate = 1/40, dr.rate = 1/100)
init <- init.dcm(s.num = 600, i.num = 1, r.num = 400)
control <- control.dcm(type = "SIR", nsteps = 600)

mod4 <- dcm(param, init, control)

plot(mod4) ## F6

param <- param.dcm(inf.prob = 0.5, act.rate = 0.25, rec.rate = 0.01, a.rate = 1/80, ds.rate = 1/100, di.rate = seq(0.01, 0.05, 0.01), dr.rate = 1/100)
init <- init.dcm(s.num = 500, i.num = 1, r.num = 500)
control <- control.dcm(type = "SIR", nsteps = 600)

mod5 <- dcm(param, init, control)

## F7
par(mfrow = c(1,1))
plot(mod5, y = "si.flow", col = "Greens", alpha = 0.8, main = "Disease Incidence")


#######################################################
####   Stochastic Models
#######################################################

param <- param.icm(inf.prob = 0.5, act.rate = 0.25, rec.rate = 0.01)
init <- init.icm(s.num = 600, i.num = 1, r.num = 400)
control <- control.icm(type = "SIR", nsims = 20, nsteps = 300)
mod6 <- icm(param, init, control)

mod6
# EpiModel Simulation
# =======================
#   Model class: icm
# 
# Simulation Summary
# -----------------------
#   Model type: SIR
# No. simulations: 20
# No. time steps: 300
# No. groups: 1
# 
# Model Parameters
# -----------------------
#   inf.prob = 0.5
# act.rate = 0.25
# rec.rate = 0.01
# 
# Model Output
# -----------------------
#   Variables: s.num i.num num r.num si.flow ir.flow

summary(mod6, at = 33)
# EpiModel Summary
# =======================
#   Model class: icm
# 
# Simulation Details
# -----------------------
#   Model type: SIR
# No. simulations: 20
# No. time steps: 300
# No. groups: 1
# 
# Model Statistics
# ------------------------------
#   Time: 33 
# ------------------------------ 
#   mean     sd    pct
# Suscept.   594.65  7.110  0.594
# Infect.      5.50  6.621  0.005
# Recov.     400.85  0.988  0.400
# Total     1001.00  0.000  1.000
# S -> I       0.25  0.550     NA
# I -> R       0.05  0.224     NA
# ------------------------------

plot(mod6) ## F8

plot(mod6, qnts = FALSE, sim.lines = TRUE, mean.line = FALSE) ## F9

plot(mod6, y = "i.num", qnts = FALSE, sim.lines = TRUE, mean.line = FALSE) ## F10




