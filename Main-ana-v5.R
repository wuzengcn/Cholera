# Choloera stochastic model
# Main script
# Author: Wu Zeng (2020)
# Set up libraries and paths ----------------------------------------------
library(readr)
library(foreach)
library(doMC)
library(lubridate)
library(magrittr)
library(coda)
library(tidyverse)
library(rootSolve)
library(mgcv)
library(ggplot2)
library(ggpubr)
library("lattice")
library(viridisLite)
#setwd("/Users/Haksoon/OneDrive - brandeis.edu/Cholera/Code")
#setwd("/Volumes/GoogleDrive/My Drive/GW/Proposal/Cholera/Code")
#setwd("/Users/haksoonahn/Google Drive/GW/Proposal/Cholera/Code")
setwd("/Users/haksoonahn/OneDrive - brandeis.edu/Res/Cholera/Code")

rm(list=ls(all=TRUE))
options(scipen = 999,digits = 9)

# Load model and plotting functions
source("Herd.R")
source("nProcess-v8.R")
source("nStochastic-v6.R")
source("nSensitivity-v4.R")

# Load model parameters
thetaR_IC <- read_csv("input/parameters-v2.csv")
theta <- c( pop=as.numeric(thetaR_IC[thetaR_IC$param=="population","value"]), # note this is only IC - SMC estimates this
            pople5=as.numeric(thetaR_IC[thetaR_IC$param=="pople5","value"]),
            pop5to15 = as.numeric(thetaR_IC[thetaR_IC$param=="pop5to15","value"]),
            popgt15 = as.numeric(thetaR_IC[thetaR_IC$param=="popgt15","value"]),
            birth.rate = as.numeric(thetaR_IC[thetaR_IC$param=="birth.rate","value"]),
            mature1.rate = as.numeric(thetaR_IC[thetaR_IC$param=="mature1.rate","value"]),
            mature2.rate=as.numeric(thetaR_IC[thetaR_IC$param=="mature2.rate","value"]),
            incidence.1=as.numeric(thetaR_IC[thetaR_IC$param=="incidence.1","value"]),
            incidence.2 = as.numeric(thetaR_IC[thetaR_IC$param=="incidence.2","value"]),
            incidence.3 = as.numeric(thetaR_IC[thetaR_IC$param=="incidence.3","value"]),
            incidence.1.sd=as.numeric(thetaR_IC[thetaR_IC$param=="incidence.1","sd"]),
            incidence.2.sd = as.numeric(thetaR_IC[thetaR_IC$param=="incidence.2","sd"]),
            incidence.3.sd = as.numeric(thetaR_IC[thetaR_IC$param=="incidence.3","sd"]),            
            ndeath.1 =as.numeric(thetaR_IC[thetaR_IC$param=="ndeath.1","value"]), 
            ndeath.2 =as.numeric(thetaR_IC[thetaR_IC$param=="ndeath.2","value"]), 
            ndeath.3 = as.numeric(thetaR_IC[thetaR_IC$param=="ndeath.3","value"]),
            vac.onedose.efficacy = as.numeric(thetaR_IC[thetaR_IC$param=="vac.onedose.efficacy","value"]),
            vac.onedose.eff.sd = as.numeric(thetaR_IC[thetaR_IC$param=="vac.onedose.efficacy","sd"]),
            vac.onedose.cov = as.numeric(thetaR_IC[thetaR_IC$param=="vac.onedose.cov","value"]),
            vac.onedose.duration = as.numeric(thetaR_IC[thetaR_IC$param=="vac.onedose.duration","value"]),
            vac.twodose.efficacy = as.numeric(thetaR_IC[thetaR_IC$param=="vac.twodose.efficacy","value"]),
            vac.twodose.efficacy.dt.1 = as.numeric(thetaR_IC[thetaR_IC$param=="vac.twodose.efficacy.dt.1","value"]),
            vac.twodose.eff.sd = as.numeric(thetaR_IC[thetaR_IC$param=="vac.twodose.efficacy","sd"]),
            vac.twodose.cov.1 = as.numeric(thetaR_IC[thetaR_IC$param=="vac.twodose.cov.1","value"]),
            vac.twodose.cov.dif = as.numeric(thetaR_IC[thetaR_IC$param=="vac.twodose.cov.dif","value"]),
            vac.twodose.cov.dif.sd = as.numeric(thetaR_IC[thetaR_IC$param=="vac.twodose.cov.dif","sd"]),
            vac.twodose.duration = as.numeric(thetaR_IC[thetaR_IC$param=="vac.twodose.duration","value"]), 
            p.casefatalityrate.1 = as.numeric(thetaR_IC[thetaR_IC$param=="p.casefatalityrate.1","value"]),
            p.casefatalityrate.2 = as.numeric(thetaR_IC[thetaR_IC$param=="p.casefatalityrate.2","value"]),  
            p.casefatalityrate.3 = as.numeric(thetaR_IC[thetaR_IC$param=="p.casefatalityrate.3","value"]),
            p.casefatalityrate.1.sd = as.numeric(thetaR_IC[thetaR_IC$param=="p.casefatalityrate.1","sd"]),            
            p.loss.immunity = as.numeric(thetaR_IC[thetaR_IC$param=="p.loss.immunity","value"]),
            c.onedose = as.numeric(thetaR_IC[thetaR_IC$param=="c.onedose.vac","value"]),
            c.twodose = as.numeric(thetaR_IC[thetaR_IC$param=="c.twodose.vac","value"]),
            c.delivery = as.numeric(thetaR_IC[thetaR_IC$param=="c.delivery","value"]),
            c.treat.1 = as.numeric(thetaR_IC[thetaR_IC$param=="c.treatment.1","value"]),
            c.treat.1.sd = as.numeric(thetaR_IC[thetaR_IC$param=="c.treatment.1","sd"]),
            c.treat.2 = as.numeric(thetaR_IC[thetaR_IC$param=="c.treatment.2","value"]),
            c.treat.2.sd = as.numeric(thetaR_IC[thetaR_IC$param=="c.treatment.2","sd"]),
            c.treat.3 = as.numeric(thetaR_IC[thetaR_IC$param=="c.treatment.3","value"]),
            c.treat.3.sd = as.numeric(thetaR_IC[thetaR_IC$param=="c.treatment.2","sd"]),
            c.travel  = as.numeric(thetaR_IC[thetaR_IC$param=="c.travel","value"]),
            c.tml  = as.numeric(thetaR_IC[thetaR_IC$param=="c.tml","value"]),
            c.tml.sd  = as.numeric(thetaR_IC[thetaR_IC$param=="c.tml","sd"]),
            c.ploss  = as.numeric(thetaR_IC[thetaR_IC$param=="c.ploss","value"]),
            c.ploss.sd  = as.numeric(thetaR_IC[thetaR_IC$param=="c.ploss","sd"]),
            p.sym  = as.numeric(thetaR_IC[thetaR_IC$param=="c.sym","value"]),
            wastage = as.numeric(thetaR_IC[thetaR_IC$param=="wastage","value"]),
            d.cholera = as.numeric(thetaR_IC[thetaR_IC$param=="d.cholera","value"]),
            d.cholera.sd = as.numeric(thetaR_IC[thetaR_IC$param=="d.cholera","sd"]),
            t.cholera = as.numeric(thetaR_IC[thetaR_IC$param=="t.cholera","value"]),
            t.cholera.sd = as.numeric(thetaR_IC[thetaR_IC$param=="t.cholera","sd"]),
            lifeexp.1 = as.numeric(thetaR_IC[thetaR_IC$param=="lifeexp.1","value"]),
            lifeexp.2 = as.numeric(thetaR_IC[thetaR_IC$param=="lifeexp.2","value"]),
            lifeexp.3 = as.numeric(thetaR_IC[thetaR_IC$param=="lifeexp.3","value"]),
            dc = as.numeric(thetaR_IC[thetaR_IC$param=="dc","value"]),
            de = as.numeric(thetaR_IC[thetaR_IC$param=="de","value"])
)

theta_Names <- c("sus","vac","inf", "recovered", "ndeath","cdeath") # also defines groups to use in model

#Effection at analysis at mean
combine <-process_model(12, theta, d.twodose.1=1, d.onedose.2.3=1, d.twodose.2.3 =0, 3,
                       theta["vac.onedose.efficacy"], theta["vac.twodose.efficacy"], 
                       theta["vac.twodose.cov.dif"], theta["c.treat.1"], theta["c.treat.2"],
                       theta["c.treat.3"], theta["c.tml"], theta["c.ploss"], theta["t.cholera"] , theta[["incidence.1"]], theta[["incidence.2"]], 
                       theta[["incidence.3"]], theta[["p.casefatalityrate.1"]] ,theta["d.cholera"], herd = 1)
#twodose <-process_model(10, theta, d.twodose.1=1, d.onedose.2.3=0, d.twodose.2.3 =1, 3,
#                        theta["vac.onedose.efficacy"], theta["vac.twodose.efficacy"], theta["vac.twodose.cov.dif"],theta["c.treat.1"], theta["c.treat.2"], theta["c.treat.3"], theta["d.cholera"])
#twodose1 <-process_model(10, theta, d.twodose.1=1, d.onedose.2.3=0, d.twodose.2.3 =0, 3,
#                         theta["vac.onedose.efficacy"], theta["vac.twodose.efficacy"], theta["vac.twodose.cov.dif"],theta["c.treat.1"], theta["c.treat.2"], theta["c.treat.3"], theta["d.cholera"])
#onedose <-process_model(10, theta, d.twodose.1=0, d.onedose.2.3=1, d.twodose.2.3 =0, 3,
#                        theta["vac.onedose.efficacy"], theta["vac.twodose.efficacy"], theta["vac.twodose.cov.dif"], theta["c.treat.1"], theta["c.treat.2"], theta["c.treat.3"], theta["d.cholera"])
#control <-process_model(10, theta, d.twodose.1=0, d.onedose.2.3=0, d.twodose.2.3 =0, 3,
#                        theta["vac.onedose.efficacy"], theta["vac.twodose.efficacy"], theta["vac.twodose.cov.dif"], theta["c.treat.1"], theta["c.treat.2"], theta["c.treat.3"], theta["d.cholera"])
#Monte carlo simuation 
run1 <- stochastic_model(n.trials=1000, n.t.end=12, theta, n.schedule=3, herd=1)
#sensility analysis on projection period and schedule of vaccine 

sen1<- sensitivity_model() 
sendta <- subset(sen1[["sentdata"]],projectionperiod>4 & projectionperiod<19)
levelplot(icer.combine.so~schedule*projectionperiod, data = sendta,col.regions = viridis(100),
          ylab = "Projection period (years)", xlab = "Schedule")

save(sen1[["sentdata"]], file = "sensitivity.RData")
  
