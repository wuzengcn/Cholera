stochastic_model <- function(n.trials=1, n.t.end=1, theta, n.schedule=1, herd=1) {
  
  #random variables
  #efficacy of two doses of vaccine
  #efficacy of one dose of vaccine
  #cost of treatment  
  vacdatastore<-array(NA, dim = c(n.t.end, 5, n.trials), dimnames = list(NULL,NULL,NULL))  #order: n.years, n.scenario,, n.simulation,  store data on coverage of vac
  datastore<- matrix(rep(NA, n.trials*90), nrow = n.trials)    
  for (ii in 1:n.trials) {
    set.seed(ii+1)
    simvaceff1 <- rbeta(1, theta["vac.onedose.efficacy"]^2*(1-theta["vac.onedose.efficacy"])/(theta["vac.onedose.eff.sd"])^2+
                          theta["vac.onedose.efficacy"]^2/(1-theta["vac.onedose.efficacy"]) - theta["vac.onedose.efficacy"]/(1-theta["vac.onedose.efficacy"]),
                        theta["vac.onedose.efficacy"]*(1-theta["vac.onedose.efficacy"])^2/(theta["vac.onedose.eff.sd"])^2 + 
                          theta["vac.onedose.efficacy"] - 1)
    simvaceff2 <- rbeta(1, theta["vac.twodose.efficacy"]^2*(1-theta["vac.twodose.efficacy"])/(theta["vac.twodose.eff.sd"])^2+
                          theta["vac.twodose.efficacy"]^2/(1-theta["vac.twodose.efficacy"]) - theta["vac.twodose.efficacy"]/(1-theta["vac.twodose.efficacy"]),
                        theta["vac.twodose.efficacy"]*(1-theta["vac.twodose.efficacy"])^2/(theta["vac.twodose.eff.sd"])^2 + 
                          theta["vac.twodose.efficacy"] - 1)
    simcovdif <- rbeta(1, theta["vac.twodose.cov.dif"]^2*(1-theta["vac.twodose.cov.dif"])/(theta["vac.twodose.cov.dif.sd"])^2+
                         theta["vac.twodose.cov.dif"]^2/(1-theta["vac.twodose.cov.dif"]) - theta["vac.twodose.cov.dif"]/(1-theta["vac.twodose.cov.dif"]),
                       theta["vac.twodose.cov.dif"]*(1-theta["vac.twodose.cov.dif"])^2/(theta["vac.twodose.cov.dif.sd"])^2 + 
                         theta["vac.twodose.cov.dif"] - 1)
  
    simctreat1 <- simctreat2 <- simctreat3 <- rgamma(1,shape = theta["c.treat.1"]^2/theta["c.treat.1.sd"]^2 , 
                                                     scale= theta["c.treat.1.sd"]^2/theta["c.treat.1"])
    
    simctml <- rgamma(1,shape = theta["c.tml"]^2/theta["c.tml.sd"]^2 , 
                                                     scale= theta["c.tml.sd"]^2/theta["c.tml"])
    simcploss <- rgamma(1,shape = theta["c.ploss"]^2/theta["c.ploss.sd"]^2 , 
                                                  scale= theta["c.ploss.sd"]^2/theta["c.ploss"])
    simtcholera <- rgamma(1,shape = theta["t.cholera"]^2/theta["t.cholera.sd"]^2 , 
                        scale= theta["t.cholera.sd"]^2/theta["t.cholera"])
    
    simdw <- rbeta(1, theta["d.cholera"]^2*(1-theta["d.cholera"])/(theta["d.cholera.sd"])^2+
                     theta["d.cholera"]^2/(1-theta["d.cholera"]) - theta["d.cholera"]/(1-theta["d.cholera"]),
                   theta["d.cholera"]*(1-theta["d.cholera"])^2/(theta["d.cholera.sd"])^2 + 
                     theta["d.cholera"] - 1)
    
    siminc1 <- rbeta(1, theta["incidence.1"]^2*(1-theta["incidence.1"])/(theta["incidence.1.sd"])^2+
                       theta["incidence.1"]^2/(1-theta["incidence.1"]) - theta["incidence.1"]/(1-theta["incidence.1"]),
                     theta["incidence.1"]*(1-theta["incidence.1"])^2/(theta["incidence.1.sd"])^2 + 
                       theta["incidence.1"] - 1)
    siminc2 <- rbeta(1, theta["incidence.2"]^2*(1-theta["incidence.2"])/(theta["incidence.2.sd"])^2+
                       theta["incidence.2"]^2/(1-theta["incidence.2"]) - theta["incidence.2"]/(1-theta["incidence.2"]),
                     theta["incidence.2"]*(1-theta["incidence.2"])^2/(theta["incidence.2.sd"])^2 + 
                       theta["incidence.2"] - 1)
    siminc3 <- rbeta(1, theta["incidence.3"]^2*(1-theta["incidence.3"])/(theta["incidence.3.sd"])^2+
                       theta["incidence.3"]^2/(1-theta["incidence.3"]) - theta["incidence.3"]/(1-theta["incidence.3"]),
                     theta["incidence.3"]*(1-theta["incidence.3"])^2/(theta["incidence.3.sd"])^2 + 
                       theta["incidence.3"] - 1)
    simfatality<- rbeta(1, theta["p.casefatalityrate.1"]^2*(1-theta["p.casefatalityrate.1"])/(theta["p.casefatalityrate.1.sd"])^2+
                          theta["p.casefatalityrate.1"]^2/(1-theta["p.casefatalityrate.1"]) - theta["p.casefatalityrate.1"]/(1-theta["p.casefatalityrate.1"]),
                        theta["p.casefatalityrate.1"]*(1-theta["p.casefatalityrate.1"])^2/(theta["p.casefatalityrate.1.sd"])^2 + 
                          theta["p.casefatalityrate.1"] - 1)
    
    combine <-process_model(n.t.end, theta, d.twodose.1=1, d.onedose.2.3=1, d.twodose.2.3 =0, n.schedule,
                            simvaceff1, simvaceff2, simcovdif, simctreat1, simctreat2, simctreat3,simctml, simcploss, simtcholera, siminc1, siminc2, 
                            siminc3, simfatality,simdw, herd)
    twodose <-process_model(n.t.end, theta, d.twodose.1=1, d.onedose.2.3=0, d.twodose.2.3 =1, n.schedule,
                            simvaceff1, simvaceff2, simcovdif,simctreat1, simctreat2, simctreat3, simctml, simcploss, simtcholera,siminc1, siminc2, 
                            siminc3, simfatality,simdw, herd)
    twodose1 <-process_model(n.t.end, theta, d.twodose.1=1, d.onedose.2.3=0, d.twodose.2.3 =0, n.schedule,
                             simvaceff1, simvaceff2, simcovdif,simctreat1, simctreat2, simctreat3, simctml, simcploss, simtcholera, siminc1, siminc2, 
                             siminc3, simfatality,simdw, herd)
    onedose <-process_model(n.t.end, theta, d.twodose.1=0, d.onedose.2.3=1, d.twodose.2.3 =0, n.schedule,
                            simvaceff1, simvaceff2, simcovdif, simctreat1, simctreat2, simctreat3, simctml, simcploss, simtcholera, siminc1, siminc2, 
                            siminc3, simfatality,simdw, herd)
    control <-process_model(n.t.end, theta, d.twodose.1=0, d.onedose.2.3=0, d.twodose.2.3 =0, n.schedule,
                            simvaceff1, simvaceff2, simcovdif, simctreat1, simctreat2, simctreat3, simctml, simcploss, simtcholera, siminc1, siminc2, 
                            siminc3, simfatality,simdw,herd)
    #Calculate number of vaccinated number and number of population in each year and vaccine coverage 
    vacdatastore[ ,1, ii]<-control[["v.vac.coverage"]]
    vacdatastore[ ,2, ii]<-onedose[["v.vac.coverage"]]
    vacdatastore[ ,3, ii]<-twodose[["v.vac.coverage"]]
    vacdatastore[ ,4, ii]<-combine[["v.vac.coverage"]]
    vacdatastore[ ,5, ii]<-twodose1[["v.vac.coverage"]]
    #calculate number of vaccine used and vaccine cost, treatment cost and productivity loss 
    
    #calculate effectiveness: number of nature deaths, number death from cholera, number of infections
    
    #Generate summary cost and effectiveness     
    dif.c.onedose.hs <- onedose[["cost.hs"]] - control[["cost.hs"]]
    dif.c.onedose.so <- onedose[["cost.so"]] - control[["cost.so"]]
    dif.e.onedose <- control[["daly"]] - onedose[["daly"]]
    icer.onedose.hs <- dif.c.onedose.hs/dif.e.onedose
    icer.onedose.so <- dif.c.onedose.so/dif.e.onedose
    
    dif.c.twodose.hs <- twodose[["cost.hs"]] - control[["cost.hs"]]
    dif.c.twodose.so <- twodose[["cost.so"]] - control[["cost.so"]]
    dif.e.twodose <- control[["daly"]] - twodose[["daly"]]
    icer.twodose.hs <- dif.c.twodose.hs/dif.e.twodose
    icer.twodose.so <- dif.c.twodose.so/dif.e.twodose
    
    dif.c.combine.hs <- combine[["cost.hs"]] - control[["cost.hs"]]
    dif.c.combine.so <- combine[["cost.so"]] - control[["cost.so"]]
    dif.e.combine <- control[["daly"]] - combine[["daly"]]
    icer.combine.hs <- dif.c.combine.hs/dif.e.combine
    icer.combine.so <- dif.c.combine.so/dif.e.combine
    
    dif.c.twodose1.hs <- twodose1[["cost.hs"]] - control[["cost.hs"]]
    dif.c.twodose1.so <- twodose1[["cost.so"]] - control[["cost.so"]]
    dif.e.twodose1 <- control[["daly"]] - twodose1[["daly"]]
    icer.twodose1.hs <- dif.c.twodose1.hs/dif.e.twodose1
    icer.twodose1.so <- dif.c.twodose1.so/dif.e.twodose1    
    
    
    tempdta <- c(control[["cost.tr.hs"]], control[["cost.vac.hs"]], control[["cost.hs"]],
                 control[["cost.tr.so"]], control[["cost.vac.so"]], control[["cost.so"]], control[["cost.death"]], control[["cost.tot.so"]], control[["daly"]], 
                 onedose[["cost.tr.hs"]], onedose[["cost.vac.hs"]], onedose[["cost.hs"]],
                 onedose[["cost.tr.so"]], onedose[["cost.vac.so"]], onedose[["cost.so"]], onedose[["cost.death"]], onedose[["cost.tot.so"]], onedose[["daly"]], 
                 twodose[["cost.tr.hs"]], twodose[["cost.vac.hs"]], twodose[["cost.hs"]],
                 twodose[["cost.tr.so"]], twodose[["cost.vac.so"]], twodose[["cost.so"]], twodose[["cost.death"]], twodose[["cost.tot.so"]], twodose[["daly"]],                        
                 combine[["cost.tr.hs"]], combine[["cost.vac.hs"]], combine[["cost.hs"]],
                 combine[["cost.tr.so"]], combine[["cost.vac.so"]], combine[["cost.so"]], combine[["cost.death"]], combine[["cost.tot.so"]], combine[["daly"]], 
                 twodose1[["cost.tr.hs"]], twodose1[["cost.vac.hs"]], twodose1[["cost.hs"]],
                 twodose1[["cost.tr.so"]], twodose1[["cost.vac.so"]], twodose1[["cost.so"]], twodose1[["cost.death"]], twodose1[["cost.tot.so"]], twodose1[["daly"]],                   
                 dif.c.onedose.hs, dif.c.onedose.so, dif.e.onedose, icer.onedose.hs, icer.onedose.so, 
                 dif.c.twodose.hs, dif.c.twodose.so, dif.e.twodose, icer.twodose.hs, icer.twodose.so, 
                 dif.c.combine.hs, dif.c.combine.so, dif.e.combine, icer.combine.hs, icer.combine.so,
                 dif.c.twodose1.hs, dif.c.twodose1.so, dif.e.twodose1, icer.twodose1.hs, icer.twodose1.so,
                 control[["v.process.ind"]], onedose[["v.process.ind"]], twodose[["v.process.ind"]], 
                 combine[["v.process.ind"]],twodose1[["v.process.ind"]])
    
    datastore[ii,] <- tempdta
  }
  simdata.vac <- vacdatastore
  temp1<-apply(simdata.vac[,1,], 1, mean)
  temp2<-apply(simdata.vac[,2,], 1, mean)
  temp3<-apply(simdata.vac[,3,], 1, mean)
  temp4<-apply(simdata.vac[,4,], 1, mean)
  temp5<-apply(simdata.vac[,5,], 1, mean)
  temp6<-apply(simdata.vac[,1,], 1, sd)
  temp7<-apply(simdata.vac[,2,], 1, sd)
  temp8<-apply(simdata.vac[,3,], 1, sd)
  temp9<-apply(simdata.vac[,4,], 1, sd)
  temp10<-apply(simdata.vac[,5,], 1, sd)
  
  tempall <- cbind(temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10)
  simdata.vac.d <- data.frame(tempall)
  colnames(simdata.vac.d) = c("control.vac.cov.m", "onedose.vac.cov.m", "twodose.vac.cov.m",
                              "combine.vac.cov.m", "twodose1.vac.cov.m", "control.vac.cov.sd",
                              "onedose.vac.cov.sd", "twodose.vac.cov.sd",
                              "combine.vac.cov.sd", "twodose1.vac.cov.sd")
  
  simdata.outcome <- data.frame(datastore)
  colnames(simdata.outcome) = c("control_tr_hs","control_vac_hs","control_cost_hs", 
                        "control_tr_so","control_vac_so","control_cost_so","control_cost_death_so","control_cost_tot_so","control_daly", 
                        "onedose_tr_hs","onedose_vac_hs","onedose_cost_hs", 
                        "onedose_tr_so","onedose_vac_so","onedose_cost_so","onedose_cost_death_so","onedose_cost_tot_so","onedose_daly",
                        "twodose_tr_hs","twodose_vac_hs","twodose_cost_hs", 
                        "twodose_tr_so","twodose_vac_so","twodose_cost_so","twodose_cost_death_so","twodose_cost_tot_so", "twodose_daly",
                        "combine_tr_hs","combine_vac_hs","combine_cost_hs", 
                        "combine_tr_so","combine_vac_so","combine_cost_so","combine_cost_death_so","combine_cost_tot_so", "combine_daly",
                        "twodose1_tr_hs","twodose1_vac_hs","twodose1_cost_hs", 
                        "twodose1_tr_so","twodose1_vac_so","twodose1_cost_so","twodose1_cost_death_so","twodose1_cost_tot_so","twodose1_daly",                      
                        "dif_cost_one_hs", "dif_cost_one_so", "dif.daly.one",
                        "icer.one.hs", "icer.one.so","dif_cost_two_hs", 
                        "dif_cost_two_so", "dif.daly.two","icer.two.hs", 
                        "icer.two.so", "dif_cost_combine_hs", "dif_cost_combine_so",
                        "dif.combine.one","icer.combine.hs", "icer.combine.so",
                        "dif_cost_two1_hs", "dif_cost_two1_so", "dif.daly.two1",
                        "icer.two1.hs", "icer.two1.so", "control_ninf", "control_ndeath", "control_cdeath",
                        "control_n.onedose", "control_n.twodose", 
                        "onedose_ninf", "onedose_ndeath", "onedose_cdeath", "onedose_n.onedose", "onedose_n.twodose",
                        "twodose_ninf","twodose_ndeath", "twodose_cdeath","twodose_n.onedose", "twodose_n.twodose",
                        "combine_ninf", "combine_ndeath", "combine_cdeath","combine_n.onedose", "combine_n.twodose",
                        "twodose1_ninf", "twodose1_ndeath", "twodose1_cdeath", "twodose1_n.onedose", "twodose1_n.twodose"
                        )
  
  simres<- list(simdata.outcome=simdata.outcome,simdata.vac.d=simdata.vac.d)
  return(simres)
}

#function to generate mean and standard deviation of variables 
desc<- function(x){
  temp1<-mean(x)
  temp2<-sd(x)
  return(c(temp1, temp2))
}

#function to combine cost and effectiveness information
cedesc <-function(x1,x2,x3,x4,x5,x6,x7,x8,y1,y2,y3,y4){
  c1 <- desc(x1)
  c2 <- desc(x2)
  c3 <- desc(x3)
  c4 <- desc(x4)
  c5 <- desc(x5)
  c6 <- desc(x6)
  c7 <- desc(x7)
  c8 <- desc(x8)
  
  c.x <- rbind(c1,c2,c3,c4,c5,c6,c7,c8)
  
  e1 <- desc(y1)
  e2 <- desc(y2)
  e3 <- desc(y3)
  e4 <- desc(y4)
  e.x <- rbind(e1,e2,e3,e4)
  return(rbind(c.x,e.x))
}

# function to generate acceptability curve
genaccurve <- function(x, n.trials=1000) {
  max1<- ceiling(max(x))
  min1<- floor(min(x))
  p<-matrix(rep(NA,((max1-min1)/1+1)*2),nrow = ((max1-min1)/1+1))  
  for (wtp in seq(from=min1, to=max1, by=1)) {
    nbelow <- sum(x < wtp)
    pbelow <- nbelow/n.trials
    pbelowr <- c(wtp, pbelow)
    p[(wtp-min1)/1+1,]<-pbelowr
  }
  acceptability <- data.frame(p)
  acceptability$X3 <- 1-acceptability$X2
  c<-subset(acceptability,X1==2784)                 #2783.6 is 1.5 GDP per capital in Bangaldesh in 2019. 
  x<- ggplot(data = acceptability, aes(x=X1)) +
    geom_line(aes(y=X2), color = "blue") +
    geom_line(aes(y=X3), color="black") +
    geom_segment(x=c[1,1], y=0, xend=c[1,1], yend=c[1,2], linetype = "dashed",
                 color= "darkred") +
    geom_text(x=c[1,1]+200, y=c[1,2]-0.1, label=c[1,2]) +
    #  geom_segment(x=0, y=c[1,2], xend=3000, yend=c[1,2]) +
    theme_classic() + 
    labs(x="Willingness to pay ($/DALY)", y="Probability of acceptance")
return(x)
}
