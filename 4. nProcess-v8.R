#process model for the simulation ---

#six compartment model: suseptible, vaccinated, infectious, recovered, nature dealth, dealth from cholera

process_model <- function(n.t.end, theta, d.twodose.1=0, d.onedose.2.3=0, 
                          d.twodose.2.3 =0, n.schedule=1, 
                          simvaceff1, simvaceff2, simcovdif, simctreat1, 
                          simctreat2, simctreat3, simctml, simcploss, simtcholera, siminc1, siminc2, 
                          siminc3, simfatality, simdw, herd=1){
  #demographic information 
  n.pop <- theta["pop"]
  s.pop.1 <- theta["pople5"]
  s.pop.2 <- theta["pop5to15"]
  s.pop.3 <- theta["popgt15"]
  n.pop.1 <- n.pop*s.pop.1
  n.pop.2 <- n.pop*s.pop.2
  n.pop.3 <- n.pop*s.pop.3
  #transption probabilities and other parameters
  r.birth <- theta["birth.rate"]
  r.mature.1 <- theta["mature1.rate"]
  r.mature.2 <- theta["mature2.rate"]
  r.loss.immunity <- theta["p.loss.immunity"]
  #use of vaccine and vaccine schedule
  d.twodose.1 <- d.twodose.1  # use two doses of shanchol for children under 5 or not (0 no, 1 yes)
  d.onedose.2.3 <- d.onedose.2.3    # use one dose of shanchol for people above 5 or not (0 no, 1 yes) 
  d.twodose.2.3 <- d.twodose.2.3 # use two dose of shanchol for people above 5 or not (0 no, 1 yes) 
  n.schedule <- n.schedule   #vaccinationa schedule
  
  #incidence of cholera in three age groups
  p.incidence.1 <- siminc1
  p.incidence.2 <- siminc2
  p.incidence.3 <- siminc3
  
  #nature mortality rate in three age agroups
  p.ndeath.1 <- theta["ndeath.1"]
  p.ndeath.2 <- theta["ndeath.2"]
  p.ndeath.3 <- theta["ndeath.3"]
  
  #vaccine efficiancy, coverage and duration of protection 
  p.vac.one.eff <- simvaceff1      #efficiency of one dose shanchol
  p.vac.one.cov <- theta["vac.onedose.cov"]      #coverage of one dose shanchol
  n.vac.one.dur <- theta["vac.onedose.duration"]      #duration of the protection of one dose shanchol 
  p.vac.two.eff <- simvaceff2                    #efficiency of two dose shanchol for aged 5+
  p.vac.two.eff.1 <- p.vac.two.eff*theta["vac.twodose.efficacy.dt.1"]    #efficiency of two dose shanchol for aged 1-4
  p.vac.two.cov.1 <- theta["vac.twodose.cov.1"]      #coverage of two dose shanchol for age under 5
  p.vac.two.cov <- p.vac.one.cov - simcovdif      #coverage of two dose shanchol for age above 5
  n.vac.two.dur <- theta["vac.twodose.duration"]      #duration of the protection of two dose shanchol 

  
  #probabilities
  #probability from v to s for age group 1
  if (d.twodose.1==1) {p.v.to.s.1 <- 1/n.vac.two.dur} else {p.v.to.s.1 <- 0.286}    #for those who do not have two dose of vaccine, the rate is 0.286 (0.2*0.5/0.7+0.5*0.2/0.7)       
  if (d.twodose.2.3==1) {
    p.v.to.s.2.3 <- (1/n.vac.two.dur)*(p.vac.two.cov/p.vac.one.cov)+(1/n.vac.one.dur)*(1-p.vac.two.cov/p.vac.one.cov)
    p.v.to.i.2 <- ((1-p.vac.one.eff)*(1-p.vac.two.cov/p.vac.one.cov)+(1-p.vac.two.eff)*(p.vac.two.cov/p.vac.one.cov))*p.incidence.2
    p.v.to.i.3 <- ((1-p.vac.one.eff)*(1-p.vac.two.cov/p.vac.one.cov)+(1-p.vac.two.eff)*(p.vac.two.cov/p.vac.one.cov))*p.incidence.3
    } else if (d.onedose.2.3==1) {
      p.v.to.s.2.3 <- 1/n.vac.one.dur
      p.v.to.i.2 <-  (1-p.vac.one.eff)*p.incidence.2
      p.v.to.i.3 <-  (1-p.vac.one.eff)*p.incidence.3
      } else {
        p.v.to.s.2.3 <- 0.35
        p.v.to.i.2 <- p.incidence.2
        p.v.to.i.3 <- p.incidence.3
        }            #probability from v to s for age group 1
  
  #case fatality rate in three age groups
  #p.i.to.d.1 <- theta["p.casefatalityrate.1"]             #case fatality rate for age group 1
  #p.i.to.d.2 <- theta["p.casefatalityrate.2"]             #case fatality rate for age group 2
  #p.i.to.d.3 <- theta["p.casefatalityrate.3"]             #case fatality rate for age group 3
  
  p.i.to.d.1 <- p.i.to.d.2 <- p.i.to.d.3 <- simfatality
  
  
  #phase in and phase out rate
  p.trout.1 <- c(r.mature.1, r.mature.1, r.mature.1, r.mature.1, 0, 0)
  p.trout.2 <- c(r.mature.2, r.mature.2, r.mature.2, r.mature.2, 0, 0)
  #intiation states in each group 
  m.s.1 <- matrix(rep(NA,n.t.end*length(theta_Names)), nrow = n.t.end)           #matrix of health status 
  m.s.2 <- matrix(rep(NA,n.t.end*length(theta_Names)), nrow = n.t.end)           
  m.s.3 <- matrix(rep(NA,n.t.end*length(theta_Names)), nrow = n.t.end)

  m.s.r.1 <- matrix(rep(NA,n.t.end*length(theta_Names)), nrow = n.t.end)         #matrix adding new population and removing matured population
  m.s.r.2 <- matrix(rep(NA,n.t.end*length(theta_Names)), nrow = n.t.end)
  m.s.r.3 <- matrix(rep(NA,n.t.end*length(theta_Names)), nrow = n.t.end) 
  
  m.si.to.v.1 <- matrix(rep(NA, n.t.end*2), nrow = n.t.end)                      #number of susceptible and infections people moving to vaccinated status  
  m.si.to.v.2 <- matrix(rep(NA, n.t.end*2), nrow = n.t.end)
  m.si.to.v.3 <- matrix(rep(NA, n.t.end*2), nrow = n.t.end)
  
  m.si.to.v.1.s <- c(0, 0)                                                       #starting value of moving from susceptible and infection to vaccinate
  m.si.to.v.2.s <- c(0, 0)            
  m.si.to.v.3.s <- c(0, 0)
  
  vacone.1 <- matrix(rep(NA, n.t.end*2), nrow = n.t.end)                         #number of peple receiving one dose of vaccine onely
  vacone.2 <- matrix(rep(NA, n.t.end*2), nrow = n.t.end)
  vacone.3 <- matrix(rep(NA, n.t.end*2), nrow = n.t.end)
  
  vacone.1.s <- c(0, 0)                                                           #starting value of peple receiving one dose of vaccine onely
  vacone.2.s <- c(0, 0)
  vacone.3.s <- c(0, 0)
  
  vactwo.1 <- matrix(rep(NA, n.t.end*2), nrow = n.t.end)                           #number of peple receiving two doses of vaccine onely
  vactwo.2 <- matrix(rep(NA, n.t.end*2), nrow = n.t.end)
  vactwo.3 <- matrix(rep(NA, n.t.end*2), nrow = n.t.end)
  
  vactwo.1.s <- c(0, 0)                                                            #starting value of peple receiving two dose of vaccine onely
  vactwo.2.s <- c(0, 0)
  vactwo.3.s <- c(0, 0)
  

  m.s.1.s <- c(n.pop.1, rep(0,(length(theta_Names)-1)))                                                #starting value of health status
  m.s.2.s <- c(n.pop.2, rep(0,(length(theta_Names)-1)))
  m.s.3.s <- c(n.pop.3, rep(0,(length(theta_Names)-1)))

  m.s.r.1.s <- c(n.pop.1, rep(0,(length(theta_Names)-1)))
  m.s.r.2.s <- c(n.pop.2, rep(0,(length(theta_Names)-1)))
  m.s.r.3.s <- c(n.pop.3, rep(0,(length(theta_Names)-1)))  
  
  
  birthin.1<- matrix(rep(NA,n.t.end*length(theta_Names)), nrow = n.t.end)
  birthin.1.s <- c(rep(0,length(theta_Names)))
  
  trout.1 <- matrix(rep(NA,n.t.end*length(theta_Names)), nrow = n.t.end)
  trout.2 <- matrix(rep(NA,n.t.end*length(theta_Names)), nrow = n.t.end)
  trin.2 <- matrix(rep(NA,n.t.end*length(theta_Names)), nrow = n.t.end)
  trin.3 <- matrix(rep(NA,n.t.end*length(theta_Names)), nrow = n.t.end)
  
  trout.1.s <- c(rep(0,length(theta_Names)))
  trin.2.s <- c(rep(0,length(theta_Names)))
  trout.2.s <- c(rep(0,length(theta_Names)))
  trin.3.s <- c(rep(0,length(theta_Names)))
  
  m.n <-matrix(rep(NA,n.t.end), nrow = n.t.end)
  m.n.s <- n.pop
  
  vac.cov <- matrix(rep(NA,n.t.end), nrow = n.t.end)
  protect<- matrix(rep(NA,n.t.end), nrow = n.t.end)
  m.tp.1 <- array(rep(NA,length(theta_Names)*length(theta_Names)*n.t.end), 
                  dim = c(length(theta_Names), c(length(theta_Names), n.t.end)))
  m.tp.2 <- array(rep(NA,length(theta_Names)*length(theta_Names)*n.t.end), 
                  dim = c(length(theta_Names), c(length(theta_Names), n.t.end)))
  m.tp.3 <- array(rep(NA,length(theta_Names)*length(theta_Names)*n.t.end), 
                  dim = c(length(theta_Names), c(length(theta_Names), n.t.end)))
  #transition probablity that does not change
  p.r.to.s.1 <- p.r.to.s.2  <- p.r.to.s.3 <- r.loss.immunity
  
 for (tt in 1:n.t.end) {

 if ((tt-1) %% n.schedule == 0 & d.twodose.1==1) {p.s.to.v.1 <- p.i.to.v.1 <- p.vac.two.cov.1
 } else {p.s.to.v.1 <- p.i.to.v.1 <- 0}                                     #probability from s to v status
 if ((tt-1) %% n.schedule == 0 & d.onedose.2.3==1) {p.s.to.v.2.3 <- p.i.to.v.2.3 <- p.vac.one.cov
 } else if((tt-1) %% n.schedule ==0 & d.twodose.2.3==1) {p.s.to.v.2.3 <- p.i.to.v.2.3 <- p.vac.one.cov  # the coverage of one dose is used, but the rate of change from v to s is different
 } else {p.s.to.v.2.3 <- p.i.to.v.2.3 <- 0}                                 #probability from s to v status

 #  p.i.to.s.1 <- theta["p.loss.immunity"]*(1-p.s.to.v.1)                    #probability of losing nature immunity , assume no nature death
#   p.i.to.s.2 <- theta["p.loss.immunity"]*(1-p.s.to.v.2.3)                #probability of losing nature immunity 
#   p.i.to.s.3 <- theta["p.loss.immunity"]*(1-p.s.to.v.2.3)                #probability of losing nature immunity    
 #transition probability for group 1
 p.v.to.i.1 <- (1-p.vac.two.eff.1*d.twodose.1)*p.incidence.1
 ifelse(tt==1,v.s.o.1 <-c(rep(NA,length(theta_Names))), 
        v.s.o.1 <- c(1-p.s.to.v.1-p.incidence.1*(1-protect[tt-1]*vac.cov[tt-1]*herd)-p.ndeath.1, p.s.to.v.1, p.incidence.1*(1-protect[tt-1]*vac.cov[tt-1]*herd), 0, p.ndeath.1, 0))
 v.s.o.1.s <- c(1-p.s.to.v.1-p.incidence.1-p.ndeath.1, p.s.to.v.1, p.incidence.1, 0, p.ndeath.1, 0)  
 v.v.o.1 <- c(p.v.to.s.1, 1-p.v.to.s.1-p.v.to.i.1-p.ndeath.1, p.v.to.i.1, 0, p.ndeath.1, 0)
 v.i.o.1 <- c(0, p.i.to.v.1, 0, 1-p.i.to.v.1-(p.ndeath.1+p.i.to.d.1),p.ndeath.1, p.i.to.d.1)
 v.r.o.1 <- c(p.r.to.s.1, 0, 0, 1-p.r.to.s.1-p.ndeath.1, p.ndeath.1, 0 )
 v.nd.o.1 <- c(0,0,0,0,1,0)
 v.cd.o.1 <- c(0,0,0,0,0,1)
 m.tp.1.s <- rbind(v.s.o.1.s, v.v.o.1,v.i.o.1,v.r.o.1, v.nd.o.1, v.cd.o.1)                      #trasition probability for age group 1
 m.tp.1[,,tt-1] <- rbind(v.s.o.1,v.v.o.1,v.i.o.1,v.r.o.1,v.nd.o.1, v.cd.o.1)                      #trasition probability for age group 1
 #transition probability for group 2  
 ifelse(tt==1,v.s.o.2 <-c(rep(NA,length(theta_Names))), 
        v.s.o.2 <- c(1-p.s.to.v.2.3 -p.incidence.2*(1-protect[tt-1]*vac.cov[tt-1]*herd)-p.ndeath.2, p.s.to.v.2.3, p.incidence.2*(1-protect[tt-1]*vac.cov[tt-1]*herd), 0, p.ndeath.2, 0))   
 v.s.o.2.s<-c(1-p.s.to.v.2.3 -p.incidence.2-p.ndeath.2, p.s.to.v.2.3, p.incidence.2, 0, p.ndeath.2, 0)
 v.v.o.2 <- c(p.v.to.s.2.3, 1-p.v.to.s.2.3-p.v.to.i.2-p.ndeath.2, p.v.to.i.2, 0, p.ndeath.2, 0)
 v.i.o.2 <- c(0, p.i.to.v.2.3, 0, 1-p.i.to.v.2.3-(p.ndeath.2+p.i.to.d.2), p.ndeath.2, p.i.to.d.2)
 v.r.o.2 <- c(p.r.to.s.2, 0, 0, 1-p.r.to.s.2-p.ndeath.2, p.ndeath.2, 0 )
 v.nd.o.2 <- c(0,0,0,0,1,0)
 v.cd.o.2 <- c(0,0,0,0,0,1)
 m.tp.2[,,tt-1] <- rbind(v.s.o.2,v.v.o.2,v.i.o.2,v.r.o.2, v.nd.o.2, v.cd.o.2)                      #trasition probability for age group 2
 m.tp.2.s <- rbind(v.s.o.2.s,v.v.o.2,v.i.o.2,v.r.o.2, v.nd.o.2, v.cd.o.2)                  #trasition probability for age group 2,starting points
 
 #transition probability for group 3    
 ifelse(tt==1,v.s.o.3 <-c(rep(NA,length(theta_Names))), 
        v.s.o.3 <- c(1-p.s.to.v.2.3-p.incidence.3*(1-protect[tt-1]*vac.cov[tt-1]*herd)-p.ndeath.3, p.s.to.v.2.3, p.incidence.3*(1-protect[tt-1]*vac.cov[tt-1]*herd), 0, p.ndeath.3,0))   
 v.s.o.3.s <- c(1-p.s.to.v.2.3-p.incidence.3-p.ndeath.3, p.s.to.v.2.3, p.incidence.3, 0, p.ndeath.3,0)   
 v.v.o.3 <- c(p.v.to.s.2.3, 1-p.v.to.s.2.3-p.v.to.i.3-p.ndeath.3, p.v.to.i.3, 0, p.ndeath.3,0)
 v.i.o.3 <- c(0, p.i.to.v.2.3, 0, 1-p.i.to.v.2.3-(p.ndeath.3+p.i.to.d.3), p.ndeath.3, p.i.to.d.3)
 v.r.o.3 <- c(p.r.to.s.3, 0, 0, 1-p.r.to.s.3-p.ndeath.3, p.ndeath.3, 0 )
 v.nd.o.3 <- c(0,0,0,0,1,0) 
 v.cd.o.3 <- c(0,0,0,0,0,1)
 m.tp.3[,,tt-1] <- rbind(v.s.o.3, v.v.o.3, v.i.o.3, v.r.o.3, v.nd.o.3, v.cd.o.3)                   #trasition probability for age group 3
 m.tp.3.s <- rbind(v.s.o.3.s, v.v.o.3, v.i.o.3, v.r.o.3, v.nd.o.3, v.cd.o.3)                   #trasition probability for age group 3
 
 
 ifelse(tt==1, m.s.1[tt,] <- m.s.r.1.s%*%m.tp.1.s, m.s.1[tt,] <- m.s.r.1[tt-1,]%*%m.tp.1[,,tt-1])   #calculation transition
 ifelse(tt==1, m.s.2[tt,] <- m.s.r.2.s%*%m.tp.2.s, m.s.2[tt,] <- m.s.r.2[tt-1,]%*%m.tp.2[,,tt-1])   #calculation transition
 ifelse(tt==1, m.s.3[tt,] <- m.s.r.3.s%*%m.tp.3.s, m.s.3[tt,] <- m.s.r.3[tt-1,]%*%m.tp.3[,,tt-1])   #calculation transition
 
#  m.s.1[tt,] <- m.s.r.1[tt-1,]%*%m.tp.1
#  m.s.2[tt,] <- m.s.r.2[tt-1,]%*%m.tp.2
#  m.s.3[tt,] <- m.s.r.3[tt-1,]%*%m.tp.3
 
  
  m.n[tt,] <- sum(m.s.1[tt,c(1,2,3,4)],m.s.2[tt,c(1,2,3,4)],m.s.3[tt,c(1,2,3,4)])   #alive population
  vac.cov[tt] <- sum(m.s.1[tt,2],m.s.2[tt,2],m.s.3[tt,2])/m.n[tt,]
  protect[tt] <- exp(model1[["coefficients"]][2]*vac.cov[tt]+model1[["coefficients"]][1])/
    (1+exp(model1[["coefficients"]][2]*vac.cov[tt]+model1[["coefficients"]][1]))
  
  if (tt==1) {
    birthin.1[1, ] <- c(m.n.s,rep(0,(length(theta_Names)-1)))*r.birth               #number of births entering into the cohort
    trout.1[1,] <- m.s.1.s*p.trout.1                #matured pop for children under 5
    trout.2[1,] <- m.s.2.s*p.trout.2                #matured pop from those 5-15
    trin.2[1,]<- 0                                                                  #number of pop entering to pop 5-15
    trin.3[1,]<- 0  
  } else {
   birthin.1[tt,]<- c(m.n[tt-1,], rep(0,(length(theta_Names)-1)))*r.birth 
   trout.1[tt,] <- m.s.1[tt-1,]*p.trout.1
   trout.2[tt,]<- m.s.2[tt-1,]*p.trout.2 
   trin.2[tt,]<- trout.1[tt-1,]                                                  #number of pop entering to pop 5-15
   trin.3[tt,]<- trout.2[tt-1,]                                                  #number of pop entering into pop 15+
   }
  
  ifelse(tt==1, m.si.to.v.1[tt,] <- m.s.r.1.s[c(1,3)]*rbind(p.s.to.v.1, p.i.to.v.1),
         m.si.to.v.1[tt,] <- m.s.r.1[tt-1,c(1,3)]*rbind(p.s.to.v.1, p.i.to.v.1))
  ifelse(tt==1, m.si.to.v.2[tt,] <- m.s.r.2.s[c(1,3)]*rbind(p.s.to.v.2.3, p.i.to.v.2.3),
         m.si.to.v.2[tt,] <- m.s.r.2[tt-1,c(1,3)]*rbind(p.s.to.v.2.3, p.i.to.v.2.3)) 
  ifelse(tt==1, m.si.to.v.3[tt,] <- m.s.r.3.s[c(1,3)]*rbind(p.s.to.v.2.3, p.i.to.v.2.3),
         m.si.to.v.3[tt,] <- m.s.r.3[tt-1,c(1,3)]*rbind(p.s.to.v.2.3, p.i.to.v.2.3))   
  
#  m.si.to.v.1[tt,] <- m.s.r.1[tt-1,c(1,3)]* rbind(p.s.to.v.1, p.i.to.v.1)
#  m.si.to.v.2[tt,] <- m.s.r.2[tt-1,c(1,3)]* rbind(p.s.to.v.2.3, p.i.to.v.2.3)
#  m.si.to.v.3[tt,] <- m.s.r.3[tt-1,c(1,3)]* rbind(p.s.to.v.2.3, p.i.to.v.2.3)
  
  if (d.twodose.1==1) {
    vacone.1[tt,] <- c(0, 0)
    vactwo.1[tt,] <- m.si.to.v.1[tt,]
  } else {vacone.1[tt,] <- vactwo.1[tt,] <- c(0,0)}   
  
  if (d.onedose.2.3==1 & d.twodose.2.3==0) {
    vacone.2[tt,] <- m.si.to.v.2[tt,]
    vactwo.2[tt,] <- c(0,0)
    vacone.3[tt,] <- m.si.to.v.3[tt,]
    vactwo.3[tt,] <- c(0,0)
  } else if(d.onedose.2.3==0 & d.twodose.2.3==1) {
    vacone.2[tt,] <- m.si.to.v.2[tt,]*(1-p.vac.two.cov/p.vac.one.cov)
    vactwo.2[tt,] <- m.si.to.v.2[tt,]*p.vac.two.cov/p.vac.one.cov
    vacone.3[tt,] <- m.si.to.v.3[tt,]*(1-p.vac.two.cov/p.vac.one.cov)
    vactwo.3[tt,] <- m.si.to.v.3[tt,]*p.vac.two.cov/p.vac.one.cov
  } else {vacone.2[tt,] <- vactwo.2[tt,] <- vacone.3[tt,] <- vactwo.3[tt,] <- 0}                
  
  m.s.r.1[tt,]<- m.s.1[tt,] + birthin.1[tt,]-trout.1[tt,]
  m.s.r.2[tt,]<- m.s.2[tt,] + trin.2[tt,]-trout.2[tt,]
  m.s.r.3[tt,]<- m.s.3[tt,] + trin.3[tt,]
 }
  dta1<- m.s.1
  dta2<- m.s.2
  dta3<- m.s.3
  dta4<-m.n
  dta5<-cbind(vacone.1, vactwo.1, 
              vacone.2, vactwo.2,
              vacone.3, vactwo.3)                             #estimation number of people in each health status
  v.vac.coverage <- (dta1[,2]+dta2[,2]+dta3[,2])/dta4          #vaccine coverage
  n.inf <- sum(dta1[,3])+sum(dta2[,3])+sum(dta3[,3])          #number of infections
  n.ndeath <- dta1[n.t.end,5]+dta2[n.t.end,5]+dta3[n.t.end,5] #number of nature death
  n.cdeath <- dta1[n.t.end,6]+dta2[n.t.end,6]+dta3[n.t.end,6] #number of death due to cholera
  n.vac.one <- sum(dta5[,c(1,2, 5,6, 9,10)])
  n.vac.two <- sum(dta5[,c(3,4,7,8,11, 12)])
  
  v.process.ind <- c(n.inf, n.ndeath, n.cdeath, n.vac.one, n.vac.two)  #combine data on number of nature death, number of vaccines
  
  #estimation vaccine cost
  c.onedose.t.hs <- (theta["c.onedose"]+theta["c.delivery"])*(1+theta["wastage"])+ theta["c.travel"]
  c.twodose.t.hs <- c.onedose.t.hs*2
  
  c.onedose.t.so <- c.onedose.t.hs 
  c.twodose.t.so <- c.onedose.t.so*2
  
  # simulation direct treatment cost
  c.treat.1.t.hs <- simctreat1
  c.treat.2.t.hs <- simctreat2
  c.treat.3.t.hs <- simctreat3
  #import indirect cost 
  c.tml <- simctml
  c.ploss <- simcploss
  
  c.treat.1.t.so <- c.treat.1.t.hs + c.tml + c.ploss   #calculation cost including productivity loss and travel
  c.treat.2.t.so <- c.treat.2.t.hs + c.tml + c.ploss
  c.treat.3.t.so <- c.treat.2.t.hs + c.tml + c.ploss
  
  dc<- theta["dc"] #discount rate for cost
  de<- theta["de"] #discount rate for effectiveness
  
  le.1 <- theta["lifeexp.1"]
  le.2 <- theta["lifeexp.2"]
  le.3 <- theta["lifeexp.3"]
  
  v.daly.1 <- c(0,0,simdw*simtcholera/365, 0, 1/de*(1-exp(-de*le.1)), 1/de*(1-exp(-de*le.1)))
  v.daly.2 <- c(0,0,simdw*simtcholera/365, 0, 1/de*(1-exp(-de*le.2)), 1/de*(1-exp(-de*le.2)))
  v.daly.3 <- c(0,0,simdw*simtcholera/365, 0, 1/de*(1-exp(-de*le.3)), 1/de*(1-exp(-de*le.3)))
  
  tempc.vac.hs <- c(c.onedose.t.hs, c.onedose.t.hs, 
                    c.twodose.t.hs, c.twodose.t.hs, 
                    c.onedose.t.hs, c.onedose.t.hs, 
                    c.twodose.t.hs, c.twodose.t.hs,
                    c.onedose.t.hs, c.onedose.t.hs, 
                    c.twodose.t.hs, c.twodose.t.hs)
  
  tempc.vac.so <- c(c.onedose.t.so, c.onedose.t.so, 
                    c.twodose.t.so, c.twodose.t.so, 
                    c.onedose.t.so, c.onedose.t.so, 
                    c.twodose.t.so, c.twodose.t.so,
                    c.onedose.t.so, c.onedose.t.so, 
                    c.twodose.t.so, c.twodose.t.so)
  
  v.c.vac.hs <- matrix(tempc.vac.hs, ncol = 1)
  v.c.vac.so <- matrix(tempc.vac.so, ncol = 1)
  
  v.cost.1.hs <- matrix(c(0, 0, c.treat.1.t.hs, 0, 0, 0), nrow=length(theta_Names))
  v.cost.2.hs <- matrix(c(0, 0, c.treat.2.t.hs, 0, 0, 0), nrow=length(theta_Names))
  v.cost.3.hs <- matrix(c(0, 0, c.treat.3.t.hs, 0, 0, 0), nrow=length(theta_Names))
  
  v.cost.1.so <- matrix(c(0, 0, c.treat.1.t.so, 0, 0, 0),nrow=length(theta_Names))
  v.cost.2.so <- matrix(c(0, 0, c.treat.2.t.so, 0, 0, 0), nrow=length(theta_Names))
  v.cost.3.so <- matrix(c(0, 0, c.treat.3.t.so, 0, 0, 0), nrow=length(theta_Names))
  
  df.c <- matrix(1/(1+theta["dc"])^(0:(n.t.end-1)), nrow = n.t.end)
  df.d <- matrix(1/(1+theta["de"])^(0:(n.t.end-1)), nrow = n.t.end)
  
  c.1.hs <- dta1%*%v.cost.1.hs
  c.2.hs <- dta2%*%v.cost.2.hs
  c.3.hs <- dta3%*%v.cost.3.hs
  c.tr.hs <- c.1.hs + c.2.hs +c.3.hs
  
  c.1.so <- dta1%*%v.cost.1.so
  c.2.so <- dta2%*%v.cost.2.so
  c.3.so <- dta3%*%v.cost.3.so
  c.tr.so <- c.1.so + c.2.so +c.3.so
  
  c.vac.hs <- dta5%*%v.c.vac.hs
  c.vac.so <- dta5%*%v.c.vac.so
  
  c.tot.hs <- c.tr.hs + c.vac.hs
  c.tot.so <- c.tr.so + c.vac.so
  
  c.dis.tr.hs <-t(c.tr.hs)%*%df.c
  c.dis.vac.hs <-t(c.vac.hs)%*%df.c
  c.dis.hs <- t(c.tot.hs)%*%df.c
  
  c.dis.tr.so <- t(c.tr.so)%*%df.c
  c.dis.vac.so <- t(c.vac.so)%*%df.c
  c.dis.so <- t(c.tot.so)%*%df.c                                                          #social cost without indirect cost from death
  c.dis.death.so.1 <- dta1[n.t.end,c(5,6)]*(1855.74/dc*(1-exp(-dc*le.1)))                  #GDP.percapita 1855.74 in 2019 as proxy (WB)
  c.dis.death.so.2 <- dta2[n.t.end,c(5,6)]*(1855.74/dc*(1-exp(-dc*le.1))) 
  c.dis.death.so.3 <- dta3[n.t.end,c(5,6)]*(1855.74/dc*(1-exp(-dc*le.1))) 
  c.dis.death.so <- sum(c(c.dis.death.so.1,c.dis.death.so.2,c.dis.death.so.3))  
  c.dis.tot.so <- c.dis.so + c.dis.death.so  
  
  d.d.1 <- dta1[,c(1,2,3,4)]%*%matrix(v.daly.1[c(1,2,3,4)], nrow = 4)               #DALYs for disability 
  d.d.2 <- dta2[,c(1,2,3,4)]%*%matrix(v.daly.2[c(1,2,3,4)], nrow = 4)
  d.d.3 <- dta3[,c(1,2,3,4)]%*%matrix(v.daly.3[c(1,2,3,4)], nrow = 4)
  d.d.total <- d.d.1+d.d.2+d.d.3
  d.d.dis <- t(d.d.total)%*%df.d
  
  d.m.1 <- dta1[n.t.end,c(5,6)]%*%v.daly.1[c(5,6)]                  #DALYs from mortality 
  d.m.2 <- dta2[n.t.end,c(5,6)]%*%v.daly.2[c(5,6)]
  d.m.3 <- dta3[n.t.end,c(5,6)]%*%v.daly.3[c(5,6)]
  d.m.dis <- d.m.1+d.m.2+d.m.3
  
  d.dis <- d.d.dis+d.m.dis
  
#  data.all.age <-cbind(dta1,dta2,dta3,dta4,dta5)
  
  res<- list(datale5=dta1, data5to15=dta2, datagt15=dta3, datan = dta4, 
             dvac = dta5, vac.cov = vac.cov, v.vac.coverage=v.vac.coverage,
             protect = protect, v.process.ind=v.process.ind,
             d.daly=d.d.dis, m.daly=d.m.dis, daly=d.dis, cost.tr.hs = c.dis.tr.hs, 
             cost.vac.hs = c.dis.vac.hs, cost.hs =c.dis.hs, cost.tr.so = c.dis.tr.so, 
             cost.vac.so = c.dis.vac.so, cost.so =c.dis.so, cost.death = c.dis.death.so,
             cost.tot.so = c.dis.tot.so)
  return(res)
}

