#sensitivity analysis 
sensitivity_model <- function(){
  i=0
  sendatastore<- matrix(rep(NA, 5*26*10), nrow = 5*26)    #5:five schedules, 26:26 projection periods; 10:10 outputs
  for (aa in 5:30) {
    for (bb in 1:4) {
      i=i+1
      temp<- stochastic_model(n.trials=1000, n.t.end=aa, theta, n.schedule=bb)
      rec<- c(aa, bb, mean(temp[["simdata.outcome"]]$icer.one.hs), mean(temp[["simdata.outcome"]]$icer.one.so),
              mean(temp[["simdata.outcome"]]$icer.two.hs), mean(temp[["simdata.outcome"]]$icer.two.so), mean(temp[["simdata.outcome"]]$icer.combine.hs),
              mean(temp[["simdata.outcome"]]$icer.combine.so),mean(temp[["simdata.outcome"]]$icer.two1.hs),mean(temp[["simdata.outcome"]]$icer.two1.so))
      sendatastore[i,]<-rec
    }
  }
  sentdata <- data.frame(sendatastore)
  colnames(sentdata) = c("projectionperiod","schedule","icer.one.hs", "icer.one.so",
                         "icer.two.hs", "icer.two.so", "icer.combine.hs", "icer.combine.so",
                         "icer.two1.hs","icer.two1.so")
  senres<-list(sentdata=sentdata)
  return(senres)
}