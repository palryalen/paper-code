library(data.table)
library(timereg)
# install.packages("ahw")
library(ahw)

load("sim_data.RData")



tms = seq(0,5,length.out = 1e3)


makeAnalysis = function(dta,BootWeights=rep(1,nrow(dta)),willPlot=F){
  
  dta$btWeights = BootWeights
  # Fitting aalen regressions for time-to-follow-up for the factual (RNA) and the hypothetical (DNA) scenarios
  faFit <- aalen(Surv(from,to,to.state == "follow-up")~1,data=dta[dta$technology == "RNA" & dta$from.state == "none",], 
                 weights = dta[dta$technology == "RNA" & dta$from.state == "none",]$btWeights )
  cfaFit <- aalen(Surv(from,to,to.state == "follow-up")~1,data=dta[dta$technology == "DNA" & dta$from.state == "none",],
                  weights = dta[dta$technology == "DNA" & dta$from.state == "none",]$btWeights)
  
  
  # Calculating weights for the RNA/Proofer group and plotting the weight trajectories. See documentation using example(makeContWeights). 
  frame <- makeContWeights(faFit,cfaFit,dta[dta$technology == "RNA",],"none","follow-up","from","to",
                           "from.state","to.state","id",0.1082,willPlotWeights = willPlot)
  
  
  # Calculating Nelson-Aalen estimates for the CIN2+ end point
  AA1 <- aalen(Surv(from,to, to.state == "cin2+")~1,data=dta[dta$technology == "RNA",], weights = dta[dta$technology == "RNA",]$btWeights)
  AA2 <- aalen(Surv(from,to, to.state == "cin2+")~1,data=dta[dta$technology == "DNA",], weights = dta[dta$technology == "DNA",]$btWeights)
  
  
  # Calculating the Kaplan-Meier estimators for CIN2+ end point
  ss1 <- cumprod(1 - diff(c(0,AA1$cum[,2])))
  ss2 <- cumprod(1 - diff(c(0,AA2$cum[,2])))
  
  
  # Weighted Nelson-Aalen estimates for the CIN2+ end point, Proofer group
  AA1_cf <- aalen(Surv(from,to, to.state == "cin2+")~1,data=frame[technology == "RNA",],weights = frame[technology == "RNA"]$weights * frame[technology == "RNA"]$btWeights)
  
  # The weighted Kaplan-Meier estimator for the Proofer group
  ss1_cf <- cumprod(1 - diff(c(0,AA1_cf$cum[,2])))
  
  
  tInd1 = sapply(tms,function(tmm)max(which(AA1$cum[,1] <= tmm)))
  tInd2 = sapply(tms,function(tmm)max(which(AA2$cum[,1] <= tmm)))
  tInd1cf = sapply(tms,function(tmm)max(which(AA1_cf$cum[,1] <= tmm)))
  
  
  # list(ss1 = ss1[tInd1], ss2 = ss2[tInd2], ss1_cf = ss1_cf[tInd1cf])
  
  # Returns marginal analyses and the difference
  c(ss1[tInd1],ss2[tInd2],ss1_cf[tInd1cf],ss2[tInd2] - ss1_cf[tInd1cf])
}


analysis = makeAnalysis(sim_data)

tLen = length(tms)

plot(tms,1-analysis[1:tLen],type="s",xaxs="i",yaxs="i",xlim=c(0,4),ylim=c(0,0.085),xlab="Years",ylab="",main="Proportion CIN2+ detected")
lines(tms,1-analysis[(tLen+1):(2*tLen)],type="s",col="gray")
lines(tms,1-analysis[(2*tLen+1):(3*tLen)],type="s",lty=2)
legend("topleft",c("Proofer-group actual follow-up","Proofer-group follow-up as DNA","DNA-group actual follow-up"),
       lty=c(1,2,1),col=c(1,1,"gray"),bty="n")


plot(tms, analysis[(3*tLen+1):(4*tLen)],main="Difference with subsequent testing regime as DNA group",type='s',ylim=c(-0.005,0.03))




getBootVec = function(dta){
  nIds = length(unique(dta$id))
  idReps = as.numeric(table(dta$id))
  
  btSamp = sample(1:nIds,replace=T)
  btVec = rep(0,nIds)
  tab = table(btSamp)
  btVec[as.numeric(names(tab))] = as.numeric(tab)
  btVec = rep(btVec,times=idReps)
  
  btVec
}



# Bootstrap 

nBoot = 20


BootMat = matrix(0,nrow=nBoot, ncol=4*tLen)
for(bb in 1:nBoot){
  
  btVec = getBootVec(sim_data)
  BootMat[bb,] = makeAnalysis(sim_data,btVec)
}


boot_var = apply(BootMat,2,var)




# Plot with bootstrapped confidence intervals

plot(tms,1-analysis[1:tLen],type="s",xaxs="i",yaxs="i",xlim=c(0,4),ylim=c(0,0.085),xlab="Years",ylab="",main="Proportion CIN2+ detected")
lines(tms,1-analysis[1:tLen] + 1.96 * sqrt(boot_var[1:tLen]),type="s",lty=2)
lines(tms,1-analysis[1:tLen] - 1.96 * sqrt(boot_var[1:tLen]),type="s",lty=2)
lines(tms,1-analysis[(tLen+1):(2*tLen)],type="s",col="gray")
lines(tms,1-analysis[(tLen+1):(2*tLen)] + 1.96 * sqrt(boot_var[(tLen+1):(2*tLen)]),type="s",col="gray",lty=2)
lines(tms,1-analysis[(tLen+1):(2*tLen)] - 1.96 * sqrt(boot_var[(tLen+1):(2*tLen)]),type="s",col="gray",lty=2)
lines(tms,1-analysis[(2*tLen+1):(3*tLen)],type="s",lty=2)
lines(tms,1-analysis[(2*tLen+1):(3*tLen)] + 1.96 * sqrt(boot_var[(2*tLen+1):(3*tLen)]),type="s",lty=2)
lines(tms,1-analysis[(2*tLen+1):(3*tLen)] - 1.96 * sqrt(boot_var[(2*tLen+1):(3*tLen)]),type="s",lty=2)
legend("topleft",c("Proofer-group actual follow-up","Proofer-group follow-up as DNA","DNA-group actual follow-up"),
       lty=c(1,2,1),col=c(1,1,"gray"),bty="n")


plot(tms, analysis[(3*tLen+1):(4*tLen)],main="Difference with subsequent testing regime as DNA group",type='s',ylim=c(-0.005,0.035))
lines(tms, analysis[(3*tLen+1):(4*tLen)] + 1.96 * sqrt(boot_var[(3*tLen+1):(4*tLen)]) ,lty=2)
lines(tms, analysis[(3*tLen+1):(4*tLen)] - 1.96 * sqrt(boot_var[(3*tLen+1):(4*tLen)]) ,lty=2)



