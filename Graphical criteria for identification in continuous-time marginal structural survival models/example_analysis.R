library(data.table)
library(timereg)
devtools::install_github("palryalen/ahw")
library(ahw)

load("sim_data.RData")

# Fitting aalen regressions for time-to-follow-up for the factual (RNA) and the hypothetical (DNA) scenarios
faFit <- aalen(Surv(from,to,to.state == "follow-up")~1,data=sim_data[sim_data$technology == "RNA" & sim_data$from.state == "none",])
cfaFit <- aalen(Surv(from,to,to.state == "follow-up")~1,data=sim_data[sim_data$technology == "DNA" & sim_data$from.state == "none",])


# Calculating weights for the RNA/Proofer group and plotting the weight trajectories. See documentation using example(makeContWeights). 
frame <- makeContWeights(faFit,cfaFit,sim_data[sim_data$technology == "RNA",],"none","follow-up","from","to",
                         "from.state","to.state","id",0.1082,willPlotWeights = T)


# Calculating Nelson-Aalen estimates for the CIN2+ end point
AA1 <- aalen(Surv(from,to, to.state == "cin2+")~1,data=sim_data[sim_data$technology == "RNA",])
AA2 <- aalen(Surv(from,to, to.state == "cin2+")~1,data=sim_data[sim_data$technology == "DNA",])


# Calculating the Kaplan-Meier estimators for CIN2+ end point
ss1 <- cumprod(1 - diff(c(0,AA1$cum[,2])))
ss2 <- cumprod(1 - diff(c(0,AA2$cum[,2])))


# Weighted Nelson-Aalen estimates for the CIN2+ end point, Proofer group
AA1_cf <- aalen(Surv(from,to, to.state == "cin2+")~1,data=frame[technology == "RNA",],weights = frame[technology == "RNA"]$weights)

# The weighted Kaplan-Meier estimator for the Proofer group
ss1_cf <- cumprod(1 - diff(c(0,AA1_cf$cum[,2])))


# plotting the results
plot(AA1$cum[,1],1-ss1,type="s",xaxs="i",yaxs="i",xlim=c(0,4),ylim=c(0,0.085),xlab="Years",ylab="",main="Proportion CIN2+ detected")
lines(AA2$cum[,1],1-ss2,type="s",col="gray")
lines(AA1_cf$cum[,1],1-ss1_cf,type="s",lty=2)
legend("topleft",c("Proofer-group actual follow-up","Proofer-group follow-up as DNA","DNA-group actual follow-up"),
       lty=c(1,2,1),col=c(1,1,"gray"),bty="n")
