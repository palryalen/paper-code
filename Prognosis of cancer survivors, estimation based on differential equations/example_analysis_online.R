

# Install and load the transform.hazards package
install_github("palryalen/transform.hazards")
library(transform.hazards)

library(survival)
library(data.table)
library(relsurv)
source("misc_online.R")


# loading sumilated example cancer registry data
load("ex_data.RData")


# loading ficticious population table data
load("pop_table.RData")



# obtaining ratetable
ratTab = transrate(men=matrix(1-pop_table$Male, length(unique(pop_table$Age)), length(unique(pop_table$Year))),
                   women=matrix(1-pop_table$Female, length(unique(pop_table$Age)), length(unique(pop_table$Year))),
                   c(1930,2014),
                   1)


lifetable_step = 0.025
t0 = sort(seq(0,100,lifetable_step))



# estimating relative survival ratio using ederer I estimator
relsurv_fit = rs.surv(Surv(time,cens)~1,rmap=list(age=age*365.241),ratetable = ratTab,data=example_data,
                        method = "ederer1")



cancer_fit = survfit(Surv(time,cens)~1,data=example_data)

# obtaining the population hazard associated with the relative survival ratio along the time scale t0
d_H_e = get_pophaz_increments(t0,cancer_fit,relsurv_fit)

# plot(relsurv_fit,ylim=c(0,1))
# lines(survfit(Surv(time,cens)~1,data=example_data))
# lines(t0*365,cumprod(1-d_H_e),type="s")


# defining Surv objects for marginal (all-cause) and cause-specific parameters
t_y_strata = Surv(time=example_data$time/365.241,event= example_data$cause != "",type='right')
t_y_cancer_spec = Surv(time=example_data$time/365.241,event= example_data$cause == "cancer",type='right')
t_y_other_spec = Surv(time=example_data$time/365.241,event= example_data$cause == "other",type='right')







Delta = 5
restrict_time <- 18



parameters = get_Delta_parameters(Delta,restrict_time,t0,t_y_strata,t_y_cancer_spec,t_y_other_spec,d_H_e)



par(mfrow=c(4,3))

plot(parameters$t_tot_large,parameters$CS_est,type="s",main="5-yr condsurv",ylab="",xlab="t")
lines(parameters$t_tot_large,parameters$CS_est + 1.96 * sqrt (parameters$CS_var),type="s",lty=2)
lines(parameters$t_tot_large,parameters$CS_est - 1.96 * sqrt (parameters$CS_var),type="s",lty=2)

plot(parameters$t_tot_large,parameters$CS_est/parameters$CS_p,type="s",main="5-yr condsurv ratio",ylab="",xlab="t")
lines(parameters$t_tot_large,parameters$CS_est/parameters$CS_p + 1.96 * sqrt (parameters$CS_var)/parameters$CS_p,type="s",lty=2)
lines(parameters$t_tot_large,parameters$CS_est/parameters$CS_p - 1.96 * sqrt (parameters$CS_var)/parameters$CS_p,type="s",lty=2)

plot(parameters$t_tot_large,parameters$CS_est-parameters$CS_p,type="s",main="5-yr condsurv diff",ylab="",xlab="t")
lines(parameters$t_tot_large,parameters$CS_est-parameters$CS_p + 1.96 * sqrt (parameters$CS_var),type="s",lty=2)
lines(parameters$t_tot_large,parameters$CS_est-parameters$CS_p - 1.96 * sqrt (parameters$CS_var),type="s",lty=2)


plot(parameters$t_tot_large,parameters$RMRL_est,type="s",main="5-yr RMRL",ylab="",xlab="t")
lines(parameters$t_tot_large,parameters$RMRL_est + 1.96 * sqrt (parameters$RMRL_var),type="s",lty=2)
lines(parameters$t_tot_large,parameters$RMRL_est - 1.96 * sqrt (parameters$RMRL_var),type="s",lty=2)

plot(parameters$t_tot_large,parameters$RMRL_est/parameters$RMRL_p,type="s",main="5-yr RMRL ratio",ylab="",xlab="t")
lines(parameters$t_tot_large,parameters$RMRL_est/parameters$RMRL_p + 1.96 * sqrt (parameters$RMRL_var)/parameters$RMRL_p,type="s",lty=2)
lines(parameters$t_tot_large,parameters$RMRL_est/parameters$RMRL_p - 1.96 * sqrt (parameters$RMRL_var)/parameters$RMRL_p,type="s",lty=2)

plot(parameters$t_tot_large,parameters$RMRL_est-parameters$RMRL_p,type="s",main="5-yr RMRL diff",ylab="",xlab="t")
lines(parameters$t_tot_large,parameters$RMRL_est-parameters$RMRL_p + 1.96 * sqrt (parameters$RMRL_var),type="s",lty=2)
lines(parameters$t_tot_large,parameters$RMRL_est-parameters$RMRL_p - 1.96 * sqrt (parameters$RMRL_var),type="s",lty=2)



plot(parameters$t_tot_large,parameters$c_spec_cond_risk,type="s",main="5-yr cond. c.spec. risk",ylab="",xlab="t")
lines(parameters$t_tot_large,parameters$c_spec_cond_risk + 1.96 * sqrt (parameters$c_spec_cond_risk_var),type="s",lty=2)
lines(parameters$t_tot_large,parameters$c_spec_cond_risk - 1.96 * sqrt (parameters$c_spec_cond_risk_var),type="s",lty=2)

plot(parameters$t_tot_large,parameters$c_spec_cond_riskRR,type="s",main="5-yr cond. rel. risk",ylab="",xlab="t")
lines(parameters$t_tot_large,parameters$c_spec_cond_riskRR + 1.96 * sqrt (parameters$c_spec_cond_riskRR_var),type="s",lty=2)
lines(parameters$t_tot_large,parameters$c_spec_cond_riskRR - 1.96 * sqrt (parameters$c_spec_cond_riskRR_var),type="s",lty=2)

plot(parameters$t_tot_large,parameters$c_spec_cond_riskRD,type="s",main="5-yr cond. risk diff",ylab="",xlab="t")
lines(parameters$t_tot_large,parameters$c_spec_cond_riskRD + 1.96 * sqrt (parameters$c_spec_cond_riskRD_var),type="s",lty=2)
lines(parameters$t_tot_large,parameters$c_spec_cond_riskRD - 1.96 * sqrt (parameters$c_spec_cond_riskRD_var),type="s",lty=2)


plot(parameters$t_tot_large,parameters$c_spec_cond_risk_proportion,type="s",main="5-yr cond. proportion",ylab="",xlab="t")
lines(parameters$t_tot_large,parameters$c_spec_cond_risk_proportion + 1.96 * sqrt (parameters$c_spec_cond_risk_proportion_var),type="s",lty=2)
lines(parameters$t_tot_large,parameters$c_spec_cond_risk_proportion - 1.96 * sqrt (parameters$c_spec_cond_risk_proportion_var),type="s",lty=2)












