
H_RMRL <- function(X,Y)(Y[2]-X[2])/X[1]
nabla_H_RMRL <- function(X,Y)matrix(c((X[2]-Y[2])/X[1]^2,-1/X[1],0,1/X[1]),nrow=1)

H_condsurv <- function(X,Y)Y[1]/X[1]
nabla_H_condsurv <- function(X,Y)matrix(c(-Y[1]/X[1]^2,1/X[1]),nrow=1)

H_cuminc <- function(X,Y)(Y[2]-X[2])/X[1]
nabla_H_cuminc <- function(X,Y)matrix(c((X[2]-Y[2])/X[1]^2,-1/X[1],0,0,1/X[1],0),nrow=1)

H_cuminc_RR <- function(X,Y)(Y[2]-X[2])/(Y[3]-X[3])
nabla_H_cuminc_RR <- function(X,Y)matrix(c( 0, -1/(Y[3]-X[3]),(Y[2]-X[2])/(Y[3]-X[3])^2,
                                            0,1/(Y[3]-X[3]),-(Y[2]-X[2])/(Y[3]-X[3])^2),nrow=1)

H_cuminc_RD <- function(X,Y)(Y[2]-X[2])/X[1] - (Y[3]-X[3])/X[1]
nabla_H_cuminc_RD <- function(X,Y)matrix(c( -(Y[2]-X[2])/X[1]^2 + (Y[3]-X[3])/X[1]^2, -1/X[1], 1/X[1],
                                            0, 1/X[1], -1/X[1]),nrow=1)


H_cuminc_proportion <- function(X,Y)(Y[2]-X[2])/(X[1]-Y[1])
nabla_H_cuminc_proportion <- function(X,Y)1/(X[1]-Y[1]) * matrix(c(-(Y[2]-X[2])/(X[1]-Y[1]),-1,0,
                                                                   (Y[2]-X[2])/(X[1]-Y[1]),1,0),nrow=1)






get_pophaz_increments = function(t0,cancer_fit,relsurv_fit){
  
  
  tmmssss = sort(unique(c(0,cancer_fit$time,relsurv_fit$time)))
  mt1=match(cancer_fit$time,tmmssss)
  mt2=match(relsurv_fit$time,tmmssss)
  
  dA_ederer1 = dA_canc = rep(0,length(tmmssss))
  dA_ederer1[mt2] = diff(c(0,-log(relsurv_fit$surv)))

  dA_canc[mt1] =diff(c(0,-log(cancer_fit$surv)))
  dA_pop = dA_canc - dA_ederer1

  
  tmmssss = tmmssss/365.241

  A_pop = cumsum(dA_pop)
  t0_iind = sapply(1:length(t0),function(i)max(which(tmmssss <= t0[i])))
  # t0_iind[t0_iind==-Inf]=1
  
  # d_H_e are the population hazard increments along the time scale t0
  d_H_e = A_pop[t0_iind]; d_H_e = diff(c(0,d_H_e))
  
  d_H_e
}





getCrossCov_Delta <- function(tms,t1_ind,Delta,X_plugin,nabla_H,H,JacobianList,hazMat){
  
 
  
  t2_ind <- max(which( tms[t1_ind] + Delta >= tms ) )
  
  
  X1 <- X_plugin$X[,t1_ind];V1 <- as.matrix(X_plugin$covariance[,,t1_ind])
  X2 <- X_plugin$X[,t2_ind];V2 <- as.matrix(X_plugin$covariance[,,t2_ind])
  
  param_t <- H(X1,X2)
  
  numJump <- t2_ind-t1_ind
  k <- ifelse(is.null(nrow(hazMat)),1,nrow(hazMat))
  if(is.na(t2_ind))
    return(NA)
  V12 <- V1
  if(numJump > 0){
    for(zz in 1:numJump){
      V_cross_temp <- matrix(0,ncol=ncol(V1),nrow=nrow(V1))
      X_last <- X_plugin$X[,zz+t1_ind-1]
      # V_last <- X_plugin$covariance[,,zz+t1_ind-1]
      dA_haz <- hazMat[,zz+t1_ind]
      for(j in 1:k){
        V_cross_temp <- V_cross_temp + JacobianList[[j]](X_last) %*% V12 * dA_haz[j]
      }
      V12 <- V12 + V_cross_temp
    }
  }
  
  BlockMatrix <- matrix(0,nrow=2*nrow(V1),ncol=2*ncol(V1))
  BlockMatrix[1:nrow(V1),1:ncol(V1)] <- V1
  BlockMatrix[((nrow(V1)+1):(2*nrow(V1))),1:ncol(V1)] <- V12
  BlockMatrix[1:nrow(V1),((ncol(V1)+1):(2*ncol(V1)))] <- t(V12)
  BlockMatrix[((nrow(V1)+1):(2*nrow(V1))),((ncol(V1)+1):(2*ncol(V1)))] <- V2
  
  # Original estimator:
  # crossCov <- nabla_H(X1,X2) %*% BlockMatrix %*% t(nabla_H(X1,X2))
  
  # Positively definite estimator:
  ee = eigen(BlockMatrix)
  SPD_matrix = ee$vectors %*% abs(diag(ee$values)) %*% solve(ee$vectors)
  crossCov <- nabla_H(X1,X2) %*% SPD_matrix %*% t(nabla_H(X1,X2))
  
  
  c(crossCov,param_t)
}






get_Delta_parameters = function(Delta,restrict_time,t0,t_y_strata,t_y_cancer_spec,t_y_other_spec,d_H_e){
  
  
  ss = survfit(t_y_strata ~ 1)
  
  t_tot <- sort(unique(c(t_y_strata[,1],t0)))
  t_tot <- t_tot[t_tot <= restrict_time]
  t_tot <- sort(unique(t_tot))
  
  
  
  t_cohort_sort <- ss$time[ss$time <= restrict_time]
  
  t_tot_large <- sort(unique(c(t_tot,t_tot - Delta)))
  
  
  
  t_tot_large_1 <- t_tot_large[t_tot_large>=0 & t_tot_large <= restrict_time]
  t_tot_large <- t_tot_large[t_tot_large>=0 & t_tot_large <= restrict_time-Delta]
  
  
  mt_event <- match(t_cohort_sort,t_tot_large)
  mt_event_Delta <- match(t_cohort_sort-Delta,t_tot_large)
  
  
  # Obtaining RMRL
  mt_event <- mt_event[!is.na(mt_event)]
  # DeltaYearRMRL and CS uncertaintly:
  
  dA <- diff(c(0,ss$cumhaz))[1:length(mt_event)]
  hazMatrix <- matrix(0,nrow = 2,ncol=length(t_tot_large_1))
  hazMatrix[1,mt_event] <- dA
  hazMatrix[2,] <- diff(c(0,t_tot_large_1))
  
  
  
  F_fun=F_restrictsurv <- function(X)matrix(c(-X[1],0,0,X[1]),nrow=2)
  X0=X_0_restrictsurv <- matrix(c(1,0),nrow=2)
  V0=V_0_restrictsurv <- matrix(0,nrow=2,ncol=2)
  JacobianList=F_restrictsurv_gradientlist <- list(function(X)matrix(c(-1,0,0,0),nrow=2),
                                                   function(X)matrix(c(0,1,0,0),nrow=2))
  
  X_restrictSurv <- pluginEstimate(max(ss$n.risk),hazMatrix,F_restrictsurv,F_restrictsurv_gradientlist,
                                   X_0_restrictsurv,V_0_restrictsurv,isLebesgue = 2)
  
  
  
  
  
  
  
  print("Obtaining RMRL...")
  RMRL_plugin_est <- sapply(2:length(t_tot_large),function(i)getCrossCov_Delta(t_tot_large_1,i,Delta,
                                                                               X_restrictSurv,nabla_H_RMRL,H_RMRL,
                                                                               F_restrictsurv_gradientlist,hazMatrix))
  
  
  RMRL_plugin_est <- cbind(c(X_restrictSurv$covariance[2,2,t_tot_large == Delta],X_restrictSurv$X[2,t_tot_large == Delta]),
                           RMRL_plugin_est)
  
  
  # Population estimates
  t_tot_pop <- unique(sort(c(seq(0.01,t_tot_large_1[length(t_tot_large_1)]-0.01,length.out = 2e4+7),
                             t0,t_cohort_sort)))
  t_tot_pop <- t_tot_pop[t_tot_pop<=restrict_time]
  
  mt_fine <- which(t_tot_pop %in% t0)
  hazMatrix_p <- matrix(0,nrow = 2,ncol=length(t_tot_pop))
  hazMatrix_p[1,mt_fine] <- d_H_e[t0<=restrict_time]
  hazMatrix_p[2,] <- diff(c(0,t_tot_pop))
  
  X_restrictSurv_p <- pluginEstimate(1,hazMatrix_p,F_restrictsurv,F_restrictsurv_gradientlist,
                                     X_0_restrictsurv,V_0_restrictsurv,isLebesgue = 2)
  
  RMRL_p <- CS_p <- rep(0,length(t_tot_pop))
  RMRL_p[1] <- X_restrictSurv_p$X[2,t_tot_pop== Delta]
  CS_p[1] <- X_restrictSurv_p$X[1,t_tot_pop== Delta]
  
  # plot(t0,cumprod(1-d_H_e),type="s",xlim=c(0,12))
  
  
  for(i in 2:length(t_tot_pop)){
    t2_ind <- max(which( t_tot_pop[i] + Delta >= t_tot_pop ) )

    
    RMRL_p[i] <- H_RMRL(X_restrictSurv_p$X[,i],X_restrictSurv_p$X[,t2_ind])
    
    CS_p[i] <- H_condsurv(X_restrictSurv_p$X[1,i],X_restrictSurv_p$X[1,t2_ind])
    # X_restrictSurv_p$X[1,t2_ind]/X_restrictSurv_p$X[1,i]
    
    
    
  }
  
  
  
  
  # plot(t_tot_pop,RMRL_p,type="s",ylim=c(4,5),xlim=c(0,restrict_time-Delta))
  
  
  iind <- sapply(t_tot_large,function(tt) min(which(t_tot_pop >= tt)) )
  RMRL_p <- RMRL_p[iind]
  
  CS_p <- CS_p[iind]
  
  # plot(t_tot_large,RMRL_p,type="s",ylim=c(4,5))
  # lines( t_tot_large,(RMS_p_Delta-RMS_p)/S_p,type="s",col=2,lty=2)
  RMRL_var <- RMRL_plugin_est[1,]
  RMRL_est <- RMRL_plugin_est[2,]
  
  
  # Conditional survival
  
  hazMat <- matrix(hazMatrix[1,],nrow=1)
  X_surv <- pluginEstimate(max(ss$n.risk),hazMat,function(X)matrix(-X,nrow=1,ncol=1),list(function(X)matrix(-1,nrow=1,ncol=1)),matrix(1),matrix(0))
  
  hazMat_p <- matrix(hazMatrix_p[1,],nrow=1)
  X_surv_p <- pluginEstimate(1,hazMat_p,function(X)matrix(-X,nrow=1,ncol=1),list(function(X)matrix(-1,nrow=1,ncol=1)),matrix(1),matrix(0))
  
  print("Obtaining conditional survival...")
  CS_plugin_est <- sapply(1:length(t_tot_large),function(i)getCrossCov_Delta(t_tot_large_1,i,Delta,
                                                                             X_surv,nabla_H_condsurv,H_condsurv,
                                                                             list(function(X)matrix(-1,nrow=1,ncol=1)),
                                                                             hazMat))
  
  CS_var <- CS_plugin_est[1,]
  CS_est <- CS_plugin_est[2,]
  
  
  
  
  
  
  
  
  # Cause-specific analysis
  
  t_loc <- survfit(t_y_strata ~ 1)$time
  # t_loc <- t_loc[1:sum(t_loc<= restrict_time)]
  dA1 = dA2 = rep(0,length(t_loc))
  mtCancer = match(survfit(t_y_other_spec ~ 1)$time,t_loc)
  mtOther = match(survfit(t_y_other_spec ~ 1)$time,t_loc)
  dA2[mtOther] <- diff(c(0,-log(survfit(t_y_other_spec ~ 1)$surv)))
  dA1[mtCancer] <- diff(c(0,-log(survfit(t_y_cancer_spec ~ 1)$surv)))
  
  
  # Setting up the cumulative incidence ODE system
  
  F_cuminc <- function(X)matrix(c(-X[1],X[1],0,-X[1],0,X[1]),nrow=3)
  X_0_cuminc <- matrix(c(1,0,0),nrow=3)
  V_0_cuminc <- matrix(0,nrow=3,ncol=3)
  F_cuminc_gradientlist <- list(function(X)matrix(c(-1,1,0,0,0,0,0,0,0),nrow=3),
                                function(X)matrix(c(-1,0,1,0,0,0,0,0,0),nrow=3))
  
  hazMatrix <- matrix(0,nrow = 2,ncol=length(t_loc))
  hazMatrix[1,] <- dA1
  hazMatrix[2,] <- dA2
  
  X_cuminc <- pluginEstimate(nrow(t_y_other_spec),hazMatrix,F_cuminc,F_cuminc_gradientlist,
                             X_0_cuminc,V_0_cuminc)
  
  
  
  
  print("Obtaining cuminc...")
  cuminc_plugin_est <- sapply(2:length(t_loc),function(i)getCrossCov_Delta(t_loc,i,Delta,
                                                                           X_cuminc,nabla_H_cuminc,H_cuminc,
                                                                           F_cuminc_gradientlist,hazMatrix))
  c_spec_cond_risk_var <- c(cuminc_plugin_est[1,1],cuminc_plugin_est[1,])
  c_spec_cond_risk <- c(cuminc_plugin_est[2,1],cuminc_plugin_est[2,])
  
  
  cuminc_pluginRR_est <- sapply(2:length(t_loc),function(i)getCrossCov_Delta(t_loc,i,Delta,
                                                                             X_cuminc,nabla_H_cuminc_RR,H_cuminc_RR,
                                                                             F_cuminc_gradientlist,hazMatrix))
  c_spec_cond_riskRR_var <- c(cuminc_pluginRR_est[1,1],cuminc_pluginRR_est[1,])
  c_spec_cond_riskRR <- c(cuminc_pluginRR_est[2,1],
                          cuminc_pluginRR_est[2,])
  
  
  cuminc_pluginRD_est <- sapply(2:length(t_loc),function(i)getCrossCov_Delta(t_loc,i,Delta,
                                                                             X_cuminc,nabla_H_cuminc_RD,H_cuminc_RD,
                                                                             F_cuminc_gradientlist,hazMatrix))
  c_spec_cond_riskRD_var <- c(cuminc_pluginRD_est[1,1],cuminc_pluginRD_est[1,])
  c_spec_cond_riskRD <- c(cuminc_pluginRD_est[2,1],
                          cuminc_pluginRD_est[2,])
  
  
  
  cuminc_plugin_proportion_est <- sapply(2:length(t_loc),function(i)getCrossCov_Delta(t_loc,i,Delta,
                                                                                      X_cuminc,nabla_H_cuminc_proportion,H_cuminc_proportion,
                                                                                      F_cuminc_gradientlist,hazMatrix))
  c_spec_cond_risk_proportion_var <- c(cuminc_plugin_proportion_est[1,1],cuminc_plugin_proportion_est[1,])
  c_spec_cond_risk_proportion <- c(cuminc_plugin_proportion_est[2,1],
                                   cuminc_plugin_proportion_est[2,])
  
  
  # convert to same time scale:
  
  t_change_inds <- sapply(t_tot_large,function(tmm)min(which(tmm <= t_loc)))
  # t_change_inds[t_change_inds == -Inf] <- 1
  
  c_spec_cond_risk <- c_spec_cond_risk[t_change_inds]
  c_spec_cond_risk_var <- c_spec_cond_risk_var[t_change_inds]
  c_spec_cond_riskRR <- c_spec_cond_riskRR[t_change_inds]
  c_spec_cond_riskRR_var <- c_spec_cond_riskRR_var[t_change_inds]
  c_spec_cond_riskRD <- c_spec_cond_riskRD[t_change_inds]
  c_spec_cond_riskRD_var <- c_spec_cond_riskRD_var[t_change_inds]
  c_spec_cond_risk_proportion <- c_spec_cond_risk_proportion[t_change_inds]
  c_spec_cond_risk_proportion_var <- c_spec_cond_risk_proportion_var[t_change_inds]
  
  # plot(t_tot_large,c_spec_cond_risk,type="s",col=3)
  # lines(t_tot_large,c_spec_cond_risk+1.96*sqrt(c_spec_cond_risk_var),type="s",col=3)
  # lines(t_tot_large,c_spec_cond_risk-1.96*sqrt(c_spec_cond_risk_var),type="s",col=3)
  ############
  
  return(list(t_tot_large=t_tot_large,CS_est=CS_est,CS_p=CS_p,CS_var=CS_var,
              RMRL_est=RMRL_est,RMRL_p=RMRL_p,RMRL_var=RMRL_var,
              c_spec_cond_risk=c_spec_cond_risk,c_spec_cond_risk_var=c_spec_cond_risk_var,
              c_spec_cond_riskRR=c_spec_cond_riskRR,c_spec_cond_riskRR_var=c_spec_cond_riskRR_var,
              c_spec_cond_riskRD=c_spec_cond_riskRD,c_spec_cond_riskRD_var=c_spec_cond_riskRD_var,
              c_spec_cond_risk_proportion=c_spec_cond_risk_proportion,
              c_spec_cond_risk_proportion_var=c_spec_cond_risk_proportion_var))
  
}



