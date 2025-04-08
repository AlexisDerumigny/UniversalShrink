


#' Ridge with target set to the identity
#' 
#' 
#' @param Y data matrix (rows are features, columns are observations).
#' TODO: transpose everything.
#' 
#' 
#' 
ridge_target_identity <- function (Y){
  
  
  ##### shrinkage Ridge
  
  hL2R<-function(u)
  {
    t<-tan(u)
    iS_Rt<-solve(S+t*Ip)
    trS1_t<-sum(diag(iS_Rt))/p
    trS2_t<-sum(diag(iS_Rt%*%iS_Rt))/p
    
    hvt<-c_n*(trS1_t-r/t)
    hvprt<--c_n*(trS2_t-r/t/t)
    ihvt<-1/hvt
    ihvt_2<-ihvt^2
    
    d0_t<-t*trS1_t
    d1_t<- -(t*trS2_t-d0_t/t)/(trS2_t-r/t/t)/c_n
    
    d0Sig_t<-ihvt/c_n-t/c_n
    d0Sig2_t<-ihvt*(q1-d0Sig_t)
    d1Sig2_t<-ihvt_2*(q1+d1_t-2*d0Sig_t)
    
    num_a_ShRt<-d0Sig_t*q2-d0Sig2_t*q1
    den_ShRt<-(d0Sig2_t+t*hvprt*d1Sig2_t)*q2-d0Sig2_t^2
    L2_ShRt<-num_a_ShRt^2/den_ShRt/q2
    return(L2_ShRt)
  }
  
  hL2R_max <- optim(1.5, hL2R,lower = eps, upper = upp, method= "L-BFGS-B", control = list(fnscale = -1))
  u_R<- hL2R_max$par
  
  t_R<-tan(u_R)
  iS_Rt1<-solve(S+t_R*Ip)
  trS1_t1<-sum(diag(iS_Rt1))/p
  trS2_t1<-sum(diag(iS_Rt1%*%iS_Rt1))/p
  
  hvt1<-c_n*(trS1_t1-r/t_R)
  hvprt1<--c_n*(trS2_t1-r/t_R/t_R)
  ihvt1<-1/hvt1
  ihvt1_2<-ihvt1^2
  
  d0_t1<-t_R*trS1_t1
  d1_t1<- -(t_R*trS2_t1-d0_t1/t_R)/(trS2_t1-r/t_R/t_R)/c_n
  
  d0Sig_t1<-ihvt1/c_n-t_R/c_n
  d0Sig2_t1<-ihvt1*(q1-d0Sig_t1)
  d1Sig2_t1<-ihvt1_2*(q1+d1_t1-2*d0Sig_t1)
  
  num_a_ShRt1<-d0Sig_t1*q2-d0Sig2_t1*q1
  num_b_ShRt1<-(d0Sig2_t1/t_R+hvprt1*d1Sig2_t1)*q1-d0Sig_t1*d0Sig2_t1/t_R
  den_ShRt1<-(d0Sig2_t1/t_R+hvprt1*d1Sig2_t1)*q2-d0Sig2_t1^2/t_R
  ha_ShRt1<-num_a_ShRt1/den_ShRt1
  hb_ShRt1<-num_b_ShRt1/den_ShRt1
  iS_ShRt1<-ha_ShRt1*iS_Rt1+hb_ShRt1*Ip
  
}

