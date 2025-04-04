

# Trace function of a matrix
tr <- function(M){
  return (sum(diag(M)))
}


#' Higher order shrinkage
#' 
#' 
#' @param Y data matrix (rows are features, columns are observations).
#' TODO: transpose everything.
#' 
#' @param Pi0 prior of the precision matrix. This a `p` by `p` matrix, used as
#' a target for the shrinkage. Default value is the identity matrix of size `p`.
#' As an advice, it should be a symmetric positive-definite matrix, but this is
#' not checked for.
#' 
#' 
higher_order_shrinkage <- function(Y, Pi0 = NULL)
{
  
  
  
  ##### Shrinkage m order
  
  M<- function(m)
  {
    v<-rep(0,2*m)
    h<- rep(0,(2*m+1))
    d<- matrix(0,nrow=2*m, ncol=2)
    s<- rep(q2,(2*m+1))
    
    h[2]<- h2
    h[3]<- h3
    
    for (i in 1:(2*m))
    {
      v[i]<- (-1)^i*factorial(i)*c_n*(1/p)*sum(diag(D_MP^(i+1)))
    }
    
    
    if (m>1)
    { 
      for (i in 3:(2*m))
      { sh<- 0
      for (k in 2:(i-1))
      {
        sh<- sh + (-1)^k*factorial(k)*h[k+1]*e_eBellPol(i,k,c(v[1:(i-k+1)],rep(0,k-1)))
      }
      h[i+1]<-(v[i]+v[1]*sh)/((v[1])^(i+1)*(-1)^(i+1)*factorial(i))
      }
    }
    
    
    for (k in 1:(2*m))
    {
      for (l in 1:2)
      {
        if (l==1) d[k,l]<- 1/c_n*(1/(hv0^(k+1))-h[k+1])
        if (l==2) ifelse(k==1, d[k,l]<-ihv0*(ihv0*(q1-ihv0*(1/c_n))-d[1,1]), d[k,l]<- ihv0*(d[k-1,2]-d[k,1]))
      }
    }
    
    
    for (j in 1:(2*m))
    {ss<-0
    for (k in 1:j)
    {
      ss<- ss+(-1)^(j+k+1)*factorial(k)/factorial(j)*d[k,2]*e_eBellPol(j,k,c(v[1:(j-k+1)],rep(0,k-1)))
    }
    s[j+1]<-ss
    }
    
    
    
    M<- s[1:(m+1)]
    for (j in 2:(m+1))
    {
      M<- cbind(M, s[j:(m+j)])
    }
    
    
    M<- apply(M, c(1,2), Re)
    return(M)
  }
  
  
  hm<- function(m)
  {
    hm<- rep(0,(m+1))
    h<- rep(0,(m+1))
    d<- rep(0,m)
    v<-rep(0,m)
    
    
    for (i in 1:m) {
      v[i]<- (-1)^i*factorial(i)*c_n*(1/p)*sum(diag(D_MP^(i+1)))
    }
    
    
    hm[1]<- q1
    
    
    h[2]<- h2
    h[3]<- h3
    
    
    if (m > 1) { 
      for (i in 3:m) {
        sh <- 0
        for (k in 2:(i-1)) {
          sh <- sh + (-1)^k * factorial(k) * h[k+1] * e_eBellPol(i,k,c(v[1:(i-k+1)],
                                                                       rep(0,k-1)))
        }
        h[i+1]<-(v[i]+v[1]*sh)/((v[1])^(i+1)*(-1)^(i+1)*factorial(i))
      }
    }
    
    for (k in 1:m)
    {
      d[k]<- 1/c_n*(1/(hv0^(k+1))-h[k+1])
    }
    
    hm[1]<-q1
    for (j in 1:m)
    {s<-0
    for (k in 1:j)
    {
      s<- s + (-1)^(j+k+1)*factorial(k)/factorial(j)*d[k]*e_eBellPol(j,k,c(v[1:(j-k+1)],rep(0,k-1)))
    }
    hm[j+1]<- s
    }
    hm<- Re(hm)
    return(hm)
  }
  
  
  
  
  
  #alpha_m1<- spdinv(M(1))%*%hm(1)
  alpha_m2<- spdinv(M(2))%*%hm(2)
  alpha_m3<-spdinv(M(3))%*%hm(3)
  alpha_m4<-spdinv(M(4))%*%hm(4)
  alpha_m5<-spdinv(M(5))%*%hm(5)
  
  
  
  
  
  
  #high_shrink1<- alpha_m1[1]*Ip+alpha_m1[2]*iS_MP
  high_shrink2<- alpha_m2[1]*Ip + alpha_m2[2]*iS_MP + alpha_m2[3]*iS_MP%*%iS_MP
  high_shrink3<-  alpha_m3[1]*Ip + alpha_m3[2]*iS_MP + alpha_m3[3]*iS_MP%*%iS_MP + alpha_m3[4]*iS_MP%*%iS_MP%*%iS_MP
  high_shrink4<- alpha_m4[1]*Ip + alpha_m4[2]*iS_MP + alpha_m4[3]*iS_MP%*%iS_MP + alpha_m4[4]*iS_MP%*%iS_MP%*%iS_MP + alpha_m4[5]*iS_MP%*%iS_MP%*%iS_MP%*%iS_MP
  high_shrink5<- alpha_m5[1]*Ip + alpha_m5[2]*iS_MP + alpha_m5[3]*iS_MP%*%iS_MP + alpha_m5[4]*iS_MP%*%iS_MP%*%iS_MP + alpha_m5[5]*iS_MP%*%iS_MP%*%iS_MP%*%iS_MP+alpha_m5[6]*iS_MP%*%iS_MP%*%iS_MP%*%iS_MP%*%iS_MP
  
  
  
  ##### shrinkage Ridge
  
  hL2R<-function(u)
  {t<-tan(u)
  iS_Rt<-spdinv(S+t*Ip)
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
  
  hL2R_max<-optim(1.5, hL2R,lower = eps, upper = upp, method= "L-BFGS-B", control = list(fnscale = -1))
  u_R<- hL2R_max$par
  
  t_R<-tan(u_R)
  iS_Rt1<-spdinv(S+t_R*Ip)
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