#setwd("~/Documents/HD_MP_inverse/Final_prog")
#setwd("/u/bodnar/ProgR/HD_MPinverse")
#setwd("C:/D/2023/HD MPInverse/figOctober2023")
#library(HDShOP)

library(Rfast)
library(Rfast2)
library(kStatistics)
library(corpcor)
#library(abind)

Sigma_sample_estimator <- function(x) {
  p <- nrow(x)
  n <- ncol(x)
  if (is.data.frame(x)) x <- as.matrix(x)
  a <- .rowMeans(x, m=p, n=n, na.rm = TRUE)
  a_x_size <- matrix(rep(a,n),nrow=p, ncol=n)
  tcrossprod(x-a_x_size)/(ncol(x)-1)
}


nonlin_shrinkLW = function(x){
  # the original version suggested that p is # of columns
  p = nrow(x)
  n = ncol(x)
  sampleC = Sigma_sample_estimator(x)
  eig = eigen(sampleC)
  #eig = eigenS
  u = eig$vectors[,p:1]
  lambda = rev(eig$values)
  
  if(p<=n){
    lambda = lambda[max(1, p-n+1):p]
    L = matrix(rep(lambda, min(p, n)), nrow = length(lambda))
    h = n^(-1/3)
    H = h * t(L)
    x <- (L - t(L)) / H # This is a different x than before
    ftilde = (3/4/sqrt(5)) * rowMeans(pmax(1-x^2/5, 0) / H)
    Hftemp = (-3/10/pi) * x + (3/4/sqrt(5)/pi) * (1 - x^2./5) * log(abs((sqrt(5) - x)/(sqrt(5) + x)))
    Hftemp[abs(x) == sqrt(5)] = (-3/10/pi) * x[abs(x) == sqrt(5)]
    Hftilde = rowMeans(Hftemp / H)
    dtilde = lambda / ((pi*(p/n)*lambda*ftilde)^2 + (1-(p/n)-pi*(p/n)*lambda*Hftilde)^2);
  }else{
    lambda = lambda[max(1, p-n+2):p]
    L = matrix(rep(lambda, min(p, n-1)), nrow = length(lambda))
    h = n^(-1/3)
    H = h * t(L)
    x <- (L - t(L)) / H # This is a different x than before
    ftilde = (3/4/sqrt(5)) * rowMeans(pmax(1-x^2/5, 0) / H)
    Hftemp = (-3/10/pi) * x + (3/4/sqrt(5)/pi) * (1 - x^2./5) * log(abs((sqrt(5) - x)/(sqrt(5) + x)))
    Hftemp[abs(x) == sqrt(5)] = (-3/10/pi) * x[abs(x) == sqrt(5)]
    Hftilde = rowMeans(Hftemp / H)
        
    Hftilde0 = (1/pi)*(3/10/h^2+3/4/sqrt(5)/h*(1-1/5/h^2) *log((1+sqrt(5)*h)/(1-sqrt(5)*h))) * mean(1/lambda)
    dtilde0=1/(pi*(p-n)/n*Hftilde0)
    dtilde1 = lambda/(pi^2*lambda^2 * (ftilde^2+Hftilde^2))
    dtilde = c(dtilde0 * rep(1, p-n+1), dtilde1)
  }
  u %*% diag(dtilde) %*% t(u)
} # analytical nonlinear shrinkage

B<-100
n<-250
c_val<-c(seq(1.5,3,0.1),seq(3.5,5,0.5))
#c_val<-c(seq(1.1,3,0.1))
p_val<- c_val*n
lc<-length(c_val)

#### for numeric optimization
eps<-1/(10^6)
upp<-pi/2-eps

MP<-matrix(0,lc,B)
SSE<-matrix(0,lc,B)
NL<-matrix(0,lc,B)
KS<-matrix(0,lc,B)
WPTZ<-matrix(0,lc,B)
ShMP<-matrix(0,lc,B)
ShRt<-matrix(0,lc,B)
ShMP2<-matrix(0,lc,B)
ShMP3<-matrix(0,lc,B)
ShMP4<-matrix(0,lc,B)
ShMP5<-matrix(0,lc,B)
NLo<-matrix(0,lc,B)


PRIAL<-matrix(0,lc,13)



for (i_c in 1:lc)
{

c_n<-c_val[i_c]
p<-p_val[i_c]
p1<-0.2*p
p2<-0.4*p
Ip<-diag(p)
t1<-p^(-1)
t2<-p^(-1/4)
c_n2<-c_n^2
r<-(c_n-1)/c_n

############################true covariance matrix
D<-diag(c(rep(1,p1),rep(3,p2),rep(10,p-p1-p2)))
sqD<-diag(sqrt(c(rep(1,p1),rep(3,p2),rep(10,p-p1-p2))))
X0<-matrnorm(p,10*p)
H<-eigen(X0%*%t(X0))$vectors

Sig<-H%*%D%*%t(H)

iSig<-H%*%spdinv(D)%*%t(H)
Sig2<-H%*%(D^2)%*%t(H)
sqSig<-H%*%sqD%*%t(H)




for (ib in 1:B)
{

dist<-"nor"


if (dist=="nor") {X<- matrnorm(p,n)}else{X<-matrix(rt(p*n,5),p,n)*sqrt(3/5)}

  Y<-sqSig%*%X
  iY<-spdinv(t(Y)%*%Y/n)
  S<-Y%*%t(Y)/n
  eigenS<-eigen(S)
  U<-eigenS$vectors
  ##### MP inverse
  iS_MP<-Y%*%iY%*%iY%*%t(Y)/n
  D_MP<- diag(eigen(iS_MP)$values)
  ##### SSE
  iS_SSE<-iS_MP*p/(n-1)
  ##### NL
  S_NL<-nonlin_shrinkLW(Y)
  iS_NL<-spdinv(S_NL)
  
  ##### LW oracle
  #IS_LWo<-U%*%diag(diag(t(U)%*%Sig2%*%U)/diag(t(U)%*%Sig%*%U))%*%t(U)
  IS_LWo<-U%*%diag(diag(t(U)%*%Sig%*%U)/diag(t(U)%*%Sig2%*%U))%*%t(U)
  
  ##### EBRT KS
  iS_KS<-p*spdinv(n*S +sum(diag(S))*Ip)

  ##### optimal ridge
  
  hL_WPTZ<- function(u)
  {t<-tan(u)
  iS_t<-spdinv(S/t+Ip)
  tr_iS_t<-sum(diag(iS_t))/p
  a1<-1-tr_iS_t
  a2<-tr_iS_t-sum(diag(iS_t%*%iS_t))/p

  hR1<-a1/(1-c_n*a1)
  hR2<-a1/((1-c_n*a1)^3)-a2/((1-c_n*a1)^4)

  hL_WPTZ<-hR1^2/hR2
  return(hL_WPTZ)
  }

  hL_WPTZ_max<-optim(1.5, hL_WPTZ,lower = eps, upper = upp, method= "L-BFGS-B", control = list(fnscale = -1))
  hL_WPTZ_bet<-tan(hL_WPTZ_max$par)

  iS_t<-spdinv(S/hL_WPTZ_bet+Ip)
  tr_iS_t<-sum(diag(iS_t))/p
  a1<-1-tr_iS_t
  a2<-tr_iS_t-sum(diag(iS_t%*%iS_t))/p

  hR1<-a1/(1-c_n*a1)
  hR2<-a1/((1-c_n*a1)^3)-a2/((1-c_n*a1)^4)
  hL_WPTZ_alp<-hR1/hR2

  iS_WPTZ<-hL_WPTZ_alp*spdinv(S +hL_WPTZ_bet*Ip)


  ##### shrinkage MP
  trS1<-sum(diag(iS_MP))/p
  trS2<-sum(diag(iS_MP%*%iS_MP))/p
  trS3<-sum(diag(iS_MP%*%iS_MP%*%iS_MP))/p

  hv0<-c_n*trS1
  ihv0<-1/hv0
  ihv0_2<-ihv0^2
  ihv0_3<-ihv0^3

  h2<-1/trS2/c_n
  h3<-trS3/(trS2^3)/c_n2

  d1<-trS1/trS2/c_n
  d2<-(trS1*trS3-trS2^2)/(trS2^3)/c_n2

  q1<-sum(diag(S))/p
  q2<-sum(diag(S%*%S))/p-c_n*(q1^2)

  d1Sig<-ihv0*(ihv0/c_n-d1)
  d1Sig2<-ihv0_2*(q1+d1-2*ihv0/c_n)
  d2Sig2<-ihv0*d1Sig2-ihv0_2*(d1Sig-d2)

  num_a_MP<-d1Sig*q2-d1Sig2*q1
  num_b_MP<--(d2Sig2-d1Sig2*h3/h2)*q1/h2-d1Sig*d1Sig2/h2
  den_MP<--(d2Sig2-d1Sig2*h3/h2)*q2/h2-d1Sig2^2/h2
  ha_MP<-num_a_MP/den_MP
  hb_MP<-num_b_MP/den_MP
  iS_ShMP<-ha_MP*iS_MP+hb_MP*Ip


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


for (i in 1:m)
{
v[i]<- (-1)^i*factorial(i)*c_n*(1/p)*sum(diag(D_MP^(i+1)))
}



hm[1]<- q1


h[2]<- h2
h[3]<- h3



if (m>1)
{ 
for (i in 3:m)
{ sh<- 0
for (k in 2:(i-1))
{
sh<- sh + (-1)^k*factorial(k)*h[k+1]*e_eBellPol(i,k,c(v[1:(i-k+1)],rep(0,k-1)))
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


  ##### Computations of PRIAL
  PRIAL[i_c,1]<-PRIAL[i_c,1]+sum(diag((iS_MP%*%Sig-Ip)%*%t(iS_MP%*%Sig-Ip)))/p/B
  
  PRIAL[i_c,2]<-PRIAL[i_c,2]+sum(diag((iS_SSE%*%Sig-Ip)%*%t(iS_SSE%*%Sig-Ip)))/p/B
  PRIAL[i_c,3]<-PRIAL[i_c,3]+sum(diag((iS_NL%*%Sig-Ip)%*%t(iS_NL%*%Sig-Ip)))/p/B
  PRIAL[i_c,4]<-PRIAL[i_c,4]+sum(diag((iS_KS%*%Sig-Ip)%*%t(iS_KS%*%Sig-Ip)))/p/B
  PRIAL[i_c,5]<-PRIAL[i_c,5]+sum(diag((iS_WPTZ%*%Sig-Ip)%*%t(iS_WPTZ%*%Sig-Ip)))/p/B
  PRIAL[i_c,6]<-PRIAL[i_c,6]+sum(diag((iS_ShMP%*%Sig-Ip)%*%t(iS_ShMP%*%Sig-Ip)))/p/B
  PRIAL[i_c,7]<-PRIAL[i_c,7]+sum(diag((iS_ShRt1%*%Sig-Ip)%*%t(iS_ShRt1%*%Sig-Ip)))/p/B
  PRIAL[i_c,8]<-PRIAL[i_c,8]+sum(diag((high_shrink2%*%Sig-Ip)%*%t(high_shrink2%*%Sig-Ip)))/p/B
  PRIAL[i_c,9]<-PRIAL[i_c,9]+sum(diag((high_shrink3%*%Sig-Ip)%*%t(high_shrink3%*%Sig-Ip)))/p/B
  PRIAL[i_c,10]<-PRIAL[i_c,10]+sum(diag((high_shrink4%*%Sig-Ip)%*%t(high_shrink4%*%Sig-Ip)))/p/B
  PRIAL[i_c,11]<-PRIAL[i_c,11]+sum(diag((high_shrink5%*%Sig-Ip)%*%t(high_shrink5%*%Sig-Ip)))/p/B
  PRIAL[i_c,12]<-PRIAL[i_c,12]+sum(diag((IS_LWo%*%Sig-Ip)%*%t(IS_LWo%*%Sig-Ip)))/p/B
  
  PRIAL[i_c,13]<-PRIAL[i_c,13]+sum(diag((Sig-Ip)%*%t(Sig-Ip)))/p/B

  MP[i_c,ib]<-sum(diag((iS_MP%*%Sig-Ip)%*%t(iS_MP%*%Sig-Ip)))/p
  SSE[i_c,ib]<-sum(diag((iS_SSE%*%Sig-Ip)%*%t(iS_SSE%*%Sig-Ip)))/p
  NL[i_c,ib]<-sum(diag((iS_NL%*%Sig-Ip)%*%t(iS_NL%*%Sig-Ip)))/p
  NLo[i_c,ib]<-sum(diag((IS_LWo%*%Sig-Ip)%*%t(IS_LWo%*%Sig-Ip)))/p/B
  KS[i_c,ib]<-sum(diag((iS_KS%*%Sig-Ip)%*%t(iS_KS%*%Sig-Ip)))/p
  WPTZ[i_c,ib]<-sum(diag((iS_WPTZ%*%Sig-Ip)%*%t(iS_WPTZ%*%Sig-Ip)))/p
  ShMP[i_c,ib]<-sum(diag((iS_ShMP%*%Sig-Ip)%*%t(iS_ShMP%*%Sig-Ip)))/p
  ShRt[i_c,ib]<-sum(diag((iS_ShRt1%*%Sig-Ip)%*%t(iS_ShRt1%*%Sig-Ip)))/p
  ShMP2[i_c,ib]<-sum(diag((high_shrink2%*%Sig-Ip)%*%t(high_shrink2%*%Sig-Ip)))/p
  ShMP3[i_c,ib]<-sum(diag((high_shrink3%*%Sig-Ip)%*%t(high_shrink3%*%Sig-Ip)))/p
  ShMP4[i_c,ib]<-sum(diag((high_shrink4%*%Sig-Ip)%*%t(high_shrink4%*%Sig-Ip)))/p
  ShMP5[i_c,ib]<-sum(diag((high_shrink5%*%Sig-Ip)%*%t(high_shrink5%*%Sig-Ip)))/p
  
}
print(i_c)
print(n)
print(dist)
}



write.table(MP,file=paste("MP_",dist,"_n",n,".csv",sep=""), sep=",", dec=".", row.names = F,col.names = F)
write.table(SSE,file=paste("SSE_",dist,"_n",n,".csv",sep=""), sep=",", dec=".", row.names = F,col.names = F)
write.table(NL,file=paste("NL_",dist,"_n",n,".csv",sep=""), sep=",", dec=".", row.names = F,col.names = F)
write.table(NLo,file=paste("NLo_",dist,"_n",n,".csv",sep=""), sep=",", dec=".", row.names = F,col.names = F)
write.table(KS,file=paste("KS_",dist,"_n",n,".csv",sep=""), sep=",", dec=".", row.names = F,col.names = F)
write.table(WPTZ,file=paste("WPTZ_",dist,"_n",n,".csv",sep=""), sep=",", dec=".", row.names = F,col.names = F)
write.table(ShMP,file=paste("ShMP_",dist,"_n",n,".csv",sep=""), sep=",", dec=".", row.names = F,col.names = F)
write.table(ShRt,file=paste("ShRt_",dist,"_n",n,".csv",sep=""), sep=",", dec=".", row.names = F,col.names = F)
write.table(ShMP2,file=paste("ShMP2_",dist,"_n",n,".csv",sep=""), sep=",", dec=".", row.names = F,col.names = F)
write.table(ShMP3,file=paste("ShMP3_",dist,"_n",n,".csv",sep=""), sep=",", dec=".", row.names = F,col.names = F)
write.table(ShMP4,file=paste("ShMP4_",dist,"_n",n,".csv",sep=""), sep=",", dec=".", row.names = F,col.names = F)
write.table(ShMP5,file=paste("ShMP5_",dist,"_n",n,".csv",sep=""), sep=",", dec=".", row.names = F,col.names = F)
write.table(PRIAL,file=paste("PRIAL",dist,"_n",n,".csv",sep=""), sep=",", dec=".", row.names = F,col.names = F)




res<-100*(1-PRIAL[,c(2:12)]/(PRIAL[,1]%*%matrix(1,1,11)))

#res<-rOSV
res<- apply(res, c(1,2), Re)
#c_val<-c(seq(1.5,3,0.1),seq(3.5,5,0.5))

#MP<- apply(MP, 1, median)/V_GMV-1
#BPS<- apply(BPS, 1, median)/V_GMV-1
#ShMP<- apply(ShMP, 1, median)/V_GMV-1
#ShRt<- apply(ShRt, 1, median)/V_GMV-1
#WPTZ<- apply(WPTZ, 1, median)/V_GMV-1
#NLo<- apply(NLo, 1, median)/V_GMV-1


#res<- cbind(NLo, BPS, ShMP, ShRt, WPTZ)
#res<- apply(res, 1:2, Re)

#ylim = range(res[,1:ncol(res)])


c_val<-c(seq(1.5,3,0.1),seq(3.5,5,0.5))
PRIAL <- read.table(paste("PRIAL", dist, "_n", n, ".csv", sep=""), sep=",", dec=".", header=FALSE)

pdf(paste("PRIAL_prec_",dist,"_n", n, ".pdf", sep=""), height=5.5, width=7)
plot(c(rep(c_val,ncol(res)) ),c(res[,1:ncol(res)]), type="n", main="", xlab=expression(c[n]), ylab="PRIAL", mgp=c(2.2,1,0), xlim = c(1.5,5), ylim = c(55,100),cex=1.5, lwd = 1.6, cex.axis=1.2,cex.lab=1.5,axes=T,frame=T,font=2)
#lines(c_val, res[,2], lty=7, col="lightblue", lwd = 1.4)
lines(c_val, res[,11], lty=1, col="blue", lwd = 1.4)
#lines(c_val, res[,6], lty=7, col="red", lwd = 1.4)
#lines(c_val, res[,3], lty=5, col="green", lwd = 1.4)
lines(c_val, res[,5], lty=1, col="orange", lwd = 1.4)
lines(c_val, res[,7], lty=2, col="orange", lwd = 1.4)
lines(c_val, res[,8], lty=3, col="orange", lwd = 1.4)
lines(c_val, res[,9], lty=4, col="orange", lwd = 1.4)
lines(c_val, res[,10], lty=5, col="orange", lwd = 1.4)
legend("topright", lty=c(1, 1,2,3,4,5), col=c("blue", rep("orange", 5)), c("NL oracle", "MP shrinkage (m=1)", "High order m=2", "High order m=3", "High order m=4", "High order m=5"), bty = "n", cex=1.0, lwd = 2.0)
dev.off()


pdf(paste("zoomedPRIAL_prec_",dist,"_n", n, ".pdf", sep=""), height=5.5, width=7)
plot(c(rep(c_val,ncol(res)) ),c(res[,1:ncol(res)]), type="n", main="", xlab=expression(c[n]), ylab="PRIAL", mgp=c(2.2,1,0), xlim = c(2.5,5), ylim = c(60,70),cex=1.5, lwd = 1.6, cex.axis=1.2,cex.lab=1.5,axes=T,frame=T,font=2)
#lines(c_val, res[,2], lty=7, col="lightblue", lwd = 1.4)
lines(c_val, res[,11], lty=1, col="blue", lwd = 1.4)
#lines(c_val, res[,6], lty=7, col="red", lwd = 1.4)
#lines(c_val, res[,3], lty=5, col="green", lwd = 1.4)
lines(c_val, res[,5], lty=1, col="orange", lwd = 1.4)
lines(c_val, res[,7], lty=2, col="orange", lwd = 1.4)
lines(c_val, res[,8], lty=3, col="orange", lwd = 1.4)
lines(c_val, res[,9], lty=4, col="orange", lwd = 1.4)
lines(c_val, res[,10], lty=5, col="orange", lwd = 1.4)
legend("topright", lty=c(1, 1,2,3,4,5), col=c("blue", rep("orange", 5)), c("NL oracle", "MP shrinkage (m=1)", "High order m=2", "High order m=3", "High order m=4", "High order m=5"), bty = "n", cex=1.0, lwd = 2.0)
dev.off()







