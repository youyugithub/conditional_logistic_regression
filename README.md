# conditional_logistic_regression

```
###################################
# conditional logistic regression #
###################################

set.seed(0)
ndim<-4
nind<-10000
nstratum<-nind/2
xxmat<-matrix(NA,nind,ndim)
yyvec<-rep(NA,nind)
beta<-seq(0.25,1,0.25)
for(ii in 1:nstratum){
  success<-F
  while(!success){
    xx_ctrl<-rnorm(ndim)
    xx_case<-rnorm(ndim)
    yy_ctrl<-rbinom(1,1,1/(1+exp(-sum(beta*xx_ctrl))))
    yy_case<-rbinom(1,1,1/(1+exp(-sum(beta*xx_case))))
    if(yy_ctrl==0&yy_case==1)success<-T
  }
  xxmat[2*ii-1,]<-xx_ctrl
  xxmat[2*ii,]<-xx_case
  yyvec[2*ii-1]<-yy_ctrl
  yyvec[2*ii]<-yy_case
}

beta_est<-rep(0,ndim)

for(iter in 1:10){
  SS<-rep(0,ndim)
  II<-matrix(0,ndim,ndim)
  expterm<-c(exp(xxmat%*%beta_est))
  for(ii in 1:nstratum){
    SS<-SS+xxmat[2*ii,]
    avec<-
      expterm[2*ii-1]/(expterm[2*ii-1]+expterm[2*ii])*xxmat[2*ii-1,]+
      expterm[2*ii]/(expterm[2*ii-1]+expterm[2*ii])*xxmat[2*ii,]
    SS<-SS-avec
    amat<-
      expterm[2*ii-1]/(expterm[2*ii-1]+expterm[2*ii])*tcrossprod(xxmat[2*ii-1,])+
      expterm[2*ii]/(expterm[2*ii-1]+expterm[2*ii])*tcrossprod(xxmat[2*ii,])
    II<-II-amat+tcrossprod(avec)
  }
  beta_est<-beta_est-solve(II,SS)
  print(beta_est)
}

beta_est<-rep(0,ndim)
for(iter in 1:10){
  expterm<-c(exp(xxmat%*%beta_est))
  for(kk in 1:ndim){
    dl<-0
    ddl<-0
    for(ii in 1:nstratum){
      dl<-dl+xxmat[2*ii,kk]
      temp1<-
        expterm[2*ii-1]/(expterm[2*ii-1]+expterm[2*ii])*xxmat[2*ii-1,kk]+
        expterm[2*ii]/(expterm[2*ii-1]+expterm[2*ii])*xxmat[2*ii,kk]
      dl<-dl-temp1
      temp2<-
        expterm[2*ii-1]/(expterm[2*ii-1]+expterm[2*ii])*xxmat[2*ii-1,kk]^2+
        expterm[2*ii]/(expterm[2*ii-1]+expterm[2*ii])*xxmat[2*ii,kk]^2
      ddl<-ddl-temp2+temp1*temp1
    }
    beta_est[kk]<-beta_est[kk]-dl/ddl
  }
  print(beta_est)
}
```
