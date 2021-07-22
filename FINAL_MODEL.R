
setwd("~/Dropbox/carStuff")

rm(list=ls())
set.seed(2021)

###
### Libraries
###

library(invgamma)
library(MASS)
library(Matrix)
library(mvtnorm)
library(raster)

###
### Functions
###

neighborhood=function(raster){
  nn=matrix(NA,length(raster[]),4)
  for(i in 1:dim(nn)[1]){
    loc=adjacent(raster,i)[,2]
    if(length(loc[which((loc+1)==i)])>0){
      ln=loc[which((loc+1)==i)]
    }else{ln=NA}
    if(length(loc[which((loc-1)==i)])>0){
      rn=loc[which((loc-1)==i)]
    }else{rn=NA}
    if(length(loc[which((loc-dim(raster)[2])==i)])>0){
      bn=loc[which((loc-dim(raster)[2])==i)]
    }else{bn=NA}
    if(length(loc[which((loc+dim(raster)[2])==i)])>0){
      tn=loc[which((loc+dim(raster)[2])==i)]
    }else{tn=NA}
    nn[i,1]=ln
    nn[i,2]=rn
    nn[i,3]=bn
    nn[i,4]=tn
  }
  nn
}

ratesLike <- function(x, prob, type, z, v, t){
  
  # Survival
  if(type == 'phi'){
    ugh <- sum(dbinom(z[[x]][,t+1],
                         v[[x]][,t],
                         prob[x], log=TRUE), na.rm=T)
    return(ugh)
  }
  
  # Harvest
  if(type == 'kappa'){
    ugh <- sum(dbinom(v[[x]][,t],
                         z[[x]][,t],
                         1-prob[x], log=T), na.rm=T)
    return(ugh)
  }
  
}


###
### Real or simulated?
###

realdata=T
pretuned=F


###
### Design
###

n=150    # individuals
S=100    # sites
T=7      # years
nr=sqrt(S)
nc=sqrt(S)
xmin=1
xmax=nc
ymin=1
ymax=nr
max=c(xmax,ymax)
min=c(xmin,ymin)

###
### Raster object
###

datad.r=datac.r=raster(,nrows=nr,ncols=nc,xmn=xmin,
                       xmx=xmax,ymn=ymin,ymx=ymax,)

###
### Neighborhood matrix
###

nn=neighborhood(datac.r)
W=matrix(0,nr*nc,nr*nc)
for(i in 1:(nr*nc)){
  W[i,nn[i,]]=1
}
D=diag(rowSums(W))

###
### True parameter values
###

##
## Correlation matrices
##

if(!realdata){
  rho.truth=0.99
  s2_kappa.truth=0.15 #0.35
  s2_phi.truth=0.25 #0.5
  Sigma_kappa.truth=s2_kappa.truth*solve(D-rho.truth*W)
  Sigma_kappa.bd.truth=kronecker(diag(T-1),Sigma_kappa.truth)
  Sigma_phi.truth=s2_phi.truth*solve(D-rho.truth*W)
  Sigma_phi.bd.truth=kronecker(diag(T-2),Sigma_phi.truth)
  
  ##
  ## means
  ##
  
  if(!realdata){
    eta.truth=matrix(NA, S, T)
    eta.truth[,1]=rnorm(S, 0, 1.5)
    for(t in 2:T){
      eta.truth[,t]=rmvnorm(1, eta.truth[,t-1],
                            Sigma_kappa.truth)
    }
    
    epsilon.truth=matrix(NA, S, T-1)
    epsilon.truth[,1]=rnorm(S, 0, 1.5)
    for(t in 2:(T-1)){
      epsilon.truth[,t]=rmvnorm(1, epsilon.truth[,t-1],
                                Sigma_phi.truth)
    }
    
    kappa.truth=exp(eta.truth)/(1+exp(eta.truth))
    phi.truth=exp(epsilon.truth)/(1+exp(epsilon.truth))
  }
  
  z=lapply(1:S, matrix, data=NA, nrow=n, ncol=T)
  v=lapply(1:S, matrix, data=NA, nrow=n, ncol=T)
  
  ##
  ## t=1
  ##
  
  t=1
  for(s in 1:S){
    z[[s]][,1]=rep(1, n)
    v[[s]][,1]=rbinom(n, z[[s]][,1],
                   1-kappa.truth[s,1])
  }
  
  ##
  ## t=2:T
  ##
  for(t in 2:T){
    for(s in 1:S){
      z[[s]][,t]=rbinom(n,size=v[[s]][,t-1],
                     prob=phi.truth[s,t-1])
      v[[s]][,t]=rbinom(n,z[[s]][,t],
                     1-kappa.truth[s,t])
    }
  }
}

###
### MCMC settings
###

n.iter=1000000

###
### Data
###
T = n.years
S = n.strat

###
### Priors
###

##
## s2_kappa ~ InvGamma(s2_kappa.prior1,s2_kappa.prior2)
##

s2_kappa.prior1=0.01
s2_kappa.prior2=0.01

##
## s2_phi ~ InvGamma(s2_phi.prior1,s2_phi.prior2)
##
s2_phi.prior1=0.01
s2_phi.prior2=0.01

###
### Starting values
###

if(!realdata){
  rho=rho.truth
  s2_kappa=s2_kappa.truth
  s2_phi=s2_phi.truth
  Sigma_kappa=s2_kappa*solve(D-rho*W)
  Sigma_kappa.bd=Sigma_kappa.bd.truth
  Sigma_phi=s2_phi*solve(D-rho*W)
  Sigma_phi.bd=Sigma_phi.bd.truth
  eta=eta.truth
  epsilon=epsilon.truth
  kappa=exp(eta)/(1+exp(eta))
  phi=exp(epsilon)/(1+exp(epsilon))
} else{
  rho=0.99
  s2_kappa=0.1
  s2_phi=0.1
  Sigma_kappa=s2_kappa*solve(D-rho*W)
  Sigma_kappa.bd=kronecker(diag(T-1),Sigma_kappa)
  Sigma_phi=s2_phi*solve(D-rho*W)
  Sigma_phi.bd=kronecker(diag(T-2),Sigma_phi)
  eta=matrix(rnorm(S*T), S, T)
  epsilon=matrix(rnorm(S*(T-1)), S, T-1)
  kappa=exp(eta)/(1+exp(eta))
  phi=exp(epsilon)/(1+exp(epsilon))
}

###
### book-keeping
###

rho.save=matrix(NA,n.iter,1)
epsilon.save=array(NA,dim=c(n.iter,S,T-1))
eta.save=array(NA,dim=c(n.iter,S,T))
s2_kappa.save=matrix(NA,n.iter,1)
s2_phi.save=matrix(NA,n.iter,1)

###
### tuning parameters
###

rho.tune=0.0001
accept.rho=0
epsilon.tune=rep(0.0001,T-1)
accept.epsilon=rep(0,T)
eta.tune=rep(0.0001,T)
accept.eta=rep(0,T)
s2_kappa.tune=0.0001
accept.s2_kappa=0
s2_phi.tune=0.0001
accept.s2_phi=0

if(pretuned){
  load("tuning.RData")
  rho.tune=tuning$rho.tune
  epsilon.tune=tuning$epsilon.tune
  eta.tune=tuning$eta.tune
  s2_kappa.tune=tuning$s2_kappa.tune
  s2_phi.tune=tuning$s2_phi.tune
}

###
### Gibbs loop
###

for(k in 1:n.iter){
  if(k%%1000==0)print(k)
  
  ##
  ## Sample rho
  ##

  rho.star=rnorm(1,rho,rho.tune)
  if(rho.star>0&rho.star<1){
    Sigma_kappa.star=s2_kappa*solve(D-rho.star*W)
    Sigma_phi.star=s2_phi*solve(D-rho.star*W)
    Sigma_kappa.bd.star=kronecker(diag(T-1),Sigma_kappa.star)
    Sigma_phi.bd.star=kronecker(diag(T-2),Sigma_phi.star)
    mh1=sum(dmvnorm(c(eta[,-1]),c(eta[,-T]),Sigma_kappa.bd.star,log=TRUE))+
      sum(dmvnorm(c(epsilon[,-1]),c(epsilon[,-((T-1):T)]),Sigma_phi.bd.star,log=TRUE))
    mh2=sum(dmvnorm(c(eta[,-1]),c(eta[,-T]),Sigma_kappa.bd,log=TRUE))+
      sum(dmvnorm(c(epsilon[,-1]),c(epsilon[,-((T-1):T)]),Sigma_phi.bd,log=TRUE))
    mh=exp(mh1-mh2)
    if(mh>runif(1) & !is.na(mh)){
      rho=rho.star
      Sigma_kappa=Sigma_kappa.star
      Sigma_kappa.bd=Sigma_kappa.bd.star
      Sigma_phi=Sigma_phi.star
      Sigma_phi.bd=Sigma_phi.bd.star
      accept.rho=accept.rho+1
    }
  }
  
  ##
  ## Sample epsilon t=1
  ##
  
  epsilon.star=rnorm(S,epsilon[,1],epsilon.tune[1])
  phi.star=exp(epsilon.star)/(1+exp(epsilon.star))
  mh1=dmvnorm(epsilon[,2],epsilon.star,Sigma_phi,log=TRUE)+
    sum(sapply(1:S, ratesLike, prob=phi.star, type='phi', z=z, v=v, t=1), na.rm=T) +
    sum(dnorm(epsilon.star,0,1.5,log=TRUE))
  mh2=dmvnorm(epsilon[,2],epsilon[,1],Sigma_phi,log=TRUE)+
    sum(sapply(1:S, ratesLike, prob=phi[,1], type='phi', z=z, v=v, t=1), na.rm=T) +
    sum(dnorm(epsilon[,1],0,1.5,log=TRUE))
  mh=exp(mh1-mh2)
  if(mh>runif(1) & !is.na(mh)){
    epsilon[,1]=epsilon.star
    phi[,1]=phi.star
    accept.epsilon[1]=accept.epsilon[1]+1
  }
  
  
  ##
  ## Sample epsilon t=2,...,T-2
  ##
  
  for(t in 2:(T-2)){
    epsilon.star=rnorm(S,epsilon[,t],epsilon.tune[t])
    phi.star=exp(epsilon.star)/(1+exp(epsilon.star))
    mh1=dmvnorm(epsilon[,t+1],epsilon.star,Sigma_phi,log=TRUE)+
      dmvnorm(epsilon.star,epsilon[,t-1],Sigma_phi,log=TRUE)+
      sum(sapply(1:S, ratesLike, prob=phi.star, type='phi', z=z, v=v, t=t), na.rm=T) 
    mh2=dmvnorm(epsilon[,t+1],epsilon[,t],Sigma_phi,log=TRUE)+
      dmvnorm(epsilon[,t],epsilon[,t-1],Sigma_phi,log=TRUE)+
      sum(sapply(1:S, ratesLike, prob=phi[,t], type='phi', z=z, v=v, t=t), na.rm=T) 
    mh=exp(mh1-mh2)
    if(mh>runif(1) & !is.na(mh)){
      epsilon[,t]=epsilon.star
      phi[,t]=phi.star
      accept.epsilon[t]=accept.epsilon[t]+1
    }
  }
  
  ##
  ## Sample epsilon t=T-1
  ##
  
  epsilon.star = rnorm(S, epsilon[,T-1], epsilon.tune[T-1])
  phi.star = exp(epsilon.star) / (1 + exp(epsilon.star))

  mh1 = dmvnorm(epsilon.star, epsilon[,T-2], Sigma_phi, log=TRUE) +
    sum(sapply(1:S, ratesLike, prob=phi.star, type='phi', z=z, v=v, t=T-1), na.rm=T)

  mh2 = dmvnorm(epsilon[,T-1], epsilon[,T-2], Sigma_phi, log=TRUE) +
    sum(sapply(1:S, ratesLike, prob=phi[,T-1], type='phi', z=z, v=v, t=T-1), na.rm=T)

  mh=exp(mh1-mh2)
  if(mh > runif(1)){
    epsilon[,T-1] = epsilon.star
    phi[,T-1] = phi.star
    accept.epsilon[T-1] = accept.epsilon[T-1]+1
  }
  
  ##
  ## Sample eta t=1
  ##
  
  eta.star=rnorm(S,eta[,1],eta.tune[1])
  kappa.star=exp(eta.star)/(1+exp(eta.star))
  mh1=dmvnorm(eta[,2],eta.star,Sigma_kappa,log=TRUE)+
    sum(sapply(1:S, ratesLike, prob=kappa.star, type='kappa', z=z, v=v, t=1), na.rm=T) +
    sum(dnorm(eta.star,0,1.5,log=TRUE))
  mh2=dmvnorm(eta[,2],eta[,1],Sigma_kappa,log=TRUE)+
    sum(sapply(1:S, ratesLike, prob=kappa[,1], type='kappa', z=z, v=v, t=1), na.rm=T) +
    sum(dnorm(eta[,1],0,1.5,log=TRUE))
  mh=exp(mh1-mh2)
  if(mh>runif(1) & !is.na(mh)){
    eta[,1]=eta.star
    kappa[,1]=kappa.star
    accept.eta[1]=accept.eta[1]+1
  }


  ##
  ## Sample eta t=2,...,T-1
  ##

  for(t in 2:(T-1)){
    eta.star=rnorm(S,eta[,t],eta.tune[t])
    kappa.star=exp(eta.star)/(1+exp(eta.star))
    mh1=dmvnorm(eta[,t+1],eta.star,Sigma_kappa,log=TRUE)+
      dmvnorm(eta.star,eta[,t-1],Sigma_kappa,log=TRUE)+
      sum(sapply(1:S, ratesLike, prob=kappa.star, type='kappa', z=z, v=v, t=t), na.rm=T)
    mh2=dmvnorm(eta[,t+1],eta[,t],Sigma_kappa,log=TRUE)+
      dmvnorm(eta[,t],eta[,t-1],Sigma_kappa,log=TRUE)+
      sum(sapply(1:S, ratesLike, prob=kappa[,t], type='kappa', z=z, v=v, t=t), na.rm=T)
    mh=exp(mh1-mh2)
    if(mh>runif(1) & !is.na(mh)){
      eta[,t]=eta.star
      kappa[,t]=kappa.star
      accept.eta[t]=accept.eta[t]+1
    }
  }

  ##
  ## Sample eta t=T
  ##

  eta.star=rnorm(S,eta[,T],eta.tune[T])
  kappa.star=exp(eta.star)/(1+exp(eta.star))
  mh1=dmvnorm(eta.star,eta[,T-1],Sigma_kappa,log=TRUE)+
    sum(sapply(1:S, ratesLike, prob=kappa.star, type='kappa', z=z, v=v, t=T), na.rm=T)
  mh2=dmvnorm(eta[,T],eta[,T-1],Sigma_kappa,log=TRUE)+
    sum(sapply(1:S, ratesLike, prob=kappa[,T], type='kappa', z=z, v=v, t=T), na.rm=T)
  mh=exp(mh1-mh2)
  if(mh>runif(1) & !is.na(mh)){
    eta[,T]=eta.star
    kappa[,T]=kappa.star
    accept.eta[T]=accept.eta[T]+1
  }


  ##
  ## Sample s2_phi
  ##

  s2_phi.star=rnorm(1,s2_phi,s2_phi.tune)
  if(s2_phi.star>0){
    Sigma_phi.star=s2_phi.star*solve(D-rho*W)
    Sigma_phi.bd.star=kronecker(diag(T-2),Sigma_phi.star)
    mh1=sum(dmvnorm(c(epsilon[,-1]),c(epsilon[,-((T-1):T)]),
                    Sigma_phi.bd.star,log=TRUE)) +
      dinvgamma(s2_phi.star,shape=s2_phi.prior1,rate=s2_phi.prior2,log=TRUE)
    mh2=sum(dmvnorm(c(epsilon[,-1]),c(epsilon[,-((T-1):T)]),
                    Sigma_phi.bd,log=TRUE))+
      dinvgamma(s2_phi,shape=s2_phi.prior1,rate=s2_phi.prior2,log=TRUE)
    mh=exp(mh1-mh2)
    if(mh>runif(1) & !is.na(mh)){
      s2_phi=s2_phi.star
      Sigma_phi=Sigma_phi.star
      Sigma_phi.bd=Sigma_phi.bd.star
      accept.s2_phi=accept.s2_phi+1
    }
  }

  ##
  ## Sample s2_kappa
  ##

  s2_kappa.star=rnorm(1,s2_kappa,s2_kappa.tune)
  if(s2_kappa.star>0){
    Sigma_kappa.star=s2_kappa.star*solve(D-rho*W)
    Sigma_kappa.bd.star=kronecker(diag(T-1),s2_kappa.star*solve(D-rho*W))
    mh1=sum(dmvnorm(c(eta[,-1]),c(eta[,-T]),Sigma_kappa.bd.star,log=TRUE))+
      dinvgamma(s2_kappa.star,s2_kappa.prior1,s2_kappa.prior2,log=TRUE)
    mh2=sum(dmvnorm(c(eta[,-1]),c(eta[,-T]),Sigma_kappa.bd,log=TRUE))+
      dinvgamma(s2_kappa,s2_kappa.prior1,s2_kappa.prior2,log=TRUE)
    mh=exp(mh1-mh2)
    if(mh>runif(1) & !is.na(mh)){
      s2_kappa=s2_kappa.star
      Sigma_kappa=Sigma_kappa.star
      Sigma_kappa.bd=Sigma_kappa.bd.star
      accept.s2_kappa=accept.s2_kappa+1
    }
  }

  ##
  ## Save output
  ##

  rho.save[k,]=rho
  epsilon.save[k,,]=epsilon
  eta.save[k,,]=eta
  s2_kappa.save[k,]=s2_kappa
  s2_phi.save[k,]=s2_phi
  
  ##
  ## Autotune
  ##
  
  if(accept.rho/k<0.3)rho.tune=rho.tune*0.9
  if(accept.rho/k>0.5)rho.tune=rho.tune*1.1
  epsilon.tune=ifelse(accept.epsilon/k<0.3,
                      epsilon.tune*0.9,
                      epsilon.tune)
  epsilon.tune=ifelse(accept.epsilon/k>0.5,
                      epsilon.tune*1.1,
                      epsilon.tune)
  eta.tune=ifelse(accept.eta/k<0.3,
                  eta.tune*0.9,
                  eta.tune)
  eta.tune=ifelse(accept.eta/k>0.5,
                  eta.tune*1.1,
                  eta.tune)
  if(accept.s2_kappa/k<0.3)s2_kappa.tune=s2_kappa.tune*0.9
  if(accept.s2_kappa/k>0.5)s2_kappa.tune=s2_kappa.tune*1.1
  if(accept.s2_phi/k<0.3)s2_phi.tune=s2_phi.tune*0.9
  if(accept.s2_phi/k>0.5)s2_phi.tune=s2_phi.tune*1.1
  
}

save.image('~/Dropbox/eps_all.RData')

###
### Save tuning parameters
###

tuning=list()
tuning$rho.tune=rho.tune
tuning$epsilon.tune=epsilon.tune
tuning$eta.tune=eta.tune
tuning$s2_kappa.tune=s2_kappa.tune
tuning$s2_phi.tune=s2_phi.tune
save(tuning,file="tuning.RData")

###
### Tuning
###

accept.rho/n.iter
accept.epsilon/n.iter
accept.eta/n.iter
accept.s2_kappa/n.iter
accept.s2_phi/n.iter


###
### Trace plots
###

burn=250000
thin=15
ind=seq(burn,n.iter,thin)

##
## rho
##

plot(rho.save[ind],type='l')
if(!realdata)abline(h=rho.truth,col=2)


##
## epsilon
##

t=1
i=1
bad.eps <- list()
for(t in 1:(T-1)){
  for(i in 1:(S)){
    plot(epsilon.save[ind,i,t],type='l',ylim=c(-5,5), main=paste(""))
    if(!realdata)abline(h=epsilon.truth[i,t],col=2)
    Sys.sleep(0.2)
  }
  q=apply(epsilon.save[ind,,t],2,quantile,c(0.025,0.975),na.rm=TRUE)
  lq=q[1,]
  uq=q[2,]
  inq=sum(ifelse(epsilon.truth[,t]>=lq&epsilon.truth[,t]<=uq,1,0))
  bad.eps[[t]] <- which(epsilon.truth[,t]<=lq|epsilon.truth[,t]>=uq)
  print((prop.in=inq/S))
}

bad.eps

##
## eta,
##

bad.eta <- list()
for(t in 1:T){
  # for(i in 1:(S)){
  #   plot(eta.save[ind,i,t],type='l', main=paste("eta: t = ", t, ", s = ",))
  #   if(!realdata)abline(h=eta.truth[i,t],col=2)
  #   Sys.sleep(0.2)
  # }
  q=apply(eta.save[ind,,t],2,quantile,c(0.025,0.975),na.rm=TRUE)
  lq=q[1,]
  uq=q[2,]
  inq=sum(ifelse(eta.truth[,t]>=lq&eta.truth[,t]<=uq,1,0))
  bad.eta[[t]] <- which(eta.truth[,t]<=lq | eta.truth[,t]>=uq)
  print((prop.in=inq/S))
}

bad.eta

##
## s2_kappa
##

plot(s2_kappa.save[ind],type='l', main="sigma kappa")
if(!realdata)abline(h=s2_kappa.truth,col=2)

##
## s2_phi
##

plot(s2_phi.save[ind],type='l', main='sigma phi')
if(!realdata)abline(h=s2_phi.truth,col=2)



##
## wonky sites
##

## eta
bad.eta

eta.truth[unique(bad.eta[[1]]),1]
mean(eta.truth[,1])

eta.truth[unique(bad.eta[[2]]),2]
mean(eta.truth[,2])

eta.truth[unique(bad.eta[[3]]),3]
mean(eta.truth[,3])

eta.truth[unique(bad.eta[[4]]),4]
mean(eta.truth[,4])


### epsilon
bad.eps

epsilon.truth[unique(bad.eps[[1]]),1]
mean(epsilon.truth[,1])

epsilon.truth[unique(bad.eps[[2]]),2]
mean(epsilon.truth[,2])

epsilon.truth[unique(bad.eps[[3]]),3]
mean(epsilon.truth[,3])

epsilon.truth[unique(bad.eps[[4]]),4]
mean(epsilon.truth[,4])
