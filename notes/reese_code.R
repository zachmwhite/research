# Basic GP code.  Stolen shamelessly form Dr. Reese at BYU. https://madison.byu.edu/bayes/OurGP.R
library(MASS)

####################
nreal<-50
np<-1000
phi10<-0.5
xx<-seq(0,1,length=np)
dd<-dist(xx,upper=T,diag=T)^2
Sig<-exp(-phi10*as.matrix(dd))
gpreal<-mvrnorm(nreal,rep(0,np),Sig)
matplot(t(gpreal),type="l",main=expression(paste(phi,'=',phi10)))

########################
mysimfunc<-function(x){
  out<-16*x^2*(1-x)^2
  out[x>=0.7]<-  out[x>=0.7] +32*(x[x>=0.7]-0.7)^3
  return(out)
}
nrep<-1000
ndat<-6
np<-100
mydatx<-runif(ndat)
#mydatx<-seq(0,1,length=ndat)
mydaty<-mysimfunc(mydatx)+rnorm(length(mydatx),0,sqrt(0.01))

phi<-10
xx<-seq(0,1,length=np)
dd<-dist(xx,upper=T,diag=T)
Sig<-exp(-phi*as.matrix(dd^2))
yy<-mvrnorm(20,rep(0,np),Sig)
matplot(t(yy),type="l")  
Sigxx<-exp(-phi*as.matrix(dist(mydatx,upper=T,diag=T)^2))
mypostmean<-rep(0,np)
signoise<-0.01
smallk1<-matrix(0,ncol=np,nrow=ndat)
smallk2<-matrix(0,ncol=ndat,nrow=np)
for(i in 1:np){
  smallk1[,i]<-exp(-phi*(xx[i]-mydatx)^2)
}
for(i in 1:ndat){
  smallk2[,i]<-exp(-phi*(mydatx[i]-xx)^2)
}
covinv<-solve(Sigxx+diag(signoise,ndat))
mumat<-t(smallk1)%*%covinv%*%mydaty
sigmat<-Sig-t(smallk1)%*%covinv%*%t(smallk2)

yfit<-mvrnorm(nrep,mumat,sigmat)
meany<-apply(yfit,2,mean)
plot(mydatx,mydaty,pch=1,cex=2.5,xlim=c(0,1),ylim=c(0,2))
for(i in 1:nrep){
  lines(xx,yfit[i,],col="lightgray")
}
lines(xx,meany,col="blue",lwd=3)
points(mydatx,mydaty,pch=1,cex=1.8)

#######################
mysimfunc<-function(x){
  out<-16*x^2*(1-x)^2
  out[x>=0.7]<-  out[x>=0.7] +32*(x[x>=0.7]-0.7)^3
  out<-out+rnorm(length(x),0,.01)
  return(out)
}
#phi<-10
np<-50
xx<-seq(0,1,length=np)
dd<-dist(xx,upper=T,diag=T)
#Sig<-exp(-phi[i-1]*as.matrix(dd^2))
#yy<-mvrnorm(20,rep(0,np),Sig)
#matplot(t(yy),type="l")  
nrep<-1000
ndat<-6
mydatx<-runif(ndat)
mydaty<-mysimfunc(mydatx)
mypostmean<-rep(0,np)
smallk1<-matrix(0,ncol=np,nrow=ndat)
smallk2<-matrix(0,ncol=ndat,nrow=np)
niter<-100
signoise<-rep(0.000001,niter)
phi<-rep(10,niter)
ar.phi<-ar.sig<-0
cs.phi<-0.001
cs.sig<-0.001
for(i in 2:niter){
  Sigxx<-exp(-phi[i-1]*as.matrix(dist(mydatx,upper=T,diag=T)^2))
  
  Sig<-exp(-phi[i-1]*as.matrix(dd^2))
  
  for(j in 1:np){
    smallk1[,j]<-exp(-phi[i-1]*(xx[j]-mydatx)^2)
  }
  for(j in 1:ndat){
    smallk2[,j]<-exp(-phi[i-1]*(mydatx[j]-xx)^2)
  }
  
  covinv<-solve(Sigxx+diag(signoise[i-1],ndat))
  mumat<-t(smallk1)%*%covinv%*%mydaty
  sigmat<-Sig-t(smallk1)%*%covinv%*%t(smallk2)
  yfit<-mvrnorm(1,mumat,sigmat)
  
  #METROPOLIS HASTINGS FOR PHI
  phi[i]<-phi[i-1]
  llo<--np/2*det(sigmat,log=T)-0.5*t(yfit-mumat)%*%solve(sigmat)%*%(yfit-mumat)+dgamma(phi[i],10,1,log=T)
  cand<-rnorm(1,phi[i-1],cs.phi)
  
  Sigxx<-exp(-cand*as.matrix(dist(mydatx,upper=T,diag=T)^2))
  
  Sig<-exp(-cand*as.matrix(dd^2))
  
  for(i in 1:np){
    smallk1[,i]<-exp(-cand*(xx[i]-mydatx)^2)
  }
  for(i in 1:ndat){
    smallk2[,i]<-exp(-cand*(mydatx[i]-xx)^2)
  }
  
  covinv<-solve(Sigxx+diag(signoise[i-1],ndat))
  mumat<-t(smallk1)%*%covinv%*%mydaty
  sigmat<-Sig-t(smallk1)%*%covinv%*%t(smallk2)
  lln<--np/2*det(sigmat,log=T)-0.5*t(yfit-mumat)%*%solve(sigmat)%*%(yfit-mumat)+dgamma(cand,10,1,log=T)
  
  if((lln-llo)>log(runif(1))){phi[i]<-cand;ar.phi<-ar.phi+1}
  #METROPOLIS HASTINGS FOR PHI
  signoise[i]<-signoise[i-1]
  llo<--np/2*det(sigmat,log=T)-0.5*t(yfit-mumat)%*%solve(sigmat)%*%(yfit-mumat)+dgamma(phi[i],10,1,log=T)
  cand<-rnorm(1,signoise[i-1],cs.sig)
  Sigxx<-exp(-phi[i]*as.matrix(dist(mydatx,upper=T,diag=T)^2))
  
  Sig<-exp(-phi[i]*as.matrix(dd^2))
  
  for(i in 1:np){
    smallk1[,i]<-exp(-phi[i]*(xx[i]-mydatx)^2)
  }
  for(i in 1:ndat){
    smallk2[,i]<-exp(-phi[i]*(mydatx[i]-xx)^2)
  }
  
  covinv<-solve(Sigxx+diag(signoise[i],ndat))
  mumat<-t(smallk1)%*%covinv%*%mydaty
  sigmat<-Sig-t(smallk1)%*%covinv%*%t(smallk2)
  lln<--np/2*det(sigmat,log=T)-0.5*t(yfit-mumat)%*%solve(sigmat)%*%(yfit-mumat)+dgamma(cand,10,1,log=T)
  
  if((lln-llo)>log(runif(1))){signoise[i]<-cand;ar.sig<-ar.sig+1}
  
}
meany<-apply(yfit,2,mean)
plot(mydatx,mydaty,pch=1,cex=2.5,xlim=c(0,1),ylim=c(0,2))
for(i in 1:nrep){
  lines(xx,yfit[i,],col="lightgray")
}
lines(xx,meany,col="blue",lwd=3)