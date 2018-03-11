library(MASS)

# Prior draws
nreal = 5
np = 1000
phi10 = .5
xx = seq(-5,5, length = np)
dd = dist(xx, upper = T, diag = T)^2
Sigma = exp(-phi10*as.matrix(dd))
gpreal = t(mvrnorm(nreal, rep(0,np), Sigma))

matplot(t(gpreal), type = "l", main = expression(paste(phi,'=', phi10)))
plot(xx,gpreal[,1], type = "l", ylim = c(-2,3))
for(i in 2:dim(gpreal)[2]){
  lines(xx,gpreal[,i], col = i)
}

#################################
# 
mysimfunc<-function(x){
  out<-3*sin(x)^3 
  return(out)
}
nrep = 1000
ndat = 6
np = 100
sigma2 = .1
nreal = 10
mydatx = runif(ndat,-5,5)
mydaty = mysimfunc(mydatx) + rnorm(ndat,0, sqrt(sigma2))

phi = 10
xx = seq(-5,5, length = np)
dd = dist(xx, upper = T, diag = T)
Sig = exp(-phi*as.matrix(dd^2))
# Prior Draws 
yy = t(mvrnorm(nreal, rep(0,np),Sig))
plot(xx,yy[,1],type = "l",ylim= c(-3,4))
for(i in 1:nreal){
  lines(xx,yy[,i], col = i)
}
# OR we can just do matplot
Sigxx = exp(-phi*as.matrix(dist(mydatx,upper = T,diag = T)^2))
mypostmean = rep(0,np)
signoise = .01
smallk1 = matrix(0, ncol = np, nrow = ndat)
smallk2 = matrix(0, ncol = ndat, nrow = np) # We could just do the transpose of the previous
for(i in 1:np){
  smallk1[,i] = exp(-phi * (xx[i] - mydatx)^2)
}
for( i in 1:ndat){
  smallk2[,i] = exp(-phi * (mydatx[i] - xx)^2)
}
covinv = solve(Sigxx + diag(signoise, ndat))
mumat = t(smallk1) %*% covinv %*% mydaty
sigmat = Sig - t(smallk1) %*% covinv %*% t(smallk2)

yfit = mvrnorm(nrep,mumat,sigmat)
meany = apply(yfit,2,mean)
plot(mydatx,mydaty,pch = 1, cex = 2.5, xlim = c(-5,5), ylim = c(-5,5))
for(i in 1:nrep){
  lines(xx, yfit[i,], col = "lightgray")
}
lines(xx,meany, col = "blue", lwd = 3)


#########################
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