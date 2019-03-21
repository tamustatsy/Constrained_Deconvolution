#setwd("~/Research/unimodal_density_measurement_error/Bayesian_codes/Kernel_deconvolution_as_comparison/")

#doOne <- function(n, df.t, sigU, lfg){
  

#-----------------------------------------------------
#Start by generating some data contaminated by noise:
#-----------------------------------------------------

# set.seed(13)
# 
# #Noise to signal ratio=varU/varX
# #NSR=1/3
# 
# #Sample size
# n=5000
# 
# #Generate data from a normal mixture
# df.t <- 5
# sigU <- 0.6
# 
# xs <- rt(n, df = df.t)
# q.zero <- 0.8
# pp <- runif(n)
# c_zero <- sum(pp < q.zero)
# sd.zero <- 0.2
# xs[pp < q.zero] <- rnorm(c_zero, mean = 0, sd = sd.zero)
# 
# u <- rnorm(n, mean = 0, sd = sigU)              
# W <- xs + u
# 


kernel_homo <- function(n, W, sigU, xs, lfg, errortype = 'norm'){

  #Grid where to estimate the true mixture density, and calculation of true density
  xx = seq(-max(xs),max(xs),length=lfg)
#  dx=xx[2]-xx[1]

#Plot the true density
#plot(xx,q.zero*dnorm(xx, mean = 0, sd = sd.zero)+(1 - q.zero)*dt(xx, df = df.t),'l',
#     xlim = c(-3,3), col='red',xlab="",ylab="")

#PI bandwidth of Delaigle and Gijbels
hPI=PI_deconvUknownth4(n,W,errortype,sigU);

#DKDE estimator without rescaling (density does not integrate exactly to 1)
y=fdecUknown(n,xx,W,hPI,errortype,sigU);

#DKDE estimator with rescaling: here xx must be equispaced and must cover the range where the estimated density is significantly non zero
y2=fdecUknown(n,xx,W,hPI,errortype,sigU,rescale=1);

#lines(xx,y2,col="green",xlab="",ylab="")
#lines(xx,y,col='black')

#legend(x="topright",legend=c( "true f","fdec, hPI", "fdec rescaled, hPI"),col=c("red","black","green"),lty=c(1,1,1,1), cex=0.55)


##CV bandwidth of Stefanski and Carroll
#hCV=CVdeconv(n,W,errortype,sigU)

##DKDE estimator without rescaling (density does not integrate exactly to 1)
#y3=fdecUknown(n,xx,W,hCV,errortype,sigU);

#lines(xx,y3,col='magenta')



#Compare with the naive KDE estimator that ignores the error (using normal reference bandwidth and standard normal kernel)
#h=1.06*sqrt(var(W))*n^(-1/5);
#xout=outerop(xx,t(W),"-");

#fnaive=apply(dnorm(xout,0,h),1,sum)/n;

#lines(xx,fnaive,col='cyan')

#-------------------------------------
#Compute the IMSE of the estimator  
#-------------------------------------

# imse_y <- sqrt(sum((y - dt(xx, df = df.t))^2*dx))
# imse_y2 <- sqrt(sum((y2 - dt(xx, df = df.t))^2*dx))
# imae_y <- sum(abs(y - dt(xx, df = df.t))*dx)
# imae_y2 <- sum(abs(y2 - dt(xx, df = df.t))*dx)
#imse_y3 <- sqrt(sum((y3 - dt(xx, df = df.t))^2*dx))

result <- list("y" = y, "y2" = y2, "xx" = xx)
return(result)

#-------------------------------------
#Example when the error is Laplace                              #from here to the end, to use the code, some modification is needed
#-------------------------------------
# windows()
# errortype="Lap"
# sigLap=sqrt(NSR*var(X)/2)
# sigU=sqrt(2)*sigLap;
# U=rlap(sigLap,1,n);
# 
# #Contaminated data
# W=as.vector(X+U);
# 
# 
# #Plot the true density
# plot(xx,truedens,'l',col='red',xlab="",ylab="")
# 
# #PI bandwidth of Delaigle and Gijbels
# hPI=PI_deconvUknownth4(n,W,errortype,sigU);
# 
# 
# #DKDE estimator without rescaling (density does not integrate exactly to 1)
# y=fdecUknown(n,xx,W,hPI,errortype,sigU);
# 
# #DKDE estimator with rescaling: here xx must be equispaced and must cover the range where the estimated density is significantly non zero
# y2=fdecUknown(n,xx,W,hPI,errortype,sigU,rescale=1);
# 
# lines(xx,y2,col="green",xlab="",ylab="")
# lines(xx,y,col='black')
# 
# 
# 
# #CV bandwidth of Stefanski and Carroll
# hCV=CVdeconv(n,W,errortype,sigU)
# 
# #DKDE estimator without rescaling (density does not integrate exactly to 1)
# y3=fdecUknown(n,xx,W,hCV,errortype,sigU);
# 
# lines(xx,y3,col='magenta')
# 
# 
# 
# #Compare with the naive KDE estimator that ignores the error (using normal reference bandwidth and standard normal kernel)
# h=1.06*sqrt(var(W))*n^(-1/5);
# xout=outerop(xx,t(W),"-");
# 
# fnaive=apply(dnorm(xout,0,h),1,sum)/n;
# 
# lines(xx,fnaive,col='cyan')
# 
# 
# legend(x="topright",legend=c( "true f","fdec, hPI", "fdec rescaled, hPI", "fdec rescaled, hCV", "naive estimator, hNR"),col=c("red","black","green","magenta","cyan"),lty=c(1,1,1,1),cex=0.73)

}#end of kernel_homo function

#}#end of the doOne function

