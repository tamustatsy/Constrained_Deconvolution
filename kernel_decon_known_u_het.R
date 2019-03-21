####This is a example for R package "fDKDEheterosc"
#library("fDKDEheterosc")

#doOne <- function(n, df.t, lfg){



#Author: Aurore Delaigle
#This code illustrates how to use the functions for computing the deconvolution kernel density estimator and its bandwidths in the case where the errors are heteroscedastic

#-----------------------------------------------------
#Start by generating some data contaminated by noise:
#-----------------------------------------------------

# xs <- rt(n, df = df.t)
# 
# 
# 
# 
# #Standard devation of the jth error, for j=1,...,n
# sigmaj <- sqrt((1 + xs/4)^2)
# u <- rnorm(n, mean = 0, sd = sigmaj) 

#Contaminated data

#W = xs + u

kernel_het <- function(n, W, sigmaj, xs, lfg, errortype = 'norm'){

  #Grid where to estimate the true mixture density, and calculation of true density
  xx = seq(-max(xs),max(xs),length=lfg)
  #dx=xx[2]-xx[1];

#Estimate the variance of X
varX=mean(W^2)-(mean(W))^2-sum(sigmaj^2)/n;

#PI bandwidth of Delaigle and Gijbels
hPI=PI_deconvUknownth4het(n,W,varX,errortype,sigmaj);

#SIMEX bandwidth of Delaigle 
#hSIMEX=hSIMEXUknown();

#DKDE estimator without rescaling (density does not integrate exactly to 1)
y=fdecUknownhet(n,xx,W,hPI,errortype,sigmaj);

#DKDE estimator with rescaling: here xx must be equispaced and must cover the range where the estimated density is significantly non zero
y2=fdecUknownhet(n,xx,W,hPI,errortype,sigmaj,rescale=1);

#Plot the true density
# plot(xx, dt(xx, df = df.t),'l',col='red',xlab="",ylab="")
# lines(xx,y,col='black')
# lines(xx,y2,col="green")
# 
# #Example of how to provide the vector of phiU_k's instead of the error type and the standard deviations
# 
# phiUkvec=c()
# for(k in 1:n)
# {	
# 	phiUk<-function(tt,k) {return(exp(-sigmaj[k]^2*tt^2/2));}
# 	phiUkvec=c(phiUkvec,phiUk)
# }
# 
# #DKDE estimator without rescaling (density does not integrate exactly to 1)
# y3=fdecUknownhet(n,xx,W,hPI,phiUkvec=phiUkvec);
# lines(xx,y3,col='magenta',lty=2)
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
# legend(x="topright",legend=c( "true f","fdec, hPI", "fdec rescaled, hPI", "fdec hPI v2", "naive estimator, hNR"),col=c("red","black","green","magenta","cyan"),lty=c(1,1,1,2,1),cex=0.73)

#-------------------------------------
#Compute the MISE of the estimator  
#-------------------------------------
# 
# imse_y <- sqrt(sum((y - dt(xx, df = df.t))^2*dx))
# imse_y2 <- sqrt(sum((y2 - dt(xx, df = df.t))^2*dx))
# imae_y <- sum(abs(y - dt(xx, df = df.t))*dx)
# imae_y2 <- sum(abs(y2 - dt(xx, df = df.t))*dx)
# 
# result <- c(imae_y, imse_y, imae_y2, imse_y2)

result <- list("y" = y, "y2" = y2, "xx" = xx)
return(result)

}#end of function kernel_het
#}#end of function doOne

