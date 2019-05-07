###########################################################################################################################
#The current file is built for batch runs. It can be run on its own by specifying lines 42-44 and comment out lines 39-40 # 
###########################################################################################################################
##################################
# Add source files and libraries #
##################################

if(!require('coda')) install.packages('coda')
library(coda)                 #required by MCMCpack
if(!require('MCMCpack')) install.packages('MCMCpack')
library(MCMCpack)
if(!require('evd')) install.packages('evd')
library(evd)                  #required by truncdist
if(!require('truncdist')) install.packages('truncdist')
library(truncdist)
if(!require('Rcpp')) install.packages('Rcpp')
library(Rcpp)
if(!require('RcppArmadillo')) install.packages('RcppArmadillo')
library(RcppArmadillo)
if(!require('VGAM')) install.packages('VGAM')
library(VGAM)                #required by Mysample.cpp
if(!require('deconvolve')) devtools::install_github("timothyhyndman/deconvolve")
library(deconvolve)
if(!require('batch')) install.packages('batch')
library(batch)               #required if conducting batch code for parallel
if(!require('transport')) install.packages('transport')
library(transport)               #required if conducting batch code for parallel

#--files containing functions for Constrained Bayes--#
source("../../../Myfunction.R")                                                  #this is used in computing the marginal density of X 
source("../../../ddsc_laplaceError.R")                                           #this contains the main function to conduct the hybrid Gibbs sampler
sourceCpp("../../../Mysample.cpp")                                               #this is used in updating the group indicator variables z
sourceCpp("../../../Mydist.cpp")                                                 #also used in updating the group indicator variables z

#####################
#Global parameters  #
#####################

parseCommandArgs()           #seed, n
set.seed(seed)

# seed <- 1100
# set.seed(seed)
# n <- 5000

#####################
# Hyper parameters  #
#####################

m <- 20.0                    #prior on mixing probabilities
K <- 8                       #number of mixture components
Xi_1 <- 1.0                  #Xi_1, Xi_2 are for prior on rate of theta
Xi_2 <- 4.0
lambda <- 2                  #prior on shape of theta
parMH <- 2                   #proposal distribution, in the paper, it is fixed at 2 and not treated as hyperparameter
tt <- 2.5                    #truncated value for alpha_k and its proposal, the larger it is, the smoother
                             #(or flatter) the peak is   

########################################
# Parameters that are used in the MCMC #
########################################

n.burnin <- 1000
n.MCMC <- 5000

#--length of grid to plot density--#
lfg <- 200                                 #indicates length for grid
x.grid <- seq(0, 3, length = lfg)


######################################################
# Generate data x using a sym and unimodal density   #
# and proxy variable w by adding a normal error to x #
######################################################

df.t <- 5
x <- rt(n, df = df.t)
q.zero <- 0.8
pp <- runif(n)
c_zero <- sum(pp < q.zero)
sd.zero <- 0.2
x[pp < q.zero] <- rnorm(c_zero, mean = 0, sd = sd.zero)
b_lap <- 0.6/sqrt(2)                                 #2b^2 = 0.6^2, the variance of laplace is same as variance of xs
u <- rlaplace(n, location = 0, scale = b_lap)              
w <- x + u

##################
# Initialization #
##################

beta <- rgamma(K, shape = Xi_1, rate = Xi_2)
alpha <- rtrunc(K, "exp", a = tt, b = Inf, rate = lambda)
p <- rdirichlet(1, rep(m/K, K))
z <- sample(1:K, size = n, prob = p, replace = TRUE)
theta <- numeric(n)
for(ii in 1:n){
  theta[ii] <- rgamma(1, shape = alpha[z[ii]], rate = beta[z[ii]])
}

###################################################################
# Two density deconvolution estimators are obtained below, one    #
# using our Constrained Bayes method and one using the Kernel     #
# density deconvolution method (Delaigle and Meister, 2008).      #
# Code for the Kernel deconvolution estimator is available at     # 
# http://researchers.ms.unimelb.edu.au/%7Eaurored/links.html#Code #
###################################################################
density.CB <- ddsc_mcmc(w, b_lap, n.burnin, n.MCMC, n, K, m, lambda, parMH, tt,  
                        Xi_1, Xi_2, z, theta, beta, alpha, p, x.grid)

res.kernel <- deconvolve(W1 = w, xx = seq(-max(x.grid),max(x.grid),length=2*lfg), 
                         errortype = 'laplace', sd_U = b_lap*sqrt(2))
density.Ker <- res.kernel$pdf
xx.Ker <- res.kernel$x

#--compute integrated absolute error and square root of integrated squared error--#
density.true <- q.zero*dnorm(x.grid, mean = 0, sd = sd.zero) + (1-q.zero)*dt(x.grid, df = df.t)
imae_CB <- sum(abs(density.CB - density.true)*(x.grid[2]-x.grid[1]))*2        #factor 2 for the other half
imse_CB <- sqrt(sum((density.CB - density.true)^2*(x.grid[2]-x.grid[1])))*2
density_true <- q.zero*dnorm(xx.Ker, mean = 0, sd = sd.zero) + (1-q.zero)*dt(xx.Ker, df = df.t)
imae_Ker <- sum(abs(density.Ker - density_true)*(xx.Ker[2]-xx.Ker[1]))
imse_Ker <- sqrt(sum((density.Ker - density_true)^2*(xx.Ker[2]-xx.Ker[1])))
#--wasserstein distance--#
dsample_CB <- sample(x = c(rev(x.grid),x.grid[-1]), prob = c(rev(density.CB),density.CB[-1]), 
                     size = 1000000, replace = TRUE)
dsample_0 <- sample(x = c(rev(x.grid),x.grid[-1]), prob = c(rev(density.true),density.true[-1]), 
                    size = 1000000, replace = TRUE)
wass_CB <- wasserstein1d(dsample_CB, dsample_0, p = 2)
dsample_Ker <- sample(x = xx.Ker, prob = density.Ker, size = 1000000, replace = TRUE)
dsample_0 <- sample(x = xx.Ker, prob = density_true, size = 1000000, replace = TRUE)
wass_Ker <- wasserstein1d(dsample_Ker, dsample_0, p = 2)
#--absolute difference of exceedance--#
x.grid.cut <- x.grid > sd.zero*3
excprob_CB <- abs(sum(density.CB[x.grid.cut]-density.true[x.grid.cut])*(x.grid[2]-x.grid[1]))*2
xx.cut <- abs(xx.Ker) > 3*sd.zero
excprob_Ker <- abs(sum(density.Ker[xx.cut]-density_true[xx.cut])*(xx.Ker[2] - xx.Ker[1]))

#--write the density estimators and other parameters in a csv file--#
write.table(data.frame(seed, imae_CB, imse_CB, wass_CB, excprob_CB, imae_Ker, 
            imse_Ker, wass_Ker, excprob_Ker, x.grid, density.CB, xx.Ker, density.Ker), 
            file=paste("./result/", n, "/seed", seed, "_HOMO.csv", sep = ""), 
            append = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")









