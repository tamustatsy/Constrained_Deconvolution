######################################################################
# ddsc_mcmc is the main function to conduct the hybrid Gibbs sampler #
#--------------------------------------------------------------------#
# Input I: data and hyper parameters                                 #
# w: vector of proxy observations                                    #
# sd_u: vector of standard deviation(s) for the measurement error,   #
#       scalar when the error has homoscedastic variance             #
# n.burnin: number of burn-in steps                                  #
# n.MCMC: total number of MCMC steps                                 #
# n: sample size, equals the length of w                             #
# K: number of components in mixtures of Gamma distributions         #
# m: the concentration parameter in Dirichlet prior                  #
# lambda: the parameter in the prior of shape parameters of Gamma    #
# parMH: the parameter in the proposal distribution of shape         #
#        parameters of Gamma, default value 2                        #
# tt: the lower bound of shape parameters of Gamma                   #
# Xi_1: the 1st parameter in the prior of rate parameters of Gamma   #
# Xi_2: the 2nd parameter in the prior of rate parameters of Gamma   #
#--------------------------------------------------------------------#
# Input II: initial values for variables in the hierarchical prior   #
#           including z, theta, beta, alpha, p. These notations      #
#           agree with those in our paper                            #
#--------------------------------------------------------------------#
# Input III: the specs associated with the output                    #
# x.grid: the grid points for the density output, necessarily across #
# an interval [0, R] for some R > 0 due to symmetry of density       #
#--------------------------------------------------------------------#
# Output: the density estimator                                      #
# density.x: the estimated density on the specified grids, equals    #
#            the average of MCMC samples after burn-in               #
######################################################################

ddsc_mcmc <- function(w, sd_u, n.burnin, n.MCMC, n, K, m, lambda, parMH = 2,
                      tt, Xi_1, Xi_2, z, theta, beta, alpha, p, x.grid){

    lfg <- length(x.grid)
    density.x <- rep(0, lfg)

    ################
    # MCMC storage #
    ################
    d.ordinates <- matrix(0, nrow = n, ncol = K)
    n.kk.z <- numeric(K)
    s.kk.theta <- numeric(K)
    s.kk.ltheta <- numeric(K)
    mh.ratio <- numeric(K)
    log.mh.ratio <- numeric(K)
    alpha.tilde <- numeric(K)
    den_i <- numeric(lfg)

    ##############
    # Start MCMC #
    ##############
    for(iii in 1:n.MCMC){

    if(iii%%100==0)
        print(iii)

    #-------updating x_1, ..., x_n from a truncated normal distribution--------#
    ############################################################################
    #X_i = W_i - U_i                                                           #
    #1. generate U_i in a vectorized way:                                      #
    #the cutoff for U_i is [W_i - \theta_i, W_i + \theta_i]                    #
    #generate uniform variable between N(W_i - \theta_i) and N(W_i + \theta_i) #
    #get the quantiles corresponding to these uniform variables                #
    ############################################################################
    uu_cut_low <- w - theta
    uu_cut_up <- w + theta
    puu_cut_low <- pnorm(uu_cut_low, mean = 0, sd = sd_u)
    puu_cut_up <- pnorm(uu_cut_up, mean = 0, sd = sd_u)

    uu <- runif(n)
    uu_ext <- (puu_cut_up - puu_cut_low)*uu + puu_cut_low
    uu_ext[uu_ext == 1] <- 1 - 1E-5
    uu_ext[uu_ext == 0] <- 1E-5
    x <- w - qnorm(uu_ext, mean = 0, sd = sd_u)
    
    #--------------------------updating z_1, ..., z_n--------------------------#
    d.ordinates = mydgamma(theta, s=alpha, r=beta)             #see Mydist.cpp
    z <- mysample(1:K, d.ordinates, p, TRUE)                 #see Mysample.cpp

    #----updating theta_1, ..., theta_n from a truncated gamma distribution----#
    sg <- alpha[z] - 1
    rg <- beta[z]
    puu_cut_low <- pgamma(abs(x), shape = sg, rate = rg)
    uu <- runif(n)
    uu_ext <- (1 - puu_cut_low)*uu + puu_cut_low
    uu_ext[uu_ext == 1] <- 1 - 1E-5
    uu_ext[uu_ext == 0] <- 1E-5
    theta <- qgamma(uu_ext, shape = sg, rate = rg)

    ltheta <- log(theta)
    for(kk in 1:K){
    index <- (z == kk)
    n.kk.z[kk] <- sum(index)                          #represents group sum of z
    s.kk.theta[kk] <- sum(theta[index])
    s.kk.ltheta[kk] <- sum(ltheta[index])
    }

    #-------------------------updating p_1, ..., p_K---------------------------#
    p = rdirichlet(1, m/K + n.kk.z)

    #----------------------updating beta_1, ..., beta_K------------------------#
    beta <- rgamma(rep(1, K), shape = alpha*n.kk.z + Xi_1,
                    rate = s.kk.theta + Xi_2)

    #---------------------updating alpha_1, ..., alpha_K-----------------------#
    for(kk in 1:K){
        alpha.tilde[kk] <- rtrunc(1, "gamma", a = tt, b = Inf, shape = parMH,
                                  rate = parMH/alpha[kk])
        tp_proposed <- 1 - pgamma(tt, parMH , rate = parMH/alpha[kk])
        tp_current <- 1 - pgamma(tt, parMH , rate = parMH/alpha.tilde[kk])
        term_1 <- parMH *(alpha[kk]/alpha.tilde[kk] - alpha.tilde[kk]/alpha[kk])
        term_2 <- (alpha.tilde[kk] - alpha[kk])*(lambda - n.kk.z[kk]*log(beta[kk])
                  - s.kk.ltheta[kk])
        log.mh.ratio[kk] <- (2*parMH  - 1)*(log(alpha[kk]) - log(alpha.tilde[kk])) +
                            n.kk.z[kk]*(lgamma(alpha[kk]) - lgamma(alpha.tilde[kk])) -
                            term_1 - term_2 + log(tp_proposed)-log(tp_current)
        mh.ratio[kk] = exp(log.mh.ratio[kk])
    }
    mh.ratio[which(is.nan(mh.ratio)==T)] = 0
    acc.prob = runif(K)
    inds.to.replace = (1:K)[acc.prob<mh.ratio]
    alpha[inds.to.replace] = alpha.tilde[inds.to.replace]

    #-----------------save the density estimator after burnin------------------#
    if(iii > n.burnin){
        den_i <- 0
        p.vec <- as.vector(p)
        for(gg in 1:(lfg - 1)){
            den_i[gg] <- integrate(fun_integ, lower=x.grid[gg], upper=x.grid[gg+1],
                                   p = p.vec, shape = alpha, rate = beta,
                                   subdivisions = 2000)$value
        }
        den_i[lfg] <- integrate(fun_integ, lower = x.grid[lfg], upper = Inf,
                                p = p.vec, shape = alpha, rate = beta,
                                subdivisions = 2000)$value
        density.x <- density.x + rev(cumsum(rev(den_i)))
    }

    }#MCMC iteration

    density.x <- density.x/(n.MCMC - n.burnin)
    return(density.x)

}#ddsc_mcmc
