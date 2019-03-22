## Nonparametric Bayesian Deconvolution of a Symmetric Unimodal Density
### Description
This repository contains the R codes to 
> 1. facilitate the application of the developed method. 

> 2. reproduce the simulation results including Figures and Tables in the titled paper. 
### Usage
We provide one R function (ddsc_mcmc) and several supporting R/Rcpp functions to deliver the proposed method. In addition, we provide the R file to generate data, perform analysis.

Two versions of ddsc_mcmc are provided in ddsc_normalError.R and ddsc_laplaceError.R corresponding to error distribution to be normal and laplace separately. The input and output variables are contained therein. 

The folders Simulation_1 and Simulation_2 are the codes for reproducing simulation results in the paper. Within each folder, there are subfolders in accordance to the homoscedastic/heteroscedastic and normal/laplace error. In the end folder, there is one R file for running analysis on one data set and another R file for carrying out parallelling analysis for various sample size and seeds.   

### License
All the codes are free to use with a proper citation of the paper. 
