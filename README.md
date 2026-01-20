# qsrr
Quantile share ratio regression

Code for the paper "Quantile Share Ratio Regression"

This project reports R code for implementing the methods Farcomeni "Quantile Share Ratio Regression"

The following libraries shall be installed in R before using code from this project:

library(Qtools)

library(modi)

library(numDeriv)

library(compiler)

library(np)

The main function is called qsrr, with syntax:

qsrr(y,x,tau1=0.8,tau2=0.2,n.sample=NULL,n.sample2=NULL,se=FALSE,weights=NULL)

the inputs are

y: a vector with the continuous response (e.g., disposable income) 

x: a design matrix or data frame with the covariates (and intercept) 

tau1: the upper quantile of interest, default: 0.8

tau2: the lower quantile of interest, default: 0.2

n.sample: number of points at which the conditional CDF is to be estimated (v in the paper), default: length(y)  

n.sample2: number of points at this the conditional CDF is to be estimated when estimating the standard error (v in the paper), default: length(y)

se: should standard errors be estimated?, default: FALSE 

weights: vector with weights (e.g., sampling weights), default rep(1,length(y))

the outputs are 

betahat: regression parameters 

se: standard errors (or rep(NA,ncol(x)) if these are not estimated) 


