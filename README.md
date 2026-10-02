# SurvDeepFIT

## Overview

SurvDeepFIT is an R package for fitting stable deep survival models and performing feature-level inference in clinical biomarker discovery. It accompanies the paper *Interpretable deep survival learning for prognostic discovery from heterogeneous health data*.

![Overview of the SurvDeepFIT framework.](paper/SurvDeepFIT_framework.png)

## Reproducibility

Code used to reproduce the paper's results is available in the `paper/` subfolder. The numerical experiments used 1,000 repetitions for each method, and the real applications used 100 repetitions for each method. Each simulation or application script processes one repetition per invocation, with the first command-line argument determining the repetition-specific random seed and output folder. We recommend running multiple repetitions in parallel through job arrays on a high-performance computing cluster. Corresponding feature selection and prediction scripts should use the same repetition ID. 

The data used in the real applications are not provided in this repository, but it can be accessed by following instructions in the `Data availability` statement in the manuscript.

## Installation

System requirements: Windows users need Rtools. No additional system tools are required for a standard macOS or Linux installation.

Install the package from GitHub in R:

``` r
devtools::install_github("BZou-lab/SurvDeepFIT")
```

Alternatively, you can install it with `pak`:

``` r
install.packages("pak")
pak::pak("BZou-lab/SurvDeepFIT")
```

For macOS users, if R reports that *libquadmath.dylib* or a similar file is missing, install gcc and gfortran first. Then copy the folder `/usr/local/gfortran/lib` to `/usr/local/lib`.

## Check Package Dependencies

Please use the following code to check whether dependencies are successfully loaded:

``` r
check_dependency()
```

## Example Simulation with a Gompertz Baseline Hazard

First, use the following code to simulate a survival dataset:

``` r
rm(list = ls())
suppressMessages(library(SurvDeepFIT))
library(MASS)
library(survival)
library(tidyverse)
library(simsurv)

# Sample size, number of features, and correlation among features
N = 1200
N_train = 1000
p_block = 10
p = 50
rho = 0.5

# Survival simulation parameters
alpha.t = 3
lambda.t = 0.00002
alpha.c = 1.5

# Simulate Continuous Features

X_Scov = data.frame(mvrnorm(n=N,mu=rep(0,p_block),
                              Sigma = (1-rho)*diag(p_block)+rho),
                      mvrnorm(n=N,mu=rep(0,p_block),
                              Sigma = (1-rho)*diag(p_block)+rho),
                      mvrnorm(n=N,mu=rep(0,p_block),
                              Sigma = (1-rho)*diag(p_block)+rho),
                      mvrnorm(n=N,mu=rep(0,p_block),
                              Sigma = (1-rho)*diag(p_block)+rho))
  colnames(X_Scov) = paste("X",1:(4*p_block),sep = "")

# Simulate Binary Features
for (z in 1:p_block) {
    if (z==1){
      Z_Scov = data.frame(rbinom(N,1,0.4))
    }else{
      Z_random = data.frame(rbinom(N,1,0.4))
      Z_Scov = cbind(Z_Scov,Z_random)
    }
}
  
colnames(Z_Scov) = paste("Z",1:(p_block),sep = "")

# Simulate Full Feature Matrix
C_Scov = cbind(X_Scov,Z_Scov)

# Significant variables: Z1, Z2, X1, X11, X21, X31
  
sims_beta = c(1,2,0.5,1,rep(0,p-1))
  
names(sims_beta) = c("Z1","X1Z2",
                    paste("X",p_block+1,"Square",sep = ""),
                    paste("X",2*p_block+1,"X",3*p_block+1,sep = ""),
                    paste("X",1:(4*p_block),sep = ""),paste("Z",2:p_block,sep = ""))
  
# Generate non-linear terms

C_Scov[,"X1Z2"] = C_Scov[,"X1"]*C_Scov[,"Z2"]
C_Scov[,paste("X",p_block+1,"Square",sep = "")] = (C_Scov[,paste("X",p_block+1,sep = "")])^2
C_Scov[,paste("X",2*p_block+1,"X",3*p_block+1,sep = "")] = C_Scov[,paste("X",2*p_block+1,sep = "")]*C_Scov[,paste("X",3*p_block+1,sep = "")]

# Generate Survival Time Using Package simsurv

Y_Scov= simsurv(dist = "gompertz",lambdas = lambda.t,gammas = alpha.t,
                betas = sims_beta,
                x= C_Scov)

C_Scov$id=1:nrow(C_Scov)

# This is the Dataset Containing True Survival Time and Features
Full_Scov = left_join(Y_Scov,C_Scov,by = "id")

# Generate censoring distribution

theta = 3.5

set.seed(1234)

Full_Scov$censor_time = rweibull(N,shape = alpha.c,scale = theta)

# Generate Survival Status and Observed Survival Time

Full_Scov$eventtime1 = apply(Full_Scov,1,function(x){
    Cens = x["censor_time"]
    Event = x["eventtime"]
    return(min(c(Cens,Event)))
})

Full_Scov$status1 = apply(Full_Scov,1,function(x){
    Cens = x["censor_time"]
    Event = x["eventtime"]
    if (Event >= Cens){
        return(0)
    } else {
        return(1)
    }
})

# Generate the dnnetSurvInput object

Train_Scov = importDnnetSurv(x = Full_Scov[1:N_train,
                                           c(paste("X",1:(4*p_block),sep = ""),
                                             paste("Z",1:(p_block),sep = ""))],
                             y = Full_Scov$eventtime1[1:N_train],
                             e = Full_Scov$status1[1:N_train])
Valid_Scov = importDnnetSurv(x = Full_Scov[(N_train+1):N,
                                           c(paste("X",1:(4*p_block),sep = ""),
                                             paste("Z",1:(p_block),sep = ""))],
                             y = Full_Scov$eventtime1[(N_train+1):N],
                             e = Full_Scov$status1[(N_train+1):N])
```

## Define the Hyperparameters Used for SurvDeepFIT

``` r
esCtrl <- list(n.hidden = c(50, 40, 30, 20), activate = "relu",
               l1.reg = 10**-4, early.stop.det = 1000, n.batch = 50,
               n.epoch = 1000, learning.rate.adaptive = "adam", plot = FALSE)
```

## PermFIT for SurvDeepFIT

``` r
PermSurvDeepFIT = permfit_survival(train = Train_Scov,n_perm =100,method = "ensemble_dnnet",
            k_fold = 5,
            n.ensemble = 100, esCtrl = esCtrl) %>% try()
View(PermSurvDeepFIT@importance)
```

## PermFIT for XGBoost

``` r
PermXGboost = permfit_survival(train = Train_Scov,n_perm =100,method = "Xgboost",
                                   k_fold = 5,nrounds = 50) %>% try()
```

## PermFIT for RSF

``` r
PermRSF = permfit_survival(train = Train_Scov,n_perm =100,method = "random_forest",
                 k_fold = 5,ntrees = 500) %>% try()
```

## PermFIT for Cox

``` r
Perm_cox_Scov = permfit_survival(train = Train_Scov,n_perm =100,method = "survival_cox",
                k_fold = 5) %>% try()
```

## Fit SurvDeepFIT Model and Predict Relative Risk Using This Model

``` r 
SurvDeepFIT_model = mod_permfit(method = "ensemble_dnnet",model.type = "survival",
                            object = Train_Scov,
                            n.ensemble = 100, esCtrl = esCtrl) %>% try()

centralized_model = model.centralize(SurvDeepFIT_model,Train_Scov) %>% try()

# Predict Relative Risk Using the Test Dataset

Predict_SurvDeepFIT = predict(centralized_model[[1]],Valid_Scov@x) %>% try() 

# Calculate the C-index

Cindex= try(1 - rcorr.cens(Predict_SurvDeepFIT,Surv(Valid_Scov@y,Valid_Scov@e))[1])
```

## Predict Survival Probability

``` r
# Predicted Survival Probability of Participants in the Test Dataset
surv_df = predict_surv_df(centralized_model,Valid_Scov)
View(surv_df)
```
