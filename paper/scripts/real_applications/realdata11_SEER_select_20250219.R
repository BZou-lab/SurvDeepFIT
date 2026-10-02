rm(list = ls())
library(survival)
library(tidyverse)
library(simsurv)
library(MASS)
library(deepTL)
library(survivalmodels)
library(randomForestSRC)
library(survivalsvm)
library(reticulate)
library(Hmisc)
library(xgboost)
library(survivalsvm)
library(pheatmap)
library(timeROC)
library(BART)

##### 0. Args Setting #####
use_condaenv("/nas/longleaf/home/{onyen}/.local/share/r-miniconda/envs/deepsurv/bin/python")
args1 = commandArgs(TRUE)
Folder_args = as.numeric(args1[1])
start.seed = as.numeric(args1[1]) *2399
Fold_Name = "/work/users/{o}/{n}/{onyen}/Permfit_Sim/Realdata11/"
#dir.create(paste(Fold_Name,Folder_args,sep = ""))
##### 1. Load RealData Set #####
Paindata_Use = readRDS("/work/users/{o}/{n}/{onyen}/Permfit_Sim/Realdata11/SEER_HCC_2018_2021.rds") %>% as.data.frame()

#Paindata_Use$GSTATUS = as.numeric(Paindata_Use$GSTATUS)
#Paindata_Use$PTIME = Paindata_Use$PTIME - Paindata_Use$LOS

#Paindata_Use = Paindata_Use %>% filter(PTIME > 0)

#colnames(Paindata_Use)[48:49] = c("OS_MONTHS","OS_STATUS")

# Make zero survival times positive before splitting, regardless of event status.
Paindata_Use$OS_MONTHS[which(Paindata_Use$OS_MONTHS == 0)] = 0.1

use_seed = 256 +start.seed
set.seed(use_seed)

n_train = floor(nrow(Paindata_Use)*0.8)
train_id = sample(nrow(Paindata_Use))[1:n_train]
valid_id = which(c(1:nrow(Paindata_Use)) %in% train_id == F)

trainN = paste(Fold_Name,Folder_args,"/traindata","_",use_seed,".rds",sep = "")
validN = paste(Fold_Name,Folder_args,"/validdata","_",use_seed,".rds",sep = "")

if (!file.exists(trainN) || !file.exists(validN)) {
  Ptrain = Paindata_Use[train_id,]
  Pvalid = Paindata_Use[valid_id,]
  # Save unstandardized splits; each script scales them in memory.
  saveRDS(Ptrain, file = trainN)
  saveRDS(Pvalid, file = validN)
} else {
  Ptrain = readRDS(trainN)
  Pvalid = readRDS(validN)
}

# Make zero survival times positive for both events and censored observations.
Ptrain$OS_MONTHS[which(Ptrain$OS_MONTHS == 0)] = 0.1
Pvalid$OS_MONTHS[which(Pvalid$OS_MONTHS == 0)] = 0.1

# Standardize both partitions using training statistics only.
for (i in c(20,25)) {
  train_mean = mean(Ptrain[,i], na.rm = TRUE)
  train_sd = sd(Ptrain[,i], na.rm = TRUE)
  if (!is.finite(train_sd) || train_sd <= 0) {
    stop(paste("Cannot standardize column", i, ": training SD must be positive and finite."))
  }
  Ptrain[,i] = (Ptrain[,i] - train_mean) / train_sd
  Pvalid[,i] = (Pvalid[,i] - train_mean) / train_sd
}

PainData_train = importDnnetSurv(x = Ptrain[,c(1:10,16:20,23:27)],
                                 y = Ptrain[,"OS_MONTHS"],
                                 e = Ptrain[,"OS_STATUS"])
PainData_valid = importDnnetSurv(x = Pvalid[,c(1:10,16:20,23:27)],
                                 y = Pvalid[,"OS_MONTHS"],
                                 e = Pvalid[,"OS_STATUS"])

########### Load Package and functions ###########
library(SurvDeepFIT)
BS_t = function(t_index,surv_df,valid = PainData_valid){
  t_i = as.numeric(colnames(surv_df)[t_index])
  death_vector = valid@e
  censor_vector = ifelse(death_vector == 1,0,1)
  T_vector = valid@y 
  km_censor = summary(survfit(Surv(T_vector,censor_vector)~1))
  bs_t = 0
  for (i in 1:nrow(surv_df)) {
    S_hat_i = as.numeric(surv_df[i,t_index])
    T_i = T_vector[i]
    D_i = death_vector[i]
    if (T_i <= t_i &D_i == 1){
      G_Hat_Ti =km_censor$surv[sum(T_i >=km_censor$time)+1]
      if (is.na(G_Hat_Ti)|G_Hat_Ti == 0){
        G_Hat_Ti =km_censor$surv[sum(T_i >=km_censor$time)]
      } else{
        G_Hat_Ti = G_Hat_Ti
      }
      bs_it = S_hat_i^2/G_Hat_Ti
    }else if (T_i > t_i){
      G_Hat_ti =km_censor$surv[sum(t_i >=km_censor$time)+1]
      if (is.na(G_Hat_ti)|G_Hat_ti == 0){
        G_Hat_ti =km_censor$surv[sum(t_i >=km_censor$time)]
      } else{
        G_Hat_ti = G_Hat_ti
      }
      bs_it = (1-S_hat_i)^2/G_Hat_ti
      if (is.infinite(bs_it) == T){
        bs_it = 0
      }
    }else{
      bs_it = 0
    }
    bs_t = bs_t+bs_it
  }
  bs_t = bs_t/nrow(surv_df)
  return(bs_t)
}

Integrated_BS = function(Brier_t,Time){
  IBS = sum((Time[2:length(Time)]-Time[1:(length(Time)-1)])*(Brier_t[1:(length(Time)-1)]+Brier_t[2:length(Time)])/2)/(Time[length(Time)]-Time[1])
  return(IBS)
}

# Perm DNNet 1 
FileN = paste(Fold_Name,Folder_args,"/PermFit_dnnet1","_",use_seed,".rds",sep = "")

if (FileN%in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F) {
  Perm_dnnet1 = permfit_survival(train = PainData_train,n_perm =100,method = "dnnet",
                                 k_fold = 5,
                                 pathway_list = list(c("race_black","race_API"),
                                                     c("grade_2","grade_3"),
                                                     c("surg_destruction","surg_resection"),
                                                     c("fibrosis_low","fibrosis_high")),
                                 n.hidden = c(50, 40, 30, 20),family = "coxph",
                                 l1.reg = 10**-4, early.stop.det = 1000, n.batch = 50,
                                 n.epoch = 1000, learning.rate.adaptive = "adam",
                                 plot = FALSE) %>% try()
  Perm_dnnet1@y_hat = NA %>% try()
  saveRDS(Perm_dnnet1,file = FileN)
  rm(Perm_dnnet1)
}
print(FileN)
print(Sys.time())

# Perm DNNet 1 
FileN = paste(Fold_Name,Folder_args,"/PermFit_dnnet1","_",use_seed,"_200.rds",sep = "")

if (FileN%in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F) {
  Perm_dnnet1 = permfit_survival(train = PainData_train,n_perm =100,method = "dnnet",
                                 k_fold = 5,
                                 pathway_list = list(c("race_black","race_API"),
                                                     c("grade_2","grade_3"),
                                                     c("surg_destruction","surg_resection"),
                                                     c("fibrosis_low","fibrosis_high")),
                                 n.hidden = c(50, 40, 30, 20),family = "coxph",
                                 l1.reg = 10**-4, early.stop.det = 1000, n.batch = 200,
                                 n.epoch = 1000, learning.rate.adaptive = "adam",
                                 plot = FALSE) %>% try()
  Perm_dnnet1@y_hat = NA %>% try()
  saveRDS(Perm_dnnet1,file = FileN)
  rm(Perm_dnnet1)
}
print(FileN)
print(Sys.time())
# Perm DNNet
FileN = paste(Fold_Name,Folder_args,"/PermFit_dnnet","_",use_seed,".rds",sep = "")

if (FileN%in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F) {
  esCtrl <- list(n.hidden = c(50, 40, 30, 20), activate = "relu",
                 l1.reg = 10**-4, early.stop.det = 1000, n.batch = 50,
                 n.epoch = 1000, learning.rate.adaptive = "adam", plot = FALSE)
  Perm_dnnet = permfit_survival(train = PainData_train,n_perm =100,method = "ensemble_dnnet",
                                k_fold = 5,
                                pathway_list = list(c("race_black","race_API"),
                                                    c("grade_2","grade_3"),
                                                    c("surg_destruction","surg_resection"),
                                                    c("fibrosis_low","fibrosis_high")),
                                n.ensemble = 100, esCtrl = esCtrl) %>% try()
  Perm_dnnet@y_hat = NA %>% try()
  saveRDS(Perm_dnnet,file = FileN)
  rm(Perm_dnnet)
}
print(FileN)
print(Sys.time())

# Perm DNNet
FileN = paste(Fold_Name,Folder_args,"/PermFit_dnnet","_",use_seed,"_200.rds",sep = "")

if (FileN%in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F) {
  esCtrl <- list(n.hidden = c(50, 40, 30, 20), activate = "relu",
                 l1.reg = 10**-4, early.stop.det = 1000, n.batch = 200,
                 n.epoch = 1000, learning.rate.adaptive = "adam", plot = FALSE)
  Perm_dnnet = permfit_survival(train = PainData_train,n_perm =100,method = "ensemble_dnnet",
                                k_fold = 5,
                                pathway_list = list(c("race_black","race_API"),
                                                    c("grade_2","grade_3"),
                                                    c("surg_destruction","surg_resection"),
                                                    c("fibrosis_low","fibrosis_high")),
                                n.ensemble = 100, esCtrl = esCtrl) %>% try()
  Perm_dnnet@y_hat = NA %>% try()
  saveRDS(Perm_dnnet,file = FileN)
  rm(Perm_dnnet)
}
print(FileN)
print(Sys.time())
# Perm Xgboost
FileN = paste(Fold_Name,Folder_args,"/PermFit_xgb","_",use_seed,".rds",sep = "")

if (FileN%in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F) {
  Perm_xgb = permfit_survival(train = PainData_train,n_perm =100,method = "Xgboost",
                              k_fold = 5,
                              pathway_list = list(c("race_black","race_API"),
                                                  c("grade_2","grade_3"),
                                                  c("surg_destruction","surg_resection"),
                                                  c("fibrosis_low","fibrosis_high")),
                              nrounds = 50) %>% try()
  Perm_xgb@y_hat = NA %>% try()
  saveRDS(Perm_xgb,file = FileN)
  rm(Perm_xgb)
}
print(FileN)
print(Sys.time())
# Perm RSF
FileN = paste(Fold_Name,Folder_args,"/PermFit_rsf","_",use_seed,".rds",sep = "")

if (FileN%in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F) {
  Perm_rsf = permfit_survival(train = PainData_train,n_perm =100,method = "random_forest",
                              k_fold = 5,
                              pathway_list = list(c("race_black","race_API"),
                                                  c("grade_2","grade_3"),
                                                  c("surg_destruction","surg_resection"),
                                                  c("fibrosis_low","fibrosis_high")),
                              ntrees = 500) %>% try()
  Perm_rsf@y_hat = NA %>% try()
  saveRDS(Perm_rsf,file = FileN)
  rm(Perm_rsf)
}
print(FileN)
print(Sys.time())
# Perm COX
FileN = paste(Fold_Name,Folder_args,"/PermFit_cox","_",use_seed,".rds",sep = "")

if (FileN%in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F) {
  Perm_cox = permfit_survival(train = PainData_train,n_perm =100,method = "survival_cox",
                              k_fold = 5,
                              pathway_list = list(c("race_black","race_API"),
                                                  c("grade_2","grade_3"),
                                                  c("surg_destruction","surg_resection"),
                                                  c("fibrosis_low","fibrosis_high"))) %>% try()
  Perm_cox@y_hat = NA %>% try()
  saveRDS(Perm_cox,file = FileN)
  rm(Perm_cox)
}
print(FileN)
print(Sys.time())
# Perm deephit
FileN = paste(Fold_Name,Folder_args,"/PermFit_dhit","_",use_seed,".rds",sep = "")

if (FileN%in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F) {
  Perm_dhit = permfit_survival(train = PainData_train,n_perm =100,method = "DeepHit",
                               k_fold = 5,
                               pathway_list = list(c("race_black","race_API"),
                                                   c("grade_2","grade_3"),
                                                   c("surg_destruction","surg_resection"),
                                                   c("fibrosis_low","fibrosis_high")),
                               activation = "relu",
                               frac = 0.2,early_stopping = T,
                               num_nodes = c(50L, 40L, 30L,20L),epochs = 1000,
                               batch_size = 50,mod_alpha = 0.2, sigma = 1, cuts = 100) %>% try()
  Perm_dhit@y_hat = NA %>% try()
  saveRDS(Perm_dhit,file = FileN)
  rm(Perm_dhit)
}
print(FileN)
print(Sys.time())

# Perm deephit
FileN = paste(Fold_Name,Folder_args,"/PermFit_dhit","_",use_seed,"_200.rds",sep = "")

if (FileN%in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F) {
  Perm_dhit = permfit_survival(train = PainData_train,n_perm =100,method = "DeepHit",
                               k_fold = 5,
                               pathway_list = list(c("race_black","race_API"),
                                                   c("grade_2","grade_3"),
                                                   c("surg_destruction","surg_resection"),
                                                   c("fibrosis_low","fibrosis_high")),
                               activation = "relu",
                               frac = 0.2,early_stopping = T,
                               num_nodes = c(50L, 40L, 30L,20L),epochs = 1000,
                               batch_size = 200,mod_alpha = 0.2, sigma = 1, cuts = 100) %>% try()
  Perm_dhit@y_hat = NA %>% try()
  saveRDS(Perm_dhit,file = FileN)
  rm(Perm_dhit)
}
print(FileN)
print(Sys.time())
# Perm deepsurv
FileN = paste(Fold_Name,Folder_args,"/PermFit_dsurv","_",use_seed,".rds",sep = "")

if (FileN%in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F) {
  Perm_dsurv = permfit_survival(train = PainData_train,n_perm =100,method = "DeepSurv",
                                k_fold = 5,
                                pathway_list = list(c("race_black","race_API"),
                                                    c("grade_2","grade_3"),
                                                    c("surg_destruction","surg_resection"),
                                                    c("fibrosis_low","fibrosis_high")),
                                activation = "relu",
                                frac = 0.2,early_stopping = T,
                                num_nodes = c(50L, 40L, 30L,20L),epochs = 1000,
                                batch_size = 50) %>% try()
  Perm_dsurv@y_hat = NA %>% try()
  saveRDS(Perm_dsurv,file = FileN)
  rm(Perm_dsurv)
}
print(FileN)
print(Sys.time())

# Perm deepsurv
FileN = paste(Fold_Name,Folder_args,"/PermFit_dsurv","_",use_seed,"_200.rds",sep = "")

if (FileN%in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F) {
  Perm_dsurv = permfit_survival(train = PainData_train,n_perm =100,method = "DeepSurv",
                                k_fold = 5,
                                pathway_list = list(c("race_black","race_API"),
                                                    c("grade_2","grade_3"),
                                                    c("surg_destruction","surg_resection"),
                                                    c("fibrosis_low","fibrosis_high")),
                                activation = "relu",
                                frac = 0.2,early_stopping = T,
                                num_nodes = c(50L, 40L, 30L,20L),epochs = 1000,
                                batch_size = 200) %>% try()
  Perm_dsurv@y_hat = NA %>% try()
  saveRDS(Perm_dsurv,file = FileN)
  rm(Perm_dsurv)
}
print(FileN)
print(Sys.time())
