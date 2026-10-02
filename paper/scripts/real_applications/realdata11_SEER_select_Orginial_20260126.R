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

###################### Orgiginal-cox ##############################

# Cox
FileN = paste(Fold_Name,Folder_args,"/Org_cox","_",use_seed,".rds",sep = "")

if (FileN%in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F) {
  Org_Cox = mod_permfit(method = "survival_cox",model.type = "survival",
                        object = PainData_train) %>% try()
  saveRDS(Org_Cox,file = FileN)
  rm(Org_Cox)
}

print(FileN)
print(Sys.time())

# AFT
FileN = paste(Fold_Name,Folder_args,"/Org_AFT","_",use_seed,".rds",sep = "")

if (FileN%in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F) {
  Org_aft = mod_permfit(method = "survival_aft",model.type = "survival",
                        object = PainData_train) %>% try()
  saveRDS(Org_aft,file = FileN)
  rm(Org_aft)
}

print(FileN)
print(Sys.time())

# RSF
FileN = paste(Fold_Name,Folder_args,"/Org_RSF","_",use_seed,".rds",sep = "")

if (FileN%in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F) {
  Org_rsf = mod_permfit(method = "random_forest",model.type = "survival",
                        object = PainData_train,ntrees = 500,importance = T) %>% try()
  subrsf = subsample(Org_rsf) %>% try()
  Sub_rsf = extract.subsample(subrsf) %>% try()
  saveRDS(Sub_rsf,file = FileN)
  rm(Sub_rsf)
}

print(FileN)
print(Sys.time())


# XGB
FileN = paste(Fold_Name,Folder_args,"/Org_XGB","_",use_seed,".rds",sep = "")

if (FileN%in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F) {
  Org_xgb = mod_permfit(method = "Xgboost",model.type = "survival",
                        object = PainData_train,nrounds = 50) %>% try()
  xgb_imp = xgb.importance(model = Org_xgb)
  saveRDS(xgb_imp,file = FileN)
  rm(xgb_imp)
}

print(FileN)
print(Sys.time())




