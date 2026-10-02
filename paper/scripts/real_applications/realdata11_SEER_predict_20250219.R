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

##### 0. Args Setting #####
use_condaenv("/nas/longleaf/home/{onyen}/.local/share/r-miniconda/envs/deepsurv/bin/python")

args1 = commandArgs(TRUE)
Folder_args = as.numeric(args1[1])
start.seed = as.numeric(args1[1]) *2399
Fold_Name = "/work/users/{o}/{n}/{onyen}/Permfit_Sim/Realdata11/"
# missing = readRDS("/work/users/{o}/{n}/{onyen}/Permfit_Sim/Realdata8/shhs1-0.21.0-all-cause-mortality-discard.rds") %>% as.data.frame()

p_thershold = 0.05

use_seed = 256 +start.seed
set.seed(use_seed)

trainN = paste(Fold_Name,Folder_args,"/traindata","_",use_seed,".rds",sep = "")
validN = paste(Fold_Name,Folder_args,"/validdata","_",use_seed,".rds",sep = "")

Ptrain = readRDS(trainN)
Pvalid = readRDS(validN)
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


# 1.Load Functions 

#source("/nas/longleaf/home/{onyen}/R-Code/source_permfit_survival_20230610.R")
library(SurvDeepFIT)
predict_surv_df_cox = function (centralized_model, test_object, time_point = NULL) 
{
  mod = centralized_model[[1]]
  train_object = centralized_model[[2]]
  fail_time = train_object@y[which(train_object@e == 1)]
  unique_time_in_train = unique(train_object@y)
  sort_time_in_train = sort(unique_time_in_train)
  unique_fail_time = unique(fail_time)
  sort_fail_time = sort(unique_fail_time)
  train_risk_score = predict(mod, data.frame(train_object@x))
  test_risk_score = predict(mod, data.frame(test_object@x))
  dH0 = c()
  for (i in 1:length(sort_time_in_train)) {
    T_i = sort_time_in_train[i]
    Death_at_T_i = which(train_object@y == T_i & train_object@e == 
                           1)
    D_i = length(Death_at_T_i)
    Riskset_i = risk.set(t_threshold = T_i, times = train_object@y)
    dH_0_Ti = D_i/sum(exp(train_risk_score)[which(train_object@y >= 
                                                    T_i)])
    dH0 = append(dH0, dH_0_Ti)
  }
  if (is.null(time_point)) {
    H0 = c()
    for (i in 1:length(sort_fail_time)) {
      t = sort_fail_time[i]
      H0_t = sum(dH0[1:last(which(sort_time_in_train <= 
                                    t))])
      H0 = append(H0, H0_t)
    }
    Surv_df_test = data.frame()
    for (i in 1:length(test_risk_score)) {
      surv_df = t(as.data.frame(exp(-H0 * exp(test_risk_score[i]))))
      colnames(surv_df) = round(sort_fail_time, 3)
      rownames(surv_df) = i
      Surv_df_test = rbind(Surv_df_test, surv_df)
    }
    return(Surv_df_test)
  }
  else {
    sort_fail_time = sort(time_point)
    H0 = c()
    for (i in 1:length(sort_fail_time)) {
      t = sort_fail_time[i]
      H0_t = sum(dH0[1:last(which(sort_time_in_train <= 
                                    t))])
      H0 = append(H0, H0_t)
    }
    Surv_df_test = data.frame()
    for (i in 1:length(test_risk_score)) {
      surv_df = t(as.data.frame(exp(-H0 * exp(test_risk_score[i]))))
      colnames(surv_df) = round(sort_fail_time, 3)
      rownames(surv_df) = i
      Surv_df_test = rbind(Surv_df_test, surv_df)
    }
    return(Surv_df_test)
  }
}

predict_surv_df_dnnet = function (centralized_model, test_object, time_point = NULL) 
{
  mod = centralized_model[[1]]
  train_object = centralized_model[[2]]
  fail_time = train_object@y[which(train_object@e == 1)]
  unique_time_in_train = unique(train_object@y)
  sort_time_in_train = sort(unique_time_in_train)
  unique_fail_time = unique(fail_time)
  sort_fail_time = sort(unique_fail_time)
  train_risk_score = predict(mod, train_object@x)
  test_risk_score = predict(mod, test_object@x)
  dH0 = c()
  for (i in 1:length(sort_time_in_train)) {
    T_i = sort_time_in_train[i]
    Death_at_T_i = which(train_object@y == T_i & train_object@e == 
                           1)
    D_i = length(Death_at_T_i)
    Riskset_i = risk.set(t_threshold = T_i, times = train_object@y)
    dH_0_Ti = D_i/sum(exp(train_risk_score)[which(train_object@y >= 
                                                    T_i)])
    dH0 = append(dH0, dH_0_Ti)
  }
  if (is.null(time_point)) {
    H0 = c()
    for (i in 1:length(sort_fail_time)) {
      t = sort_fail_time[i]
      H0_t = sum(dH0[1:last(which(sort_time_in_train <= 
                                    t))])
      H0 = append(H0, H0_t)
    }
    Surv_df_test = data.frame()
    for (i in 1:length(test_risk_score)) {
      surv_df = t(as.data.frame(exp(-H0 * exp(test_risk_score[i]))))
      colnames(surv_df) = round(sort_fail_time, 3)
      rownames(surv_df) = i
      Surv_df_test = rbind(Surv_df_test, surv_df)
    }
    return(Surv_df_test)
  }
  else {
    sort_fail_time = sort(time_point)
    H0 = c()
    for (i in 1:length(sort_fail_time)) {
      t = sort_fail_time[i]
      H0_t = sum(dH0[1:last(which(sort_time_in_train <= 
                                    t))])
      H0 = append(H0, H0_t)
    }
    Surv_df_test = data.frame()
    for (i in 1:length(test_risk_score)) {
      surv_df = t(as.data.frame(exp(-H0 * exp(test_risk_score[i]))))
      colnames(surv_df) = round(sort_fail_time, 3)
      rownames(surv_df) = i
      Surv_df_test = rbind(Surv_df_test, surv_df)
    }
    return(Surv_df_test)
  }
}

predict_surv_df_xgb = function (centralized_model, test_object, time_point = NULL) 
{
  mod = centralized_model[[1]]
  train_object = centralized_model[[2]]
  fail_time = train_object@y[which(train_object@e == 1)]
  unique_time_in_train = unique(train_object@y)
  sort_time_in_train = sort(unique_time_in_train)
  unique_fail_time = unique(fail_time)
  sort_fail_time = sort(unique_fail_time)
  train_risk_score = predict(mod, train_object@x) %>% log()
  test_risk_score = predict(mod, test_object@x) %>% log()
  dH0 = c()
  for (i in 1:length(sort_time_in_train)) {
    T_i = sort_time_in_train[i]
    Death_at_T_i = which(train_object@y == T_i & train_object@e == 
                           1)
    D_i = length(Death_at_T_i)
    Riskset_i = risk.set(t_threshold = T_i, times = train_object@y)
    dH_0_Ti = D_i/sum(exp(train_risk_score)[which(train_object@y >= 
                                                    T_i)])
    dH0 = append(dH0, dH_0_Ti)
  }
  if (is.null(time_point)) {
    H0 = c()
    for (i in 1:length(sort_fail_time)) {
      t = sort_fail_time[i]
      H0_t = sum(dH0[1:last(which(sort_time_in_train <= 
                                    t))])
      H0 = append(H0, H0_t)
    }
    Surv_df_test = data.frame()
    for (i in 1:length(test_risk_score)) {
      surv_df = t(as.data.frame(exp(-H0 * exp(test_risk_score[i]))))
      colnames(surv_df) = round(sort_fail_time, 3)
      rownames(surv_df) = i
      Surv_df_test = rbind(Surv_df_test, surv_df)
    }
    return(Surv_df_test)
  }
  else {
    sort_fail_time = sort(time_point)
    H0 = c()
    for (i in 1:length(sort_fail_time)) {
      t = sort_fail_time[i]
      H0_t = sum(dH0[1:last(which(sort_time_in_train <= 
                                    t))])
      H0 = append(H0, H0_t)
    }
    Surv_df_test = data.frame()
    for (i in 1:length(test_risk_score)) {
      surv_df = t(as.data.frame(exp(-H0 * exp(test_risk_score[i]))))
      colnames(surv_df) = round(sort_fail_time, 3)
      rownames(surv_df) = i
      Surv_df_test = rbind(Surv_df_test, surv_df)
    }
    return(Surv_df_test)
  }
}
####################################


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

# PermFit Cox

FileN = paste(Fold_Name,Folder_args,"/PermFit_cox","_",use_seed,".rds",sep = "")

result = readRDS(file = FileN) %>% try()

if(class(result) != "try-error") {
  imp = result@importance
  
  imp_block = result@block_importance
  
  p.value.race = imp_block[1,"importance_pval_x"]
  p.value.grade = imp_block[2,"importance_pval_x"]
  p.value.surg = imp_block[3,"importance_pval_x"]
  p.value.fibrosis = imp_block[4,"importance_pval_x"]
  
  if (p.value.race > p_thershold) {imp = imp[-which(imp$var_name %in% c("race_black","race_API")),]}
  if (p.value.grade > p_thershold) {imp = imp[-which(imp$var_name %in% c("grade_2","grade_3")),]}
  if (p.value.surg > p_thershold) {imp = imp[-which(imp$var_name %in% c("surg_destruction","surg_resection")),]}
  if (p.value.fibrosis > p_thershold) {imp = imp[-which(imp$var_name %in% c("fibrosis_low","fibrosis_high")),]}
  
  
  sig_var = imp$var_name[which(imp$importance_pval_x < p_thershold)]
  
  if(length(sig_var) < 5){
    arrange = imp %>% arrange(importance_pval_x)
    sig_var = arrange$var_name[1:5]
  }
  
  selected_train = importDnnetSurv(x = Ptrain[,sig_var],
                                   y = Ptrain[,"OS_MONTHS"],
                                   e = Ptrain[,"OS_STATUS"])
  selected_valid = importDnnetSurv(x = Pvalid[,sig_var],
                                   y = Pvalid[,"OS_MONTHS"],
                                   e = Pvalid[,"OS_STATUS"])
  
  
  org_model = mod_permfit(method = "survival_cox",model.type = "survival",
                          object = PainData_train) %>% try()
  
  selected_model = mod_permfit(method = "survival_cox",model.type = "survival",
                               object = selected_train) %>% try()
  
  pred_org = predict_mod_permfit(mod = org_model,
                                 method = "survival_cox",
                                 model.type = "survival",
                                 object = PainData_valid) %>% try()
  
  pred_selected = predict_mod_permfit(mod = selected_model,
                                      method = "survival_cox",
                                      model.type = "survival",
                                      object = selected_valid) %>% try()
  
  
  Ccox_org = try(1 - rcorr.cens(pred_org[[1]],Surv(PainData_valid@y,PainData_valid@e))[1])
  Ccox_selected = try(1 - rcorr.cens(pred_selected[[1]],Surv(PainData_valid@y,PainData_valid@e))[1])
  
  ##### IBS
  # org on test
  pred_org_surv = predict_surv_df_cox(centralized_model = list(org_model,PainData_train),
                                      test_object = PainData_valid) %>% try()
  # selected on test
  pred_selected_surv = predict_surv_df_cox(centralized_model = list(selected_model,selected_train),
                                           test_object = selected_valid) %>% try()
  
  bss = c()
  ts=c()
  for (i in 1:ncol(pred_org_surv)) {
    ts = append(ts,as.numeric(colnames(pred_org_surv)[i])) %>% try()
    bss = append(bss,BS_t(t_index = i,surv_df = pred_org_surv,valid = PainData_valid))%>% try()
  }
  IBS_cox_org = Integrated_BS(bss,ts) %>% try()
  
  bss = c()
  ts=c()
  for (i in 1:ncol(pred_selected_surv)) {
    ts = append(ts,as.numeric(colnames(pred_selected_surv)[i])) %>% try()
    bss = append(bss,BS_t(t_index = i,surv_df = pred_selected_surv,valid = selected_valid))%>% try()
  }
  IBS_cox_selected = Integrated_BS(bss,ts) %>% try()
  
  
}else {
  Ccox_org = NA
  Ccox_selected = NA
  IBS_cox_org = NA
  IBS_cox_selected = NA
}

# PermFit-dnnet and dnnet

FileN = paste(Fold_Name,Folder_args,"/PermFit_dnnet","_",use_seed,".rds",sep = "")

esCtrl <- list(n.hidden = c(50, 40, 30, 20), activate = "relu",
               l1.reg = 10**-4, early.stop.det = 1000, n.batch = 50,
               n.epoch = 1000, learning.rate.adaptive = "adam", plot = FALSE) %>% try()

result = readRDS(file = FileN) %>% try()

if(class(result) != "try-error"){
  imp = result@importance
  
  imp_block = result@block_importance
  
  p.value.race = imp_block[1,"importance_pval_x"]
  p.value.grade = imp_block[2,"importance_pval_x"]
  p.value.surg = imp_block[3,"importance_pval_x"]
  p.value.fibrosis = imp_block[4,"importance_pval_x"]
  
  if (p.value.race > p_thershold) {imp = imp[-which(imp$var_name %in% c("race_black","race_API")),]}
  if (p.value.grade > p_thershold) {imp = imp[-which(imp$var_name %in% c("grade_2","grade_3")),]}
  if (p.value.surg > p_thershold) {imp = imp[-which(imp$var_name %in% c("surg_destruction","surg_resection")),]}
  if (p.value.fibrosis > p_thershold) {imp = imp[-which(imp$var_name %in% c("fibrosis_low","fibrosis_high")),]}
  
  sig_var = imp$var_name[which(imp$importance_pval_x < p_thershold)]
  
  if(length(sig_var) < 5){
    arrange = imp %>% arrange(importance_pval_x)
    sig_var = arrange$var_name[1:5]
  }
  
  selected_train = importDnnetSurv(x = Ptrain[,sig_var],
                                   y = Ptrain[,"OS_MONTHS"],
                                   e = Ptrain[,"OS_STATUS"])
  selected_valid = importDnnetSurv(x = Pvalid[,sig_var],
                                   y = Pvalid[,"OS_MONTHS"],
                                   e = Pvalid[,"OS_STATUS"])
  
  org_model = mod_permfit(method = "ensemble_dnnet",model.type = "survival",
                          object = PainData_train,
                          n.ensemble = 100, esCtrl = esCtrl) %>% try()
  
  selected_model = mod_permfit(method = "ensemble_dnnet",model.type = "survival",
                               object = selected_train,
                               n.ensemble = 100, esCtrl = esCtrl) %>% try()
  
  
  pred_org = predict_mod_permfit(mod = org_model,
                                 method = "ensemble_dnnet",
                                 model.type = "survival",
                                 object = PainData_valid) %>% try()
  
  pred_selected = predict_mod_permfit(mod = selected_model,
                                      method = "ensemble_dnnet",
                                      model.type = "survival",
                                      object = selected_valid) %>% try()
  
  
  Cdnnet_org = try(1 - rcorr.cens(pred_org[[1]],Surv(PainData_valid@y,PainData_valid@e))[1])
  Cdnnet_selected = try(1 - rcorr.cens(pred_selected[[1]],Surv(PainData_valid@y,PainData_valid@e))[1])
  
  centralized_model_org = model.centralize(org_model,PainData_train) %>% try()
  centralized_model_selected = model.centralize(selected_model,selected_train) %>% try()
  
  pred_org_surv = predict_surv_df_dnnet(centralized_model = centralized_model_org,
                                        test_object = PainData_valid) %>% try()
  pred_selected_surv = predict_surv_df_dnnet(centralized_model = centralized_model_selected,
                                             test_object = selected_valid) %>% try()
  
  bss = c()
  ts=c()
  for (i in 1:ncol(pred_org_surv)) {
    ts = append(ts,as.numeric(colnames(pred_org_surv)[i])) %>% try()
    bss = append(bss,BS_t(t_index = i,surv_df = pred_org_surv,valid = PainData_valid))%>% try()
  }
  IBS_dnnet_org = Integrated_BS(bss,ts) %>% try()
  
  bss = c()
  ts=c()
  for (i in 1:ncol(pred_selected_surv)) {
    ts = append(ts,as.numeric(colnames(pred_selected_surv)[i])) %>% try()
    bss = append(bss,BS_t(t_index = i,surv_df = pred_selected_surv,valid = selected_valid))%>% try()
  }
  IBS_dnnet_selected = Integrated_BS(bss,ts) %>% try()
  
  
}else {
  Cdnnet_org = NA
  Cdnnet_selected = NA
  IBS_dnnet_org = NA
  IBS_dnnet_selected = NA
}


# PermFit-dnnet and dnnet batch = 200

FileN = paste(Fold_Name,Folder_args,"/PermFit_dnnet","_",use_seed,"_200.rds",sep = "")

esCtrl <- list(n.hidden = c(50, 40, 30, 20), activate = "relu",
               l1.reg = 10**-4, early.stop.det = 1000, n.batch = 200,
               n.epoch = 1000, learning.rate.adaptive = "adam", plot = FALSE) %>% try()

result = readRDS(file = FileN) %>% try()

if(class(result) != "try-error"){
  imp = result@importance
  
  imp_block = result@block_importance
  
  p.value.race = imp_block[1,"importance_pval_x"]
  p.value.grade = imp_block[2,"importance_pval_x"]
  p.value.surg = imp_block[3,"importance_pval_x"]
  p.value.fibrosis = imp_block[4,"importance_pval_x"]
  
  if (p.value.race > p_thershold) {imp = imp[-which(imp$var_name %in% c("race_black","race_API")),]}
  if (p.value.grade > p_thershold) {imp = imp[-which(imp$var_name %in% c("grade_2","grade_3")),]}
  if (p.value.surg > p_thershold) {imp = imp[-which(imp$var_name %in% c("surg_destruction","surg_resection")),]}
  if (p.value.fibrosis > p_thershold) {imp = imp[-which(imp$var_name %in% c("fibrosis_low","fibrosis_high")),]}
  
  sig_var = imp$var_name[which(imp$importance_pval_x < p_thershold)]
  
  if(length(sig_var) < 5){
    arrange = imp %>% arrange(importance_pval_x)
    sig_var = arrange$var_name[1:5]
  }
  
  selected_train = importDnnetSurv(x = Ptrain[,sig_var],
                                   y = Ptrain[,"OS_MONTHS"],
                                   e = Ptrain[,"OS_STATUS"])
  selected_valid = importDnnetSurv(x = Pvalid[,sig_var],
                                   y = Pvalid[,"OS_MONTHS"],
                                   e = Pvalid[,"OS_STATUS"])
  
  org_model = mod_permfit(method = "ensemble_dnnet",model.type = "survival",
                          object = PainData_train,
                          n.ensemble = 100, esCtrl = esCtrl) %>% try()
  
  selected_model = mod_permfit(method = "ensemble_dnnet",model.type = "survival",
                               object = selected_train,
                               n.ensemble = 100, esCtrl = esCtrl) %>% try()
  
  
  pred_org = predict_mod_permfit(mod = org_model,
                                 method = "ensemble_dnnet",
                                 model.type = "survival",
                                 object = PainData_valid) %>% try()
  
  pred_selected = predict_mod_permfit(mod = selected_model,
                                      method = "ensemble_dnnet",
                                      model.type = "survival",
                                      object = selected_valid) %>% try()
  
  
  Cdnnet_org200 = try(1 - rcorr.cens(pred_org[[1]],Surv(PainData_valid@y,PainData_valid@e))[1])
  Cdnnet_selected200 = try(1 - rcorr.cens(pred_selected[[1]],Surv(PainData_valid@y,PainData_valid@e))[1])
  
  centralized_model_org = model.centralize(org_model,PainData_train) %>% try()
  centralized_model_selected = model.centralize(selected_model,selected_train) %>% try()
  
  pred_org_surv = predict_surv_df_dnnet(centralized_model = centralized_model_org,
                                        test_object = PainData_valid) %>% try()
  pred_selected_surv = predict_surv_df_dnnet(centralized_model = centralized_model_selected,
                                             test_object = selected_valid) %>% try()
  
  bss = c()
  ts=c()
  for (i in 1:ncol(pred_org_surv)) {
    ts = append(ts,as.numeric(colnames(pred_org_surv)[i])) %>% try()
    bss = append(bss,BS_t(t_index = i,surv_df = pred_org_surv,valid = PainData_valid))%>% try()
  }
  IBS_dnnet_org200 = Integrated_BS(bss,ts) %>% try()
  
  bss = c()
  ts=c()
  for (i in 1:ncol(pred_selected_surv)) {
    ts = append(ts,as.numeric(colnames(pred_selected_surv)[i])) %>% try()
    bss = append(bss,BS_t(t_index = i,surv_df = pred_selected_surv,valid = selected_valid))%>% try()
  }
  IBS_dnnet_selected200 = Integrated_BS(bss,ts) %>% try()
  
  
}else {
  Cdnnet_org200 = NA
  Cdnnet_selected200 = NA
  IBS_dnnet_org200 = NA
  IBS_dnnet_selected200 = NA
}

# PermFit RSF and RSF

FileN = paste(Fold_Name,Folder_args,"/PermFit_rsf","_",use_seed,".rds",sep = "")

result = readRDS(file = FileN) %>% try()

if(class(result) != "try-error"){
  imp = result@importance
  
  imp_block = result@block_importance
  
  p.value.race = imp_block[1,"importance_pval_x"]
  p.value.grade = imp_block[2,"importance_pval_x"]
  p.value.surg = imp_block[3,"importance_pval_x"]
  p.value.fibrosis = imp_block[4,"importance_pval_x"]
  
  if (p.value.race > p_thershold) {imp = imp[-which(imp$var_name %in% c("race_black","race_API")),]}
  if (p.value.grade > p_thershold) {imp = imp[-which(imp$var_name %in% c("grade_2","grade_3")),]}
  if (p.value.surg > p_thershold) {imp = imp[-which(imp$var_name %in% c("surg_destruction","surg_resection")),]}
  if (p.value.fibrosis > p_thershold) {imp = imp[-which(imp$var_name %in% c("fibrosis_low","fibrosis_high")),]}
  
  sig_var = imp$var_name[which(imp$importance_pval_x < p_thershold)]
  
  if(length(sig_var) < 5){
    arrange = imp %>% arrange(importance_pval_x)
    sig_var = arrange$var_name[1:5]
  }
  
  selected_train = importDnnetSurv(x = Ptrain[,sig_var],
                                   y = Ptrain[,"OS_MONTHS"],
                                   e = Ptrain[,"OS_STATUS"])
  selected_valid = importDnnetSurv(x = Pvalid[,sig_var],
                                   y = Pvalid[,"OS_MONTHS"],
                                   e = Pvalid[,"OS_STATUS"])
  
  org_model = mod_permfit(method = "random_forest",model.type = "survival",
                          object = PainData_train,ntrees = 500) %>% try()
  
  selected_model = mod_permfit(method = "random_forest",model.type = "survival",
                               object = selected_train, ntrees = 500) %>% try()
  
  
  pred_org = predict_mod_permfit(mod = org_model,
                                 method = "random_forest",
                                 model.type = "survival",
                                 object = PainData_valid) %>% try()
  
  pred_selected = predict_mod_permfit(mod = selected_model,
                                      method = "random_forest",
                                      model.type = "survival",
                                      object = selected_valid) %>% try()
  
  
  Crsf_org = try(1 - rcorr.cens(pred_org[[1]],Surv(PainData_valid@y,PainData_valid@e))[1])
  Crsf_selected = try(1 - rcorr.cens(pred_selected[[1]],Surv(PainData_valid@y,PainData_valid@e))[1])
  
  
  org_rsf = predict(org_model,data.frame(x = PainData_valid@x)) %>% try()
  pred_org_surv = org_rsf$survival %>% try()
  colnames(pred_org_surv) = org_rsf$time.interest %>% try()
  
  selected_rsf = predict(selected_model,data.frame(x = selected_valid@x)) %>% try()
  pred_selected_surv = selected_rsf$survival %>% try()
  colnames(pred_selected_surv) = selected_rsf$time.interest %>% try()
  
  
  bss = c()
  ts=c()
  for (i in 1:ncol(pred_org_surv)) {
    ts = append(ts,as.numeric(colnames(pred_org_surv)[i])) %>% try()
    bss = append(bss,BS_t(t_index = i,surv_df = pred_org_surv,valid = PainData_valid))%>% try()
  }
  IBS_rsf_org = Integrated_BS(bss,ts) %>% try()
  
  bss = c()
  ts=c()
  for (i in 1:ncol(pred_selected_surv)) {
    ts = append(ts,as.numeric(colnames(pred_selected_surv)[i])) %>% try()
    bss = append(bss,BS_t(t_index = i,surv_df = pred_selected_surv,valid = selected_valid))%>% try()
  }
  IBS_rsf_selected = Integrated_BS(bss,ts) %>% try()
  
  
}else {
  Crsf_org = NA
  Crsf_selected = NA
  IBS_rsf_org = NA
  IBS_rsf_selected = NA
}

# PermFit XGboost and XGboost

FileN = paste(Fold_Name,Folder_args,"/PermFit_xgb","_",use_seed,".rds",sep = "")

result = readRDS(file = FileN) %>% try()

if(class(result) != "try-error"){
  imp = result@importance
  
  imp_block = result@block_importance
  
  p.value.race = imp_block[1,"importance_pval_x"]
  p.value.grade = imp_block[2,"importance_pval_x"]
  p.value.surg = imp_block[3,"importance_pval_x"]
  p.value.fibrosis = imp_block[4,"importance_pval_x"]
  
  if (p.value.race > p_thershold) {imp = imp[-which(imp$var_name %in% c("race_black","race_API")),]}
  if (p.value.grade > p_thershold) {imp = imp[-which(imp$var_name %in% c("grade_2","grade_3")),]}
  if (p.value.surg > p_thershold) {imp = imp[-which(imp$var_name %in% c("surg_destruction","surg_resection")),]}
  if (p.value.fibrosis > p_thershold) {imp = imp[-which(imp$var_name %in% c("fibrosis_low","fibrosis_high")),]}
  
  sig_var = imp$var_name[which(imp$importance_pval_x < p_thershold)]
  
  if(length(sig_var) < 5){
    arrange = imp %>% arrange(importance_pval_x)
    sig_var = arrange$var_name[1:5]
  }
  
  selected_train = importDnnetSurv(x = Ptrain[,sig_var],
                                   y = Ptrain[,"OS_MONTHS"],
                                   e = Ptrain[,"OS_STATUS"])
  selected_valid = importDnnetSurv(x = Pvalid[,sig_var],
                                   y = Pvalid[,"OS_MONTHS"],
                                   e = Pvalid[,"OS_STATUS"])
  
  
  org_model = mod_permfit(method = "Xgboost",model.type = "survival",
                          object = PainData_train,nrounds = 50) %>% try()
  
  selected_model = mod_permfit(method = "Xgboost",model.type = "survival",
                               object = selected_train,nrounds = 50) %>% try()
  
  pred_org = predict_mod_permfit(mod = org_model,
                                 method = "Xgboost",
                                 model.type = "survival",
                                 object = PainData_valid) %>% try()
  
  pred_selected = predict_mod_permfit(mod = selected_model,
                                      method = "Xgboost",
                                      model.type = "survival",
                                      object = selected_valid) %>% try()
  
  
  Cxgb_org = try(1 - rcorr.cens(pred_org[[1]],Surv(PainData_valid@y,PainData_valid@e))[1])
  Cxgb_selected = try(1 - rcorr.cens(pred_selected[[1]],Surv(selected_valid@y,selected_valid@e))[1])
  
  ##### IBS
  # org on test
  pred_org_surv = predict_surv_df_xgb(centralized_model = list(org_model,PainData_train),
                                      test_object = PainData_valid) %>% try()
  # selected on test
  pred_selected_surv = predict_surv_df_xgb(centralized_model = list(selected_model,selected_train),
                                           test_object = selected_valid) %>% try()
  
  bss = c()
  ts=c()
  for (i in 1:ncol(pred_org_surv)) {
    ts = append(ts,as.numeric(colnames(pred_org_surv)[i])) %>% try()
    bss = append(bss,BS_t(t_index = i,surv_df = pred_org_surv,valid = PainData_valid))%>% try()
  }
  IBS_xgb_org = Integrated_BS(bss,ts) %>% try()
  
  bss = c()
  ts=c()
  for (i in 1:ncol(pred_selected_surv)) {
    ts = append(ts,as.numeric(colnames(pred_selected_surv)[i])) %>% try()
    bss = append(bss,BS_t(t_index = i,surv_df = pred_selected_surv,valid = selected_valid))%>% try()
  }
  IBS_xgb_selected = Integrated_BS(bss,ts) %>% try()
  
  
}else {
  Cxgb_org = NA
  Cxgb_selected = NA
  IBS_xgb_org = NA
  IBS_xgb_selected = NA
}

# PermFit DeepSurv and DeepSurv

FileN = paste(Fold_Name,Folder_args,"/PermFit_dsurv","_",use_seed,".rds",sep = "")

result = readRDS(file = FileN) %>% try()

if(class(result) != "try-error"){
  imp = result@importance
  
  imp_block = result@block_importance
  
  p.value.race = imp_block[1,"importance_pval_x"]
  p.value.grade = imp_block[2,"importance_pval_x"]
  p.value.surg = imp_block[3,"importance_pval_x"]
  p.value.fibrosis = imp_block[4,"importance_pval_x"]
  
  if (p.value.race > p_thershold) {imp = imp[-which(imp$var_name %in% c("race_black","race_API")),]}
  if (p.value.grade > p_thershold) {imp = imp[-which(imp$var_name %in% c("grade_2","grade_3")),]}
  if (p.value.surg > p_thershold) {imp = imp[-which(imp$var_name %in% c("surg_destruction","surg_resection")),]}
  if (p.value.fibrosis > p_thershold) {imp = imp[-which(imp$var_name %in% c("fibrosis_low","fibrosis_high")),]}
  
  
  sig_var = imp$var_name[which(imp$importance_pval_x < p_thershold)]
  
  if(length(sig_var) < 2){
    arrange = imp %>% arrange(importance_pval_x)
    sig_var = arrange$var_name[1:5]
  }
  
  selected_train = importDnnetSurv(x = Ptrain[,sig_var],
                                   y = Ptrain[,"OS_MONTHS"],
                                   e = Ptrain[,"OS_STATUS"])
  selected_valid = importDnnetSurv(x = Pvalid[,sig_var],
                                   y = Pvalid[,"OS_MONTHS"],
                                   e = Pvalid[,"OS_STATUS"])
  
  
  org_model = mod_permfit(method = "DeepSurv",model.type = "survival",
                          object = PainData_train,activation = "relu",
                          frac = 0.2,early_stopping = T,
                          num_nodes = c(50L, 40L, 30L,20L),epochs = 1000,
                          batch_size = 50) %>% try()
  
  selected_model = mod_permfit(method = "DeepSurv",model.type = "survival",
                               object = selected_train,activation = "relu",
                               frac = 0.2,early_stopping = T,
                               num_nodes = c(50L, 40L, 30L,20L),epochs = 1000,
                               batch_size = 50) %>% try()
  
  
  pred_org = predict_mod_permfit(mod = org_model,
                                 method = "DeepSurv",
                                 model.type = "survival",
                                 object = PainData_valid) %>% try()
  
  pred_selected = predict_mod_permfit(mod = selected_model,
                                      method = "DeepSurv",
                                      model.type = "survival",
                                      object = selected_valid) %>% try()
  
  
  Cdsurv_org = try(1 - rcorr.cens(pred_org[[1]],Surv(PainData_valid@y,PainData_valid@e))[1])
  Cdsurv_selected = try(1 - rcorr.cens(pred_selected[[1]],Surv(selected_valid@y,selected_valid@e))[1])
  
  pred_org_surv = predict(org_model,data.frame(x = PainData_valid@x),type = "survival") %>% try()
  pred_selected_surv = predict(selected_model,data.frame(x = selected_valid@x),type = "survival") %>% try()
  
  bss = c()
  ts=c()
  for (i in 1:ncol(pred_org_surv)) {
    ts = append(ts,as.numeric(colnames(pred_org_surv)[i])) %>% try()
    bss = append(bss,BS_t(t_index = i,surv_df = pred_org_surv,valid = PainData_valid))%>% try()
  }
  IBS_dsurv_org = Integrated_BS(bss,ts) %>% try()
  
  bss = c()
  ts=c()
  for (i in 1:ncol(pred_selected_surv)) {
    ts = append(ts,as.numeric(colnames(pred_selected_surv)[i])) %>% try()
    bss = append(bss,BS_t(t_index = i,surv_df = pred_selected_surv,valid = selected_valid))%>% try()
  }
  IBS_dsurv_selected = Integrated_BS(bss,ts) %>% try()
  
  
}else {
  Cdsurv_org = NA
  Cdsurv_selected = NA
  IBS_dsurv_org = NA
  IBS_dsurv_selected = NA
}


# PermFit DeepSurv and DeepSurv batch = 200

FileN = paste(Fold_Name,Folder_args,"/PermFit_dsurv","_",use_seed,"_200.rds",sep = "")

result = readRDS(file = FileN) %>% try()

if(class(result) != "try-error"){
  imp = result@importance
  
  imp_block = result@block_importance
  
  p.value.race = imp_block[1,"importance_pval_x"]
  p.value.grade = imp_block[2,"importance_pval_x"]
  p.value.surg = imp_block[3,"importance_pval_x"]
  p.value.fibrosis = imp_block[4,"importance_pval_x"]
  
  if (p.value.race > p_thershold) {imp = imp[-which(imp$var_name %in% c("race_black","race_API")),]}
  if (p.value.grade > p_thershold) {imp = imp[-which(imp$var_name %in% c("grade_2","grade_3")),]}
  if (p.value.surg > p_thershold) {imp = imp[-which(imp$var_name %in% c("surg_destruction","surg_resection")),]}
  if (p.value.fibrosis > p_thershold) {imp = imp[-which(imp$var_name %in% c("fibrosis_low","fibrosis_high")),]}
  
  
  sig_var = imp$var_name[which(imp$importance_pval_x < p_thershold)]
  
  if(length(sig_var) < 2){
    arrange = imp %>% arrange(importance_pval_x)
    sig_var = arrange$var_name[1:5]
  }
  
  selected_train = importDnnetSurv(x = Ptrain[,sig_var],
                                   y = Ptrain[,"OS_MONTHS"],
                                   e = Ptrain[,"OS_STATUS"])
  selected_valid = importDnnetSurv(x = Pvalid[,sig_var],
                                   y = Pvalid[,"OS_MONTHS"],
                                   e = Pvalid[,"OS_STATUS"])
  
  
  org_model = mod_permfit(method = "DeepSurv",model.type = "survival",
                          object = PainData_train,activation = "relu",
                          frac = 0.2,early_stopping = T,
                          num_nodes = c(50L, 40L, 30L,20L),epochs = 1000,
                          batch_size = 200) %>% try()
  
  selected_model = mod_permfit(method = "DeepSurv",model.type = "survival",
                               object = selected_train,activation = "relu",
                               frac = 0.2,early_stopping = T,
                               num_nodes = c(50L, 40L, 30L,20L),epochs = 1000,
                               batch_size = 200) %>% try()
  
  
  pred_org = predict_mod_permfit(mod = org_model,
                                 method = "DeepSurv",
                                 model.type = "survival",
                                 object = PainData_valid) %>% try()
  
  pred_selected = predict_mod_permfit(mod = selected_model,
                                      method = "DeepSurv",
                                      model.type = "survival",
                                      object = selected_valid) %>% try()
  
  
  Cdsurv_org200 = try(1 - rcorr.cens(pred_org[[1]],Surv(PainData_valid@y,PainData_valid@e))[1])
  Cdsurv_selected200 = try(1 - rcorr.cens(pred_selected[[1]],Surv(selected_valid@y,selected_valid@e))[1])
  
  pred_org_surv = predict(org_model,data.frame(x = PainData_valid@x),type = "survival") %>% try()
  pred_selected_surv = predict(selected_model,data.frame(x = selected_valid@x),type = "survival") %>% try()
  
  bss = c()
  ts=c()
  for (i in 1:ncol(pred_org_surv)) {
    ts = append(ts,as.numeric(colnames(pred_org_surv)[i])) %>% try()
    bss = append(bss,BS_t(t_index = i,surv_df = pred_org_surv,valid = PainData_valid))%>% try()
  }
  IBS_dsurv_org200 = Integrated_BS(bss,ts) %>% try()
  
  bss = c()
  ts=c()
  for (i in 1:ncol(pred_selected_surv)) {
    ts = append(ts,as.numeric(colnames(pred_selected_surv)[i])) %>% try()
    bss = append(bss,BS_t(t_index = i,surv_df = pred_selected_surv,valid = selected_valid))%>% try()
  }
  IBS_dsurv_selected200 = Integrated_BS(bss,ts) %>% try()
  
  
}else {
  Cdsurv_org200 = NA
  Cdsurv_selected200 = NA
  IBS_dsurv_org200 = NA
  IBS_dsurv_selected200 = NA
}

# PermFit DeepHit

FileN = paste(Fold_Name,Folder_args,"/PermFit_dhit","_",use_seed,".rds",sep = "")

result = readRDS(file = FileN) %>% try()

if(class(result) != "try-error"){
  imp = result@importance
  
  imp_block = result@block_importance
  
  p.value.race = imp_block[1,"importance_pval_x"]
  p.value.grade = imp_block[2,"importance_pval_x"]
  p.value.surg = imp_block[3,"importance_pval_x"]
  p.value.fibrosis = imp_block[4,"importance_pval_x"]
  
  if (p.value.race > p_thershold) {imp = imp[-which(imp$var_name %in% c("race_black","race_API")),]}
  if (p.value.grade > p_thershold) {imp = imp[-which(imp$var_name %in% c("grade_2","grade_3")),]}
  if (p.value.surg > p_thershold) {imp = imp[-which(imp$var_name %in% c("surg_destruction","surg_resection")),]}
  if (p.value.fibrosis > p_thershold) {imp = imp[-which(imp$var_name %in% c("fibrosis_low","fibrosis_high")),]}
  
  
  sig_var = imp$var_name[which(imp$importance_pval_x < p_thershold)]
  
  if(length(sig_var) < 5){
    arrange = imp %>% arrange(importance_pval_x)
    sig_var = arrange$var_name[1:5]
  }
  
  
  org_model = mod_permfit(method = "DeepHit",model.type = "survival",
                          object = PainData_train,activation = "relu",
                          frac = 0.2,early_stopping = T,
                          num_nodes = c(50L, 40L, 30L,20L),epochs = 1000,
                          batch_size = 50,mod_alpha = 1,sigma = 5,cuts = 10) %>% try()
  
  selected_model = mod_permfit(method = "DeepHit",model.type = "survival",
                               object = selected_train,activation = "relu",
                               frac = 0.2,early_stopping = T,
                               num_nodes = c(50L, 40L, 30L,20L),epochs = 1000,
                               batch_size = 50,mod_alpha = 1,sigma = 5,cuts = 10) %>% try()
  
  pred_org = predict_mod_permfit(mod = org_model,
                                 method = "DeepHit",
                                 model.type = "survival",
                                 object = PainData_valid) %>% try()
  
  pred_selected = predict_mod_permfit(mod = selected_model,
                                      method = "DeepHit",
                                      model.type = "survival",
                                      object = selected_valid) %>% try()
  
  
  Cdhit_org = try(1 - rcorr.cens(pred_org[[1]],Surv(PainData_valid@y,PainData_valid@e))[1])
  Cdhit_selected = try(1 - rcorr.cens(pred_selected[[1]],Surv(selected_valid@y,selected_valid@e))[1])
  
  pred_org_surv = predict(org_model,data.frame(x = PainData_valid@x),type = "survival") %>% try()
  pred_selected_surv = predict(selected_model,data.frame(x = selected_valid@x),type = "survival") %>% try()
  
  bss = c()
  ts=c()
  for (i in 1:ncol(pred_org_surv)) {
    ts = append(ts,as.numeric(colnames(pred_org_surv)[i])) %>% try()
    bss = append(bss,BS_t(t_index = i,surv_df = pred_org_surv,valid = PainData_valid))%>% try()
  }
  IBS_dhit_org = Integrated_BS(bss,ts) %>% try()
  
  bss = c()
  ts=c()
  for (i in 1:ncol(pred_selected_surv)) {
    ts = append(ts,as.numeric(colnames(pred_selected_surv)[i])) %>% try()
    bss = append(bss,BS_t(t_index = i,surv_df = pred_selected_surv,valid = selected_valid))%>% try()
  }
  IBS_dhit_selected = Integrated_BS(bss,ts) %>% try()
  
  
}else {
  Cdhit_org = NA
  Cdhit_selected = NA
  IBS_dhit_org = NA
  IBS_dhit_selected = NA
  
}

# PermFit DeepHit batch = 200

FileN = paste(Fold_Name,Folder_args,"/PermFit_dhit","_",use_seed,"_200.rds",sep = "")

result = readRDS(file = FileN) %>% try()

if(class(result) != "try-error"){
  imp = result@importance
  
  imp_block = result@block_importance
  
  p.value.race = imp_block[1,"importance_pval_x"]
  p.value.grade = imp_block[2,"importance_pval_x"]
  p.value.surg = imp_block[3,"importance_pval_x"]
  p.value.fibrosis = imp_block[4,"importance_pval_x"]
  
  if (p.value.race > p_thershold) {imp = imp[-which(imp$var_name %in% c("race_black","race_API")),]}
  if (p.value.grade > p_thershold) {imp = imp[-which(imp$var_name %in% c("grade_2","grade_3")),]}
  if (p.value.surg > p_thershold) {imp = imp[-which(imp$var_name %in% c("surg_destruction","surg_resection")),]}
  if (p.value.fibrosis > p_thershold) {imp = imp[-which(imp$var_name %in% c("fibrosis_low","fibrosis_high")),]}
  
  
  sig_var = imp$var_name[which(imp$importance_pval_x < p_thershold)]
  
  if(length(sig_var) < 5){
    arrange = imp %>% arrange(importance_pval_x)
    sig_var = arrange$var_name[1:5]
  }
  
  
  org_model = mod_permfit(method = "DeepHit",model.type = "survival",
                          object = PainData_train,activation = "relu",
                          frac = 0.2,early_stopping = T,
                          num_nodes = c(50L, 40L, 30L,20L),epochs = 1000,
                          batch_size = 200,mod_alpha = 1,sigma = 5,cuts = 10) %>% try()
  
  selected_model = mod_permfit(method = "DeepHit",model.type = "survival",
                               object = selected_train,activation = "relu",
                               frac = 0.2,early_stopping = T,
                               num_nodes = c(50L, 40L, 30L,20L),epochs = 1000,
                               batch_size = 200,mod_alpha = 1,sigma = 5,cuts = 10) %>% try()
  
  pred_org = predict_mod_permfit(mod = org_model,
                                 method = "DeepHit",
                                 model.type = "survival",
                                 object = PainData_valid) %>% try()
  
  pred_selected = predict_mod_permfit(mod = selected_model,
                                      method = "DeepHit",
                                      model.type = "survival",
                                      object = selected_valid) %>% try()
  
  
  Cdhit_org200 = try(1 - rcorr.cens(pred_org[[1]],Surv(PainData_valid@y,PainData_valid@e))[1])
  Cdhit_selected200 = try(1 - rcorr.cens(pred_selected[[1]],Surv(selected_valid@y,selected_valid@e))[1])
  
  pred_org_surv = predict(org_model,data.frame(x = PainData_valid@x),type = "survival") %>% try()
  pred_selected_surv = predict(selected_model,data.frame(x = selected_valid@x),type = "survival") %>% try()
  
  bss = c()
  ts=c()
  for (i in 1:ncol(pred_org_surv)) {
    ts = append(ts,as.numeric(colnames(pred_org_surv)[i])) %>% try()
    bss = append(bss,BS_t(t_index = i,surv_df = pred_org_surv,valid = PainData_valid))%>% try()
  }
  IBS_dhit_org200 = Integrated_BS(bss,ts) %>% try()
  
  bss = c()
  ts=c()
  for (i in 1:ncol(pred_selected_surv)) {
    ts = append(ts,as.numeric(colnames(pred_selected_surv)[i])) %>% try()
    bss = append(bss,BS_t(t_index = i,surv_df = pred_selected_surv,valid = selected_valid))%>% try()
  }
  IBS_dhit_selected200 = Integrated_BS(bss,ts) %>% try()
  
  
}else {
  Cdhit_org200 = NA
  Cdhit_selected200 = NA
  IBS_dhit_org200 = NA
  IBS_dhit_selected200 = NA
  
}

c_all = list(Ccox_org,Ccox_selected,IBS_cox_org,IBS_cox_selected,
             Cdnnet_org,Cdnnet_selected,IBS_dnnet_org,IBS_dnnet_selected,
             Cdnnet_org200,Cdnnet_selected200,IBS_dnnet_org200,IBS_dnnet_selected200,
             Crsf_org,Crsf_selected,IBS_rsf_org,IBS_rsf_selected,
             Cxgb_org,Cxgb_selected,IBS_xgb_org,IBS_xgb_selected,
             #Csvm_org,Csvm_selected,
             Cdsurv_org,Cdsurv_selected,IBS_dsurv_org,IBS_dsurv_selected,
             Cdsurv_org200,Cdsurv_selected200,IBS_dsurv_org200,IBS_dsurv_selected200,
             Cdhit_org,Cdhit_selected,IBS_dhit_org,IBS_dhit_selected,
             Cdhit_org200,Cdhit_selected200,IBS_dhit_org200,IBS_dhit_selected200)

saveRDS(c_all,file = paste(Fold_Name,Folder_args,"/Summary_Survival_SEER_final_005_CIBS.rds",sep = ""))







