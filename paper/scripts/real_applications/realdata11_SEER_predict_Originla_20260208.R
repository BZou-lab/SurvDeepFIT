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

############### Original Cox ###########################

FileN = paste(Fold_Name,Folder_args,"/Org_cox","_",use_seed,".rds",sep = "")

result = readRDS(file = FileN) %>% try()

if(class(result) != "try-error"){
  vimp = summary(result)$coef
  sig_var = rownames(vimp)[which(vimp[,"Pr(>|z|)"] < 0.05)]
  
  if(length(sig_var) < 5){
    arrange = vimp[order(vimp[,"Pr(>|z|)"]),]
    sig_var = rownames(arrange)[1:5]
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
}

################ Original RSF ###########################
FileN = paste(Fold_Name,Folder_args,"/Org_RSF","_",use_seed,".rds",sep = "")

result = readRDS(file = FileN) %>% try()

if(class(result) != "try-error"){
  vimp = result[["var.sel"]]
  
  sig_var = rownames(vimp)[which(vimp[,"pvalue"] < 0.05)]
  sig_var <- sub("^x\\.", "", sig_var)
  
  if(length(sig_var) < 5){
    arrange = vimp[order(vimp[,"pvalue"]),]
    sig_var = rownames(arrange)[1:5]
    sig_var <- sub("^x\\.", "", sig_var)
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
}


########## Original XGB

FileN = paste(Fold_Name,Folder_args,"/Org_XGB","_",use_seed,".rds",sep = "")

result = readRDS(file = FileN) %>% try()

if(class(result)[1] != "try-error"){
  vimp = result
  arrange = vimp[order(vimp[,"Gain"],decreasing = T),]
  sig_var = arrange$Feature[1:10]
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
}


########## Original AFT

FileN = paste(Fold_Name,Folder_args,"/Org_AFT","_",use_seed,".rds",sep = "")

result = readRDS(file = FileN) %>% try()


if(class(result) != "try-error"){
  vimp = summary(result)$table
  vimp = vimp[2:(nrow(vimp)-1),]
  sig_var = rownames(vimp)[which(vimp[,"p"] < 0.05)]
  
  if(length(sig_var) < 5){
    arrange = vimp[order(vimp[,"p"]),]
    sig_var = rownames(arrange)[1:5]
  }
  
  
  selected_train = importDnnetSurv(x = Ptrain[,sig_var],
                                   y = Ptrain[,"OS_MONTHS"],
                                   e = Ptrain[,"OS_STATUS"])
  selected_valid = importDnnetSurv(x = Pvalid[,sig_var],
                                   y = Pvalid[,"OS_MONTHS"],
                                   e = Pvalid[,"OS_STATUS"])
  
  mod_permfit_logisticaft = function(object){
    df_train = cbind(object@y, object@e, object@x) %>% 
      as.data.frame()
    colnames(df_train) = c("EventTime", "EventStatus", 
                           colnames(object@x))
    mod <- survival::survreg(Surv(EventTime, EventStatus) ~ 
                               ., data = df_train, dist = "logistic")
  }
  
  org_model = mod_permfit_logisticaft(object = PainData_train) %>% try()
  
  selected_model = mod_permfit_logisticaft(object = selected_train)
  
  # ---- Function 1: C-index ----
  calc_cindex <- function(model, object) {
    df_valid <- cbind(object@y, object@e, object@x) %>% as.data.frame()
    colnames(df_valid) <- c("EventTime", "EventStatus", colnames(object@x))
    
    lp <- predict(model, newdata = df_valid, type = "lp")
    
    conc <- survival::concordance(Surv(EventTime, EventStatus) ~ I(-lp), data = df_valid)
    
    return(1 - conc$concordance)
  }
  
  # ---- Function 2: Integrated Brier Score ----
  calc_ibs <- function(model, object) {
    df_valid <- cbind(object@y, object@e, object@x) %>% as.data.frame()
    colnames(df_valid) <- c("EventTime", "EventStatus", colnames(object@x))
    df_valid <- df_valid[df_valid$EventTime > 0, ]
    
    lp <- predict(model, newdata = df_valid, type = "lp")
    scale <- model$scale
    pfun <- switch(model$dist,
                   "logistic" = plogis, "lognormal" = pnorm, "gaussian" = pnorm,
                   "weibull" = function(q) 1 - exp(-exp(q)),
                   "exponential" = function(q) 1 - exp(-exp(q)))
    
    # Restrict to 10th-90th percentile of event times
    event_times <- df_valid$EventTime[df_valid$EventStatus == 1]
    t_low <- quantile(event_times, 0.10)
    t_high <- quantile(event_times, 0.90)
    times_eval <- seq(t_low, t_high, length.out = 100)
    
    sp_matrix <- sapply(times_eval, function(t) pmax(0, pmin(1, 1 - pfun((log(t) - lp) / scale))))
    
    n <- nrow(df_valid)
    cens_fit <- survfit(Surv(EventTime, 1 - EventStatus) ~ 1, data = df_valid)
    G <- stepfun(cens_fit$time, c(1, cens_fit$surv))
    
    bs_vec <- sapply(seq_along(times_eval), function(j) {
      t <- times_eval[j]
      s <- sp_matrix[, j]
      w1 <- (df_valid$EventTime <= t & df_valid$EventStatus == 1)
      w2 <- (df_valid$EventTime > t)
      sum(s[w1]^2 / pmax(G(df_valid$EventTime[w1]), 1e-6)) / n +
        sum((1 - s[w2])^2 / pmax(G(t), 1e-6)) / n
    })
    
    dt <- diff(times_eval)
    sum((bs_vec[-length(bs_vec)] + bs_vec[-1]) / 2 * dt) / (t_high - t_low)
  }
  
  
  # For selected_model
  Caft_org = calc_cindex(model = selected_model, object = selected_valid)
  IBS_aft_selected = calc_ibs(model = selected_model, object = selected_valid)
  
  # For org_model — define org_valid similarly
  Caft_selected = calc_cindex(model = org_model, object = PainData_valid)
  IBS_aft_org = calc_ibs(model = org_model, object = PainData_valid)
  
}

c_all = list(Ccox_org,Ccox_selected,IBS_cox_org,IBS_cox_selected,
             Crsf_org,Crsf_selected,IBS_rsf_org,IBS_rsf_selected,
             Cxgb_org,Cxgb_selected,IBS_xgb_org,IBS_xgb_selected,
             #Csvm_org,Csvm_selected,
             Caft_org,Caft_selected,IBS_aft_org,IBS_aft_selected)

saveRDS(c_all,file = paste(Fold_Name,Folder_args,"/Summary_Survival_SEER_ORG_005_CIBS.rds",sep = ""))

