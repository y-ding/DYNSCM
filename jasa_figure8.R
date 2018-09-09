rm(list=ls())
# Loading R packages
library(rEDM)
library(Synth)
library(astsa)
library(forecast)
library(ggplot2)
library(gridExtra)

source("ccm_function.R")

# prepare covariates for artificial states: permutation index
artcov_prep <-function(){
  state_name = unemploy$state
  new_state = paste0("f_",state_name)
  unemploy$state = new_state
  unemploy$unit.num = unemploy$unit.num + 39
  
  ## permute index for lnincome
  lnincome_index = sample(39)
  lnincome_mat_v = matrix(unemploy$lnincome,nrow = 31,ncol = 39)
  lnincome_mat = lnincome_mat_v[, lnincome_index]
  unemploy$lnincome = c(lnincome_mat)
  
  ## permute index for beer
  beer_index = sample(39)
  beer_mat_v = matrix(unemploy$beer,nrow = 31,ncol = 39)
  beer_mat = beer_mat_v[, beer_index]
  unemploy$beer = c(beer_mat)
  
  ## permute index for age15to24
  age15to24_index = sample(39)
  age15to24_mat_v = matrix(unemploy$age15to24,nrow = 31,ncol = 39)
  age15to24_mat = age15to24_mat_v[, age15to24_index]
  unemploy$age15to24 = c(age15to24_mat)
  
  ## permute index for retprice
  retprice_index = sample(39)
  retprice_mat_v = matrix(unemploy$retprice,nrow = 31,ncol = 39)
  retprice_mat = retprice_mat_v[, retprice_index]
  unemploy$retprice = c(retprice_mat)
  
  return(unemploy)
}

# prepare cigsale data for artificial states: using unemployment data
artdata_prep <- function(unemploy_mat, unemploy){
  ## add extra term to make unempolyment data scale fit cigsale data
  year = c(1:31)
  unemploy_cig = unemploy_mat[year,]
  unemploy$cigsale = c(unemploy_cig)
  smoking_unemploy = rbind(smoking,unemploy)
  return(smoking_unemploy)
}

synth_ccm_M0<-function(data){
  tr.in = 19
  num_units = 78
  ## Data normalization for all data instead of preintervention
  smo_mat = matrix(smoking_unemploy$cigsale,nrow = 31,ncol = num_units)
  smo_mat_n = scale(smo_mat)
  cntl_index = c(1:2,4:num_units)
  Y0 = smo_mat_n
  Y1 = smo_mat_n[,3]
  l=length(Y1)
  
  index_pre=c(1:tr.in)
  index_post=c((tr.in+1):31)
  # preintervention Y00 : 20 X 38 matrix (20 years of smoking data for 38 control states)
  Y00 = Y0[index_pre,]
  # postintervention Y01 : 11 X 38 matrix (20 years of smoking data for 38 control states)
  Y01 = Y0[index_post,]
  # preintervention Y10 : 20 X 1 matrix (11 years of smoking data for 1 treated state)
  Y10 = Y1[index_pre]
  # postintervention Y11 : 11 X 1 matrix (11 years of smoking data for 1 treated state)
  Y11 = Y1[index_post]
  
  ## Data normalization for preintervention data
  Y10_n = Y10 #((Y10-mean(Y10))/sd(Y10))
  Y00_n = Y00
  
  index_cut=16
  k = 78
  top_k_index = setdiff(c(1:78),3)
  non_selected = setdiff(c(1:78),top_k_index)
  
  # create matrices from panel data that provide inputs for synth()
  dataprep.out<-
    dataprep(
      foo = data,
      predictors = c("lnincome", "beer", "age15to24","retprice"),
      predictors.op = "mean",
      dependent = "cigsale",
      unit.variable = "unit.num",
      time.variable = "year",
      special.predictors = list(
        list("cigsale", 1988, "mean"),
        list("cigsale", 1980, "mean"),
        list("cigsale", 1975, "mean")
      ),
      treatment.identifier = 3,
      controls.identifier = top_k_index,
      time.predictors.prior = c(1970:1988),
      time.optimize.ssr = c(1970:1988),
      unit.names.variable = c("state"),
      time.plot = c(1970:2000)
    )
  
  synth.out <- synth(dataprep.out)
  ## there are two ways to summarize the results
  ## we can either access the output from synth.out directly
  round(synth.out$solution.w,2)
  # contains the unit weights or
  synth.out$solution.v
  
  gaps<- dataprep.out$Y1plot-(dataprep.out$Y0plot%*%synth.out$solution.w) 
  synth.tables <- synth.tab(dataprep.res = dataprep.out,synth.res = synth.out)
  print(synth.tables)
  
  # ccm plots for control with non-zero weights
  a=synth.tables$tab.w
  nonZeroIndex = a$unit.numbers[which(a$w.weights>0.002 & a$unit.numbers>39)]
  nonZeroWeightAll = a$w.weights[which(a$w.weights>0.002)]
  nonZeroWeightAdv = a$w.weights[which(a$w.weights>0.002 & a$unit.numbers>39)]
  nonZeroState = a$unit.names[which(a$w.weights>0.002 & a$unit.numbers>39)]
  # N_yc: # artificial states selected
  N_yc = length(nonZeroIndex)
  # ATE = synthCA - realCA
  ATE = -sum(gaps[(tr.in+1):l])/(l+1-tr.in)
  # AWALL: average weights for all states
  AWALL = sum(nonZeroWeightAll)/length(nonZeroWeightAll)
  # AWADV: average weights for all adversarial states
  AWADV = sum(nonZeroWeightAdv)/length(nonZeroWeightAdv)
  
  return(c(N_yc, ATE, AWALL, AWADV))
}


synth_ccm_M1<-function(data){
  tr.in = 19
  num_units = 78
  ## Data normalization for all data instead of preintervention
  smo_mat = matrix(smoking_unemploy$cigsale,nrow = 31,ncol = num_units)
  smo_mat_n = scale(smo_mat)
  
  cntl_index = c(1:2,4:num_units)
  
  Y0 = smo_mat_n
  Y1 = smo_mat_n[,3]
  l=length(Y1)
  
  index_pre=c(1:tr.in)
  index_post=c((tr.in+1):31)
  # preintervention Y00 : 20 X 38 matrix (20 years of smoking data for 38 control states)
  Y00 = Y0[index_pre,]
  # postintervention Y01 : 11 X 38 matrix (20 years of smoking data for 38 control states)
  Y01 = Y0[index_post,]
  # preintervention Y10 : 20 X 1 matrix (11 years of smoking data for 1 treated state)
  Y10 = Y1[index_pre]
  # postintervention Y11 : 11 X 1 matrix (11 years of smoking data for 1 treated state)
  Y11 = Y1[index_post]
  
  ## Data normalization for preintervention data
  Y10_n = Y10 #((Y10-mean(Y10))/sd(Y10))
  Y00_n = Y00
  index_cut=16
  k = 39
  data_name = "smoking_income_M1"
  ## Select the most informative control units by CCM. The ranking is robust.
  rank_index = ccm_all(data_name, Y10_n, Y00_n, index_cut, 3)
  print(rank_index)
  ccm_index = setdiff(rank_index, c(3)) ## fix.
  
  # create matrices from panel data that provide inputs for synth()
  dataprep.out<-
    dataprep(
      foo = smoking_unemploy,
      predictors = c("lnincome", "beer", "age15to24","retprice"),
      predictors.op = "mean",
      dependent = "cigsale",
      unit.variable = "unit.num",
      time.variable = "year",
      special.predictors = list(
        list("cigsale", 1988, "mean"),
        list("cigsale", 1980, "mean"),
        list("cigsale", 1975, "mean")
      ),
      treatment.identifier = 3,
      controls.identifier = ccm_index,
      time.predictors.prior = c(1970:1988),
      time.optimize.ssr = c(1970:1988),
      unit.names.variable = c("state"),
      time.plot = c(1970:2000)
    )
  
  synth.out <- synth(dataprep.out)
  ## there are two ways to summarize the results
  ## we can either access the output from synth.out directly
  round(synth.out$solution.w,2)
  # contains the unit weights or
  synth.out$solution.v
  
  gaps<- dataprep.out$Y1plot-(dataprep.out$Y0plot%*%synth.out$solution.w) 
  synth.tables <- synth.tab(dataprep.res = dataprep.out,synth.res = synth.out)
  print(synth.tables)
  
  # ccm plots for control with non-zero weights
  a=synth.tables$tab.w
  nonZeroIndex = a$unit.numbers[which(a$w.weights>0.002 & a$unit.numbers>39)]
  nonZeroWeightAll = a$w.weights[which(a$w.weights>0.002)]
  nonZeroWeightAdv = a$w.weights[which(a$w.weights>0.002 & a$unit.numbers>39)]
  nonZeroState = a$unit.names[which(a$w.weights>0.002 & a$unit.numbers>39)]
  # N_yc: # artificial states selected
  N_yc = length(nonZeroIndex)
  # ATE = synthCA - realCA
  ATE = -sum(gaps[(tr.in+1):l])/(l+1-tr.in)
  # AWALL: average weights for all states
  AWALL = sum(nonZeroWeightAll)/length(nonZeroWeightAll)
  # AWADV: average weights for all adversarial states
  AWADV = sum(nonZeroWeightAdv)/length(nonZeroWeightAdv)
  
  return(c(N_yc, ATE, AWALL, AWADV))
}


synth_ccm_M2<-function(data,r){
  tr.in = 19
  num_units = 78
  ## Data normalization for all data instead of preintervention
  smo_mat = matrix(smoking_unemploy$cigsale,nrow = 31,ncol = num_units)
  smo_mat_n = scale(smo_mat)
  
  cntl_index = c(1:2,4:num_units)
  
  Y0 = smo_mat_n
  Y1 = smo_mat_n[,3]
  l=length(Y1)
  
  index_pre=c(1:tr.in)
  index_post=c((tr.in+1):31)
  # preintervention Y00 : 20 X 38 matrix (20 years of smoking data for 38 control states)
  Y00 = Y0[index_pre,]
  # postintervention Y01 : 11 X 38 matrix (20 years of smoking data for 38 control states)
  Y01 = Y0[index_post,]
  # preintervention Y10 : 20 X 1 matrix (11 years of smoking data for 1 treated state)
  Y10 = Y1[index_pre]
  # postintervention Y11 : 11 X 1 matrix (11 years of smoking data for 1 treated state)
  Y11 = Y1[index_post]
  
  ## Data normalization for preintervention data
  Y10_n = Y10 #((Y10-mean(Y10))/sd(Y10))
  Y00_n = Y00
  
  index_cut=16
  data_name = paste0("smoking_income_",r)
  ## Select the most informative control units by CCM. 
  rank_index = ccm_all(data_name, Y10_n, Y00_n, index_cut, 3)
  print(rank_index)
  
  ## Step 1: get the selected index by CCM
  ccm_index = setdiff(rank_index, c(3)) ## fix.
  
  # create matrices from panel data that provide inputs for synth()
  dataprep.out<-
    dataprep(
      foo = smoking_unemploy,
      predictors = c("lnincome", "beer", "age15to24","retprice"),
      predictors.op = "mean",
      dependent = "cigsale",
      unit.variable = "unit.num",
      time.variable = "year",
      special.predictors = list(
        list("cigsale", 1988, "mean"),
        list("cigsale", 1980, "mean"),
        list("cigsale", 1975, "mean")
      ),
      treatment.identifier = 3,
      controls.identifier = c(1:2,4:78),
      time.predictors.prior = c(1970:1988),
      time.optimize.ssr = c(1970:1988),
      unit.names.variable = c("state"),
      time.plot = c(1970:2000)
    )
  synth.out <- synth(dataprep.out)
  ## there are two ways to summarize the results
  ## we can either access the output from synth.out directly
  round(synth.out$solution.w,2)
  # contains the unit weights or
  synth.out$solution.v
  gaps<- dataprep.out$Y1plot-(dataprep.out$Y0plot%*%synth.out$solution.w) 
  synth.tables <- synth.tab(dataprep.res = dataprep.out,synth.res = synth.out)
  # ccm plots for control with non-zero weights
  a=synth.tables$tab.w
  
  ## Step 2: get the selected index by SCM
  scm_index = a$unit.numbers[which(a$w.weights>0.000)]
  print(paste("scm_index",scm_index))
  print(paste("ccm_index",ccm_index))
  ## Step 3: get the new index by intersecting ccm_index and scm_index
  inter_index = intersect(ccm_index,scm_index)
  print(paste("inter_index",inter_index))
  
  if(length(inter_index)<2){
    return(c(0,0,0,0))
  }else{
    ## Step 4: run SCM again by using inter_index as control units
    dataprep.out<-
      dataprep(
        foo = smoking_unemploy,
        predictors = c("lnincome", "beer", "age15to24","retprice"),
        predictors.op = "mean",
        dependent = "cigsale",
        unit.variable = "unit.num",
        time.variable = "year",
        special.predictors = list(
          list("cigsale", 1988, "mean"),
          list("cigsale", 1980, "mean"),
          list("cigsale", 1975, "mean")
        ),
        treatment.identifier = 3,
        controls.identifier = inter_index,
        time.predictors.prior = c(1970:1988),
        time.optimize.ssr = c(1970:1988),
        unit.names.variable = c("state"),
        time.plot = c(1970:2000)
      )
    synth.out <- synth(dataprep.out)
    synth.out <- synth(dataprep.out)
    ## there are two ways to summarize the results
    ## we can either access the output from synth.out directly
    round(synth.out$solution.w,2)
    # contains the unit weights or
    synth.out$solution.v
    gaps<- dataprep.out$Y1plot-(dataprep.out$Y0plot%*%synth.out$solution.w) 
    synth.tables <- synth.tab(dataprep.res = dataprep.out,synth.res = synth.out)
    print(synth.tables)
    
    ## Do CCM plot
    a=synth.tables$tab.w
    nonZeroIndex = a$unit.numbers[which(a$w.weights>0.002 & a$unit.numbers>39)]
    nonZeroWeightAll = a$w.weights[which(a$w.weights>0.002)]
    nonZeroWeightAdv = a$w.weights[which(a$w.weights>0.002 & a$unit.numbers>39)]
    nonZeroState = a$unit.names[which(a$w.weights>0.002 & a$unit.numbers>39)]
    # N_yc: # artificial states selected
    N_yc = length(nonZeroIndex)
    # ATE = synthCA - realCA
    ATE = -sum(gaps[(tr.in+1):l])/(l+1-tr.in)
    # AWALL: average weights for all states
    AWALL = sum(nonZeroWeightAll)/length(nonZeroWeightAll)
    # AWADV: average weights for all adversarial states
    AWADV = sum(nonZeroWeightAdv)/length(nonZeroWeightAdv)
    
    return(c(N_yc, ATE, AWALL, AWADV))
  }
}
#### Generate adversarial data ####

dat01 = read.csv("data/income01.csv",header = F)
data01 = 70-(dat01$V1)*300
data1 = data01[c(TRUE, FALSE)]
x = data1[20:50]
l = length(x)
tsx = ts(x)
plot(x,type = "l",col = "blue")

## Fit time series with AR(p) model
x_fit <- arima(x, order=c(2,0,0),method="ML")
x_fitted = x - x_fit$residuals
lines(x_fitted, col="red")

n_inc = 31
incmat = matrix(0, nrow = n_inc, ncol = 39)

len=31
# pdf(file="fig/scdata_inc.pdf", width = 15, height = 5)
# par(mfrow=c(4,5),mar=c(4,4,1,1), mgp = c(2.5, 1, 0))
# for(i in c(1:39)){
#   plot(1:len, incmat[,i], main = paste0(i),type = "l", col = "blue",lwd=2, xlab = "Time Index", ylab = "Value")
# }
# dev.off()

load("data/smoking.RData")
# unemploy = smoking
# unemploy=artcov_prep()
# unemploy_mat = incmat
# save(unemploy,unemploy_mat, file = "data/income_base.RData")
load("data/income_base.RData")

r_l=seq(2,20,2)

sc=1
M0_res = list()
M1_res = list()
M2_res = list()
for(r in r_l){
  for(i in c(1:39)){
    incmat[,i] = x_fitted[1:n_inc] + r*rnorm(n_inc)
  }
  unemploy_mat = incmat
  smoking_unemploy = artdata_prep(unemploy_mat, unemploy)
  
  result0 = synth_ccm_M0(smoking_unemploy)
  M0_res = c(M0_res, result0)
  
  result1 = synth_ccm_M1(smoking_unemploy)
  M1_res = c(M1_res, result1)
  
  result2 = synth_ccm_M2(smoking_unemploy,r)
  if(sum(result2[1:2])!=0){
    M2_res = c(M2_res, result2)
  }
  
}
M0_res_mat = matrix(unlist(M0_res), ncol = length(r_l), byrow = F)
M1_res_mat = matrix(unlist(M1_res), ncol = length(r_l), byrow = F)
M2_res_mat = matrix(unlist(M2_res), ncol = length(M2_res)/4, byrow = F)


save(M0_res_mat,file = paste0("res/M0_inc_norm.RData"))
save(M1_res_mat,file = paste0("res/M1_inc_norm.RData"))
save(M2_res_mat,file = paste0("res/M2_inc_norm.RData"))


