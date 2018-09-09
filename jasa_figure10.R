rm(list=ls())
# Loading R packages
library(rEDM)
library(Synth)
library(astsa)
source("ccm_function.R")

# prepare covariates for artificial states: permutation index
artcov_prep <-function(){
  country_name = advers$country
  new_country = paste0("f_",country_name)
  advers$country = new_country
  advers$unit.num = advers$unit.num + 31
  
  ## permute index for consumption
  consumption_index = sample(31)
  consumption_mat_v = matrix(advers$consumption,nrow = 96,ncol = 31)
  consumption_mat = consumption_mat_v[, consumption_index]
  advers$consumption = c(consumption_mat)
  
  ## permute index for investment
  investment_index = sample(31)
  investment_mat_v = matrix(advers$investment,nrow = 96,ncol = 31)
  investment_mat = investment_mat_v[, investment_index]
  advers$investment = c(investment_mat)
  
  ## permute index for export
  export_index = sample(31)
  export_mat_v = matrix(advers$export,nrow = 96,ncol = 31)
  export_mat = export_mat_v[, export_index]
  advers$export = c(export_mat)
  
  ## permute index for short_interest
  short_interest_index = sample(31)
  short_interest_mat_v = matrix(advers$short_interest,nrow = 96,ncol = 31)
  short_interest_mat = short_interest_mat_v[, short_interest_index]
  advers$short_interest = c(short_interest_mat)
  
  ## permute index for exchange
  exchange_index = sample(31)
  exchange_mat_v = matrix(advers$exchange,nrow = 96,ncol = 31)
  exchange_mat = exchange_mat_v[, exchange_index]
  advers$exchange = c(exchange_mat)
  
  ## permute index for inflation
  inflation_index = sample(31)
  inflation_mat_v = matrix(advers$inflation,nrow = 96,ncol = 31)
  inflation_mat = inflation_mat_v[, inflation_index]
  advers$inflation = c(inflation_mat)
  
  return(advers)
}

# prepare cigsale data for artificial countries: using daily call data
artdata_prep <- function(advers_mat, advers, sc, k){
  ## add extra term to make daily call data scale fit gdp data
  year = c(1:96)+k
  advers_gdp = advers_mat[year,]*sc
  advers$gdp = c(advers_gdp)
  gdp_advers = rbind(gdp, advers)
  return(gdp_advers)
}

synth_ccm_M0<-function(data){
  ## 2005Q1: 40 ; 2008Q1: 52 ; 2009Q3: 58 ; 2012Q1: 68
  tr.in = 85
  l=96
  controls = c(1:29,31:62)
  
  # create matrices from panel data that provide inputs for synth()
  dataprep.out<-
    dataprep(
      foo = gdp_advers,
      predictors = c("consumption","investment","export","short_interest","exchange","inflation"),
      predictors.op = "mean",
      dependent = "gdp",
      unit.variable = "unit.num",
      time.variable = "quarter_index",
      special.predictors = list(
        list("gdp", 1, "mean"),
        list("gdp", 20, "mean"),
        list("gdp", 40, "mean"),
        list("gdp", 60, "mean"),
        list("gdp", 65, "mean"),
        list("gdp", 70, "mean"),
        list("gdp", 80, "mean")
      ),
      treatment.identifier = 30,
      controls.identifier = controls,
      time.predictors.prior = c(1:tr.in),
      time.optimize.ssr = c(1:tr.in),
      unit.names.variable = c("country"),
      time.plot = c(1:96)
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
  nonZeroIndex = a$unit.numbers[which(a$w.weights>0.002 & a$unit.numbers>31)]
  nonZeroWeightAll = a$w.weights[which(a$w.weights>0.002)]
  nonZeroWeightAdv = a$w.weights[which(a$w.weights>0.002 & a$unit.numbers>31)]
  nonZeroState = a$unit.names[which(a$w.weights>0.002 & a$unit.numbers>31)]
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
  ## 2005Q1: 40 ; 2008Q1: 52 ; 2009Q3: 58 ; 2012Q1: 68
  tr.in = 85
  num_units = 62
  
  ## Data normalization for all data instead of preintervention
  gdp_mat = matrix(gdp_advers$gdp,nrow = 96,ncol = num_units)
  gdp_mat_n = scale(gdp_mat)
  
  cntl_index = c(1:29,31:num_units)
  
  Y0 = gdp_mat_n
  Y1 = gdp_mat_n[,30]
  l=length(Y1)
  
  index_pre=c(1:tr.in)
  index_post=c((tr.in+1):96)
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
  
  ## Simplex projection for control and treatment to select optimal E
  index_cut = 50
  data_name = "gdp_advers_M1"
  ## Select the most informative control units for MTV. The ranking is robust.
  rank_index = ccm_all(data_name, Y10_n, Y00_n, index_cut, 3)
  print(rank_index)
  ccm_index = setdiff(rank_index, c(30)) ## fix.
  
  # create matrices from panel data that provide inputs for synth()
  dataprep.out<-
    dataprep(
      foo = gdp_advers,
      predictors = c("consumption","investment","export","short_interest","exchange","inflation"),
      predictors.op = "mean",
      dependent = "gdp",
      unit.variable = "unit.num",
      time.variable = "quarter_index",
      special.predictors = list(
        list("gdp", 1, "mean"),
        list("gdp", 20, "mean"),
        list("gdp", 40, "mean"),
        list("gdp", 60, "mean"),
        list("gdp", 65, "mean"),
        list("gdp", 70, "mean"),
        list("gdp", 80, "mean")
      ),
      treatment.identifier = 30,
      controls.identifier = ccm_index,
      time.predictors.prior = c(1:tr.in),
      time.optimize.ssr = c(1:tr.in),
      unit.names.variable = c("country"),
      time.plot = c(1:96)
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
  nonZeroIndex = a$unit.numbers[which(a$w.weights>0.002 & a$unit.numbers>31)]
  nonZeroWeightAll = a$w.weights[which(a$w.weights>0.002)]
  nonZeroWeightAdv = a$w.weights[which(a$w.weights>0.002 & a$unit.numbers>31)]
  nonZeroState = a$unit.names[which(a$w.weights>0.002 & a$unit.numbers>31)]
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

synth_ccm_M2<-function(data){
  ## 2005Q1: 40 ; 2008Q1: 52 ; 2009Q3: 58 ; 2012Q1: 68
  tr.in = 85
  num_units = 62
  controls = c(1:29,31:62)
  ## Data normalization for all data instead of preintervention
  gdp_mat = matrix(gdp_advers$gdp,nrow = 96,ncol = num_units)
  gdp_mat_n = scale(gdp_mat)
  
  cntl_index = c(1:29,31:num_units)
  
  Y0 = gdp_mat_n
  Y1 = gdp_mat_n[,30]
  l=length(Y1)
  
  index_pre=c(1:tr.in)
  index_post=c((tr.in+1):96)
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
  
  ## Simplex projection for control and treatment to select optimal E
  index_cut = 50
  data_name = "gdp_advers_M1"
  ## Select the most informative control units for MTV. The ranking is robust.
  rank_index = ccm_all(data_name, Y10_n, Y00_n, index_cut, 3)
  print(rank_index)
  
  ## Step 1: get the selected index by CCM
  ccm_index = setdiff(rank_index, c(30)) ## fix.
  
  # create matrices from panel data that provide inputs for synth()
  dataprep.out<-
    dataprep(
      foo = gdp_advers,
      predictors = c("consumption","investment","export","short_interest","exchange","inflation"),
      predictors.op = "mean",
      dependent = "gdp",
      unit.variable = "unit.num",
      time.variable = "quarter_index",
      special.predictors = list(
        list("gdp", 1, "mean"),
        list("gdp", 20, "mean"),
        list("gdp", 40, "mean"),
        list("gdp", 60, "mean"),
        list("gdp", 65, "mean"),
        list("gdp", 70, "mean"),
        list("gdp", 80, "mean")
      ),
      treatment.identifier = 30,
      controls.identifier = controls,
      time.predictors.prior = c(1:tr.in),
      time.optimize.ssr = c(1:tr.in),
      unit.names.variable = c("country"),
      time.plot = c(1:96)
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
        foo = gdp_advers,
        predictors = c("consumption","investment","export","short_interest","exchange","inflation"),
        predictors.op = "mean",
        dependent = "gdp",
        unit.variable = "unit.num",
        time.variable = "quarter_index",
        special.predictors = list(
          list("gdp", 1, "mean"),
          list("gdp", 20, "mean"),
          list("gdp", 40, "mean"),
          list("gdp", 60, "mean"),
          list("gdp", 65, "mean"),
          list("gdp", 70, "mean"),
          list("gdp", 80, "mean")
        ),
        treatment.identifier = 30,
        controls.identifier = inter_index,
        time.predictors.prior = c(1:tr.in),
        time.optimize.ssr = c(1:tr.in),
        unit.names.variable = c("country"),
        time.plot = c(1:96)
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
    nonZeroIndex = a$unit.numbers[which(a$w.weights>0.002 & a$unit.numbers>31)]
    nonZeroWeightAll = a$w.weights[which(a$w.weights>0.002)]
    nonZeroWeightAdv = a$w.weights[which(a$w.weights>0.002 & a$unit.numbers>31)]
    nonZeroState = a$unit.names[which(a$w.weights>0.002 & a$unit.numbers>31)]
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


load("data/brexit.RData")
advers = gdp
advers=artcov_prep()
# dat = read.csv("data/daily-calls.csv",header = F)
# dat = dat$V1
# # dat = dat[25:131]/dat[1]
# dat = dat[25:130]
# dat = dat/dat[1]
# log_dat = log(dat) * 50
# log_dat[55:73] = log_dat[55:73] + 4
# # plot(log_dat, type = "l", col="red")
# l = length(log_dat)
# x = log_dat
# tsx = ts(x)
# ## Fit time series with AR(p) model
# x_fit <- sarima(x, p = 2, d = 1, q = 1)
# x_fitted = x - x_fit$fit$residuals
# # plot(tsx,col="green")
# # lines(c(1:l),x_fitted, col="red")
# n_art = 106
# advers_mat = matrix(0, nrow = n_art, ncol = 31)
# for(i in c(1:31)){
#   advers_mat[,i] = x_fitted[1:n_art] + 2*rnorm(n_art)
# }
# len=106
# pdf(file="fig/scdata_art.pdf", width = 15, height = 5)
# par(mfrow=c(4,5),mar=c(4,4,1,1), mgp = c(2.5, 1, 0))
# for(i in c(1:31)){
#   plot(1:len, advers_mat[,i], main = paste0(i),type = "l", col = "blue",lwd=2, xlab = "Time Index", ylab = "Value")
# }
# dev.off()

# save(advers,advers_mat, file = "data/call_base.RData")

load("data/call_base.RData")

y_l = c(0:10)
sc_l = c(1)

for(sc in sc_l){
  M0_res = list()
  for(y in y_l){
    gdp_advers = artdata_prep(advers_mat, advers, sc, y)
    result = synth_ccm_M0(gdp_advers)
    M0_res = c(M0_res, result)
  }
  M0_res_mat = matrix(unlist(M0_res), ncol = length(y_l), byrow = F)
  save(M0_res_mat,file = paste0("res/M0_",sc,".RData"))
}

for(sc in sc_l){
  M1_res = list()
  for(y in y_l){
    gdp_advers = artdata_prep(advers_mat, advers, sc, y)
    result = synth_ccm_M1(gdp_advers)
    M1_res = c(M1_res, result)
  }
  M1_res_mat = matrix(unlist(M1_res), ncol = length(y_l), byrow = F)
  save(M1_res_mat,file = paste0("res/M1_",sc,".RData"))
}

for(sc in sc_l){
  M2_res = list()
  for(y in y_l){
    gdp_advers = artdata_prep(advers_mat, advers, sc, y)
    result = synth_ccm_M2(gdp_advers)
    if(sum(result[1:2])!=0){
      M2_res = c(M2_res, result)
    }
  }
  M2_res_mat = matrix(unlist(M2_res), ncol = length(y_l), byrow = F)
  save(M2_res_mat,file = paste0("res/M2_",sc,".RData"))
}


