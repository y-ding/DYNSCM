## Results of synthetic control methods on tobacco and unemployed data

rm(list=ls())

# Loading R packages
library(rEDM)
library(Synth)
library(ggplot2)

source("ccm_function.R")

data_name_v = "smoking_unemploy"
suffix = "_40_1976"

data_name = paste0(data_name_v,suffix)
tr.in = 19
l = 31
num_units = 78
## Data loading

load(paste("data/",data_name,".RData", sep=""))

## Data normalization for all data instead of preintervention
smo_mat = matrix(smoking_unemploy$cigsale,nrow = 31,ncol = num_units)
smo_mat_n = scale(smo_mat)

cntl_index = c(1:2,4:num_units)
Y0 = smo_mat_n
Y1 = smo_mat_n[,3]
l=length(Y1)

k = 48
index_col = c(1:k)
index_pre=c(1:tr.in)
index_post=c((tr.in+1):31)
# preintervention Y00 : 20 X 38 matrix (20 years of smoking data for 38 control states)
Y00 = Y0[index_pre,index_col]
# postintervention Y01 : 11 X 38 matrix (20 years of smoking data for 38 control states)
Y01 = Y0[index_post,index_col]
# preintervention Y10 : 20 X 1 matrix (11 years of smoking data for 1 treated state)
Y10 = Y1[index_pre]
# postintervention Y11 : 11 X 1 matrix (11 years of smoking data for 1 treated state)
Y11 = Y1[index_post]

top_k_index = setdiff(c(1:k),3)

smoking_unemploy$cigsale[smoking_unemploy$unit.num>39]=smoking_unemploy$cigsale[smoking_unemploy$unit.num>39]
smoking_unemploy$cigsale[smoking_unemploy$year>1988 & smoking_unemploy$unit.num>39]=-50

# create matrices from panel data that provide inputs for synth()
dataprep2.out<-
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
    controls.identifier = top_k_index,
    time.predictors.prior = c(1970:1988),
    time.optimize.ssr = c(1970:1988),
    unit.names.variable = c("state"),
    time.plot = c(1970:2000)
  )

synth2.out <- synth(dataprep2.out)
## there are two ways to summarize the results
## we can either access the output from synth.out directly
round(synth2.out$solution.w,2)
# contains the unit weights or
synth2.out$solution.v

gaps2<- dataprep2.out$Y1plot-(dataprep2.out$Y0plot%*%synth2.out$solution.w)
synth2.tables <- synth.tab(dataprep.res = dataprep2.out,synth.res = synth2.out)
print(synth2.tables)
print(paste0("ATE2: ", sum(-gaps2[(tr.in+1):31])/(l-tr.in)))
print(paste0("CCM+SCM raw-all  gap norm: ", norm(gaps2,"2")))
print(paste0("CCM+SCM raw-pre  gap norm: ", norm(gaps2[1:tr.in],"2")))
print(paste0("CCM+SCM raw-post gap norm: ", norm(gaps2[(tr.in+1):31],"2")))

pdf(file="fig/jasa_figure3.pdf", width = 6, height = 5)
path.plot(synth.res = synth2.out,
          dataprep.res = dataprep2.out,
          tr.intake = 1988,
          # Main = "CCM+SCM",
          Ylab = c("per-capita cigarette sales (in packs)"),
          Xlab = c("year"),
          Legend = c("California","synthetic California"),
)
dev.off()