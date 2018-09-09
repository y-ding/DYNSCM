## Reproduce synthetic control methods on tobacco data from Abardie et al 

rm(list=ls())

# Loading R packages
library(rEDM)
library(Synth)
library(ggplot2)

source("ccm_function.R")
## Data loading
data_name = "smoking"
tr.in = 19
l = 31
## Data loading
load(paste("data/",data_name,".RData", sep="")) 

# create matrices from panel data that provide inputs for synth()
dataprep.out<-
  dataprep(
    foo = smoking,
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
    controls.identifier = c(1:2,4:39),
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
print(paste0("ATE1: ", sum(-gaps[(tr.in+1):31])/(l-tr.in)))
print(paste0("SCM raw-all  gap norm: ", norm(gaps,"2")))
print(paste0("SCM raw-pre  gap norm: ", norm(gaps[1:tr.in],"2")))
print(paste0("SCM raw-post gap norm: ", norm(gaps[(tr.in+1):31],"2")))

pdf(file="fig/jasa_figure2.pdf", width = 6, height = 5)
path.plot(synth.res = synth.out,
          dataprep.res = dataprep.out,
          tr.intake = 1988,
          Ylab = c("per-capita cigarette sales (in packs)"),
          Xlab = c("year"),
          Legend = c("California","synthetic California"),
)
dev.off()