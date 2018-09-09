rm(list=ls())

# Loading R packages
library(rEDM)
library(Synth)
library(latex2exp)
library(ggplot2)
library(gridExtra)
library(grid)
source("grid_function.R")
## Data loading
data_name = "smoking"

## Data loading
load(paste("data/",data_name,".RData", sep="")) 

## Data normalization for all data instead of preintervention
smo_mat = matrix(smoking$cigsale,nrow = 31,ncol = 39)

smo_mat_new = matrix(0,nrow = 31,ncol = 40)
smo_mat_new[,c(1:39)]=smo_mat
CT_new = smo_mat[,5]
CT_new[c(3:6)] = CT_new[c(3:6)] + 8
CT_new[c(7)] = CT_new[c(7)] + 6
smo_mat_new[,40] = CT_new

smo_mat_n = scale(smo_mat_new)

tr.in = 19
cntl_index = c(1:2,4:40)

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

index_cut = 16
em = 3

n_pre = length(Y10_n)
n_ctl = dim(Y00_n)[2]
idx_set = c(4,5,19,21,34,40)
state_set = c('Colorado','Connecticut','Montana','Nevada','Utah','Smooth Connecticut')
names(state_set) = state_set

j=1
plot_list = list()
theme_set(theme_bw()) 
pdf(file=paste("fig/jasa_figure5.pdf"), width = 8, height = 6)
op = par(mfrow=c(2,3),mar=c(3.5,3,1,1), mgp = c(2, 0.5, 0),oma=c(4,0,0,0))
for(i in idx_set){
  ctl = Y00_n[,i]
  cmp_pre = data.frame(Y10_n,ctl)
  lib_seq = c(seq(2, index_cut, by = 2))
  # N1 cross-mapping N2 (i.e. testing N2 as a cause of N1)
  y00_xmap_y10 <- ccm(cmp_pre, E=em, lib_column = "ctl", target_column = "Y10_n", 
                      lib_sizes = lib_seq, num_samples = 100,
                      random_libs = F, replace = TRUE, silent = TRUE)
  y10_xmap_y00 <- ccm(cmp_pre, E=em, lib_column = "Y10_n", target_column = "ctl",
                      lib_sizes = lib_seq, num_samples = 100, 
                      random_libs = F, replace = TRUE, silent = TRUE)
  y00_xmap_y10_means <- ccm_means(y00_xmap_y10)
  y10_xmap_y00_means <- ccm_means(y10_xmap_y00)
  rho1 <- pmax(0, y00_xmap_y10_means$rho)
  rho2 <- pmax(0, y10_xmap_y00_means$rho)
  
  L= y00_xmap_y10_means$lib_size
  df=data.frame(L, rho1, rho2) 
  plot_list[[j]] = ggplot(df, aes(L)) + geom_line(aes(y = rho1, colour = "CCM(CA|~)")) +
    geom_line(aes(y = rho2, colour = "CCM(~|CA)"))+
    theme(legend.position = "right",legend.title=element_blank())+
    ylab(TeX('$\\rho$'))+ggtitle(paste0(state_set[j]))
  j=j+1
}

grid_arrange_shared_legend(plot_list, nrow=3,ncol = 2)
dev.off()