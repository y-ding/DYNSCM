# Timeseries data - Internet connections and traffic at the MPI for Intelligent Systems.
# X(t) - bytes sent at minute t.
# Y(t) - open http connections during that minute 
# ground truth:
#   Y --> X

rm(list=ls())
# Loading R packages
library(rEDM)
library(latex2exp)
library(ggplot2)
library(gridExtra)
source("grid_function.R")

idx=68

data = read.csv(paste0("data/t",idx,".csv"),header = F)
x = as.numeric(data$V1)
y = as.numeric(data$V2)
x = scale(x)
y = scale(y)
x_log = as.numeric(log(data$V1))
y_log = as.numeric(log(data$V2))
x_log = scale(x_log)
y_log = scale(y_log)

df = data.frame(x,y)
df_log = data.frame(x_log,y_log)

n = length(x)

lib_seq = c(seq(20, 340, by = 40))

y_xmap_x <- ccm(df, E=3, lib_column = "y", target_column = "x", lib_sizes = lib_seq, num_samples = 200,
                random_libs = T, replace = TRUE, silent = TRUE)
x_xmap_y <- ccm(df, E=3, lib_column = "x", target_column = "y",lib_sizes = lib_seq, num_samples = 200, 
                random_libs = T, replace = TRUE, silent = TRUE)

y_xmap_x_means <- ccm_means(y_xmap_x)
x_xmap_y_means <- ccm_means(x_xmap_y)
rho1_3 <- pmax(0, y_xmap_x_means$rho)
rho2_3 <- pmax(0, x_xmap_y_means$rho)

y_xmap_x_log <- ccm(df_log, E=3, lib_column = "y_log", target_column = "x_log", lib_sizes = lib_seq, num_samples = 200,
                random_libs = T, replace = TRUE, silent = TRUE)
x_xmap_y_log <- ccm(df_log, E=3, lib_column = "x_log", target_column = "y_log", lib_sizes = lib_seq, num_samples = 200, 
                random_libs = T, replace = TRUE, silent = TRUE)

y_xmap_x_log_means <- ccm_means(y_xmap_x_log)
x_xmap_y_log_means <- ccm_means(x_xmap_y_log)
rho1_log_3 <- pmax(0, y_xmap_x_log_means$rho)
rho2_log_3 <- pmax(0, x_xmap_y_log_means$rho)

y_xmap_x <- ccm(df, E=4, lib_column = "y", target_column = "x", lib_sizes = lib_seq, num_samples = 200,
                random_libs = T, replace = TRUE, silent = TRUE)
x_xmap_y <- ccm(df, E=4, lib_column = "x", target_column = "y",lib_sizes = lib_seq, num_samples = 200, 
                random_libs = T, replace = TRUE, silent = TRUE)

y_xmap_x_means <- ccm_means(y_xmap_x)
x_xmap_y_means <- ccm_means(x_xmap_y)
rho1_4 <- pmax(0, y_xmap_x_means$rho)
rho2_4 <- pmax(0, x_xmap_y_means$rho)

y_xmap_x_log <- ccm(df_log, E=4, lib_column = "y_log", target_column = "x_log", lib_sizes = lib_seq, num_samples = 200,
                    random_libs = T, replace = TRUE, silent = TRUE)
x_xmap_y_log <- ccm(df_log, E=4, lib_column = "x_log", target_column = "y_log", lib_sizes = lib_seq, num_samples = 200, 
                    random_libs = T, replace = TRUE, silent = TRUE)

y_xmap_x_log_means <- ccm_means(y_xmap_x_log)
x_xmap_y_log_means <- ccm_means(x_xmap_y_log)
rho1_log_4 <- pmax(0, y_xmap_x_log_means$rho)
rho2_log_4 <- pmax(0, x_xmap_y_log_means$rho)

L = y_xmap_x_means$lib_size

save(rho1_3,rho2_3,rho1_4,rho2_4,rho1_log_3,rho2_log_3,rho1_log_4,rho2_log_4,file = paste0("res/r",idx,".RData"))
load(paste0("res/r",idx,".RData"))
## Plot rho (accuracy) vs Lib size
pdf(file=paste0("fig/jasa_figure12.pdf"), width = 8, height = 7)
op = par(mfrow=c(2,2),mar=c(3.5,3,1,1), mgp = c(2, 0.5, 0),oma=c(4,0,0,0))
L= y_xmap_x_means$lib_size
plot_list = list()
df = data.frame(L,rho1_3,rho2_3,rho1_4,rho2_4,rho1_log_3,rho2_log_3,rho1_log_4,rho2_log_4)
plot_list[[1]] = ggplot(df, aes(L)) + geom_line(aes(y = rho1_3, colour = "CCM(X|Y)")) +
  geom_line(aes(y = rho2_3, colour = "CCM(Y|X)"))+theme(legend.position = "bottom",legend.title=element_blank())+
  ylab(TeX('$\\rho$'))+ggtitle("Original data, d=3")+ ylim(0.45, 0.9)

plot_list[[2]] = ggplot(df, aes(L)) + geom_line(aes(y = rho1_4, colour = "CCM(X|Y)")) +
  geom_line(aes(y = rho2_4, colour = "CCM(Y|X)"))+theme(legend.position = "bottom",legend.title=element_blank())+
  ylab(TeX('$\\rho$'))+ggtitle("Original data, d=4")+ ylim(0.45, 0.9)

plot_list[[3]] = ggplot(df, aes(L)) + geom_line(aes(y = rho1_log_3, colour = "CCM(X|Y)")) +
  geom_line(aes(y = rho2_log_3, colour = "CCM(Y|X)"))+theme(legend.position = "bottom",legend.title=element_blank())+
  ylab(TeX('$\\rho$'))+ggtitle("Transformed data, d=4")+ ylim(0.45, 0.9)

plot_list[[4]] = ggplot(df, aes(L)) + geom_line(aes(y = rho1_log_4, colour = "CCM(X|Y)")) +
  geom_line(aes(y = rho2_log_4, colour = "CCM(Y|X)"))+theme(legend.position = "bottom",legend.title=element_blank())+
  ylab(TeX('$\\rho$'))+ggtitle("Transformed data, d=4")+ ylim(0.45, 0.9)

grid_arrange_shared_legend(plot_list,nrow=2,ncol=2)
dev.off()
