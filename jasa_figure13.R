library(latex2exp)
library(ggplot2)
library(gridExtra)
library(grid)

source("grid_function.R")
theme_set(theme_bw()) 

pdf(file=paste0("fig/jasa_figure13.pdf"), width = 8, height = 7)
plot_list = list()

load(paste0("res/sim_beta0.RData"))
L = outX$L
pbeta0 = data.frame(L, outX$val, outY$val)
plot_list[[1]] = ggplot(pbeta0, aes(L)) +geom_line(aes(y = pbeta0$outX.val, colour = "CCM(X|Y)")) +
  geom_line(aes(y = pbeta0$outY.val, colour = "CCM(Y|X)"))+
  theme(legend.position = "right",legend.title=element_blank())+
  ylab(TeX('MAE'))+ggtitle(TeX('$\\beta$=0'))
# print(p_beta0)

load(paste0("res/sim_beta05.RData"))
L = outX$L
pbeta05 = data.frame(L, outX$val, outY$val)
plot_list[[2]] = ggplot(pbeta05, aes(L)) +geom_line(aes(y = pbeta05$outX.val, colour = "CCM(X|Y)")) +
  geom_line(aes(y = pbeta05$outY.val, colour = "CCM(Y|X)"))+
  theme(legend.position = "right",legend.title=element_blank())+
  ylab(TeX('MAE'))+ggtitle(TeX('$\\beta$=0.5'))
# print(p_beta05)

load(paste0("res/sim_beta100.RData"))
L = outX$L
pbeta100 = data.frame(L, outX$val, outY$val)
plot_list[[3]] = ggplot(pbeta100, aes(L)) +geom_line(aes(y = pbeta100$outX.val, colour = "CCM(X|Y)")) +
  geom_line(aes(y =pbeta100$outY.val, colour = "CCM(Y|X)"))+
  theme(legend.position = "right",legend.title=element_blank())+
  ylab(TeX('MAE'))+ggtitle(TeX('$\\beta$=100'))
# print(p_beta100)

load(paste0("res/sim_sigma0.RData"))
L = outX$L
psigma0 = data.frame(L, outX$val, outY$val)
plot_list[[4]] = ggplot(psigma0, aes(L)) +geom_line(aes(y = psigma0$outX.val, colour = "CCM(X|Y)")) +
  geom_line(aes(y = psigma0$outY.val, colour = "CCM(Y|X)"))+
  theme(legend.position = "right",legend.title=element_blank())+
  ylab(TeX('MAE'))+ggtitle(TeX('$\\sigma_X$=0'))
# print(p_sigma0)

load(paste0("res/sim_sigmaX100.RData"))
L = outX$L
psigmaX100 = data.frame(L, outX$val, outY$val)
plot_list[[5]] = ggplot(psigmaX100, aes(L)) +geom_line(aes(y = psigmaX100$outX.val, colour = "CCM(X|Y)")) +
  geom_line(aes(y = psigmaX100$outY.val, colour = "CCM(Y|X)"))+
  theme(legend.position = "right",legend.title=element_blank())+
  ylab(TeX('MAE'))+ggtitle(TeX('$\\sigma_X$=100'))
# print(p_sigmaX100)

load(paste0("res/sim_sigmaY0.RData"))
L = outX$L
psigmaY0 = data.frame(L, outX$val, outY$val)
plot_list[[6]] = ggplot(psigmaY0, aes(L)) +geom_line(aes(y = psigmaY0$outX.val, colour = "CCM(X|Y)")) +
  geom_line(aes(y = psigmaY0$outY.val, colour = "CCM(Y|X)"))+
  theme(legend.position = "right",legend.title=element_blank())+
  ylab(TeX('MAE'))+ggtitle(TeX('$\\sigma_Y$=0'))
# print(p_sigmaY0)

grid_arrange_shared_legend(plot_list,nrow=3,ncol=2)
dev.off()
