rm(list=ls())
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(cowplot)
source("grid_function.R")
# plot_list = list()

# pal_style = "Set1"
pal_style = "Dark2"
wid = 6
hei = 3

########### Unemployment ########### 
i = 6
suffix0 = paste0(0,"_",i)
load(paste0("res/M",suffix0,".RData"))
nstates0 = M0_res_mat[1,]
ate0 = M0_res_mat[2,]
suffix1 = paste0(1,"_",i)
load(paste0("res/M",suffix1,".RData"))
nstates1 = M1_res_mat[1,]
ate1 = M1_res_mat[2,]

nstates = c(nstates0,nstates1)
ATE = c(ate0, ate1)
Unemployment = c(rep(0, length(nstates0)),rep(1, length(nstates1)))

df_unemploy=data.frame(nstates, ATE, Unemployment)

pdf(file=paste("fig/box_unemploy.pdf"), width = wid, height = hei)

df_unemploy$Unemployment <- factor(df_unemploy$Unemployment, labels = c("SCM", "CCM+SCM"))
p_nstates <- ggplot(df_unemploy, aes(x = Unemployment, y = nstates, color=Unemployment)) + geom_boxplot()+
  scale_color_brewer(palette=pal_style)

df_unemploy$Unemployment <- factor(df_unemploy$Unemployment, labels = c("SCM", "CCM+SCM"))
p_ate <- ggplot(df_unemploy, aes(x = Unemployment, y = ATE, color=Unemployment)) + geom_boxplot()+
  scale_color_brewer(palette=pal_style)

# arrange the three plots in a single row
p <- plot_grid(p_nstates + theme(legend.position="none",axis.title.x = element_blank()),
                      p_ate + theme(legend.position="none",axis.title.x = element_blank()),
                      align = 'vh',
                      label_size=18,
                      nrow = 1
)
print(p)
dev.off()

########### Income ########### 
load(paste0("res/M0_inc_norm.RData"))
nstates0 = M0_res_mat[1,]
ate0 = M0_res_mat[2,]
load(paste0("res/M1_inc_norm.RData"))
nstates1 = M1_res_mat[1,]
ate1 = M1_res_mat[2,]

nstates = c(nstates0,nstates1)
ATE = c(ate0, ate1)
Income = c(rep(0, length(nstates0)),rep(1, length(nstates1)))

df_income=data.frame(nstates, ATE, Income)

pdf(file=paste("fig/box_income.pdf"), width = wid, height = hei)

df_income$Income <- factor(df_income$Income, labels = c("SCM", "CCM+SCM"))
p_nstates <- ggplot(df_income, aes(x = Income, y = nstates, color=Income)) + geom_boxplot()+
  scale_color_brewer(palette=pal_style)

df_income$Income <- factor(df_income$Income, labels = c("SCM", "CCM+SCM"))
p_ate <- ggplot(df_income, aes(x = Income, y = ATE, color=Income)) + geom_boxplot()+
  scale_color_brewer(palette=pal_style)

# arrange the three plots in a single row
p <- plot_grid(p_nstates + theme(legend.position="none",axis.title.x = element_blank()),
               p_ate + theme(legend.position="none",axis.title.x = element_blank()),
               align = 'vh',
               # label_size=18,
               nrow = 1
)
print(p)
dev.off()

########### Downshift ########### 
load(paste0("res/M0_art_44.RData"))
nstates0 = M0_res_mat[1,]
ate0 = M0_res_mat[2,]
load(paste0("res/M1_art_44.RData"))
nstates1 = M1_res_mat[1,]
ate1 = M1_res_mat[2,]

nstates = c(nstates0,nstates1)
ATE = c(ate0, ate1)
Downshift = c(rep(0, length(nstates0)),rep(1, length(nstates1)))

df_downshift=data.frame(nstates, ATE, Downshift)

pdf(file=paste("fig/box_downshift.pdf"), width = wid, height = hei)

df_downshift$Downshift <- factor(df_downshift$Downshift, labels = c("SCM", "CCM+SCM"))
p_nstates <- ggplot(df_downshift, aes(x = Downshift, y = nstates, color=Downshift)) + geom_boxplot()+
  scale_color_brewer(palette=pal_style)

df_downshift$Downshift <- factor(df_downshift$Downshift, labels = c("SCM", "CCM+SCM"))
p_ate <- ggplot(df_downshift, aes(x = Downshift, y = ATE, color=Downshift)) + geom_boxplot()+
  scale_color_brewer(palette=pal_style)

# arrange the three plots in a single row
p <- plot_grid(p_nstates + theme(legend.position="none",axis.title.x = element_blank()),
               p_ate + theme(legend.position="none",axis.title.x = element_blank()),
               align = 'vh',
               # label_size=18,
               nrow = 1
)
print(p)
dev.off()

