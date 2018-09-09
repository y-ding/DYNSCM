rm(list=ls())
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(cowplot)

# pal_style = "Set1"
pal_style = "Dark2"
wid = 6
hei = 3

########### Brexit data ########### 

load(paste0("res/M0_1.RData"))
ncountries0 = M0_res_mat[1,]
ate0 = M0_res_mat[2,]

load(paste0("res/M1_1.RData"))
ncountries1 = M1_res_mat[1,]
ate1 = M1_res_mat[2,]

ncountries = c(ncountries0,ncountries1)
ATE = c(ate0, ate1)
Call = c(rep(0, length(ncountries0)),rep(1, length(ncountries1)))

df_call=data.frame(ncountries, ATE, Call)

pdf(file=paste("fig/box_call.pdf"), width = wid, height = hei)

df_call$Call <- factor(df_call$Call, labels = c("SCM", "CCM+SCM"))
p_ncountries <- ggplot(df_call, aes(x = Call, y = ncountries, color=Call)) + geom_boxplot()+
  scale_color_brewer(palette=pal_style)

df_call$Call <- factor(df_call$Call, labels = c("SCM", "CCM+SCM"))
p_ate <- ggplot(df_call, aes(x = Call, y = ATE, color=Call)) + geom_boxplot()+
  scale_color_brewer(palette=pal_style)

# arrange the three plots in a single row
p <- plot_grid(p_ncountries + theme(legend.position="none",axis.title.x = element_blank()),
               p_ate + theme(legend.position="none",axis.title.x = element_blank()),
               align = 'vh',
               label_size=18,
               nrow = 1
)
print(p)
dev.off()