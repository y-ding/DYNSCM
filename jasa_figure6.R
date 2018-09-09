rm(list=ls())
library(latex2exp)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(grid)
source("grid_function.R")

load("data/smoking.RData")

# get old CT 
smo_mat = matrix(smoking$cigsale,nrow = 31,ncol = 39)
CT_old = smo_mat[,5]
smo_mat_new = matrix(0,nrow = 31,ncol = 42)
smo_mat_new[,c(1:39)]=smo_mat
smo_mat_n = scale(smo_mat_new)

CT_new = CT_old
CT_new[c(3:6)] = CT_old[c(3:6)] + 8
CT_new[c(7)] = CT_old[c(7)] + 6
smo_mat_new[,40] = CT_new

ind_set = c(3,5,40)
name_set = c("California","Connecticut","Smooth Connecticut")

theme_set(theme_bw()) 
pdf(file="fig/jasa_figure6.pdf", width = 9, height = 3)
par(mfrow=c(1,3),mar=c(4,4,1,1), mgp = c(2.5, 1, 0))
j=1
plot_list = list()
L = c(1970:2000)
for(i in ind_set){
  smoo = smo_mat_new[,i]
  df=data.frame(L, smoo) 
  plot_list[[j]] = ggplot(df, aes(L)) + geom_line(aes(y = smoo)) +
    theme(legend.position = "right",legend.title=element_blank())+ xlab("year") +
    ylab("per-capita cigarette sales (in packs)")+ggtitle(paste0(name_set[j]))
  j = j+1
}

p1 = plot_list[[1]]; p2 = plot_list[[2]]; p3 = plot_list[[3]]; 
p_all <- plot_grid(p1+ theme(legend.position="none"),
                   p2 + theme(legend.position="none"),
                   p3+ theme(legend.position="none"),
                   align = 'vh',
                   label_size=18,
                   nrow = 1
)
print(p_all)
dev.off()