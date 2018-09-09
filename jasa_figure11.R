library(latex2exp)
library(ggplot2)
library(gridExtra)
library(grid)

source("grid_function.R")

idx = c(67,69)
ylim_set = cbind(rbind(0.5,0.7),rbind(0.1,0.5))
plot_list = list()
pdf(file=paste0("fig/jasa_figure11",".pdf"), width = 8, height = 4)
par(mfrow=c(1,2),mar=c(4,4,2,2), mgp = c(2.5, 1, 0))  # set up margins for plotting
j=1

# yaxis = range(min(min(rho1),min(rho2)),max(max(rho1),max(rho2)))
yaxis = range(0.2, 0.6)
for(i in idx){
  res = load(paste0("res/r",i,".RData"))
  rho1 <- pmax(0, y_xmap_x_means$rho)
  rho2 <- pmax(0, x_xmap_y_means$rho)
  # L = y_xmap_x_means$lib_size
  L = x_xmap_y_means$lib_size
  
  df=data.frame(L, rho1, rho2) 
  plot_list[[j]] = ggplot(df, aes(L))+geom_line(aes(y = rho1, colour = "CCM(X|Y)"))+
    geom_line(aes(y = rho2, colour = "CCM(Y|X)"))+ theme(legend.position = "bottom",legend.title=element_blank())+ 
    ylab(TeX('$\\rho$')) + ggtitle(paste0(i))+ ylim(0.2, 0.6)
  
  j=j+1
}

grid_arrange_shared_legend(plot_list, nrow=1,ncol = 2)


# grid.arrange(plot_list[[1]], plot_list[[2]],ncol = 2,nrow=1)
# p1 = plot_list[[1]]; p2 = plot_list[[2]]; 
# p_all <- plot_grid(p1+ theme(legend.position="none"),
#                    p2 + theme(legend.position="none"),
#                    align = 'vh',
#                    label_size=18,
#                    nrow = 1)
# print(p_all)
dev.off()

