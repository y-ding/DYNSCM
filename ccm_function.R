ccm_all <- function(data_name, Y10_n, Y00_n, index_cut, em){
  ## Simplex projection for control and treatment to select optimal E
  n_pre = length(Y10_n)
  n_ctl = dim(Y00_n)[2]
  # n_ctl = 20
  rho1_vec = rep(0, n_ctl)
  rho2_vec = rep(0, n_ctl)
  rho_vec = rep(0, n_ctl)
  gap_vec = rep(0, n_ctl)
  filtered_index = list()

  pdf(file=paste("fig/ccm_vs_L_",data_name,".pdf", sep=""), width = 15, height = 5)
  par(mfrow=c(4,5),mar=c(4,4,1,1), mgp = c(2.5, 1, 0))  # set up margins for plotting
  
  for(i in c(1:n_ctl)){
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
    
    ## Plot rho (accuracy) vs Lib size
    plot(y00_xmap_y10_means$lib_size, main=paste(i), 
         rho1, type = "l", col = "red",lwd=2, xlab = "Library Size", ylab = "Cross Map Skill", ylim = c(0, 1))
    lines(y10_xmap_y00_means$lib_size, rho2, col = "blue",lwd=2)
    legend(x = "topleft", legend = c("control_xmap_treatment", "treatment_xmap_control"), col = c("red", "blue"), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)
    rho1[is.na(rho1)] = 0
    rho2[is.na(rho2)] = 0
    rho1_vec[i] = max(rho1)
    rho2_vec[i] = max(rho2)
    rho_vec[i] = max(rho1+rho2)
    gap_vec[i] = max(abs(rho1-rho2))
    # gap_vec[i]<0.18 for most cases
    if(rho1_vec[i]>0.5 & rho2_vec[i]>0.5 & gap_vec[i]<0.18){
      filtered_index=c(filtered_index, i)
    }
  }
  dev.off()
  final_index = unlist(filtered_index, use.names=FALSE) 
  print(final_index)
  return(final_index)
}

ccm_plot <- function(data_name, Y10_n, Y00_n, index_cut, em, nonZeroWeight, nonZeroIndex){
  ## Simplex projection for control and treatment to select optimal E
  n_pre = length(Y10_n)
  n_ctl = dim(Y00_n)[2]
  
  pdf(file=paste("fig/ccm_vs_L_",data_name,".pdf", sep=""), width = 10, height = 2)
  par(mfrow=c(1,6),mar=c(4,4,1,1), mgp = c(2.5, 1, 0))  # set up margins for plotting
  # par(xpd=T, mar=par()$mar+c(0,0,0,7))
  par(oma = c(1, 1, 1, 1))
  j=1
  for(i in nonZeroIndex){
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
    
    ## Plot rho (accuracy) vs Lib size
    plot(y00_xmap_y10_means$lib_size, font.main = 1, main=paste0(i, "(", nonZeroWeight[j],")"), 
         rho1, type = "l", col = "red",lwd=2, xlab = "Library Size", ylab = "Cross Map Skill", ylim = c(0, 1))
    lines(y10_xmap_y00_means$lib_size, rho2, col = "blue",lwd=2)
    legend(x = "topleft", legend = c("c_xmap_t", "t_xmap_c"),
           col = c("red", "blue"), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)
    j=j+1
  }

  dev.off()
}

