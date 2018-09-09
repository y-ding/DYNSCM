# Compute CCM(X_t | Y_t) using autoregressive model
# X_t = alpha * X_{t-1} + mu + ex
# Y_t = beta * X_{t-1} + mu + ey
# X causes Y, need to predict X_hat from Y

rm(list=ls())
library(latex2exp)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(dplyr)

### Change parameters accordingly for different simulations
sigmax = 1
sigmay = 1
beta = 0.1
alpha = 0.9
mu = 0
X0 = 1
T = 1000 # max time.


nn = 5 # number of nearest neighbors
num_embed = nn-1 # embedding dimension


gen_data <- function()  {
  # generate time series
  X = vector(mode="numeric", length=T)
  Y = vector(mode="numeric", length=T)
  # errors
  ex =  rnorm(T, 0, sd=sigmax)
  ey = rnorm(T, 0, sd=sigmay)
  #
  for(t in 1:T) {
    X[t] = alpha * ifelse(t == 1, X0, X[t-1]) + mu + ex[t]
    Y[t] = beta * ifelse(t == 1, X0, X[t-1]) + mu + ey[t]
  }
  X = scale(X)
  Y = scale(Y)
  # Embedding (d=3)
  Y_emb = embed(Y, num_embed)
  X_emb = embed(X, num_embed)
  return(list(X=X, Y=Y, Y_emb=Y_emb, X_emb=X_emb))
}

L_vals = as.integer(seq(10, 200, by=10))
# evaluation of CCM (either MAE or Cor)
evalsX = matrix(0, nrow=0, ncol=2)
evalsY = matrix(0, nrow=0, ncol=2)
adhoc = matrix(0, nrow=0, ncol=2)
colnames(evalsX) <- c("L", "val")
colnames(evalsY) <- c("L", "val")
colnames(adhoc) <- c("L", "val")

weights = matrix(0, nrow=0, ncol=nn+1)
times = matrix(0, nrow=0, ncol=nn+2)
colnames(weights) <- c("type", paste(1:nn))
colnames(times) <- c("type", "L", paste(1:nn))

for(L in L_vals) {
  t_vals = c(T - 2* num_embed) # as.integer(seq(L+2, T-num_embed, length.out=20))
  # t_vals = seq(L+1, T-num_embed)
  nreps = 200
  
  for(irep in 1:nreps) {
    # For every L we create nreps replications.
    #
    D = gen_data()
    Y = D$Y
    X = D$X
    Y_emb = D$Y_emb
    X_emb = D$X_emb
    
    XL = X[1:L]
    YL = Y[1:L]
    Y_embL = Y_emb[1:L, ]
    X_embL = X_emb[1:L, ]
    XL_plus = X[t_vals]
    YL_plus = Y[t_vals]
    
    X_hat = c()
    Y_hat = c()
    adhoc_hat = c()
    # Compute X_hat by weights from Y
    # Do both ways X=> Y and Y=> X
    for(t in t_vals) {
      yt = Y_emb[t, ]
      xt = X_emb[t, ]
      Y_norm = apply(Y_embL, 1, function(row) sqrt(sum((row - yt)^2)))
      X_norm = apply(X_embL, 1, function(row) sqrt(sum((row - xt)^2)))
      
      indY = order(Y_norm)[1:nn] ## indices for nn
      # indY = sample(1:L, size=nn, replace = F)  ## uniform (uninformative sampling.)
      
      indX  = order(X_norm)[1:nn] ## indices for nn
      # indX  = sample(1:L, size=nn, replace = F)  ## uniform (uninformative sampling.)
      
      ## compute weights
      dy = Y_norm[indY]
      dx = X_norm[indX]
      
      wy = exp(-dy / dy[1]); wy = wy / sum(wy)
      wx = exp(-dx / dx[1]); wx = wx / sum(wx)
      
      weights = rbind(weights, c(0, wy))
      weights = rbind(weights, c(1, wx))
      times = rbind(times, c(0, L, indY))
      times = rbind(times, c(1, L, indX))
      
      # X_hat = weighted_ave(w, Xl)
      X_hat = c(X_hat, sum(wy * XL[indY]))
      Y_hat = c(Y_hat, sum(wx * YL[indX]))
      
      # adhoc estimation of Xt
      distHoc = abs(X[1:L] - X[t-1])
      k = order(distHoc)[1:3]
      # print(w)
      
      adhoc_hat = c(adhoc_hat, mean(X[k+1]))
      # print(adhoc_hat)
      
    }
    # evals X = CCM(X | Y)  and evalsY = CCM(Y | X)
    # evalsX = rbind(evalsX, c(L, cor(XL_plus, X_hat)))
    # evalsY = rbind(evalsY, c(L, cor(YL_plus, Y_hat)))
    evalsX = rbind(evalsX, c(L, mean(abs(XL_plus - X_hat))))
    evalsY = rbind(evalsY, c(L, mean(abs(YL_plus - Y_hat))))
    adhoc = rbind(adhoc, c(L, cor(XL_plus, adhoc_hat)))
    # evals = rbind(evals, c(L, mean(abs(XL_plus - X_hat))))
  }
  
  outX = as.data.frame(evalsX) %>% group_by(L) %>% summarise_all(mean)
  outY = as.data.frame(evalsY) %>% group_by(L) %>% summarise_all(mean)
  out_ad = as.data.frame(adhoc) %>% group_by(L) %>% summarise_all(mean)
  
  print("X|Y")
  print(outX[nrow(outX), ])
  print("Y|X")
  print(outY[nrow(outY), ])
  print("Ad hoc")
  print(out_ad[nrow(outY), ])
  
  outW = as.data.frame(weights) %>% group_by(type) %>% summarise_all(mean)
  print(outW)
  # outT = as.data.frame(times) %>% group_by(L, type) %>% summarise_all(mean)
  # print(outT)
  # blue = cor(Xt, Xt^)  and red = cor(Yt, Yt^)
  plot(outX$L, outX$val, type = "b",pch=18, col = "black", ylim=c(0, max(c(outX$val, outY$val))),xlab="L",  ylab="MAE",lwd=2)
  lines(outY$L, outY$val, type = "b",pch=17, col = "blue", lwd=2)
  lines(out_ad$L, out_ad$val, col="green", lwd=2)
  legend("topright", legend=c("CCM(X|Y)", "CCM(Y|X)"), col=c("black","blue"), lty=c(2,2),
         cex=0.6, lwd=c(2, 2))
}

save(outX,outY,file = paste0("res/sim_beta01.RData"))
