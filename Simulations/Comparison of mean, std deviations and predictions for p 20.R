rm(list=ls())

library(BayesLogit)
library(ggplot2)
library(gridExtra)
library(rstan)

p = 20
set.seed(9)
sim_mix <- function(n, p, m, v)
{
  k <- sample(1:length(m), n, replace = TRUE, prob = p)
  return(rnorm(n, m[k], v[k]))
}
beta_p20 <- sim_mix(p, c(0.3,0.5,0.4,0.4), c(0.6,-0.1,2.5,-2.4), c(0.7,0.7,0.7,0.7))

seeds = 1:50

for(seed in seeds){
  n = 1000
  set.seed(seed)
  X <- matrix(rep(NA, n*p), n, p)
  for(j in 1:4)
  {
    X[,j] <- rnorm(n, 0, 1)
    X[,j+4] <- runif(n, -2, 1)
    X[,j+8] <- rgamma(n, 2,2)
    X[,j+12] <- rnorm(n, 2,2)
  }
  X[, 17] = rbinom(n, 1, 0.4)
  X[, 18] = rbinom(n, 1, 0.7)
  X[,19:20] = rexp(2*n, 2) + rpois(n*2, 1)
  X[,1] = 1
  head(X)
  
  for(j in c(2:16,19:20))
  {
    X[,j] = (X[,j] - mean(X[,j]))/sd(X[,j])
  }
  
  ix = ( exp(X%*%beta_p20) > 1 ) & ( exp(X%*%beta_p20) < 200 )
  if(sum(ix) < 200) stop()
  X = X[ix,]
  
  #-------------#   n=25   #-------------#
  n = 25
  ixn = sample.int(sum(ix),n)
  X_p20_n25 <- X[ixn,]
  
  set.seed(seed)
  y_p20_n25 <- rpois(n, exp(X_p20_n25%*%beta_p20))
  
  
  #-------------#   n=50   #-------------#
  n = 50
  ixn = sample.int(sum(ix),n)
  X_p20_n50 <- X[ixn,]
  
  set.seed(seed)
  y_p20_n50 <- rpois(n, exp(X_p20_n50%*%beta_p20))
  
  
  #-------------#   n=100   #-------------#
  n = 100
  ixn = sample.int(sum(ix),n)
  X_p20_n100 <- X[ixn,]
  
  set.seed(seed)
  y_p20_n100 <- rpois(n, exp(X_p20_n100%*%beta_p20))
  
  #-------------#   n=200   #-------------#
  n = 200
  ixn = sample.int(sum(ix),n)
  X_p20_n200 <- X[ixn,]
  
  set.seed(seed)
  y_p20_n200 <- rpois(n, exp(X_p20_n200%*%beta_p20))
  
  
}


b <- rep(0,p)
B <- diag(10,p)
max_r <- 100
trunc_lambda <- 1000


########### n = 25 ###########
train_idx <- sample(1:25, size = 0.75*25)
X_p20_n25_train <- X_p20_n25[train_idx,]
y_p20_n25_train <- y_p20_n25[train_idx]

## POlya-Gamma Gibbs Sampler
gs_p20_n25 <- Gibbs_fixed_r(R = 3000, burn_in = 3000, y = y_p20_n25_train, X = X_p20_n25_train,
                            B = B, b = b,max_r = max_r, trunc_lambda = trunc_lambda)

means_gs <- colMeans(gs_p20_n25[[1]])
sd_gs <- apply(as.matrix(gs_p20_n25[[1]]), 2, sd)
preds_gs <- rowMeans(exp(X_p20_n25[-train_idx,]%*%t(as.matrix(gs_p20_n25[[1]]))))

## Expectation-Propagation
EP_p20_n25 <- getParamsEP(X = X_p20_n25_train, y = y_p20_n25_train,
                          family = "poisson", Omega0 = B, beta0 = b)

m_pred <- X_p20_n25[-train_idx,]%*%EP_p20_n25$meanBeta
s2_pred <- rowSums(X_p20_n25[-train_idx,]*(X_p20_n25[-train_idx,]%*%diag(EP_p20_n25$diagOmega)))
predMeanEP <- exp(m_pred+s2_pred/2)

## Hamiltonian Monte Carlo
data_p20_n25 <- list(
  N = 18,
  K = p,
  X = X_p20_n25_train,
  y = y_p20_n25_train
)

tmp <- stan(
  file = "C:/Users/almar/Downloads/poisson_regression.stan",
  data = data_p20_n25,
  iter = 2000,
  chains = 5,
  seed = 123
)

betaHMC_p20_n25 <- t(extract(tmp)$beta)
meanBetaHMC_p20_n25 = apply(betaHMC_p20_n25,1,mean)
sdBetaHMC_p20_n25 = apply(betaHMC_p20_n25,1,sd)
predMeanHMC_p20_n25 = rowMeans(exp(X_p20_n25[-train_idx,]%*%betaHMC_p20_n25))

meanBetaEP = EP_p20_n25$meanBeta
sdBetaEP = sqrt(EP_p20_n25$diagOmega)
EP_par <- list(meanBetaEP, sdBetaEP)
save(EP_par, file = "EP_p20_n25.RData")

## Plots on posterior means, standrd deviations and predictions
library(ggplot2)
library(reshape2)
library(scales)

### EP
colnames(meanBetaEP) = c("EP")
meanData = melt(meanBetaEP)[,-1]
meanData$group1 = c("Expectations")
meanData$y = c(meanBetaHMC_p20_n25)

sdBetaEP = as.matrix(sdBetaEP, ncols = 1)
colnames(sdBetaEP) = c("EP")
sdData = melt(sdBetaEP)[,-1]
sdData$group1 = c("Standard Deviations")
sdData$y = c(sdBetaHMC_p20_n25)

colnames(predMeanEP) = c("EP")
predData = melt(predMeanEP)[,-1]
predData$group1 = c("Predictive Means")
predData$y = c(predMeanHMC_p20_n25)

data_points = rbind(meanData,sdData,predData)
data_points$group1 = factor(data_points$group1, levels=c("Expectations", "Standard Deviations","Predictive Means" ))
P_Pois = ggplot(data_points, aes(y=value, x=y,color=c("red"))) +
  geom_point(shape=21,fill="red",alpha=0.5)+facet_wrap(group1~., scales="free")+
  geom_abline(slope = 1, intercept = 0, alpha = 0.7, lty = 2)+theme_bw()+
  labs(x="Quantities computed from the exact posterior", y = "Quantities computed from the EP posterior")+
  theme(strip.text.x = element_text(size = 15),axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15), axis.text.y = element_text(size=15),
        plot.margin = margin(0.1, 0.1, 0.05, 0.2, "cm"),legend.position = "none")+
  scale_color_manual(values=c("red"))+
  scale_shape_manual(values=c(1))+
  scale_fill_manual(values=c("red"))
windows()
show(P_Pois)

## GS
means_gs <- matrix(means_gs, ncol = 1)
colnames(means_gs) = c("GS")
meanData = melt(means_gs)[,-1]
meanData$group1 = c("Expectations")
meanData$y = c(meanBetaHMC_p20_n25)

sdBetaGS = as.matrix(sd_gs, ncols = 1)
colnames(sdBetaGS) = c("GS")
sdData = melt(sdBetaGS)[,-1]
sdData$group1 = c("Standard Deviations")
sdData$y = c(sdBetaHMC_p20_n25)

preds_gs <- as.matrix(preds_gs, ncols = 1)
colnames(preds_gs) = c("GS")
predData = melt(preds_gs)[,-1]
predData$group1 = c("Predictive Means")
predData$y = c(predMeanHMC_p20_n25)

data_points = rbind(meanData,sdData,predData)
data_points$group1 = factor(data_points$group1, levels=c("Expectations", "Standard Deviations","Predictive Means" ))
P_Pois = ggplot(data_points, aes(y=value, x=y,color=c("red"))) +
  geom_point(shape=21,fill="red",alpha=0.5)+facet_wrap(group1~., scales="free")+
  geom_abline(slope = 1, intercept = 0, alpha = 0.7, lty = 2)+theme_bw()+
  labs(x="Quantities computed from the exact posterior", y = "Quantities computed from the GS posterior")+
  theme(strip.text.x = element_text(size = 15),axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15), axis.text.y = element_text(size=15),
        plot.margin = margin(0.1, 0.1, 0.05, 0.2, "cm"),legend.position = "none")+
  scale_color_manual(values=c("red"))+
  scale_shape_manual(values=c(1))+
  scale_fill_manual(values=c("red"))
windows()
show(P_Pois)


########### n = 50 ###########
train_idx <- sample(1:50, size = 0.75*50)
X_p20_n50_train <- X_p20_n50[train_idx,]
y_p20_n50_train <- y_p20_n50[train_idx]

## Polya-Gamma Gibbs Sampler
gs_p20_n50 <- Gibbs_fixed_r(R = 3000, burn_in = 3000, y = y_p20_n50_train, X = X_p20_n50_train,
                            B = B, b = b,max_r = max_r, trunc_lambda = trunc_lambda)

means_gs <- colMeans(gs_p20_n50[[1]])
sd_gs <- apply(as.matrix(gs_p20_n50[[1]]), 2, sd)
preds_gs <- rowMeans(exp(X_p20_n50[-train_idx,]%*%t(as.matrix(gs_p20_n50[[1]]))))

## Expectation-Propagation
EP_p20_n50 <- getParamsEP(X = X_p20_n50_train, y = y_p20_n50_train,
                          family = "poisson", Omega0 = B, beta0 = b)

m_pred <- X_p20_n50[-train_idx,]%*%EP_p20_n50$meanBeta
s2_pred <- rowSums(X_p20_n50[-train_idx,]*(X_p20_n50[-train_idx,]%*%diag(EP_p20_n50$diagOmega)))
predMeanEP <- exp(m_pred+s2_pred/2)

## Hamiltonian Monte Carlo
data_p20_n50 <- list(
  N = 37,
  K = p,
  X = X_p20_n50_train,
  y = y_p20_n50_train
)

tmp <- stan(
  file = "C:/Users/almar/Downloads/poisson_regression.stan",
  data = data_p20_n50,
  iter = 2000,
  chains = 5,
  seed = 123
)

betaHMC_p20_n50 <- t(extract(tmp)$beta)
meanBetaHMC_p20_n50 = apply(betaHMC_p20_n50,1,mean)
sdBetaHMC_p20_n50 = apply(betaHMC_p20_n50,1,sd)
predMeanHMC_p20_n50 = rowMeans(exp(X_p20_n50[-train_idx,]%*%betaHMC_p20_n50))

meanBetaEP = EP_p20_n50$meanBeta
sdBetaEP = sqrt(EP_p20_n50$diagOmega)
EP_par <- list(meanBetaEP, sdBetaEP)
save(EP_par, file = "EP_p20_n50.RData")

## Plots on posterior means, standaard deviations and predictions
library(ggplot2)
library(reshape2)
library(scales)

### EP
colnames(meanBetaEP) = c("EP")
meanData = melt(meanBetaEP)[,-1]
meanData$group1 = c("Expectations")
meanData$y = c(meanBetaHMC_p20_n50)

sdBetaEP = as.matrix(sdBetaEP, ncols = 1)
colnames(sdBetaEP) = c("EP")
sdData = melt(sdBetaEP)[,-1]
sdData$group1 = c("Standard Deviations")
sdData$y = c(sdBetaHMC_p20_n50)

colnames(predMeanEP) = c("EP")
predData = melt(predMeanEP)[,-1]
predData$group1 = c("Predictive Means")
predData$y = c(predMeanHMC_p20_n50)

data_points = rbind(meanData,sdData,predData)
data_points$group1 = factor(data_points$group1, levels=c("Expectations", "Standard Deviations","Predictive Means" ))
P_Pois = ggplot(data_points, aes(y=value, x=y,color=c("red"))) +
  geom_point(shape=21,fill="red",alpha=0.5)+facet_wrap(group1~., scales="free")+
  geom_abline(slope = 1, intercept = 0, alpha = 0.7, lty = 2)+theme_bw()+
  labs(x="Quantities computed from the exact posterior", y = "Quantities computed from the EP posterior")+
  theme(strip.text.x = element_text(size = 15),axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15), axis.text.y = element_text(size=15),
        plot.margin = margin(0.1, 0.1, 0.05, 0.2, "cm"),legend.position = "none")+
  scale_color_manual(values=c("red"))+
  scale_shape_manual(values=c(1))+
  scale_fill_manual(values=c("red"))
windows()
show(P_Pois)

## GS
means_gs <- matrix(means_gs, ncol = 1)
colnames(means_gs) = c("GS")
meanData = melt(means_gs)[,-1]
meanData$group1 = c("Expectations")
meanData$y = c(meanBetaHMC_p20_n50)

sdBetaGS = as.matrix(sd_gs, ncols = 1)
colnames(sdBetaGS) = c("GS")
sdData = melt(sdBetaGS)[,-1]
sdData$group1 = c("Standard Deviations")
sdData$y = c(sdBetaHMC_p20_n50)

preds_gs <- as.matrix(preds_gs, ncols = 1)
colnames(preds_gs) = c("GS")
predData = melt(preds_gs)[,-1]
predData$group1 = c("Predictive Means")
predData$y = c(predMeanHMC_p20_n50)

data_points = rbind(meanData,sdData,predData)
data_points$group1 = factor(data_points$group1, levels=c("Expectations", "Standard Deviations","Predictive Means" ))
P_Pois = ggplot(data_points, aes(y=value, x=y,color=c("red"))) +
  geom_point(shape=21,fill="red",alpha=0.5)+facet_wrap(group1~., scales="free")+
  geom_abline(slope = 1, intercept = 0, alpha = 0.7, lty = 2)+theme_bw()+
  labs(x="Quantities computed from the exact posterior", y = "Quantities computed from the GS posterior")+
  theme(strip.text.x = element_text(size = 15),axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15), axis.text.y = element_text(size=15),
        plot.margin = margin(0.1, 0.1, 0.05, 0.2, "cm"),legend.position = "none")+
  scale_color_manual(values=c("red"))+
  scale_shape_manual(values=c(1))+
  scale_fill_manual(values=c("red"))
windows()
show(P_Pois)



########### n = 100 ###########
train_idx <- sample(1:100, size = 0.75*100)
X_p20_n100_train <- X_p20_n100[train_idx,]
y_p20_n100_train <- y_p20_n100[train_idx]

## Polya-Gamma Gibbs Sampler
gs_p20_n100 <- Gibbs_fixed_r(R = 3000, burn_in = 3000, y = y_p20_n100_train, X = X_p20_n100_train,
                             B = B, b = b,max_r = max_r, trunc_lambda = trunc_lambda)

means_gs <- colMeans(gs_p20_n100[[1]])
sd_gs <- apply(as.matrix(gs_p20_n100[[1]]), 2, sd)
preds_gs <- rowMeans(exp(X_p20_n100[-train_idx,]%*%t(as.matrix(gs_p20_n100[[1]]))))

## Expectation-Propagation
EP_p20_n100 <- getParamsEP(X = X_p20_n100_train, y = y_p20_n100_train,
                           family = "poisson", Omega0 = B, beta0 = b)

m_pred <- X_p20_n100[-train_idx,]%*%EP_p20_n100$meanBeta
s2_pred <- rowSums(X_p20_n100[-train_idx,]*(X_p20_n100[-train_idx,]%*%diag(EP_p20_n100$diagOmega)))
predMeanEP <- exp(m_pred+s2_pred/2)

## Hamiltonian Monte Carlo
data_p20_n100 <- list(
  N = 75,
  K = p,
  X = X_p20_n100_train,
  y = y_p20_n100_train
)

tmp <- stan(
  file = "C:/Users/almar/Downloads/poisson_regression.stan",
  data = data_p20_n100,
  iter = 2000,
  chains = 5,
  seed = 123
)

betaHMC_p20_n100 <- t(extract(tmp)$beta)
meanBetaHMC_p20_n100 = apply(betaHMC_p20_n100,1,mean)
sdBetaHMC_p20_n100 = apply(betaHMC_p20_n100,1,sd)
predMeanHMC_p20_n100 = rowMeans(exp(X_p20_n100[-train_idx,]%*%betaHMC_p20_n100))

meanBetaEP = EP_p20_n100$meanBeta
sdBetaEP = sqrt(EP_p20_n100$diagOmega)
EP_par <- list(meanBetaEP, sdBetaEP)
save(EP_par, file = "EP_p20_n100.RData")

## Plots on posterior means, standard deviations and predictions
library(ggplot2)
library(reshape2)
library(scales)

#### EP
colnames(meanBetaEP) = c("EP")
meanData = melt(meanBetaEP)[,-1]
meanData$group1 = c("Expectations")
meanData$y = c(meanBetaHMC_p20_n100)

sdBetaEP = as.matrix(sdBetaEP, ncols = 1)
colnames(sdBetaEP) = c("EP")
sdData = melt(sdBetaEP)[,-1]
sdData$group1 = c("Standard Deviations")
sdData$y = c(sdBetaHMC_p20_n100)

colnames(predMeanEP) = c("EP")
predData = melt(predMeanEP)[,-1]
predData$group1 = c("Predictive Means")
predData$y = c(predMeanHMC_p20_n100)

data_points = rbind(meanData,sdData,predData)
data_points$group1 = factor(data_points$group1, levels=c("Expectations", "Standard Deviations","Predictive Means" ))
P_Pois = ggplot(data_points, aes(y=value, x=y,color=c("red"))) +
  geom_point(shape=21,fill="red",alpha=0.5)+facet_wrap(group1~., scales="free")+
  geom_abline(slope = 1, intercept = 0, alpha = 0.7, lty = 2)+theme_bw()+
  labs(x="Quantities computed from the exact posterior", y = "Quantities computed from the EP posterior")+
  theme(strip.text.x = element_text(size = 15),axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15), axis.text.y = element_text(size=15),
        plot.margin = margin(0.1, 0.1, 0.05, 0.2, "cm"),legend.position = "none")+
  scale_color_manual(values=c("red"))+
  scale_shape_manual(values=c(1))+
  scale_fill_manual(values=c("red"))
windows()
show(P_Pois)

## GS
means_gs <- matrix(means_gs, ncol = 1)
colnames(means_gs) = c("GS")
meanData = melt(means_gs)[,-1]
meanData$group1 = c("Expectations")
meanData$y = c(meanBetaHMC_p20_n100)

sdBetaGS = as.matrix(sd_gs, ncols = 1)
colnames(sdBetaGS) = c("GS")
sdData = melt(sdBetaGS)[,-1]
sdData$group1 = c("Standard Deviations")
sdData$y = c(sdBetaHMC_p20_n100)

preds_gs <- as.matrix(preds_gs, ncols = 1)
colnames(preds_gs) = c("GS")
predData = melt(preds_gs)[,-1]
predData$group1 = c("Predictive Means")
predData$y = c(predMeanHMC_p20_n100)

data_points = rbind(meanData,sdData,predData)
data_points$group1 = factor(data_points$group1, levels=c("Expectations", "Standard Deviations","Predictive Means" ))
P_Pois = ggplot(data_points, aes(y=value, x=y,color=c("red"))) +
  geom_point(shape=21,fill="red",alpha=0.5)+facet_wrap(group1~., scales="free")+
  geom_abline(slope = 1, intercept = 0, alpha = 0.7, lty = 2)+theme_bw()+
  labs(x="Quantities computed from the exact posterior", y = "Quantities computed from the GS posterior")+
  theme(strip.text.x = element_text(size = 15),axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15), axis.text.y = element_text(size=15),
        plot.margin = margin(0.1, 0.1, 0.05, 0.2, "cm"),legend.position = "none")+
  scale_color_manual(values=c("red"))+
  scale_shape_manual(values=c(1))+
  scale_fill_manual(values=c("red"))
windows()
show(P_Pois)



########### n = 200 ###########
train_idx <- sample(1:200, size = 0.75*200)
X_p20_n200_train <- X_p20_n200[train_idx,]
y_p20_n200_train <- y_p20_n200[train_idx]

## Polya-Gamma Gibbs Sampler
gs_p20_n200 <- Gibbs_fixed_r(R = 3000, burn_in = 3000, y = y_p20_n200_train, X = X_p20_n200_train,
                             B = B, b = b,max_r = max_r, trunc_lambda = trunc_lambda)

means_gs <- colMeans(gs_p20_n200[[1]])
sd_gs <- apply(as.matrix(gs_p20_n200[[1]]), 2, sd)
preds_gs <- rowMeans(exp(X_p20_n200[-train_idx,]%*%t(as.matrix(gs_p20_n200[[1]]))))

## Expectation-Propagation
EP_p20_n200 <- getParamsEP(X = X_p20_n200_train, y = y_p20_n200_train,
                           family = "poisson", Omega0 = B, beta0 = b)

m_pred <- X_p20_n200[-train_idx,]%*%EP_p20_n200$meanBeta
s2_pred <- rowSums(X_p20_n200[-train_idx,]*(X_p20_n200[-train_idx,]%*%diag(EP_p20_n200$diagOmega)))
predMeanEP <- exp(m_pred+s2_pred/2)

## Hamiltonian Monte Carlo
data_p20_n200 <- list(
  N = 150,
  K = p,
  X = X_p20_n200_train,
  y = y_p20_n200_train
)

tmp <- stan(
  file = "C:/Users/almar/Downloads/poisson_regression.stan",
  data = data_p20_n200,
  iter = 2000,
  chains = 5,
  seed = 123
)

betaHMC_p20_n200 <- t(extract(tmp)$beta)
meanBetaHMC_p20_n200 = apply(betaHMC_p20_n200,1,mean)
sdBetaHMC_p20_n200 = apply(betaHMC_p20_n200,1,sd)
predMeanHMC_p20_n200 = rowMeans(exp(X_p20_n200[-train_idx,]%*%betaHMC_p20_n200))

meanBetaEP = EP_p20_n200$meanBeta
sdBetaEP = sqrt(EP_p20_n200$diagOmega)
EP_par <- list(meanBetaEP, sdBetaEP)
save(EP_par, file = "EP_p20_n200.RData")

## Plots on posterior means, standrd deviations and predictions
library(ggplot2)
library(reshape2)
library(scales)

#### EP
colnames(meanBetaEP) = c("EP")
meanData = melt(meanBetaEP)[,-1]
meanData$group1 = c("Expectations")
meanData$y = c(meanBetaHMC_p20_n200)

sdBetaEP = as.matrix(sdBetaEP, ncols = 1)
colnames(sdBetaEP) = c("EP")
sdData = melt(sdBetaEP)[,-1]
sdData$group1 = c("Standard Deviations")
sdData$y = c(sdBetaHMC_p20_n200)

colnames(predMeanEP) = c("EP")
predData = melt(predMeanEP)[,-1]
predData$group1 = c("Predictive Means")
predData$y = c(predMeanHMC_p20_n200)

data_points = rbind(meanData,sdData,predData)
data_points$group1 = factor(data_points$group1, levels=c("Expectations", "Standard Deviations","Predictive Means" ))
P_Pois = ggplot(data_points, aes(y=value, x=y,color=c("red"))) +
  geom_point(shape=21,fill="red",alpha=0.5)+facet_wrap(group1~., scales="free")+
  geom_abline(slope = 1, intercept = 0, alpha = 0.7, lty = 2)+theme_bw()+
  labs(x="Quantities computed from the exact posterior", y = "Quantities computed from the EP posterior")+
  theme(strip.text.x = element_text(size = 15),axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15), axis.text.y = element_text(size=15),
        plot.margin = margin(0.1, 0.1, 0.05, 0.2, "cm"),legend.position = "none")+
  scale_color_manual(values=c("red"))+
  scale_shape_manual(values=c(1))+
  scale_fill_manual(values=c("red"))
windows()
show(P_Pois)


## GS
means_gs <- matrix(means_gs, ncol = 1)
colnames(means_gs) = c("GS")
meanData = melt(means_gs)[,-1]
meanData$group1 = c("Expectations")
meanData$y = c(meanBetaHMC_p20_n200)

sdBetaGS = as.matrix(sd_gs, ncols = 1)
colnames(sdBetaGS) = c("GS")
sdData = melt(sdBetaGS)[,-1]
sdData$group1 = c("Standard Deviations")
sdData$y = c(sdBetaHMC_p20_n200)

preds_gs <- as.matrix(preds_gs, ncols = 1)
colnames(preds_gs) = c("GS")
predData = melt(preds_gs)[,-1]
predData$group1 = c("Predictive Means")
predData$y = c(predMeanHMC_p20_n200)

data_points = rbind(meanData,sdData,predData)
data_points$group1 = factor(data_points$group1, levels=c("Expectations", "Standard Deviations","Predictive Means" ))
P_Pois = ggplot(data_points, aes(y=value, x=y,color=c("red"))) +
  geom_point(shape=21,fill="red",alpha=0.5)+facet_wrap(group1~., scales="free")+
  geom_abline(slope = 1, intercept = 0, alpha = 0.7, lty = 2)+theme_bw()+
  labs(x="Quantities computed from the exact posterior", y = "Quantities computed from the GS posterior")+
  theme(strip.text.x = element_text(size = 15),axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15), axis.text.y = element_text(size=15),
        plot.margin = margin(0.1, 0.1, 0.05, 0.2, "cm"),legend.position = "none")+
  scale_color_manual(values=c("red"))+
  scale_shape_manual(values=c(1))+
  scale_fill_manual(values=c("red"))
windows()
show(P_Pois)
