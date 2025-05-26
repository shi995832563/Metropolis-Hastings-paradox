#Pois-lognormal

#sample lambda's and sigma2
#weak informative prior--IG(epsilon)
set.seed(123)
rrr <- 1000
mse_estimates_map <- matrix(0,rrr,3)
mse_mean1 <- mse_mean2 <- mse_mean3 <- matrix(0,rrr,3)
mse_mode1 <- mse_mode2 <- mse_mode3 <- matrix(0,rrr,3)

library(mvtnorm)
library(MASS)
library(invgamma)
library(LaplacesDemon)

Modes <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}
Iter = 10000

#<1----------------have a look at the data

logLike_Lambda <- function(lambda, mu, sigma2){
  return(log(dnorm(lambda, mu, sqrt(sigma2))))
}

logLike_Y <- function(Y, lambda, tau2){
  ll <- sum(log(dnorm(Y, lambda,sqrt(tau2))))  
  return(ll)
}

logLike_L <- function(Y, mu, sigma2, tau2){
  ll <- sum(log(dnorm(Y, mu, sqrt(tau2+sigma2))))  
  return(ll)
}

tau2 <- 0.1
#sim function: simulated dataset Y, X, Z, Lambda
sim <- function(N=100,k=100){
  Y <- NULL
  mu <- 2
  sigma2 <- 0.1      #sigma^2 = 0.1  
  tau2 <- 0.1
  m <- N/k
  
  Lambda_group <- rnorm(k, mu , sqrt(sigma2))
  Lambda <- NULL
  for (i in 1:k) {
    Lambda <- c(Lambda,rep(Lambda_group[i],m))
    Y[(m*(i-1)+1):(m*(i-1)+m)] <- rnorm(m,Lambda_group[i],sqrt(tau2))
  }
  return(list(Y = Y, Lambda = Lambda, parameters = list(mu = mu, sigma2 = sigma2)))
}

##########################################################################
for (ii in 1:rrr) {
#simulated data
N = 200
k <- N#200#2
m <- N/k
res <- sim(N,k)

Y <- res$Y

#----------------2>
#prior for mu: normal distribution
#mu_0 = 2.0
#mu_0 = res$parameters$mu
#sd_mu = 0.50

#Prior for sigma2: inverse gamma distribution
#mean: scale/(shape-1), var: scale^2/(shape-1)^2(shape-2)
gam_shape <- 0.01
gam_scale <- 0.01

#storage and Initial values
mu <- sigma2 <- NULL
mu[1] <- res$parameters$mu;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
sigma2[1] <- 0.005#0.004, 
lambda <- matrix(NA, N, Iter);
lambda[,1] <- res$parameters$mu#+runif(1,0,0.01)#+1/2*res$parameters$sigma2
#lambda[2,,1] <- res$parameters$mu#+runif(1,0,0.01)#exp(res$parameters$mu+1/2*res$parameters$sigma2)
initial_set <- seq(1,N-m+1,m)
log_posterior <- NULL

#record the log-posterior, the first element will keep as NA
sd1 <- 0.1   #for lambda
set_run <- seq(1,N-m+1,m)
  for(t in 2:Iter){
    #record previous logP
    for(i in 1:k){
      lambda_new <- rnorm(1, lambda[set_run[i],t-1], sd1)#-1/2*sd1^2
      l_p <- lambda[set_run[i],t-1]
      tilde_l <- (sum(Y[set_run[i]:(set_run[i]+m-1)])/tau2+mu[t-1]/sigma2[t-1])/(m/tau2+1/sigma2[t-1])
      tilde_v <- 1/(m/tau2+1/sigma2[t-1])
      r <- prod(dnorm(Y[(m*(i-1)+1):(m*(i-1)+m)], lambda_new,sqrt(tau2))) * (dnorm(lambda_new, mu[t-1], sqrt(sigma2[t-1])))^1 /
        prod(dnorm(Y[(m*(i-1)+1):(m*(i-1)+m)], l_p,sqrt(tau2))) / dnorm(l_p, mu[t-1], sqrt(sigma2[t-1])) 
     
      random_number <- runif(1)
      lambda[(m*(i-1)+1):(m*(i-1)+m),t] <- ifelse(random_number < r, lambda_new, lambda[m*(i-1)+1,t-1])
    }
    
    #To sample mu
    mu_proposal <- rnorm(1,mu[t-1],sd1)
    r <- dnorm(mu_proposal,mean(lambda[,t]),sqrt(sigma2[t-1]/k))/
      dnorm(mu[t-1],mean(lambda[,t]),sqrt(sigma2[t-1]/k))
    random_number <- runif(1)
    mu[t] <- ifelse(random_number<r,mu_proposal,mu[t-1])
    
    #sample sigma2
    sigma2_proposal <- rlnorm(1, log(sigma2[t-1])-1/2*sd1^2, sd1)
    r <- prod(dnorm(lambda[,t],mu[t],sqrt(sigma2_proposal)))*
      dinvgamma(sigma2_proposal,gam_shape,gam_scale)/
      prod(dnorm(lambda[,t],mu[t],sqrt(sigma2[t-1])))/
      dinvgamma(sigma2[t-1],gam_shape,gam_scale)*
      dlnorm(sigma2[t-1], log(sigma2_proposal)-1/2*sd1^2, sd1)/
      dlnorm(sigma2_proposal, log(sigma2[t-1])-1/2*sd1^2, sd1)
    sigma2[t] <- ifelse(runif(1) < r, sigma2_proposal, sigma2[t-1])
    ###log-posterior
    log_posterior[t] <- logLike_Y(Y, lambda[,t],tau2) + 
      sum(logLike_Lambda(lambda[set_run,t], mu[t], sigma2[t])) +
      log(dinvgamma(sigma2[t], gam_shape, gam_scale))
  }

effective_samples1 <- cbind(mu[1001:2000],sigma2[1001:2000],t(lambda[,1001:2000]))
effective_samples2 <- cbind(mu[3001:6000],sigma2[3001:6000],t(lambda[,3001:6000]))
effective_samples3 <- cbind(mu[5001:10000],sigma2[5001:10000],t(lambda[,5001:10000]))

posterior <- function(theta)
{
  mmu <- theta[1]
  ssigma2 <- theta[2]
  mui <- theta[3:202]
  logL <- logLike_Y(Y, mui,tau2) + sum(logLike_Lambda(mui, mmu, ssigma2)) +
    log(dinvgamma(ssigma2, gam_shape, gam_scale))
  return(logL)
}

estimates_mode1 <- effective_samples1[which.max(apply(effective_samples1,1,posterior)),]
estimates_mode2 <- effective_samples2[which.max(apply(effective_samples2,1,posterior)),]
estimates_mode3 <- effective_samples3[which.max(apply(effective_samples3,1,posterior)),]

estimates_mean1 <- apply(effective_samples1, 2, mean)
estimates_mean2 <- apply(effective_samples2, 2, mean)
estimates_mean3 <- apply(effective_samples3, 2, mean)

norm.fun <- function(theta,x)
{
  mmu <- theta[1]
  ssigma2 <- theta[2]
  #mui <- theta[3:202]
  logL <- sum(log(dnorm(x,mmu,sqrt(tau2+ssigma2))))+
    #(logLike_Y(x, mui,tau2) + sum(logLike_Lambda(mui, mmu, ssigma2)) +
    log(dinvgamma(ssigma2, gam_shape, gam_scale))
  return(logL)
}

#theta0 <- c(2,0.005)#,res$Lambda
#result <- optim(theta0,norm.fun,x=Y)
#map <- result$par

theta <- theta_old <- c(2,0.1,res$Lambda)
t <- 0
d <- 100
threshold <- 0.01
l <- NULL
while (abs(d)>threshold & t<500) {
  theta[1] <- mean(theta[3:202])
  theta[2] <- (sum((theta[3:202]-theta[1])^2)/2 + gam_scale)/(k/2+gam_shape+1)
  theta[3:202] <- (Y/tau2+theta[1]/theta[2])/(1/tau2+1/theta[2])
  l[t] <- norm.fun(theta,Y)
  d <- norm.fun(theta_old,Y)-norm.fun(theta,Y)
  theta_old <- theta
  t <- t+1
}

#beta_estimates_mode_marginal <- apply(effective_samples, 2, Modes)

#Evaluate results
mse_mui_mode1 <- mean((estimates_mode1[3:202]-res$Lambda)^2)#mean((beta_estimates_mode-beta_true)^2)
mse_mui_mode2 <- mean((estimates_mode2[3:202]-res$Lambda)^2)
mse_mui_mode3 <- mean((estimates_mode3[3:202]-res$Lambda)^2)

mse_par_mode1 <- (estimates_mode1[1:2]-c(res$parameters$mu,res$parameters$sigma2))^2#mean((beta_estimates_mode-beta_true)^2)
mse_par_mode2 <- (estimates_mode2[1:2]-c(res$parameters$mu,res$parameters$sigma2))^2
mse_par_mode3 <- (estimates_mode3[1:2]-c(res$parameters$mu,res$parameters$sigma2))^2

mse_mui_mean1 <- mean((estimates_mean1[3:202]-res$Lambda)^2)#mean((beta_estimates_mean-beta_true)^2)
mse_mui_mean2 <- mean((estimates_mean2[3:202]-res$Lambda)^2)
mse_mui_mean3 <- mean((estimates_mean3[3:202]-res$Lambda)^2)

mse_par_mean1 <- (estimates_mean1[1:2]-c(res$parameters$mu,res$parameters$sigma2))^2#mean((beta_estimates_mode-beta_true)^2)
mse_par_mean2 <- (estimates_mean2[1:2]-c(res$parameters$mu,res$parameters$sigma2))^2
mse_par_mean3 <- (estimates_mean3[1:2]-c(res$parameters$mu,res$parameters$sigma2))^2


mse_estimates_map[ii,1:2] <- (theta[1:2]-c(res$parameters$mu,res$parameters$sigma2))^2
mse_estimates_map[ii,3] <- mean((theta[3:202]-res$Lambda)^2)

mse_mode1[ii,] <- c(mse_par_mode1,mse_mui_mode1)
mse_mode2[ii,] <- c(mse_par_mode2,mse_mui_mode2)
mse_mode3[ii,] <- c(mse_par_mode3,mse_mui_mode3)

mse_mean1[ii,] <- c(mse_par_mean1,mse_mui_mean1)
mse_mean2[ii,] <- c(mse_par_mean2,mse_mui_mean2)
mse_mean3[ii,] <- c(mse_par_mean3,mse_mui_mean3)
}


mo1 <- apply(mse_mode1, 2, mean)
mo2 <- apply(mse_mode2, 2, mean)
mo3 <- apply(mse_mode3, 2, mean)

me1 <- apply(mse_mean1, 2, mean)
me2 <- apply(mse_mean2, 2, mean)
me3 <- apply(mse_mean3, 2, mean)

ma <- apply(mse_estimates_map,2,mean)
#ma <- c(ma,0)
res <- rbind(ma,mo1,mo2,mo3,me1,me2,me3)

setwd("C:/Users/visitor01/OneDrive - Hong Kong Baptist University/Documents/cuhk/MCMC_diagnosis/Diagnosis_summary/Diagnosis_summary/R")

write.csv(res, file ="res0206.csv")
