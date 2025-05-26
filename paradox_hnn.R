set.seed(4)
library(mvtnorm)
library(MASS)
library(invgamma)
library(LaplacesDemon)

#simulated data
tau2 <- 0.1
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


N = 200
k <- N#200#2
m <- N/k
res <- sim(N,k)

Y <- res$Y
Iter = 2000
nRun= 4

#-----------------------------------1>
start_time = proc.time()

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
mu <- matrix(NA, nRun, Iter);
mu[1:nRun,1] <- res$parameters$mu;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
sigma2 <- matrix(NA, nRun, Iter); 
sigma2[1:nRun,1] <- c(0.005, 0.5, 0.1, 0.2)#0.004, 
lambda <- array(NA, c(nRun, N, Iter));
lambda[1,,1] <- res$parameters$mu#+runif(1,0,0.01)#+1/2*res$parameters$sigma2
lambda[2,,1] <- res$parameters$mu#+runif(1,0,0.01)#exp(res$parameters$mu+1/2*res$parameters$sigma2)
initial_set <- seq(1,N-m+1,m)
for(iRun in 3:nRun){
  for (i in 1:k) 
  {
    lambda[iRun,(m*(i-1)+1):(m*(i-1)+m),1] <- ifelse(Y[(m*(i-1)+1)]==0, 0.2*iRun, Y[(m*(i-1)+1)])+runif(1,0,0.0001)
  }
}


#record the log-posterior, the first element will keep as NA
logP <- logL <- matrix(NA, nRun, Iter);  
log_posterior <- array(0, c(nRun, k+1, Iter));
p_increasing <- p_decreasing <- array(0, c(nRun, k, Iter));
rr <- matrix(NA,nRun,Iter);
v_logP <- matrix(NA,nRun,Iter);
logP_lambda1 <- logP_lambda2 <- NULL
r_matrix <- rr_matrix <- array(NA, c(nRun, N, Iter));
sd1 <- sqrt(0.03)   #for lambda
set_run <- seq(1,N-m+1,m)
for(iRun in 1:nRun){
  cat("\nstart run ", iRun, ":")
  #log_p_matrix2 <- NULL
  logL[iRun,1] <- logLike_L(Y, mu[iRun,1],sigma2[iRun,1],tau2)
  logP[iRun,1]= logLike_Y(Y, lambda[iRun,,1],tau2) + 
    sum(logLike_Lambda(lambda[iRun,set_run,1], mu[iRun, 1], sigma2[iRun,1])) +
    log(dinvgamma(sigma2[iRun,1], gam_shape, gam_scale))
  log_posterior[iRun,,1] <- logP[iRun,1]
  for(t in 2:Iter){
    #record previous logP
    logL[iRun,t-1] <- logLike_L(Y, mu[iRun,t-1],sigma2[iRun,t-1],tau2)
    logP[iRun,t-1]= logLike_Y(Y, lambda[iRun,,t-1],tau2) + 
      sum(logLike_Lambda(lambda[iRun,set_run,t-1], mu[iRun, t-1], sigma2[iRun,t-1])) +
      log(dinvgamma(sigma2[iRun,t-1], gam_shape, gam_scale))
    lambda[iRun,,t] <- lambda[iRun,,t-1]
    for(i in 1:k){
      lambda_new <- rnorm(1, lambda[iRun,set_run[i],t-1], sd1)#-1/2*sd1^2
      l_p <- lambda[iRun,set_run[i],t-1]
      tilde_l <- (sum(Y[set_run[i]:(set_run[i]+m-1)])/tau2+mu[iRun,t-1]/sigma2[iRun,t-1])/(m/tau2+1/sigma2[iRun,t-1])
      tilde_v <- 1/(m/tau2+1/sigma2[iRun,t-1])
      r <- prod(dnorm(Y[(m*(i-1)+1):(m*(i-1)+m)], lambda_new,sqrt(tau2))) * (dnorm(lambda_new, mu[iRun,t-1], sqrt(sigma2[iRun,t-1])))^1 /
        prod(dnorm(Y[(m*(i-1)+1):(m*(i-1)+m)], l_p,sqrt(tau2))) / dnorm(l_p, mu[iRun,t-1], sqrt(sigma2[iRun,t-1])) 
      p_increasing[iRun,i,t] <- abs(pnorm(2*abs(l_p-tilde_l)/sd1))-0.5
      
      mm1 <- (tilde_l/tilde_v+l_p/sd1^2)/(1/tilde_v+1/sd1^2)
      ss1 <- 1/sqrt(1/tilde_v+1/sd1^2)
      p1 <- pnorm((l_p-mm1)/ss1)
      p2<- pnorm((2*tilde_l-l_p-mm1)/ss1)
      p_decreasing[iRun,i,t] <- dnorm((tilde_l-l_p)/sqrt(tilde_v+sd1^2))*
        ifelse(l_p<tilde_l,p1+1-p2,1-p1+p2)/
        dnorm((l_p-tilde_l)/sqrt(tilde_v))
      
      #p_decreasing[iRun,i,t] <- (1-p_increasing[iRun,i,t])*0.95*exp(-1.96*sd1*(2*abs(tilde_l-l_p)+1.96*sd1)/2/tilde_v)
      
      random_number <- runif(1)
      lambda[iRun,(m*(i-1)+1):(m*(i-1)+m),t] <- ifelse(random_number < r, lambda_new, lambda[iRun,m*(i-1)+1,t-1])
      
      log_posterior[iRun,i,t] <- logLike_Y(Y, lambda[iRun,,t],tau2) + 
        sum(logLike_Lambda(lambda[iRun,,t], mu[iRun, t-1], sigma2[iRun,t-1])) +
        log(dinvgamma(sigma2[iRun,t-1], gam_shape, gam_scale))
      #log_p_matrix2 <- c(log_p_matrix2,log_posterior[iRun,i,t])
    }
    
    #To sample mu
    mu[iRun,t] <- mu[iRun,t-1]
    
    #sample sigma2
    sigma2_proposal <- rlnorm(1, log(sigma2[iRun,t-1])-1/2*sd1^2, sd1)
    r <- prod(dnorm(lambda[iRun,,t],mu[iRun,t],sqrt(sigma2_proposal)))*
      dinvgamma(sigma2_proposal,gam_shape,gam_scale)/
      prod(dnorm(lambda[iRun,,t],mu[iRun,t],sqrt(sigma2[iRun,t-1])))/
      dinvgamma(sigma2[iRun,t-1],gam_shape,gam_scale)*
      dlnorm(sigma2[iRun,t-1], log(sigma2_proposal)-1/2*sd1^2, sd1)/
      dlnorm(sigma2_proposal, log(sigma2[iRun,t-1])-1/2*sd1^2, sd1)
    sigma2[iRun,t] <- ifelse(runif(1) < r, sigma2_proposal, sigma2[iRun,t-1])
    
    log_posterior[iRun,k+1,t] <- logLike_Y(Y, lambda[iRun,,t],tau2) + 
      sum(logLike_Lambda(lambda[iRun,set_run,t], mu[iRun, t], sigma2[iRun,t])) +
      log(dinvgamma(sigma2[iRun,t], gam_shape, gam_scale))
    #progress hint
    if(0 == t%%200)
      cat(t,",")  
  }
  #calculate logP for the last iteration
  logP[iRun,t]= logLike_Y(Y, lambda[iRun,,t],tau2) + 
    sum(logLike_Lambda(lambda[iRun,set_run,t], mu[iRun,t], sigma2[iRun,t])) +
    log(dinvgamma(sigma2[iRun,t], gam_shape, gam_scale))
  logL[iRun,t] <- logLike_L(Y, mu[iRun,t], sigma2[iRun,t], tau2)
  
}

par(mfrow=c(2,1),mar=c(3,3,1,0.5),oma=c(0,0,0,0.5),mgp=c(1.5,0.5,0)) 
Iter_set <- 1:1000
#win.graph()
plot(Iter_set, sigma2[1,Iter_set], ylim=range(sigma2[c(1,2),]), xlab = "", ylab = expression(math(tau)^2), type="l", main="", col=1)
lines(Iter_set, sigma2[2,Iter_set],col=2)
abline(h = res$parameters$sigma2, lty=2, col=2)

legend(700 ,0.45 ,c("Setting (a)","Setting (b)"),
       lty=c(1,1),col=c(2,1),
       lwd=1,cex=1, box.lty=0 , bg="transparent")

true_ll <- logLike_Y(Y, res$Lambda,tau2) + 
  sum(logLike_Lambda(res$Lambda[set_run], res$parameters$mu, res$parameters$sigma2))+
  log(dinvgamma(res$parameters$sigma2, gam_shape, gam_scale))

plot(Iter_set, logP[1,Iter_set], ylim=range(logP), type="l", xlab = "Iteration",ylab ="log-posterior distribution", main="", col=1)
#for(i in 2:nRun)
lines(Iter_set, logP[2,Iter_set],col=2)
abline(h=true_ll, lty=2, col=2)


p_inc <- p_dec <- NULL
for (i in 1:Iter) 
{
  p_inc[i] <- mean(p_increasing[1,,i])
  p_dec[i] <- mean(p_decreasing[1,,i])
}

par(mfrow=c(1,1),mar=c(3,3,1,0.5),oma=c(0,0,0,0.5),mgp=c(1.5,0.5,0)) 
Iter_set <- 2:1000

plot(Iter_set, p_inc[Iter_set], ylim=c(0,1), type="l", xlab = "Iteration",ylab ="Probability", main="", col=2)
lines(Iter_set, p_dec[Iter_set], col=3)
lines(Iter_set, 1-p_inc[Iter_set]-p_dec[Iter_set], col=1)
abline(h=0.5)

legend(800 ,0.9 ,c(expression(math(p)[h]^"(t)"),expression(math(p)[h]^"(t)"),expression(math(p)[h]^"(t)")),
       lty=c(1,1,1),col=c(2,3,1),
       lwd=1,cex=1, box.lty=0 , bg="transparent")




