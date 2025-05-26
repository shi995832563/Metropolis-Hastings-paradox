set.seed(3)
library(mvtnorm)
library(MASS)
library(invgamma)
library(LaplacesDemon)

p <- 200
n <- 1
y0 <- rmvnorm(1,rep(0,p),diag(p))
mu <- 0
sigma <- 1

logL <- function(y,mu,sigma)
{
  l <- sum(log(dnorm(y,mu,sigma)))
}

Iter = 2000
nRun= 2
y <- array(NA, c(nRun, p, Iter));
y[1,,1] <- 0
y[2,,1] <- 2

#record the log-posterior, the first element will keep as NA
logP <- matrix(NA, nRun, Iter);  
log_posterior <- array(0, c(nRun, p, Iter));
p_increasing <- p_decreasing <- array(0, c(nRun, p, Iter));
rr <- matrix(NA,nRun,Iter);
v_logP <- matrix(NA,nRun,Iter);
logP_lambda1 <- logP_lambda2 <- NULL
r_matrix <- rr_matrix <- array(NA, c(nRun, n, Iter));
sd1 <- sqrt(0.03)  #for lambda
set_run <- 1:n
for(iRun in 1:nRun){
  cat("\nstart run ", iRun, ":")
  #log_p_matrix2 <- NULL
  logP[iRun,1]= logL(y[iRun,,1], mu,sigma)   ###########这里不是要估计参数 而就是要抽到这个分布的数据而已
  log_posterior[iRun,,1] <- logP[iRun,1]
  for(t in 2:Iter){
    #record previous logP
    logP[iRun,t-1]= logL(y[iRun,,t-1], mu,sigma)
    y[iRun,,t] <- y[iRun,,t-1]
    for(i in 1:p){
      y_new <- rnorm(1, y[iRun,i,t-1], sd1)#-1/2*sd1^2
      l_p <- y[iRun,i,t-1]
      tilde_l <- mu
      tilde_v <- sigma^2
      r <- dnorm(y_new,mu,sigma)/dnorm(y[iRun,i,t-1],mu,sigma)
      
      p_increasing[iRun,i,t] <- abs(pnorm(2*abs(l_p-tilde_l)/sd1))-0.5
      
      mm1 <- (tilde_l/tilde_v+l_p/sd1^2)/(1/tilde_v+1/sd1^2)
      ss1 <- 1/sqrt(1/tilde_v+1/sd1^2)
      p1 <- pnorm((l_p-mm1)/ss1)
      p2<- pnorm((2*tilde_l-l_p-mm1)/ss1)
      p_decreasing[iRun,i,t] <- dnorm((tilde_l-l_p)/sqrt(tilde_v+sd1^2))*
        ifelse(l_p<tilde_l,p1+1-p2,1-p1+p2)/
        dnorm((l_p-tilde_l)/sqrt(tilde_v))
      
      p_decreasing[iRun,i,t] <- (1-p_increasing[iRun,i,t])*0.95*exp(-1.96*sd1*(2*abs(tilde_l-l_p)+1.96*sd1)/2/tilde_v)
      
      random_number <- runif(1)
      y[iRun,i,t] <- ifelse(random_number < r, y_new, y[iRun,i,t-1])
      
      log_posterior[iRun,i,t] <- logL(y[iRun,,t], mu,sigma)
      #log_p_matrix2 <- c(log_p_matrix2,log_posterior[iRun,i,t])
    }
    #progress hint
    if(0 == t%%200)
      cat(t,",")  
  }
  #calculate logP for the last iteration
  logP[iRun,t]= logL(y[iRun,,t], mu,sigma)
}


logl <- v <- matrix(NA,2,Iter)
for (i in 1:2) {
  for (j in 1:Iter) 
  {
    yy <- rmvnorm(1,rep(0,p),diag(p))
    logl[i,j] <- logL(yy,mu,sigma)
    v[i,j] <- var(c(yy))
  }}

par(mfrow=c(1,1),mar=c(3,3,1,0.5),oma=c(0,0,0,0.5),mgp=c(1.5,0.5,0)) 
Iter_set <- 1:1000

true_ll <- logL(0, rep(0,p),1) 

#win.graph()
plot(Iter_set, logl[2,Iter_set], ylim=range(c(logP,logl[2,Iter_set])), type="l", xlab = "Iteration",ylab ="log pdf", main="", col=2)
lines(Iter_set, logP[2,Iter_set], col=3)
lines(Iter_set, logP[1,Iter_set], col=1)

legend(700 ,-400 ,c("Direct sampling","MH Setting I","MH Setting II"),
       lty=c(1,1,1),col=c(2,3,1),
       lwd=1,cex=1, box.lty=0 , bg="transparent")


####################################################
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

legend(800 ,0.9 ,c(expression(math(p)[h]^"(t)"),expression(math(p)[l]^"(t)"),expression(math(p)[s]^"(t)")),
       lty=c(1,1,1),col=c(2,3,1),
       lwd=1,cex=1, box.lty=0 , bg="transparent")

