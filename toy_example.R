mh_stn <- function(x)
{
  #proposal procedure
  p <- rnorm(1,x,0.1)
  #acceptance ratio
  r <- dnorm(p)/dnorm(x)
  #decision
  y <- ifelse(r>runif(1),p,x)
  return(y)
}
set.seed(78)
sequence1 <- sequence2 <- NULL
sequence1[1] <- 0
sequence2[1] <- 10
for (i in 2:2000) 
{
  sequence1[i] <- mh_stn(sequence1[i-1])
  sequence2[i] <- mh_stn(sequence2[i-1])
}

plot(1:2000,sequence1,ylim=c(-3,12),type="l", xlab = "Iteration",ylab ="log target distribution", main="", col=2)
lines(1:2000,sequence2,col="red")

par(mfrow=c(1,1),mar=c(3,3,1,0.5),oma=c(0,0,0,0.5),mgp=c(1.5,0.5,0)) 
plot(1:2000,log(dnorm(sequence2)),type="l", xlab = "Iteration",ylab ="log target distribution", main="", col=1)
lines(1:2000,log(dnorm(sequence1)),col=2)

legend(1300 ,-40 ,c("1) the initial value as 0","2) the initial value as 10"),
       lty=c(1,1),col=c(2,1),
       lwd=1,cex=1, box.lty=0 , bg="transparent")
