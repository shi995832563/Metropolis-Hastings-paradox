simulation_function_change2<-function(p){
  source("linear_regression.R")
  n <- c(110,200,500,1000)
  results <- matrix(0,7,length(n))
  for (k in 1:4) 
  {
    ni <- n[k]
    results[,k] <- linear_regression(ni,p)
  }
  
  file_name<-paste("o_",p,".csv",sep = "")
  write.csv(results, file =file_name)
}