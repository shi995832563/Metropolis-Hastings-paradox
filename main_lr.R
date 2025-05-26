library(parallel)
setwd("C:/Users/visitor01/OneDrive - The Chinese University of Hong Kong/New folder/Documents/cuhk/MCMC_diagnosis/Diagnosis_summary/Diagnosis_summary/R")


source(file="simulation_function_change2.R")

cores <- detectCores(logical = FALSE)
cl <- makeCluster(2)


p <- c(50,100)

system.time(
  res1 <- parLapply(cl, p , simulation_function_change2 )
)

stopCluster(cl)

