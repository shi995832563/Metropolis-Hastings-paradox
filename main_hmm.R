library(parallel)
setwd("C:/Users/visitor01/OneDrive - The Chinese University of Hong Kong/New folder/Documents/cuhk/MCMC_diagnosis/Diagnosis_summary/Diagnosis_summary/R")


source(file="hmm.R")

cores <- detectCores(logical = FALSE)
cl <- makeCluster(3)


ssd <- c(0.1,0.5,1)

system.time(
  res1 <- parLapply(cl, ssd , hmm )
)

stopCluster(cl)

