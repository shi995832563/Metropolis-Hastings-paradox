hmm <- function(ssd)
{
# Viterbi algorithm for HMM with Gaussian emissions
viterbi <- function(observations, A, means, initial_probs) {
  # observations: vector of observed data
  # A: 4x4 transition matrix
  # means: emission means for states 1-4
  # initial_probs: initial state probabilities (length 4)
  
  n_states <- 4  #number of states
  T <- length(observations)  #length of the sequence
  delta <- matrix(0, nrow = T, ncol = n_states) #claim the matrix for the best foward probability
  psi <- matrix(0, nrow = T, ncol = n_states) #claim the matrix for the parameter 
  best_path <- rep(0, T) #claim the best path
  
  # Initialize (t = 1)
  delta[1, ] <- initial_probs * dnorm(observations[1], mean = means, sd = ssd) #initialize the forward probability in the first step
  
  # Recursion (t = 2 to T)
  for (t in 2:T) {
    for (s in 1:n_states) {
      prev_probs <- delta[t-1, ] * A[, s] #Each update for t and s, which is to exhaust all the probabilities only with transition
      psi[t, s] <- which.max(prev_probs) #Selet the state with the largest forward probability
      delta[t, s] <- prev_probs[psi[t, s]] * dnorm(observations[t], mean = means[s], sd = ssd) #The best forward probability matrix with emission
    }
  }
  
  # Termination (t = T)
  best_path[T] <- which.max(delta[T, ]) #Start from the very end to select the state with the largest probability
  
  # Backtracking (t = T-1 to 1)
  for (t in (T-1):1) {
    best_path[t] <- psi[t+1, best_path[t+1]] #Select the state with the largest probability for each element of the sequence
  }
  
  return(best_path)
}

posterior_density <- function(x)
{
  l <- dim(x)
  p <- 0.25
  for (i in 1:(l[1])) 
  {
    for (j in 1:(l[2]-1)) {
      p[i] <- p[i]*dnorm(observations[j],x[i,j],ssd)*A[x[i,j],x[i,j+1]]
    }
    p[i] <- p[i]*dnorm(observations[n],x[i,n],ssd)
  }
  return(x[which.max(p),])
}
# Metropolis-Hastings sampler for latent states
mh_hmm <- function(observations, A, means, initial_probs, n_iter = 5000, burn_in = 0) {
  # observations: vector of observed data
  # A: 4x4 transition matrix
  # means: emission means for states 1-4
  # initial_probs: initial state probabilities
  # n_iter: total MCMC iterations
  # burn_in: number of burn-in samples
  
  T <- length(observations)
  states <- matrix(0, nrow = n_iter, ncol = T)
  current_states <- true_states#sample(1:4, T, replace = TRUE)  # Initialize randomly
  
  for (iter in 1:n_iter) {
    current_s <- current_states[1]
    proposed_s <- sample(1:4,1)
    # Calculate emission likelihood ratio
    likelihood_ratio <- dnorm(observations[1], mean = means[proposed_s], sd = ssd) /
      dnorm(observations[1], mean = means[current_s], sd = ssd)
    
    # Calculate transition probabilities
    trans_ratio_prev <- 1
    trans_ratio_next <- 1
    
    # Next state (if t < T)
      trans_ratio_next <- A[proposed_s, current_states[2]] / 
        A[current_s, current_states[2]]
    
    # Acceptance probability
    alpha <- min(1, likelihood_ratio * trans_ratio_next)
    
    # Accept/reject
    if (runif(1) < alpha) {
      current_states[1] <- proposed_s
    }
    for (t in 2:T) {
      current_s <- current_states[t]
      sample(1:4,1)
      #proposed_s <- sample(1:4, 1, prob = A[current_states[t-1],])  # Propose a new state
      
      # Calculate emission likelihood ratio
      likelihood_ratio <- dnorm(observations[t], mean = means[proposed_s], sd = ssd) /
        dnorm(observations[t], mean = means[current_s], sd = ssd)
      
      # Calculate transition probabilities
      trans_ratio_prev <- 1
      trans_ratio_next <- 1
      
      # Previous state (if t > 1)
      if(t<T){
        alpha <- 1
        if((A[current_s, current_states[t+1]]*A[current_states[t-1], current_s])>0){
          trans_ratio_prev <- A[current_states[t-1], proposed_s] / 
            A[current_states[t-1], current_s]
          trans_ratio_next <- A[proposed_s, current_states[t+1]] / 
            A[current_s, current_states[t+1]]
          alpha <- min(1, likelihood_ratio * trans_ratio_prev * trans_ratio_next)
        }
      }else{
        trans_ratio_prev <- A[current_states[t-1], proposed_s] / 
          A[current_states[t-1], current_s]
        alpha <- min(1, likelihood_ratio * trans_ratio_prev * trans_ratio_next)
      }
      # Acceptance probability
      
      # Accept/reject
      if (runif(1) < alpha) {
        current_states[t] <- proposed_s      }
    }
    states[iter, ] <- current_states
  }
  
  # Discard burn-in and return posterior samples
  return(states[(burn_in + 1):n_iter, ])
}


# Define model parameters
A <- matrix(c(
  0.8, 0.4/3, 0, 0.2/3,
  0.2/3, 0.8, 0.2/3, 0.2/3,
  0.2/3, 0.2/3, 0.8, 0.2/3,
  0.4/3, 0, 0.2/3, 0.8
), nrow = 4, byrow = TRUE)

means <- c(1, 2, 3, 4)  # Emission means for states 1-4
initial_probs <- rep(0.25, 4)

repeti <- 1000
# Simulate synthetic data
#set.seed(123)
p_map <- p_mode1 <- p_mode2 <- p_mode3 <- p_mean1 <- p_mean2 <- p_mean3 <- rep(0,repeti)
accuracy_map <- accuracy_mean1 <- accuracy_mean2 <- accuracy_mean3 <- accuracy_mode1 <- accuracy_mode2 <- accuracy_mode3 <- NULL
n <- 100
for (rrr in 1:repeti) 
{
true_states <- NULL
true_states[1] <- sample(1:4,1)
for (i in 2:100) {
  true_states[i] <- sample(1:4,1,prob=A[true_states[i-1],])
}

observations <- rnorm(n, mean = true_states, sd = ssd)
###########################################################


# Run Viterbi
decoded_states <- viterbi(observations, A, means, initial_probs)

# Calculate accuracy (if true states are known)
accuracy_map[rrr] <- mean(decoded_states == true_states)
mse_map <- mean((decoded_states-true_states)^2)
#cat("Viterbi Accuracy:", round(accuracy * 100, 2), "%\n")

##################################################
# Run MH sampler
posterior_samples <- mh_hmm(observations, A, means, initial_probs, n_iter = 10000)

# Posterior mode estimate (most frequent state at each time point)
posterior_mode1 <- posterior_density(posterior_samples[1001:2000,])
posterior_mode2 <- posterior_density(posterior_samples[3001:6000,])
posterior_mode3 <- posterior_density(posterior_samples[5001:10000,])

posterior_mean1 <- round(apply(posterior_samples[1001:2000,], 2, mean))
posterior_mean2 <- round(apply(posterior_samples[3001:6000,], 2, mean))
posterior_mean3 <- round(apply(posterior_samples[5001:10000,], 2, mean))


# Compare with true states
accuracy_mode1[rrr] <- mean(posterior_mode1 == true_states)
accuracy_mode2[rrr] <- mean(posterior_mode2 == true_states)
accuracy_mode3[rrr] <- mean(posterior_mode3 == true_states)

accuracy_mean1[rrr] <- mean(posterior_mean1 == true_states)
accuracy_mean2[rrr] <- mean(posterior_mean2 == true_states)
accuracy_mean3[rrr] <- mean(posterior_mean3 == true_states)

mse_mode1 <- mean((posterior_mode1-true_states)^2)
mse_mode2 <- mean((posterior_mode2-true_states)^2)
mse_mode3 <- mean((posterior_mode3-true_states)^2)

mse_mean1 <- mean((posterior_mean1-true_states)^2)
mse_mean2 <- mean((posterior_mean2-true_states)^2)
mse_mean3 <- mean((posterior_mean3-true_states)^2)

p_map[rrr] <- ((sum(decoded_states[2:n]-decoded_states[1:(n-1)]==2 & decoded_states[2:n]*decoded_states[1:(n-1)]==3)
                              +sum(decoded_states[2:n]-decoded_states[1:(n-1)]==-2 & decoded_states[2:n]*decoded_states[1:(n-1)]==8))>0)

p_mode1[rrr] <- ((sum(posterior_mode1[2:n]-posterior_mode1[1:(n-1)]==2 & posterior_mode1[2:n]*posterior_mode1[1:(n-1)]==3)
                               +sum(posterior_mode1[2:n]-posterior_mode1[1:(n-1)]==-2 & posterior_mode1[2:n]*posterior_mode1[1:(n-1)]==8))>0)
p_mode2[rrr] <- ((sum(posterior_mode2[2:n]-posterior_mode2[1:(n-1)]==2 & posterior_mode2[2:n]*posterior_mode2[1:(n-1)]==3)
                               +sum(posterior_mode2[2:n]-posterior_mode2[1:(n-1)]==-2 & posterior_mode2[2:n]*posterior_mode2[1:(n-1)]==8))>0)
p_mode3[rrr] <- ((sum(posterior_mode3[2:n]-posterior_mode3[1:(n-1)]==2 & posterior_mode3[2:n]*posterior_mode3[1:(n-1)]==3)
                               +sum(posterior_mode3[2:n]-posterior_mode3[1:(n-1)]==-2 & posterior_mode3[2:n]*posterior_mode3[1:(n-1)]==8))>0)

p_mean1[rrr] <- ((sum(posterior_mean1[2:n]-posterior_mean1[1:(n-1)]==2 & posterior_mean1[2:n]*posterior_mean1[1:(n-1)]==3)
                               +sum(posterior_mean1[2:n]-posterior_mean1[1:(n-1)]==-2 & posterior_mean1[2:n]*posterior_mean1[1:(n-1)]==8))>0)
p_mean2[rrr] <- ((sum(posterior_mean2[2:n]-posterior_mean2[1:(n-1)]==2 & posterior_mean2[2:n]*posterior_mean2[1:(n-1)]==3)
                               +sum(posterior_mean2[2:n]-posterior_mean2[1:(n-1)]==-2 & posterior_mean2[2:n]*posterior_mean2[1:(n-1)]==8))>0)
p_mean3[rrr] <- ((sum(posterior_mean3[2:n]-posterior_mean3[1:(n-1)]==2 & posterior_mean3[2:n]*posterior_mean3[1:(n-1)]==3)
                               +sum(posterior_mean3[2:n]-posterior_mean3[1:(n-1)]==-2 & posterior_mean3[2:n]*posterior_mean3[1:(n-1)]==8))>0)
}

accu <- c(mean(accuracy_map),mean(accuracy_mode1),mean(accuracy_mode2),mean(accuracy_mode3),mean(accuracy_mean1),mean(accuracy_mean2),mean(accuracy_mean3))
p_mis <- c(mean(p_map),mean(p_mode1),mean(p_mode2),mean(p_mode3),mean(p_mean1),mean(p_mean2),mean(p_mean3))
res <- cbind(accu,p_mis)

file_name<-paste("ssd_",ssd,".csv",sep = "")
write.csv(res, file =file_name)
}