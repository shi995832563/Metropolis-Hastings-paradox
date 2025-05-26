# Set seed for reproducibility
#set.seed(123)
linear_regression <- function(n,p)
{
Modes <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}


# Simulate data
#n <- 100  # Number of observations
#p <- 50   # Number of predictors
rep1 <- 1000
mse_estimates_mode1 <- mse_estimates_mean1 <- mse_estimates_mode2 <- mse_estimates_mean2 <- mse_estimates_mode3 <- mse_estimates_mean3 <- NULL
mse_estimates_map <- mse_estimates_mode_marginal <- NULL
#mse_estimates_mode_matrix <- mse_estimates_mean_matrix <- mse_estimates_map_matrix <- mse_estimates_mode_marginal_matrix <- matrix(0,rep1,p)
beta_true <- rep(1,p)
for (j in 1:rep1) 
{
# True coefficients
  # True coefficients
X <- matrix(rnorm(n * p), nrow = n, ncol = p)  # Design matrix
epsilon <- rnorm(n, mean = 0, sd = 1)  # Gaussian error
y <- X %*% beta_true + epsilon  # Response variable

# Prior parameters
prior_mean <- rep(0, p)
prior_sd <- 100  # Large prior variance for coefficients

# Likelihood function (Gaussian)
likelihood <- function(beta) {
  mu <- X %*% beta
  sum(dnorm(y, mean = mu, sd = 1, log = TRUE))  # Assume known error variance
}

# Prior distribution for coefficients (Gaussian)
prior <- function(beta) {
  sum(dnorm(beta, mean = prior_mean, sd = prior_sd, log = TRUE))
}

# Posterior function
posterior <- function(beta) {
  likelihood(beta) #+ prior(beta)
}

fit <- lm(y~X+0)
fit
# Metropolis-Hastings algorithm
metropolis_hastings <- function(initial, iterations) {
  beta_samples <- matrix(NA, nrow = iterations, ncol = p)
  beta_current <- initial
  acceptance_rate <- 0
  log_posterior <- NULL
  
  for (i in 1:iterations) {
    # Propose new beta
    beta_proposed <- rnorm(p, mean = beta_current, sd = 0.01)  # Propose step size
    
    # Calculate acceptance probability
    log_acceptance_ratio <- posterior(beta_proposed) - posterior(beta_current)
    acceptance_prob <- exp(log_acceptance_ratio)
    
    # Accept or reject
    if (runif(1) < acceptance_prob) {
      beta_current <- beta_proposed
      acceptance_rate <- acceptance_rate + 1
    }
    
    beta_samples[i, ] <- beta_current  # Store the sample
    log_posterior[i] <- posterior(beta_current)
  }
  
  list(samples = beta_samples, log_posterior = log_posterior, acceptance_rate = acceptance_rate / iterations)
}

# Run MCMC
initial_beta <- fit$coefficients#rep(0, p)  # Initial guess for coefficients
iterations <- 10000  # Number of iterations
results <- metropolis_hastings(initial_beta, iterations)

# Summarize results
effective_samples1 <- results$samples[1001:2000,]
effective_samples2 <- results$samples[3001:6000,]
effective_samples3 <- results$samples[5001:10000,]

beta_estimates_mode1 <- effective_samples1[which.max(apply(effective_samples1,1,posterior)),]
beta_estimates_mode2 <- effective_samples2[which.max(apply(effective_samples2,1,posterior)),]
beta_estimates_mode3 <- effective_samples3[which.max(apply(effective_samples3,1,posterior)),]

beta_estimates_mean1 <- apply(effective_samples1, 2, mean)
beta_estimates_mean2 <- apply(effective_samples2, 2, mean)
beta_estimates_mean3 <- apply(effective_samples3, 2, mean)

#beta_estimates_mode_marginal <- apply(effective_samples, 2, Modes)

#Evaluate results
mse_estimates_mode1[j] <- sum((beta_estimates_mode1-beta_true)^2)#mean((beta_estimates_mode-beta_true)^2)
mse_estimates_mode2[j] <- sum((beta_estimates_mode2-beta_true)^2)
mse_estimates_mode3[j] <- sum((beta_estimates_mode3-beta_true)^2)

mse_estimates_mean1[j] <- sum((beta_estimates_mean1-beta_true)^2)#mean((beta_estimates_mean-beta_true)^2)
mse_estimates_mean2[j] <- sum((beta_estimates_mean2-beta_true)^2)
mse_estimates_mean3[j] <- sum((beta_estimates_mean3-beta_true)^2)

mse_estimates_map[j] <- sum((fit$coefficients-beta_true)^2)#mean((fit$coefficients-beta_true)^2)
#mse_estimates_mode_marginal_matrix[j] <- #mean((beta_estimates_mode_marginal-beta_true)^2)
}


#mse_estimates_mode
#mse_estimates_mean
#mse_estimates_map
return(y <- c(mean(mse_estimates_map),mean(mse_estimates_mode1),mean(mse_estimates_mode2),mean(mse_estimates_mode3),mean(mse_estimates_mean1),mean(mse_estimates_mean2),mean(mse_estimates_mean3)))

}
# Plotting the trace of the samples
#par(mfrow = c(1, 1))  # Set up the plotting area for 10 coefficients
#plot(results$log_posterior[1:10000], type = 'l', main = "log-posterior", xlab = "Iteration", ylab = "Value")
#abline(h=posterior(beta_true*0))
