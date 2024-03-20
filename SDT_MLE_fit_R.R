# This code is a branch from the original package. Ported to R by 
# Trevor Caruso. Please contact trevorcaruso@ufl.edu for any questions or comments.
# For more information about metacognitive and Type 2 Signal Detection Theory (SDT), 
# visit the following website at http://www.columbia.edu/~bsm2105/type2sdt/ and refer 
# to the following publications:

# Maniscalco, B., & Lau, H. (2012). "A signal detection theoretic approach for 
# estimating metacognitive sensitivity from confidence ratings." Consciousness and 
# Cognition, 21(1), 422–430. doi:10.1016/j.concog.2011.09.021

# Maniscalco, B., & Lau, H. (2014). "Signal detection theory analysis of type 1 and 
# type 2 data: meta-d’, response-specific meta-d’, and the unequal variance SDT mode."  
# In S. M. Fleming & C. D. Frith (Eds.), The Cognitive Neuroscience of Metacognition 
# (pp.25-66). Springer.

# If you utilize these functions, please cite the aforementioned papers and scripts 
# upon which they are based.

# The same input-format as M & L is expected. As per the original code:

# INPUTS

# * nR_S1, nR_S2
# these are vectors containing the total number of responses in
# each response category, conditional on presentation of S1 and S2.

# e.g. if nR_S1 = [100 50 20 10 5 1], then when stimulus S1 was
# presented, the subject had the following response counts:
# responded S1, rating=3 : 100 times
# responded S1, rating=2 : 50 times
# responded S1, rating=1 : 20 times
# responded S2, rating=1 : 10 times
# responded S2, rating=2 : 5 times
# responded S2, rating=3 : 1 time

# The ordering of response / rating counts for S2 should be the same as it
# is for S1. e.g. if nR_S2 = [3 7 8 12 27 89], then when stimulus S2 was
# presented, the subject had the following response counts:
# responded S1, rating=3 : 3 times
# responded S1, rating=2 : 7 times
# responded S1, rating=1 : 8 times
# responded S2, rating=1 : 12 times
# responded S2, rating=2 : 27 times
# responded S2, rating=3 : 89 times

# The function's input should consist of counts for each of these responses separately for each stimulus type.


SDT_logL <- function(parameters) {
        
  # Initialize parameters
  S1mu <- -parameters[1] / 2
  S2mu <- parameters[1] / 2
  S1sd <- 1
  S2sd <- 1 / parameters[2]
  c1 <- parameters[3:length(parameters)]
  
  # Calculate response probabilities
  ci <- c(-Inf, c1)
  cj <- c(c1, Inf)
  pC_S1 <- pnorm(cj, mean = S1mu, sd = S1sd) - pnorm(ci, mean = S1mu, sd = S1sd)
  pC_S2 <- pnorm(cj, mean = S2mu, sd = S2sd) - pnorm(ci, mean = S2mu, sd = S2sd)

  # Ensure probabilities are within bounds
  pC_S1 <- ifelse(pC_S1 <= 0, 1e-10, pC_S1)
  pC_S2 <- ifelse(pC_S2 <= 0, 1e-10, pC_S2)

  # Calculate likelihood
  logL <- sum(nR_S1 * log(pC_S1)) + sum(nR_S2 * log(pC_S2))
  if (is.nan(logL) || is.infinite(logL)) {
    logL <- -Inf
  }
  return(-logL)  # Return negative log likelihood

}


SDT_MLE_fit <- function(nR_S1, nR_S2) {
  
  # Preprocess data and set up initial conditions
  nRatings <- length(nR_S1) / 2
  nCriteria <- 2 * nRatings - 1
  
  ratingHR <- numeric()
  ratingFAR <- numeric()
  for (c in 2:(nRatings * 2)) {
    ratingHR <- c(ratingHR, sum(nR_S2[c:length(nR_S2)]) / sum(nR_S2))
    ratingFAR <- c(ratingFAR, sum(nR_S1[c:length(nR_S1)]) / sum(nR_S1))
  }
  
  # Use linear model fitting as an equivalent to MATLAB's polyfit
  if (length(ratingFAR) == 0 || all(ratingFAR == ratingFAR[1])) {
    s <- 1
  } else {
    fit <- lm(qnorm(ratingHR) ~ qnorm(ratingFAR) + 0)
    s <- coef(fit)[1]
  }
  
  
  
  d1 <- (1/s) * qnorm(ratingHR) - qnorm(ratingFAR)
  d1 <- mean(d1)
  c1 <- (-1/(1+s)) * (qnorm(ratingHR) + qnorm(ratingFAR))
  guess <- c(d1, s, c1)
  

  # Optimizing using 'optim'
  optim_result <- optim(par = guess, fn = SDT_logL, method = "L-BFGS-B",
                        lower = c(-Inf, 1e-6, rep(-Inf, nCriteria)),
                        upper = c(Inf, Inf, rep(Inf, nCriteria)))
  
  # Analysis of the result
  analysis <- list(
    logL = -optim_result$value,
    k = length(guess),
    n = sum(nR_S1) + sum(nR_S2),
    AIC = -2 * -optim_result$value + 2 * length(guess),
    AICc = (-2 * -optim_result$value + 2 * length(guess) * sum(nR_S1 + nR_S2)) / (sum(nR_S1 + nR_S2) - length(guess) - 1),
    BIC = -2 * -optim_result$value + length(guess) * log(sum(nR_S1 + nR_S2))
  )

  params <- list(
    d1 = optim_result$par[1],
    s = optim_result$par[2],
    c1 = optim_result$par[3:length(optim_result$par)],
    guess = guess
  )
  
  # Printing s value
  cat("s value:", params$s, "\n")
  
  # Returning results
  return(list(analysis = analysis, params = params))
}


##################
# Example usage with the provided nR_S1 and nR_S2 (sanity check; result should be ~1)
nR_S1 = c(100, 50, 20, 10, 5, 1)
nR_S2 = c(1, 5, 10, 20, 50, 100)

# Run the SDT_MLE_fit function
result <- SDT_MLE_fit(nR_S1, nR_S2)
result

# the estimated s-value is the ratio of standard deviations for type 1 distributions. Use-case: unequal variance
# s is the slope of z(HR) vs z(FAR)
# estimated_s = sd(S1) / sd(S2)

# If using the code to extract the estimated s-value:
# Extracting the s-value from the result and storing it in the variable s_estim
s_value <- unname(result$params$s)
s_value

