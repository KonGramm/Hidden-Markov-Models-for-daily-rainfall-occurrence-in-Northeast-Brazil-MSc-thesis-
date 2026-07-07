

compute_state_probs = function(m, n, X, theta)
{
  transition_probs = list(1:(n-1))
  for (k in 1:(n-1))
  {
    s = theta$gamma + matrix(rep(X[k+1,] %*% t(theta$rho), m), m, m, byrow = TRUE)
    s_max = max(s)
    transition_probs[[k]] = exp(s - s_max)
    #transition_probs = transition_probs / matrixrep(rowSums(transition_probs), m), m, m)
    transition_probs[[k]] = transition_probs[[k]] / apply(transition_probs[[k]], 1, sum)
  }
  
  s = theta$delta + X[1,] %*% t(theta$rho)
  s_max = max(s)
  initial_probs = exp(s - s_max)
  initial_probs = initial_probs / sum(initial_probs)
  
  return (list(initial = initial_probs, transition = transition_probs))
}





pois_NHMM_lalphabeta = function(R, m, D, state_probs, lambda)
{
  #if (is.null(delta)) 
  #    delta = solve(t(diag(m) - gamma + 1), rep(1,m))
  n = length(R)
  lalpha = lbeta = matrix(NA, m, n)
  allprobs = outer(R, lambda[1:m], dpois)
  
  foo = state_probs$initial * allprobs[1,]    # foo is alpha_1(i) for i = 1..m
  sumfoo = sum(foo)
  lscale = log(sumfoo)
  foo = foo / sumfoo
  lalpha[,1] = log(foo) + lscale
  
  for (k in 2:n)
  {
    foo = foo %*% state_probs$transition[[k-1]] * allprobs[k,]
    sumfoo = sum(foo)
    lscale = lscale + log(sumfoo)
    foo = foo / sumfoo
    lalpha[,k] = log(foo) + lscale
  }
  
  lbeta[,n] = rep(0,m)
  foo = rep(1/m, m)
  lscale = log(m)
  for (k in (n-1):1)
  {
    foo = state_probs$transition[[k]] %*% (allprobs[k+1,] * foo)
    lbeta[,k] = log(foo) + lscale
    sumfoo = sum(foo)
    foo = foo / sumfoo
    lscale = lscale + log(sumfoo)
  }
  
  return (list(la = lalpha, lb = lbeta))
}






pois_NHMM_state_probs = function(R, X, m, D, lambda, theta)
{
  n = length(R)
  state_probs = compute_state_probs(m, n, X, theta)
  fb = pois_NHMM_lalphabeta(R, m, D, state_probs, lambda)
  la = fb$la
  lb = fb$lb
  c = max(la[,n])
  llk = c + log(sum(exp(la[,n] - c)))
  A = exp(la + lb - llk)
  A = A / apply(A, 2, sum)
  return (A)
  
  #     stateprobs = matrix(NA , ncol = n, nrow = m)
  #     for (i in 1:n) 
  #         stateprobs[,i] = exp(la[,i] + lb[,i] - llk)
  #     stateprobs
}



poisson_NHMM_EM = function(R, X, m, D, lambda, theta, maxiter = 30, tol = 1e-2, ...)
{
  # R is observed sequence
  # X[t,] is d-dimensional time covariate sequence
  # m is number of states
  # D is dimension of time covariate Xt
  # lambda_i is state-dependent Poisson distribution rate for modeling P(Rt | St = i)
  # gamma_ji is logistic function parameter for modeling P(St = i | St-1 = j, Xt = x)
  # delta_i is logistic function parameter for modeling P(S1 = i | X1 = x)
  # rho[i,] is d-dimensional state-dependent time covariate coefficients in 
  # .. logistic function for modeling P(St = i | St-1 = j, Xt = x) and P(S1 = i | X1 = x)
  
  cat(sprintf("******* fitting HMM with %s states ******", m), "\n")
  n = length(R)
  np = m*m + m + m*D + m
  
  lambda.next = lambda
  theta.next = theta
  
  for (iter in 1:maxiter)
  {
    cat("******* iter", iter)
    lallprobs = outer(R, lambda[1:m], dpois, log = TRUE)
    state_probs = compute_state_probs(m, n, X, theta)
    
    # forward-backward algorithm
    fb = pois_NHMM_lalphabeta(R, m, D, state_probs, lambda)
    la = fb$la
    lb = fb$lb
    c = max(la[,n])
    llk = c + log(sum(exp(la[,n] - c)))
    
    B = list(1:(n-1))
    for (k in 1:(n-1))
    {
      B[[k]] = state_probs$transition[[k]] * 
        exp(matrix(rep(la[,k], m), m, m) - llk +
              matrix(rep(lb[,k+1], m), m, m, byrow = TRUE) +
              matrix(rep(lallprobs[k+1,], m), m, m, byrow = TRUE) )
      B[[k]] = B[[k]] / sum(B[[k]])
    }
    
    A = exp(la + lb - llk)
    lambda.next = A %*% R / apply(A, 1, sum)
    A = A / apply(A, 2, sum)
    
    theta.next = NHMM_conjugate_gradient(X, m, D, n, A, B, state_probs, lambda, theta)
    #         theta_vector_optim = optim(par = c(theta$delta, theta$gamma, theta$rho), 
    #                                    fn = compute_Qs_theta_vector, 
    #                                    gr = gradient_Qs_theta_vector,
    #                                    m = m, D = D, n = n, A = A, B = B, method = "BFGS")
    #         theta.next = list(delta = theta_vector_optim$par[1:m],
    #                           gamma = matrix(theta_vector_optim$par[m + 1:(m*m)], m, m),
    #                           rho = matrix(theta_vector_optim$par[m + m*m + 1:(D*m)], m, D) )
    
    crit = sum(abs(lambda - lambda.next)) + 
      sum(abs(theta$gamma - theta.next$gamma)) + 
      sum(abs(theta$delta - theta.next$delta)) + 
      sum(abs(theta$rho - theta.next$rho))
    
    # for debugging
    if (is.nan(crit)) {
      print('poisson_NHMM_EM crit is nan')
      cat('lambda.next', lambda.next, '\n')
      cat('la', la, '\n')
      cat('lb', lb, '\n')
      cat('delta.next', theta.next$delta, '\n')
      cat('gamma.next', theta.next$gamma, '\n')
      cat('rho.next', theta.next$rho, '\n')
      return (NA)
    }
    
    AIC = -2*(llk - np)
    BIC = -2*llk + np * log(n)
    print(BIC)
    if (crit < tol)
    {
      return (list(lambda = lambda, delta = theta$delta, 
                   gamma = theta$gamma, rho = theta$rho,
                   mllk = -llk, AIC = AIC, BIC = BIC))
    }
    
    lambda = lambda.next
    theta = theta.next
  }
  
  print(paste("No convergence after ", maxiter, " iterations"))
  return (list(lambda = lambda, delta = theta$delta, 
               gamma = theta$gamma, rho = theta$rho,
               mllk = -llk, AIC = AIC, BIC = BIC))
}



compute_Qs = function(n, A, B, state_probs)
{
  Qs_last_term = list(1:(n-1))
  for (k in 1:(n-1))
  {
    #         if ( any(is.na(state_probs$transition[[k]])) | any(state_probs$transition[[k]] <= 0 ) )
    #         {
    #             print(state_probs$transition[[k]])
    #         }
    Qs_last_term[[k]] = B[[k]] * log(state_probs$transition[[k]])
  }
  
  #     if ( any(is.na(state_probs$transition[[k]])) | any(state_probs$initial <= 0 ) )
  #     {
  #         print(state_probs$initial)
  #     }
  Qs_first_term = A[,1] * log(state_probs$initial)
  Qs = sum(Qs_first_term) + sum(sapply( 1:(n-1), function(k) sum(Qs_last_term[[k]]) ))
  
  return (Qs)
}

compute_Qs_nu = function(nu, m, n, X, A, B, theta, phi)
{
  theta_nu = list(delta = theta$delta + nu * phi$delta,
                  gamma = theta$gamma + nu * phi$gamma,
                  rho = theta$rho + nu * phi$rho)
  
  state_probs = compute_state_probs(m, n, X, theta_nu)
  
  return (-compute_Qs(n, A, B, state_probs))
}



compute_dQs_dnu_nu = function(nu, m, n, X, A, B, theta, phi)
{
  theta_nu = list(delta = theta$delta + nu * phi$delta,
                  gamma = theta$gamma + nu * phi$gamma,
                  rho = theta$rho + nu * phi$rho)
  
  state_probs = compute_state_probs(m, n, X, theta_nu)
  
  return (compute_dQs_dnu(m, X, n, A, B, state_probs, theta))
}


gradient_Qs = function(m, D, n, A, B, state_probs)
{
  Qs_last_term = list(1:(n-1))
  for (k in 1:(n-1))
  {
    Qs_last_term[[k]] = B[[k]] * (1 - state_probs$transition[[k]])
  }
  
  grad_gamma = matrix(1:(m*m), m, m)
  grad_rho = matrix(rep(0, (m*D)), m, D)
  for (j in 1:m)
  {
    for (i in 1:m)
    {
      grad_gamma[j,i] = sum(sapply(1:(n-1), function(k) Qs_last_term[[k]][j,i]) )
      for (d in 1:D)
        grad_rho[i,d] = grad_rho[i,d] + 
          sum(sapply(1:(n-1), function(k) Qs_last_term[[k]][j,i] * X[k+1,d]) )
    }
  }
  
  grad_delta = A[,1] * (1 - state_probs$initial)
  
  return (list(delta = grad_delta, gamma = grad_gamma, rho = grad_rho))
}


compute_dQs_dnu = function(m, X, n, A, B, state_probs, theta)        
{
  # calculate first derivative d_Qs / d_nu
  Qs_last_term = list(1:(n-1))
  for (k in 1:(n-1))
  {
    Qs_last_term[[k]] = B[[k]] * (theta$gamma + 
                                    matrix(rep(X[k+1,] %*% t(theta$rho), m), m, m, byrow = TRUE)) *
      (1 - state_probs$transition[[k]])
  }
  
  Qs_first_term = A[,1] * (theta$delta + X[1,] %*% t(theta$rho) ) *
    (1 - state_probs$initial)
  
  Qs_first_derivative = sum(Qs_first_term) +
    sum(sapply( 1:(n-1), function(k) sum(Qs_last_term[[k]]) ))
  
  return (Qs_first_derivative)
}



NHMM_conjugate_gradient = function(X, m, D, n, A, B, state_probs, lambda, theta0,
                                   maxiter_conjugate = 30, tol = 1e-4, ...)
{
  theta = theta0
  gradQ = gradient_Qs(m, D, n, A, B, state_probs)
  phi = list(delta = -gradQ$delta, gamma = -gradQ$gamma, rho = -gradQ$rho)
  Qs = compute_Qs(n, A, B, state_probs)
  
  for (iter in 1:maxiter_conjugate)
  {
    #nu = newton_raphson(m, n, X, A, B, state_probs, theta, phi)
    #nu = uniroot(compute_dQs_dnu_nu, interval = c(-1, 1), 
    #             m = m, n = n, X = X, A = A, B = B, theta = theta, phi = phi)
    #nu = multiroot(compute_dQs_dnu_nu, start = 1, 
    #               m = m, n = n, X = X, A = A, B = B, theta = theta, phi = phi)
    #nl_min = nlm(compute_Qs_nu, p = 1,
    #             m = m, n = n, X = X, A = A, B = B, theta = theta, phi = phi)
    #nu = nl_min$estimate
    
    nu_optim = optim(1, compute_Qs_nu, gr=NULL, 
                     m = m, n = n, X = X, A = A, B = B, theta = theta, phi = phi, 
                     method = "L-BFGS-B")
    nu = nu_optim$par
    
    theta.next = list(delta = theta$delta + nu * phi$delta,
                      gamma = theta$gamma + nu * phi$gamma,
                      rho = theta$rho + nu * phi$rho)
    
    state_probs = compute_state_probs(m, n, X, theta.next)
    
    gradQ.next = gradient_Qs(m, D, n, A, B, state_probs)
    eta = sum( (gradQ.next$delta - gradQ$delta) * gradQ.next$delta,
               (gradQ.next$gamma - gradQ$gamma) * gradQ.next$gamma,
               (gradQ.next$rho - gradQ$rho) * gradQ.next$rho ) /
      sum( gradQ$delta * gradQ$delta, 
           gradQ$gamma * gradQ$gamma,
           gradQ$rho * gradQ$rho )
    phi.next = list(delta = gradQ$delta - eta * phi$delta, 
                    gamma = gradQ$gamma - eta * phi$gamma, 
                    rho = gradQ$rho - eta * phi$rho)
    
    Qs.next = compute_Qs(n, A, B, state_probs)
    crit = abs(Qs.next - Qs)
    
    # for debugging
    if (is.nan(crit)) {
      print('NHMM_conjugate_gradient crit is nan')
      cat('theta.next$delta', theta.next$delta, '\n')
      cat('theta.next$gamma', theta.next$gamma, '\n')
      cat('theta.next$rho', theta.next$rho, '\n')
      cat('gradQ.next$delta', gradQ.next$delta, '\n')
      cat('gradQ.next$gamma', gradQ.next$gamma, '\n')
      cat('gradQ.next$rho', gradQ.next$rho, '\n')
      cat('phi.next$delta', phi.next$delta, '\n')
      cat('phi.next$gamma', phi.next$gamma, '\n')
      cat('phi.next$rho', phi.next$rho, '\n')
      cat('Qs.next', Qs.next, '\n')
      return (NA)
    }
    if (crit < tol)
      return (theta)
    
    theta = theta.next
    gradQ = gradQ.next
    phi = phi.next
    Qs = Qs.next
  }
  print(paste("Conjugate gradient no convergence after ", maxiter_conjugate, " iterations"))
  return (theta)
}






pois_NHMM_forecast = function(R, X, m, D, lambda, theta, xrange = NULL, H = 1, ...)
{
  if (is.null(xrange))
    xrange <- qpois(0.001, min(lambda)) : qpois(0.999, max(lambda))
    xrange <- round(xrange)
  n = length(R)
  allprobs = outer(R, lambda, dpois)
  allprobs = ifelse (!is.na(allprobs), allprobs, 1)
  state_probs = compute_state_probs(m, n+H, X, theta)
  
  foo = state_probs$initial * allprobs[1,]
  sumfoo = sum(foo)
  lscale = log(sumfoo)
  foo = foo / sumfoo
  
  for (k in 2:n)
  {
    foo = foo %*% state_probs$transition[[k-1]] * allprobs[k,]
    sumfoo = sum(foo)
    lscale = lscale + log(sumfoo)
    foo = foo / sumfoo
  }
  
  xi = matrix(NA, nrow = m, ncol = H)
  for (k in 1:H)
  {
    foo = foo %*% state_probs$transition[[n+k-1]]
    xi[,k] = foo
  }
  
  allprobs = outer(xrange, lambda, dpois)
  fdists = allprobs %*% xi[, 1:H]
  
  list(xrange = xrange, fdists = fdists)
}






library(rootSolve)

# constants definitions
min_m = 2
max_m = 4



y <- read.table("C:\\Users\\kosti\\OneDrive\\Desktop\\Msc Thesis - Hidden Markov Models\\Rainfal_Brazil_Paper_Data\\precip_data_filled")
y <- y[,2]
y <- round(y)
y
sim_rainfall <- read.table("C:\\Users\\kosti\\OneDrive\\Desktop\\Msc Thesis - Hidden Markov Models\\Rainfal_Brazil_Paper_Data\\full28_sim_rainfall")

#remove the omitted years 1976,1978,1984,1986
rows_omit<- c(91:180,271:360,811:900,991:1080)
X<- sim_rainfall[-rows_omit, ]

X <- as.matrix(X)  



y


# Define parameters
num_truncate <- 2160# Number of observations
num_predict <- 30        # Forecast horizon
D <- 1             # Number of covariates (here, D = 1)

# Model selection
model_result <- list()
theta <- list()
for (m in min_m:max_m) {
  start_time <- Sys.time()
  cat('Number of states:', m, '\n')
  
  # Initialize parameters
  lambda0 <- round(runif(m, min = max(1, min(y[y > 0])), max = max(y)))
  
  delta0 <- rep(1 / m, m)  # Equal probability for each state initially
  gamma0 <- matrix(runif(m * m, min = 0.1, max = 0.9), m, m)
  gamma0 <- gamma0 / rowSums(gamma0)  # Normalize each row to sum to 1
  
  # Initialize rho0 with random values between 0.1 and 0.9
  rho0 <- matrix(runif(m * ncol(X), min = 0.1, max = 0.9), m, ncol(X))
  theta0 <- list(delta = delta0, gamma = gamma0, rho = rho0)
  
  # Fit the model
  p_nhmm <- poisson_NHMM_EM(y,X[1:num_truncate, , drop = FALSE], m, D, lambda0, theta0)
  
  # Check for errors
  if (any(is.na(p_nhmm))) {
    if (length(model_result$num_state) == 0) next else break
  }
  
  # Store results
  model_result$num_state <- c(model_result$num_state, m)
  model_result$aic_list <- c(model_result$aic_list, p_nhmm$AIC)
  model_result$bic_list <- c(model_result$bic_list, p_nhmm$BIC)
  model_result$lambda <- c(model_result$lambda, p_nhmm$lambda)
  model_result$delta <- c(model_result$delta, p_nhmm$delta)
  model_result$gamma <- c(model_result$gamma, p_nhmm$gamma)
  model_result$rho <- c(model_result$rho, p_nhmm$rho)
  model_result$mllk <- c(model_result$mllk, p_nhmm$mllk)
  
  finish_time <- Sys.time()
  print(finish_time - start_time)
}

# Extract optimal model
optimal_index <- which.min(model_result$bic_list)
m <- model_result$num_state[optimal_index]
lambda <- model_result$lambda[(m * (m - 1) / 2 - min_m * (min_m - 1) / 2) + 1:m]
theta <- list(
  delta = model_result$delta[(m * (m - 1) / 2 - min_m * (min_m - 1) / 2) + 1:m],
  gamma = matrix(model_result$gamma[
    (m * (m - 1) * (2 * m - 1) / 6 - min_m * (min_m - 1) * (2 * min_m - 1) / 6) + 1:(m * m)
  ], m, m),
  rho = matrix(model_result$rho[
    D * (m * (m - 1) / 2 - min_m * (min_m - 1) / 2) + 1:(D * m)
  ], m, D)
)

# Forecast
forecasts <- pois_NHMM_forecast(y, x, m, D, lambda, theta, H = num_predict)

# Print results
cat("\n")
cat("Number of HMM states fitted =", model_result$num_state, "\n\n")
cat("AIC =", model_result$aic_list, "\n")
cat("BIC =", model_result$bic_list, "\n\n")
cat("Chosen HMM Poisson lambda =", lambda, "\n\n")
print(theta)



