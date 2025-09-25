library(dplyr)
library(tidyr)

simulateData <- function(
    params = list(
      data = list(
        alpha_lamb = c(1:10, seq(10, 120, by = 4), seq(120, 3, by = -6)),
        beta_lamb  = 0.5,
        T       = 30,
        date_start = as.Date("2024-01-01"),
        D = 20
      ),
      q_model = list(
        method       = "exp",  # "nonparam", "exp", "exp_rw", "exp_ou"
        method_params= list(
          # nonparam: b, phi
          # exp: b, phi
          # exp_rw: mu_logb, sigma_logb, mu_logitphi, sigma_logitphi
          # exp_ou: init_logb, mu_logb, sigma_logb, theta_logb,
          #     theta_logitphi, init_logitphi, mu_logitphi, sigma_logitphi
          b = 0.5, phi = 0.2
        )
      )
    )
){
  # (A) data
  alpha_lamb  <- params$data$alpha_lamb
  beta_lamb   <- params$data$beta_lamb
  T       <- params$data$T
  date_start  <- params$data$date_start
  D           <- params$data$D

  # (B) model for qd
  method          <- params$q_model$method
  method_params   <- params$q_model$method_params

  if (!(length(alpha_lamb) == T)) {
    stop("The length of `alpha_lamb` should be equal to the length of `T`！")
  }

  simsQ_out <- generateQ(method = method, params = method_params, T = T, D = D)

  simulation_result <- runSimulation(
    alpha_lamb = alpha_lamb, beta_lamb  = beta_lamb, T = T, date_start = date_start,
    simsQ_out  = simsQ_out, D = D)

  return(simulation_result)
}

runSimulation <- function(
    alpha_lamb,
    beta_lamb,
    T,
    date_start,
    simsQ_out,
    D
) {
  # Generate the date sequence
  date_seq <- seq.Date(from = date_start, by = "day", length.out = T)

  # Initialize variables
  lambda     <- numeric(T)        # Disease intensity
  case_true    <- integer(T)        # True number of cases
  case_reported<- matrix(0, nrow = T, ncol = D + 1)  # Reported cases
  rownames(case_reported) <- as.character(date_seq)

  # Simulation process
  for (tt in seq_len(T)){

    # 1) λ_t ~ Gamma
    lambda[tt] <- rgamma(1, shape = alpha_lamb[tt], rate = beta_lamb)

    # 2) True number of cases
    case_true[tt] <- rpois(1, lambda = lambda[tt])

    # 3) Get the reporting proportion for the current time point
    if (is.vector(simsQ_out$q)){
      # If q is a vector (constant model)
      prob_temp <- simsQ_out$q
    } else {
      # If q is a matrix (time-varying model, etc.)
      prob_temp <- simsQ_out$q[tt, ]
    }

    # 4) Calculate single-day reporting proportions
    p_temp <- c(prob_temp[1], diff(prob_temp), 1 - prob_temp[D+1])

    # 5) Distribute the true cases to each delay day
    case_reported[tt, ] <- rmultinom(n = 1, size = case_true[tt], prob = p_temp)[1:(D+1)]

  }

  qd_out <- if(is.vector(simsQ_out$q)){
    simsQ_out$q[1:(D+1)]
  } else {
    simsQ_out$q[, 1:(D+1)]
  }

  # Convert true cases to matrix and set row names
  case_true = as.matrix(case_true)
  rownames(case_true) = as.character(date_seq)

  # cumulated reported cases
  case_reported_cumulated <- t(apply(case_reported, 1, cumsum))

  # Return the final result list
  return(list(
    # Reported cases
    case_reported = case_reported[, c(1:(D+1))],
    case_reported_cumulated = case_reported_cumulated[, c(1:(D+1))],
    # True cases
    case_true  = case_true,
    # Disease intensity
    lambda   = round(lambda, 4),
    # b(t) parameter
    b        = round(simsQ_out$b, 4),
    # intercept
    phi = round(simsQ_out$phi, 4),
    # q(d) reporting proportion
    q         = round(qd_out, 4),
    # Date sequence
    date_seq   = date_seq,
    # Used D value
    D     = D
  ))
}

generateQ <- function(method, params, T, D) {
  #----------------------------------------------------------------
  # method:         "nonparam", "exp", "exp_rw", "exp_ou"
  # params:  a list of parameters required by each method
  # T:          number of observations (time points)
  # D:              maximum delay used in certain methods
  #
  # Returns a list:
  #   $q    : either a matrix (T x (max_delay+1)) or a vector
  #   $b_t   : vector (length T) of b(t)
  #   $phi   : vector (length T) of phi(t)
  #----------------------------------------------------------------

  if (method %in% c("nonparam", "exp")) {

    if (!all(c("b", "phi") %in% names(params))) {
      stop("method=nonparam,exp requires 'b' and 'phi' in params!")
    }
    b <- params$b
    phi <- params$phi
    q <- 1 - phi * exp(-b * (0:D))

  } else if (method == "exp_rw") {

    if (!all(c("mu_logb", "sigma_logb", "mu_logitphi", "sigma_logitphi") %in% names(params))) {
      stop("method=exp_rw requires mu_logb, sigma_logb, mu_logitphi, sigma_logitphi!")
    }

    logb <- arima.sim(model = list(order = c(0, 1, 0)), sd = params$sigma_logb,
            n = T, n.start = 100) |> as.numeric()
    logitphi <- arima.sim(model = list(order = c(0, 1, 0)), sd = params$sigma_logitphi,
            n = T, n.start = 100) |> as.numeric()

    b <- exp(params$mu_logb + logb - mean(logb))
    phi <- plogis(params$mu_logitphi + logitphi - mean(logitphi))

    q <- matrix(NA, nrow = T, ncol = D + 1)
    for (i in seq_len(T)) {
      q[i, ] <- 1 - phi[i] * exp(-b[i] * (0:D))
    }

  } else if (method == "exp_ou") {

    if (!all(c("init_logb", "mu_logb", "sigma_logb", "theta_logitphi",
            "init_logitphi", "mu_logitphi", "sigma_logitphi", "theta_logb") %in% names(params))) {
      stop("method=exp_ou requires init_logb, mu_logb, sigma_logb, theta_logb, init_logitphi, mu_logitphi, sigma_logitphi, theta_logitphi!")
    }

    logb <- logitphi <- numeric(T)
    logb[1] <- params$init_logb
    logitphi[1] <- params$init_logitphi

    for (i in 2:T) {
      drift_b_star <- params$theta_logb * (params$mu_logb - logb[i-1])
      logb[i] <- logb[i-1] + drift_b_star + rnorm(1, 0, params$sigma_logb)

      drift_phi_star <- params$theta_logitphi * (params$mu_logitphi - logitphi[i-1])
      logitphi[i] <- logitphi[i-1] + drift_phi_star + rnorm(1, 0, params$sigma_logitphi)
    }

    b <- exp(logb)
    phi <- plogis(logitphi)

    q <- matrix(NA, nrow = T, ncol = D + 1)
    for (i in seq_len(T)) {
      q[i, ] <- 1 - phi[i] * exp(-b[i] * (0:D))
    }

  } else {
    stop("method must be one of: 'nonparam', 'exp', 'exp_rw', 'exp_ou'!")
  }

  return(list(q = q, b = b, phi = phi))
}
