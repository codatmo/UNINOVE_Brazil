if (!require("DirichletReg")) {
    install.packages("DirichletReg")
    library(DirichletReg)
}

# Integrate state over time
integrateODE <- function(initial_state, initial_time, times, params) {
    h <- times[1] - initial_time
    d_state_dt_initial_time <- seeiittdRate(initial_time, initial_state, params)
    k <- h * d_state_dt_initial_time

    state_estimate <- matrix(rep(0, length(times) * length(initial_state)), length(times), length(initial_state))
    state_estimate[1, ] <- initial_state + h *
        (d_state_dt_initial_time + seeiittdRate(times[1], initial_state + k, params)) / 2

    for (t in 2:length(times)) {
        h <- times[t] - times[t - 1]

        d_state_dt <- seeiittdRate(times[t - 1], state_estimate[t - 1, ], params)
        k <- h * d_state_dt

        tmpRate <- seeiittdRate(times[t], state_estimate[t - 1, ] + k, params)
        if (length(tmpRate) != length(d_state_dt)) {
            browser() # What's going on?
            tmpRate <- seeiittdRate(times[t], state_estimate[t - 1, ] + k, params)
        }

        state_estimate[t, ] <- state_estimate[t - 1, ] + h * (d_state_dt + tmpRate) / 2
    }

    state_estimate
}

# Calculate the derivatives
seeiittdRate <- function(time, state, params) {
    s <- state[1]
    e1 <- state[2]
    e2 <- state[3]
    i1 <- state[4]
    i2 <- state[5]
    t1 <- state[6]
    t2 <- state[7]
    d <- state[8]

    betaPieceIndex <- which(time >= params$beta_left_t & time < params$beta_right_t)
    beta <- params$grad_beta[betaPieceIndex] *
        (time - params$beta_left_t[betaPieceIndex]) + params$beta_left[betaPieceIndex]

    dsdt <- -beta * (i1 + i2) / params$population * s
    de1dt <- beta * (i1 + i2) / params$population * s - 2 * e1 / params$dL
    de2dt <- 2 * (e1 - e2) / params$dL
    di1dt <- 2 / params$dL * e2 - 2 / params$dI * i1
    di2dt <- 2 / params$dI * (i1 - i2)
    drdt <- 2 / params$dI * i2 * (1 - params$omega)
    dt1dt <- 2 / params$dI * i2 * params$omega - 2 / params$dT * t1
    dt2dt <- 2 / params$dT * (t1 - t2)
    dddt <- 2 / params$dT * t2

    rate <- c(dsdt, de1dt, de2dt, di1dt, di2dt, dt1dt, dt2dt, dddt)

    if (sum(is.na(rate)) != 0) {
        browser()
    }

    rate
}


# Generate some simulated data with true values sampled from the prior
generateSimulatedData <- function(seed = 0,
                                  maxTime = 100,
                                  population = 1000000,
                                  n_beta_pieces = 20,
                                  output_tweets = TRUE) {

    set.seed(seed)
    initial_time <- 0
    n_disease_states <- 8

    beta_pieceLength <- ceiling(maxTime / n_beta_pieces)
    beta_left_t <- seq(0, maxTime - 1, by = beta_pieceLength)
    beta_right_t <- c(beta_left_t[2:n_beta_pieces], maxTime + 1)

    # Deaths are reported weekly, reporting periods are every seven days
    deathStarts <- seq(1, maxTime, by = 7)
    deathStops <- seq(7, maxTime, by = 7)
    nDeathCounts <- min(c(length(deathStarts), length(deathStops)))
    deathStarts <- deathStarts[1:nDeathCounts]
    deathStops <- deathStops[1:nDeathCounts]

    # Prior parameters for dL, dI, dT
    mu_dL <- 4.00
    sigma_dL <- 0.20
    mu_dI <- 3.06
    sigma_dI <- 0.21
    mu_dT <- 16.00
    sigma_dT <- 0.71

    max_lag <- 13

    # Draw from priors from model
    initial_state_raw <- c(rbeta(1, 5.0, 0.5), rbeta(1, 1.1, 1.1))

    # beta_left and beta_right are half-normal
    beta_left <- abs(rnorm(n_beta_pieces, 0, 0.5))
    beta_right <- abs(rnorm(n_beta_pieces, 0, 0.5))
    dL <- rnorm(1, mu_dL, sigma_dL)
    dI <- rnorm(1, mu_dI, sigma_dI)
    dT <- rnorm(1, mu_dT, sigma_dT)
    omega <- rbeta(1, 100, 9803)

    # rho_tweets are half-normal
    rho_tweets <- abs(rnorm(n_beta_pieces, 0, 0.5))

    reciprocal_phi_deaths <- rexp(1, 5)
    reciprocal_phi_tweets <- rexp(1, 5)


    lag_weights_tweets <- rdirichlet(1, rep(0.1, max_lag + 1))

    # Calculate initial state
    initial_state <- rep(0.0, 8)
    initial_state[1] <- (population - 5.0) * initial_state_raw[1] + 1.0
    initial_state[2] <- (population - 5.0) * (1 - initial_state_raw[1]) * initial_state_raw[2] / 2.0 + 1.0
    initial_state[3] <- (population - 5.0) * (1 - initial_state_raw[1]) * initial_state_raw[2] / 2.0 + 1.0
    initial_state[4] <- (population - 5.0) * (1 - initial_state_raw[1]) * (1 - initial_state_raw[2]) / 2.0 + 1.0
    initial_state[5] <- (population - 5.0) * (1 - initial_state_raw[1]) * (1 - initial_state_raw[2]) / 2.0 + 1.0

    grad_beta <- (beta_right - beta_left) / (beta_right_t - beta_left_t)
    phi_deaths <- 1 / reciprocal_phi_deaths
    phi_tweets <- 1 / reciprocal_phi_tweets

    params <- list(
        n_beta_pieces = n_beta_pieces,
        beta_left = beta_left,
        beta_right = beta_right,
        beta_left_t = beta_left_t,
        beta_right_t = beta_right_t,
        grad_beta = grad_beta,
        dL = dL,
        dI = dI,
        dT = dT,
        omega = omega,
        population = population
    )

    times <- seq(initial_time + 1, maxTime)

    # Simulate
    state_estimate <- integrateODE(initial_state, initial_time, times, params)
    state_estimate <- rbind(initial_state, state_estimate)

    # Calculate daily numbers of infections and numbers of deaths
    daily_infections <- as.numeric(state_estimate[1:nrow(state_estimate) - 1, 1]) -
        as.numeric(state_estimate[2:nrow(state_estimate), 1])
    daily_deaths <- as.numeric(state_estimate[2:nrow(state_estimate), 8]) -
        as.numeric(state_estimate[1:maxTime, 8])


    # Calculate daily calls given lagged infections
    tweets_lagged_daily_infections <- lag_weights_tweets[1] * daily_infections
    for (i in 1:max_lag) {
        tweets_lagged_daily_infections <- tweets_lagged_daily_infections +
            lag_weights_tweets[i + 1] * c(rep(0, i), daily_infections[1:(length(daily_infections) - i)])
    }


    daily_tweets <- rep(0, length(daily_infections))
    for (i in 1:n_beta_pieces) {
        daily_tweets[beta_left_t[i]:(beta_right_t[i] - 1)] <-
            tweets_lagged_daily_infections[beta_left_t[i]:(beta_right_t[i] - 1)] *
                rho_tweets[i]
    }

    # Sample weekly death reports
    weekly_deaths <- sapply(seq(1, nDeathCounts), function(i) {
        rnbinom(1, mu = sum(daily_deaths[deathStarts[i]:deathStops[i]]), size = phi_deaths)
    })

    # Sample numbers tweets weekly
    weekly_tweets <- sapply(seq(1, nDeathCounts), function(i) {
        rnbinom(1, mu = daily_tweets[i], size = phi_tweets)
    })

    # Sample numbers tweets weekly
     tweets <- sapply(seq_len(maxTime), function(i) {
         rnbinom(1, mu = daily_tweets[i], size = phi_tweets)
     })


    # Data for Stan model
    stan_data <- list(
        initial_time = initial_time,
        n_beta_pieces = n_beta_pieces,
        beta_left_t = beta_left_t,
        beta_right_t = beta_right_t,
        T = maxTime,
        times = times,
        n_disease_states = n_disease_states,
        population = population,
        deaths_length = length(weekly_deaths),
        deaths_starts = deathStarts,
        deaths_stops = deathStops,
        deaths = weekly_deaths,
        real_data_length = length(beta_left_t) + length(beta_right_t) + 1,
        real_data = c(beta_left_t, beta_right_t, population),
        integer_data_length = 5,
        integer_data = c(maxTime, length(beta_left_t), length(beta_right_t), length(beta_left_t), n_disease_states),
        tweets = tweets,
        weekly_tweets = weekly_tweets
    )

    # Return the ground truth values to compare against
    ground_truth <- list(
        initial_state_raw = initial_state_raw,
        beta_left = beta_left,
        beta_right = beta_right,
        dL = dL,
        dI = dI,
        dT = dT,
        omega = omega,
        reciprocal_phi_deaths = reciprocal_phi_deaths,
        phi_deaths = phi_deaths,
        initial_state = initial_state,
        grad_beta = grad_beta,
        state_estimate = state_estimate,
        daily_infections = daily_infections,
        daily_deaths = daily_deaths,
        tweets_lagged_daily_infections = tweets_lagged_daily_infections,
        daily_tweets = daily_tweets
    )

    list(stan_data = stan_data, ground_truth = ground_truth)
}
