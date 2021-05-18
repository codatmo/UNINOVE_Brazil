functions {
  vector seeiittd(real time,
                  vector state,
                  vector params,
                  real beta_left_t,
                  real beta_left,
                  real grad_beta,
                  real population) {
    real nu = params[1];
    real gamma = params[2];
    real kappa = params[3];
    real omega = params[4];

    // Unpack state
    real S = state[1];
    real E1 = state[2];
    real E2 = state[3];
    real I1 = state[4];
    real I2 = state[5];
    real T1 = state[6];
    real T2 = state[7];
    real D = state[8];

    real infection_rate;
    real nuE1 = nu * E1;
    real nuE2 = nu * E2;
    real gammaI1 = gamma * I1;
    real gammaI2 = gamma * I2;
    real kappaT1 = kappa * T1;
    real kappaT2 = kappa * T2;

    real dS_dt;
    real dE1_dt;
    real dE2_dt;
    real dI1_dt;
    real dI2_dt;
    real dT1_dt;
    real dT2_dt;
    real dD_dt;

    real beta = grad_beta * (time - beta_left_t) + beta_left;
    infection_rate = beta * (I1 + I2) * S / population;

    dS_dt = -infection_rate;
    dE1_dt = infection_rate - nuE1;
    dE2_dt = nuE1 - nuE2;
    dI1_dt = nuE2 - gammaI1;
    dI2_dt = gammaI1 - gammaI2;
    dT1_dt = gammaI2 * omega - kappaT1;
    dT2_dt = kappaT1 - kappaT2;
    dD_dt = kappaT2;

    return [dS_dt, dE1_dt, dE2_dt, dI1_dt, dI2_dt, dT1_dt, dT2_dt, dD_dt]';
  }

// https://github.com/spinkney/helpful_stan_functions/blob/main/functions/ode/odeint_rk4.stan
vector[] odeint_rk4(real t0, vector y0, real hh, int num_steps, int num_sub_steps,
    vector params, vector beta_left_t, vector beta_right_t,
    vector beta_left, vector grad_beta, real population
){

  int d = num_elements(y0);
  vector[d] y[num_steps+1];
  vector[d] k1; vector[d] k2; vector[d] k3; vector[d] k4;
  real t = t0;
  int k = 1;
  real h = hh/num_sub_steps;
  y[1] = y0;

  // Integrate at grid of time points
  for(i in 1:num_steps){
    y[i+1] = y[i];
    for(j in 1:num_sub_steps){
      if (t >= beta_right_t[k])
        k += 1;
      k1 = h * seeiittd(t          , y[i+1]           , params, beta_left_t[k], beta_left[k], grad_beta[k], population);
      k2 = h * seeiittd(t + 0.5 * h, y[i+1] + 0.5 * k1, params, beta_left_t[k], beta_left[k], grad_beta[k], population);
      k3 = h * seeiittd(t + 0.5 * h, y[i+1] + 0.5 * k2, params, beta_left_t[k], beta_left[k], grad_beta[k], population);
      k4 = h * seeiittd(t + h      , y[i+1] + k3      , params, beta_left_t[k], beta_left[k], grad_beta[k], population);
      y[i+1] = y[i+1] + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
      t = t + h;
    }
  }

  return(y);
}
}
data {
  real initial_time;
  int<lower=1> n_beta_pieces;
  vector<lower=0>[n_beta_pieces] beta_left_t;
  vector<lower=0>[n_beta_pieces] beta_right_t;
  int<lower=1> T;
  real times[T];
  int<lower=1> n_disease_states;
  real<lower=0> population;
  int<lower=1> deaths_length;
  int<lower=1> deaths_starts[deaths_length];
  int<lower=1> deaths_stops[deaths_length];
  int<lower=0> deaths[deaths_length];
  int real_data_length;
  real real_data[real_data_length];
  int integer_data_length;
  int integer_data[integer_data_length];
  int<lower=0, upper=1> compute_likelihood;
  int beta_linear_discontinuous;
  int beta_linear_continuous;
  int beta_constant_discontinuous;
  real beta_regularization;
  int num_sub_steps;
}
transformed data {
  real mu_dL = 4.00;
  real sigma_dL = 0.20;
  real mu_dI = 3.06;
  real sigma_dI = 0.21;
  real mu_dT = 16.00;
  real sigma_dT = 0.71;
  int max_lag = 13;
  real h = times[1] - initial_time;
  int num_steps = T;
}
parameters {
  real<lower=0, upper=1> initial_state_raw[2];
  vector<lower=0>[n_beta_pieces] beta_left;
  vector<lower=0>[beta_linear_discontinuous ? n_beta_pieces : 0] beta_right;
  real<lower=0> dL;
  real<lower=0> dI;
  real<lower=0> dT;
  real<lower=0, upper=1> omega;
  real<lower=0> reciprocal_phi_deaths;
  vector<lower=0>[beta_regularization < 0 ? 1 : 0] beta_reg;
}
transformed parameters {
  vector<lower=0>[n_beta_pieces] dummy_beta_right;
  real dummy_beta_reg = beta_regularization < 0 ? beta_reg[1] : beta_regularization;
  vector[n_disease_states] initial_state;
  vector[n_beta_pieces] grad_beta;
  real nu;
  real gamma;
  real kappa;
  real phi_deaths;
  vector[n_disease_states] state_estimate[T+1];
  vector[T+1] S;
  vector[T+1] E1;
  vector[T+1] E2;
  vector[T+1] I1;
  vector[T+1] I2;
  vector[T+1] T1;
  vector[T+1] T2;
  vector[T+1] D;
  vector[T] daily_infections;
  vector[T] daily_deaths;
  vector[T] effective_reproduction_number;
  if(beta_linear_discontinuous){
    dummy_beta_right = beta_right;
  } else if(beta_linear_continuous){
    if(n_beta_pieces > 1){
      dummy_beta_right[:n_beta_pieces-1] = beta_left[2:];
    }
    dummy_beta_right[n_beta_pieces] = beta_left[n_beta_pieces];
  } else if(beta_constant_discontinuous){
    dummy_beta_right = beta_left;
  } else{
    reject("Select an ansatz!");
  }

  initial_state[1] = (population-5.0)*initial_state_raw[1] + 1.0;
  initial_state[2] = (population-5.0)*(1.0-initial_state_raw[1])*initial_state_raw[2]/2.0 + 1.0;
  initial_state[3] = (population-5.0)*(1.0-initial_state_raw[1])*initial_state_raw[2]/2.0 + 1.0;
  initial_state[4] = (population-5.0)*(1.0-initial_state_raw[1])*(1.0-initial_state_raw[2])/2.0 + 1.0;
  initial_state[5] = (population-5.0)*(1.0-initial_state_raw[1])*(1.0-initial_state_raw[2])/2.0 + 1.0;
  initial_state[6] = 0.0;
  initial_state[7] = 0.0;
  initial_state[8] = 0.0;
  grad_beta = (dummy_beta_right - beta_left) ./ (beta_right_t - beta_left_t);
  nu = 2.0/dL;
  gamma = 2.0/dI;
  kappa = 2.0/dT;
  phi_deaths = 1.0 / reciprocal_phi_deaths;

  {
    vector[4] params = [nu, gamma, kappa, omega]';
    state_estimate = odeint_rk4(
        initial_time, initial_state, h, num_steps, num_sub_steps,
        params, beta_left_t, beta_right_t, beta_left, grad_beta, population
    );
  }

  S = to_vector(state_estimate[, 1]);
  E1 = to_vector(state_estimate[, 2]);
  E2 = to_vector(state_estimate[, 3]);
  I1 = to_vector(state_estimate[, 4]);
  I2 = to_vector(state_estimate[, 5]);
  T1 = to_vector(state_estimate[, 6]);
  T2 = to_vector(state_estimate[, 7]);
  D = to_vector(state_estimate[, 8]);

  daily_infections = S[:T] - S[2:] + machine_precision();
  daily_deaths = D[2:] - D[:T]     + machine_precision();

  {
    vector[T+1] I = I1 + I2;
    effective_reproduction_number= (daily_infections ./ I[:T])*dI;
  }
}

model {
  initial_state_raw[1] ~ beta(5.0, 0.5);
  initial_state_raw[2] ~ beta(1.1, 1.1);
  beta_left ~ normal(0, 0.5);
  if(beta_regularization < 0){
    dummy_beta_reg ~ lognormal(log(-beta_regularization), 1);
  }
  if(beta_regularization && n_beta_pieces > 1){
    if (beta_linear_discontinuous){
      //Not clear what is appropriate here
      dummy_beta_right ~ normal(beta_left, dummy_beta_reg);
      beta_left[2:] ~ normal(dummy_beta_right[:n_beta_pieces-1], dummy_beta_reg);
    }else{
      beta_left[2:] ~ normal(beta_left[:n_beta_pieces-1], dummy_beta_reg);
    }
  }
  beta_right ~ normal(0, 0.5);
  dL ~ normal(mu_dL, sigma_dL);
  dI ~ normal(mu_dI, sigma_dI);
  dT ~ normal(mu_dT, sigma_dT);
  omega ~ beta(100, 9803);
  reciprocal_phi_deaths ~ exponential(5);

  if (compute_likelihood == 1) {
    for (i in 1:deaths_length) {
      target += neg_binomial_2_lpmf(deaths[i] |
                sum(daily_deaths[deaths_starts[i]:deaths_stops[i]]), phi_deaths);
    }
  }
}
generated quantities {
  vector[T-1] growth_rate = (log(daily_infections[2:]) - log(daily_infections[:T-1]))*100;

  int pred_deaths[deaths_length];

  for (i in 1:deaths_length) {
    pred_deaths[i] = neg_binomial_2_rng(sum(daily_deaths[deaths_starts[i]:deaths_stops[i]]),
                     phi_deaths);
  }
}