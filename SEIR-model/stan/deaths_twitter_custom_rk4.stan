functions {
  vector seeiittd(real time,
                  vector state,
                  vector params,
                  real[] real_data,
                  int[] integer_data) {

    // Unpack integer data values
    int T = integer_data[1];
    int n_beta_pieces = integer_data[2];
    int n_disease_states = integer_data[3];

    // Unpack real data values
    vector[n_beta_pieces] beta_left_t = to_vector(real_data[1:n_beta_pieces]);
    vector[n_beta_pieces] beta_right_t = to_vector(real_data[n_beta_pieces+1:2*n_beta_pieces]);
    real population = real_data[2*n_beta_pieces+1];

    // Unpack parameter values
    vector[n_beta_pieces] beta_left = params[1:n_beta_pieces];
    vector[n_beta_pieces] grad_beta = params[n_beta_pieces+1:2*n_beta_pieces];
    real nu = params[2*n_beta_pieces+1];
    real gamma = params[2*n_beta_pieces+2];
    real kappa = params[2*n_beta_pieces+3];
    real omega = params[2*n_beta_pieces+4];

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

    for (i in 1:n_beta_pieces) {
      if(time >= beta_left_t[i] && time < beta_right_t[i]) {
        real beta = grad_beta[i] * (time - beta_left_t[i]) + beta_left[i];
        infection_rate = beta * (I1 + I2) * S / population;
      }
    }

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
  
  /** @addtogroup fixed_step Explicit fixed-step methods
  *
  * \ingroup odeint
  *  @{ */
  
  /** 
  * **Fourth-order Runge-Kutta method**
  *
  * Info: https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
  * Taken from: https://github.com/spinkney/helpful_stan_functions/blob/main/functions/ode/odeint_rk4.stan
  * @author Juho Timonen / adapted by Jose Storopoli
  *
  * @param t0 initial time
  * @param y0 initial state, D-vector
  * @param h step size (positive)
  * @param num_steps number of steps to take
  * @param a0 array of integer inputs given to `seeiittd()`
  * @param theta parameter vector given to `seeiittd()`
  * @return array of D-vectors, length equal to `num_steps + 1`
  */
vector[] odeint_rk4(real t0, vector y0, real h, int num_steps,
    real[] a0, vector theta, int[] integer_data){

  int d = num_elements(y0);
  vector[d] y[num_steps+1];
  vector[d] k1; vector[d] k2; vector[d] k3; vector[d] k4;
  real t = t0;
  y[1] = y0;
  
  // Integrate at grid of time points
  for(i in 1:num_steps){
    k1 = h * seeiittd(t          , y[i]           , theta, a0, integer_data);
    k2 = h * seeiittd(t + 0.5 * h, y[i] + 0.5 * k1, theta, a0, integer_data);
    k3 = h * seeiittd(t + 0.5 * h, y[i] + 0.5 * k2, theta, a0, integer_data);
    k4 = h * seeiittd(t + h      , y[i] + k3      , theta, a0, integer_data);
    y[i+1] = y[i] + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    t = t + h;
  }
  
  return(y);
}
}
data {
  real initial_time;
  int<lower=1> n_beta_pieces;
  real<lower=0> beta_left_t[n_beta_pieces];
  real<lower=0> beta_right_t[n_beta_pieces];
  int<lower=1> n_rho_twitter_symptons_pieces;
  int<lower=0> rho_twitter_symptons_left_t[n_rho_twitter_symptons_pieces];
  int<lower=0> rho_twitter_symptons_right_t[n_rho_twitter_symptons_pieces];
  int<lower=1> T;
  real times[T];
  int<lower=1> n_disease_states;
  real<lower=0> population;
  int<lower=1> deaths_length;
  int<lower=1> deaths_starts[deaths_length];
  int<lower=1> deaths_stops[deaths_length];
  int<lower=0> deaths[deaths_length];
  int<lower=1> twitter_symptons_length;
  int<lower=1> twitter_symptons_start;
  int<lower=0> twitter_symptons[twitter_symptons_length];
  int real_data_length;
  real real_data[real_data_length];
  int integer_data_length;
  int integer_data[integer_data_length];
  int<lower=0, upper=1> compute_likelihood;
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
  real<lower=0> beta_right[n_beta_pieces];
  real<lower=0> dL;
  real<lower=0> dI;
  real<lower=0> dT;
  real<lower=0, upper=1> omega;
  real<lower=0> reciprocal_phi_deaths;
}
transformed parameters {
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

  initial_state[1] = (population-5.0)*initial_state_raw[1] + 1.0;
  initial_state[2] = (population-5.0)*(1.0-initial_state_raw[1])*initial_state_raw[2]/2.0 + 1.0;
  initial_state[3] = (population-5.0)*(1.0-initial_state_raw[1])*initial_state_raw[2]/2.0 + 1.0;
  initial_state[4] = (population-5.0)*(1.0-initial_state_raw[1])*(1.0-initial_state_raw[2])/2.0 + 1.0;
  initial_state[5] = (population-5.0)*(1.0-initial_state_raw[1])*(1.0-initial_state_raw[2])/2.0 + 1.0;
  initial_state[6] = 0.0;
  initial_state[7] = 0.0;
  initial_state[8] = 0.0;
  grad_beta = (to_vector(beta_right) - to_vector(beta_left))./(to_vector(beta_right_t) -
              to_vector(beta_left_t));
  nu = 2.0/dL;
  gamma = 2.0/dI;
  kappa = 2.0/dT;
  phi_deaths = 1.0 / reciprocal_phi_deaths;

  {
    vector[2*n_beta_pieces+4] params;
    params[1:n_beta_pieces] = beta_left;
    params[n_beta_pieces+1:2*n_beta_pieces] = grad_beta;
    params[2*n_beta_pieces+1] = nu;
    params[2*n_beta_pieces+2] = gamma;
    params[2*n_beta_pieces+3] = kappa;
    params[2*n_beta_pieces+4] = omega;

    
    state_estimate = odeint_rk4(initial_time, initial_state, h, num_steps, real_data, params, integer_data);
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
  
  twitter_symptons_lagged_daily_infections = lag_weights_twitter_symptons[1]*daily_infections;

  for (i in 1:max_lag) {
    twitter_symptons_lagged_daily_infections += lag_weights_twitter_symptons[i+1]*
                                          append_row(rep_vector(0.0, i), daily_infections[:T-i]);
  }

  daily_twitter_symptons = rep_vector(0.0, T);

  for (i in 1:n_rho_twitter_symptons_pieces) {
    daily_twitter_symptons[rho_twitter_symptons_left_t[i]:rho_twitter_symptons_right_t[i]-1] =
    twitter_symptons_lagged_daily_infections[rho_twitter_symptons_left_t[i]:rho_twitter_symptons_right_t[i]-1] *
    rho_twitter_symptons[i];
  }
}

model {
  initial_state_raw[1] ~ beta(5.0, 0.5);
  initial_state_raw[2] ~ beta(1.1, 1.1);
  beta_left ~ normal(0, 0.5);
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
  int pred_twitter_symptons[twitter_symptons_length];

  for (i in 1:deaths_length) {
    pred_deaths[i] = neg_binomial_2_rng(sum(daily_deaths[deaths_starts[i]:deaths_stops[i]]),
                     phi_deaths);
  }

  for (i in 1:twitter_symptons_length) {
    pred_twitter_symptons[i] = neg_binomial_2_rng(daily_twitter_symptons[twitter_symptons_start - 1 + i],
                        phi_twitter_symptons);
  }
}
