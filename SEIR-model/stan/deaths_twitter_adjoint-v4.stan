
functions {

/* version which uses bisectioning search
 */
int find_interval_elem(real x, vector sorted, int start_ind) {
  int res;
  int N = num_elements(sorted);
  int max_iter = 100 * N;
  int left_ind = start_ind;
  int right_ind = N;
  real left = sorted[left_ind ] - x;
  real right = sorted[right_ind] - x;
  int iter = 1;

  if(N == 0) return(0);
  if(0 == left) return(left_ind);
  if(0 < left) return(left_ind-1);
  if(0 >= right) return(right_ind);

  while((right_ind - left_ind) > 1  && iter != max_iter) {
    int mid_ind;
    real mid;
    // is there a controlled way without being yelled at with a
    // warning?
    mid_ind = (left_ind + right_ind) / 2;
    mid = sorted[mid_ind] - x;
    if (mid == 0) return(mid_ind-1);
    if (left  * mid < 0) { right = mid; right_ind = mid_ind; }
    if (right * mid < 0) { left  = mid; left_ind  = mid_ind; }
    iter = iter + 1;
  }
  if(iter == max_iter)
    print("Maximum number of iterations reached.");
  return(left_ind);
}

  
  vector seeiittd(real time,
                  vector state,
                  vector beta_left,
                  vector grad_beta,
                  real nu,
                  real gamma,
                  real kappa,
                  real omega,
                  //vector beta_t,
                  vector beta_left_t,
                  vector beta_right_t,
                  real population,
                  int T) {

    // Unpack integer data values
    /*
    int T = integer_data[1];
    int n_beta_pieces = integer_data[2];
    int n_disease_states = integer_data[3];
    */

    int n_beta_pieces = num_elements(beta_left_t);
    

    // Unpack real data values
    /*
    vector[n_beta_pieces] beta_left_t = to_vector(real_data[1:n_beta_pieces]);
    vector[n_beta_pieces] beta_right_t = to_vector(real_data[n_beta_pieces+1:2*n_beta_pieces]);
    real population = real_data[2*n_beta_pieces+1];
    */
    
    // Unpack parameter values
    /*
    vector[n_beta_pieces] beta_left = params[1:n_beta_pieces];
    vector[n_beta_pieces] grad_beta = params[n_beta_pieces+1:2*n_beta_pieces];
    real nu = params[2*n_beta_pieces+1];
    real gamma = params[2*n_beta_pieces+2];
    real kappa = params[2*n_beta_pieces+3];
    real omega = params[2*n_beta_pieces+4];
    */

    // Unpack state
    real S = state[1];
    real E1 = state[2];
    real E2 = state[3];
    real I1 = state[4];
    real I2 = state[5];
    real T1 = state[6];
    real T2 = state[7];
    real D = state[8];

    //int idx = find_interval_elem(time, beta_t, 1); 

    real infection_rate;// = (grad_beta[idx] * (time - beta_t[idx]) + beta_left[idx]) * (I1 + I2) * S / population;
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

    // this should be coded using a binary search
    for (i in 1:n_beta_pieces) {
      if(time >= beta_left_t[i] && time < beta_right_t[i]) {
        real beta = grad_beta[i] * (time - beta_left_t[i]) + beta_left[i];
        infection_rate = beta * (I1 + I2) * S / population;
        // maybe we can do a break here? Once you found the right bin,
        // we can escape from the loop, no?
        break;
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
  vector[n_beta_pieces+1] beta_t;
  vector[n_beta_pieces] vbeta_left_t = to_vector(beta_left_t);
  vector[n_beta_pieces] vbeta_right_t = to_vector(beta_right_t);
  beta_t[1:n_beta_pieces] = to_vector(beta_left_t);
  beta_t[n_beta_pieces+1] = beta_right_t[n_beta_pieces];
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
  real<lower=0> reciprocal_phi_twitter_symptons;
  real<lower=0> rho_twitter_symptons[n_rho_twitter_symptons_pieces];
  simplex[max_lag+1] lag_weights_twitter_symptons;
}
transformed parameters {
  vector[n_disease_states] initial_state;
  vector[n_beta_pieces] grad_beta;
  real nu;
  real gamma;
  real kappa;
  real phi_deaths;
  real phi_twitter_symptons;
  vector[n_disease_states] state_estimate[T];
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
  vector[T] twitter_symptons_lagged_daily_infections;
  vector[T] daily_twitter_symptons;

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
  phi_twitter_symptons = 1.0 / reciprocal_phi_twitter_symptons;

  /*
  {
    vector[2*n_beta_pieces+4] params;
    params[1:n_beta_pieces] = beta_left;
    params[n_beta_pieces+1:2*n_beta_pieces] = grad_beta;
    params[2*n_beta_pieces+1] = nu;
    params[2*n_beta_pieces+2] = gamma;
    params[2*n_beta_pieces+3] = kappa;
    params[2*n_beta_pieces+4] = omega;

    state_estimate = ode_rk45(seeiittd, initial_state, initial_time, times, params, real_data, integer_data);
  }
  */

  /*
  state_estimate = ode_rk45_tol(seeiittd,
                                initial_state, initial_time, times,
                                1E-6, 1E-4, 10000,
                                beta_left,
                                grad_beta,
                                nu,
                                gamma,
                                kappa,
                                omega,
                                beta_t,
  */
                                /*
                                  beta_left_t,
                                  beta_right_t,
                                */
  /*
                                population,
                                T);
  */
  
  state_estimate 
            = ode_adjoint_tol_ctl(seeiittd,
                                  initial_state, initial_time, times,
                                  1E-6, rep_vector(1E-8, n_disease_states),
                                  1E-6, rep_vector(1E-5, n_disease_states),
                                  1E-6, 1E-5,
                                  10000000,
                                  150,
                                  1,
                                  2, 1,
                                  beta_left,
                                  grad_beta,
                                  nu,
                                  gamma,
                                  kappa,
                                  omega,
                                  //beta_t,
                                  vbeta_left_t,
                                  vbeta_right_t,
                                  population,
                                  T);


  S = append_row(initial_state[1], to_vector(state_estimate[, 1]));
  E1 = append_row(initial_state[2], to_vector(state_estimate[, 2]));
  E2 = append_row(initial_state[3], to_vector(state_estimate[, 3]));
  I1 = append_row(initial_state[4], to_vector(state_estimate[, 4]));
  I2 = append_row(initial_state[5], to_vector(state_estimate[, 5]));
  T1 = append_row(initial_state[6], to_vector(state_estimate[, 6]));
  T2 = append_row(initial_state[7], to_vector(state_estimate[, 7]));
  D = append_row(initial_state[8], to_vector(state_estimate[, 8]));

  //daily_infections = S[:T] - S[2:] + machine_precision();
  daily_infections = S[:T] - S[2:] + 1E-4;
  daily_deaths = D[2:] - D[:T] + 1E-4; // I had problems with negativity

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
  reciprocal_phi_twitter_symptons ~ exponential(5);
  rho_twitter_symptons ~ normal(0, 0.5);
  lag_weights_twitter_symptons ~ dirichlet(rep_vector(0.1, max_lag+1));

  if (compute_likelihood == 1) {
    for (i in 1:deaths_length) {
      target += neg_binomial_2_lpmf(deaths[i] |
                sum(daily_deaths[deaths_starts[i]:deaths_stops[i]]), phi_deaths);
    }
    
    target += neg_binomial_2_lpmf(twitter_symptons |
            daily_twitter_symptons[twitter_symptons_start:(twitter_symptons_start-1)+twitter_symptons_length], phi_twitter_symptons);
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
