functions {

  vector seeiittd(real time,
                  vector state,
                  real beta_left,
                  real grad_beta,
                  real nu,
                  real gamma,
                  real kappa,
                  real omega,
                  real beta_t,
                  /*
                  real[] beta_left_t,
                  real[] beta_right_t,
                  */
                  real population) {

    // Unpack integer data values
    /*
    int T = integer_data[1];
    int n_beta_pieces = integer_data[2];
    int n_disease_states = integer_data[3];
    */

    //int n_beta_pieces = num_elements(beta_left_t);
    

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

    //real infection_rate = (grad_beta[idx] * (time - beta_t[idx]) + beta_left[idx]) * (I1 + I2) * S / population;
    real infection_rate = (grad_beta * (time - beta_t) + beta_left) * (I1 + I2) * S / population;
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
    /*
    for (i in 1:n_beta_pieces) {
      if(time >= beta_left_t[i] && time < beta_right_t[i]) {
        real beta = grad_beta[i] * (time - beta_left_t[i]) + beta_left[i];
        infection_rate = beta * (I1 + I2) * S / population;
        // maybe we can do a break here? Once you found the right bin,
        // we can escape from the loop, no?
        // break
      }
    }
    */

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
}
transformed data {
  real mu_dL = 4.00;
  real sigma_dL = 0.20;
  real mu_dI = 3.06;
  real sigma_dI = 0.21;
  real mu_dT = 16.00;
  real sigma_dT = 0.71;
  int max_lag = 13;
  real all_times[T+1];
  vector[n_beta_pieces+1] beta_t;
  beta_t[1:n_beta_pieces] = to_vector(beta_left_t);
  beta_t[n_beta_pieces+1] = beta_right_t[n_beta_pieces];
  all_times[1] = initial_time;
  all_times[2:(T+1)] = times;
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

  {
    vector[n_disease_states] last_state = initial_state;
    int j = 1;
    for(i in 1:n_beta_pieces) {
      int k = 0;
      if (times[j] == beta_t[i]) {
        state_estimate[j] = last_state;
        j += 1;
      }
      while(j+k <= T && times[j + k] <= beta_t[i+1]) {
        k += 1;
      }
      /*
      print("i = ", i);
      print("k = ", k);
      print("j = ", j);
      print("beta_t[i] = ", beta_t[i]);
      print("times = ", times[j:(j+k-1)]);
      print("beta_t[i+1] = ", beta_t[i+1]);
      */
      {
        /*
        vector[n_disease_states] state_piece[k+1] = ode_bdf_tol(seeiittd,
                          last_state, beta_t[i],
                          append_array(times[j:(j+k-1)], to_array_1d(beta_t[i+1:i+1])),
                          1E-6, 1E-6, 10000000,
                          beta_left[i],
                          grad_beta[i],
                          nu,
                          gamma,
                          kappa,
                          omega,
                          beta_t[i],
                          population);
        */
        /**/
        vector[n_disease_states] state_piece[k+1]
            = ode_adjoint_tol_ctl(seeiittd,
                                  last_state, beta_t[i],
                                  append_array(times[j:(j+k-1)], to_array_1d(beta_t[i+1:i+1])),
                                  1E-6, rep_vector(1E-8, n_disease_states),
                                  1E-6, rep_vector(1E-6, n_disease_states),
                                  1E-6, 1E-6,
                                  10000000,
                                  150,
                                  1,
                                  2, 2,
                                  beta_left[i],
                                  grad_beta[i],
                                  nu,
                                  gamma,
                                  kappa,
                                  omega,
                                  beta_t[i],
                                  population);

        /**/
        for(l in 1:k)
          state_estimate[j + l - 1] = state_piece[l];
        last_state = state_piece[k+1];
        j += k;
      }
    }
  }
  //print("state_estimate = ", state_estimate[T-1:T]);


  S = append_row(initial_state[1], to_vector(state_estimate[, 1]));
  E1 = append_row(initial_state[2], to_vector(state_estimate[, 2]));
  E2 = append_row(initial_state[3], to_vector(state_estimate[, 3]));
  I1 = append_row(initial_state[4], to_vector(state_estimate[, 4]));
  I2 = append_row(initial_state[5], to_vector(state_estimate[, 5]));
  T1 = append_row(initial_state[6], to_vector(state_estimate[, 6]));
  T2 = append_row(initial_state[7], to_vector(state_estimate[, 7]));
  D = append_row(initial_state[8], to_vector(state_estimate[, 8]));

  daily_infections = S[:T] - S[2:] + machine_precision();
  daily_deaths = D[2:] - D[:T];

  {
    vector[T+1] I = I1 + I2;
    effective_reproduction_number= (daily_infections ./ I[:T])*dI;
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

  for (i in 1:deaths_length) {
    pred_deaths[i] = neg_binomial_2_rng(sum(daily_deaths[deaths_starts[i]:deaths_stops[i]]),
                     phi_deaths);
  }
}
