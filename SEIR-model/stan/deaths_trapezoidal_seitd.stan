functions {
  real[] seitd(real time,
                  real[] state,
                  real[] params,
                  real[] real_data,
                  int[] integer_data) {

    // Unpack integer data values
    int T = integer_data[1];
    int n_beta_pieces = integer_data[2];
    int n_disease_states = integer_data[3];

    // Unpack real data values
    real beta_left_t[n_beta_pieces] = real_data[1:n_beta_pieces];
    real beta_right_t[n_beta_pieces] = real_data[n_beta_pieces+1:2*n_beta_pieces];
    real population = real_data[2*n_beta_pieces+1];

    // Unpack parameter values
    real beta_left[n_beta_pieces] = params[1:n_beta_pieces];
    real grad_beta[n_beta_pieces] = params[n_beta_pieces+1:2*n_beta_pieces];
    real nu = params[2*n_beta_pieces+1];
    real gamma = params[2*n_beta_pieces+2];
    real kappa = params[2*n_beta_pieces+3];
    real omega = params[2*n_beta_pieces+4];

    // Unpack state
    real S = state[1];
    real E = state[2];
    real I = state[3];
    real T_ = state[4];
    real D = state[5];

    real infection_rate;
    real nuE = nu * E;
    real gammaI = gamma * I;
    real kappaT = kappa * T_;

    real dS_dt;
    real dE_dt;
    real dI_dt;
    real dT_dt;
    real dD_dt;

    for (i in 1:n_beta_pieces) {
      if(time >= beta_left_t[i] && time < beta_right_t[i]) {
        real beta = grad_beta[i] * (time - beta_left_t[i]) + beta_left[i];
        infection_rate = beta * I * S / population;
      }
    }

    dS_dt = -infection_rate;
    dE_dt = infection_rate - nuE;
    dI_dt = nuE - gammaI;
    dT_dt = gammaI * omega - kappaT;
    dD_dt = kappaT;

    return {dS_dt, dE_dt, dI_dt, dT_dt, dD_dt};
  }

  real[ , ] integrate_ode_explicit_trapezoidal(real[] initial_state, real initial_time, real[] times, real[] params, real[] real_data, int[] integer_data) {
    real h;
    vector[size(initial_state)] dstate_dt_initial_time;
    vector[size(initial_state)] dstate_dt_tidx;
    vector[size(initial_state)] k;
    real state_estimate[size(times),size(initial_state)];

    h = times[1] - initial_time;
    dstate_dt_initial_time = to_vector(seitd(initial_time, initial_state, params, real_data, integer_data));
    k = h*dstate_dt_initial_time;
    state_estimate[1,] = to_array_1d(to_vector(initial_state) + h*(dstate_dt_initial_time + to_vector(seitd(times[1], to_array_1d(to_vector(initial_state)+k), params, real_data, integer_data)))/2);

    for (tidx in 1:size(times)-1) {
      h = (times[tidx+1] - times[tidx]);
      dstate_dt_tidx = to_vector(seitd(times[tidx], state_estimate[tidx], params, real_data, integer_data));
      k = h*dstate_dt_tidx;
      state_estimate[tidx+1,] = to_array_1d(to_vector(state_estimate[tidx,]) + h*(dstate_dt_tidx + to_vector(seitd(times[tidx+1], to_array_1d(to_vector(state_estimate[tidx,])+k), params, real_data, integer_data)))/2);
    }

    return state_estimate;
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
}
parameters {
  real<lower=0, upper=1> initial_state_raw[2];
  real<lower=0> beta_left[n_beta_pieces];
  real<lower=0> beta_right[n_beta_pieces];
  real<lower=0> dL;
  real<lower=0> dI;
  real<lower=0> dT;
  real<lower=0, upper=1> omega;
  real<lower=0> reciprocal_phi_deaths;
}
transformed parameters {
  real initial_state[n_disease_states];
  real grad_beta[n_beta_pieces];
  real nu;
  real gamma;
  real kappa;
  real phi_deaths;
  real state_estimate[T,n_disease_states];
  vector[T+1] S;
  vector[T+1] E;
  vector[T+1] I;
  vector[T+1] T_;
  vector[T+1] D;
  vector[T] daily_infections;
  vector[T] daily_deaths;
  vector[T] effective_reproduction_number;

  initial_state[1] = (population-5.0)*initial_state_raw[1] + 1.0;
  initial_state[2] = (population-5.0)*(1.0-initial_state_raw[1])*initial_state_raw[2] + 1.0;
  initial_state[3] = (population-5.0)*(1.0-initial_state_raw[1])*(1.0-initial_state_raw[2]) + 1.0;
  initial_state[4] = 0.0;
  initial_state[5] = 0.0;
  grad_beta = to_array_1d((to_vector(beta_right) - to_vector(beta_left))./(to_vector(beta_right_t) -
              to_vector(beta_left_t)));
  nu = 1.0/dL;
  gamma = 1.0/dI;
  kappa = 1.0/dT;
  phi_deaths = 1.0 / reciprocal_phi_deaths;

  {
    real params[2*n_beta_pieces+4];
    params[1:n_beta_pieces] = beta_left;
    params[n_beta_pieces+1:2*n_beta_pieces] = grad_beta;
    params[2*n_beta_pieces+1] = nu;
    params[2*n_beta_pieces+2] = gamma;
    params[2*n_beta_pieces+3] = kappa;
    params[2*n_beta_pieces+4] = omega;

    state_estimate = integrate_ode_explicit_trapezoidal(initial_state, initial_time, times, params, real_data, integer_data);
  }

  S = append_row(initial_state[1], to_vector(state_estimate[, 1]));
  E = append_row(initial_state[2], to_vector(state_estimate[, 2]));
  I = append_row(initial_state[3], to_vector(state_estimate[, 3]));
  T_ = append_row(initial_state[4], to_vector(state_estimate[, 4]));
  D = append_row(initial_state[5], to_vector(state_estimate[, 5]));

  daily_infections = S[:T] - S[2:] + machine_precision();
  daily_deaths = D[2:] - D[:T];
  effective_reproduction_number = (daily_infections ./ I[:T])*dI;
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
  vector[deaths_length] log_lik;
  for (i in 1:deaths_length) {
    log_lik[i] = neg_binomial_2_lpmf(deaths[i] |
                 sum(daily_deaths[deaths_starts[i]:deaths_stops[i]]), phi_deaths);
  }
}
