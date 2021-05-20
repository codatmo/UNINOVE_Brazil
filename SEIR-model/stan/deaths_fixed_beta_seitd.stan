functions {
  vector seitd(real time,
               vector state,
               real nu,
               real gamma,
               real kappa,
               real omega,
               real beta,
               real population,
               int T) {

    // Unpack state
    real S = state[1];
    real E = state[2];
    real I = state[3];
    real T_ = state[4];
    real D = state[5];

    real infection_rate = beta * I * S / population;
    real nuE = nu * E;
    real gammaI = gamma * I;
    real kappaT = kappa * T_;

    real dS_dt;
    real dE_dt;
    real dI_dt;
    real dT_dt;
    real dD_dt;

    dS_dt = -infection_rate;
    dE_dt = infection_rate - nuE;
    dI_dt = nuE - gammaI;
    dT_dt = gammaI * omega - kappaT;
    dD_dt = kappaT;

    return [dS_dt, dE_dt, dI_dt, dT_dt, dD_dt]';
  }
}
data {
  real initial_time;
  int<lower=1> T;
  real times[T];
  int<lower=1> n_disease_states;
  real<lower=0> population;
  int<lower=1> deaths_length;
  int<lower=1> deaths_starts[deaths_length];
  int<lower=1> deaths_stops[deaths_length];
  int<lower=0> deaths[deaths_length];
}
transformed data {
  real mu_dL = 4.00;
  real sigma_dL = 0.20;
  real mu_dI = 3.06;
  real sigma_dI = 0.21;
  real mu_dT = 16.00;
  real sigma_dT = 0.71;
}
parameters {
  real<lower=0, upper=1> initial_state_raw[2];
  real<lower=0> beta;
  real<lower=0> dL;
  real<lower=0> dI;
  real<lower=0> dT;
  real<lower=0, upper=1> omega;
  real<lower=0> reciprocal_phi_deaths;
}
transformed parameters {
  vector[n_disease_states] initial_state;
  real nu;
  real gamma;
  real kappa;
  real phi_deaths;
  vector[n_disease_states] state_estimate[T];
  vector[T+1] S;
  vector[T+1] E;
  vector[T+1] I;
  vector[T+1] T_;
  vector[T+1] D;
  vector[T] daily_infections;
  vector[T] daily_deaths;
  vector[T] effective_reproduction_number;

  initial_state[1] = (population-5.0)*initial_state_raw[1] + 1.0;
  initial_state[2] = (population-5.0)*(1.0-initial_state_raw[1])*initial_state_raw[2]/1.0 + 1.0;
  initial_state[3] = (population-5.0)*(1.0-initial_state_raw[1])*(1.0-initial_state_raw[2])/1.0 + 1.0;
  initial_state[4] = 0.0;
  initial_state[5] = 0.0;

  nu = 1.0/dL;
  gamma = 1.0/dI;
  kappa = 1.0/dT;
  phi_deaths = 1.0 / reciprocal_phi_deaths;

  state_estimate = ode_rk45(seitd,
                            initial_state, initial_time, times,
                            nu, gamma, kappa, omega, beta,
                            population, T);

  S = append_row(initial_state[1], to_vector(state_estimate[, 1]));
  E = append_row(initial_state[2], to_vector(state_estimate[, 2]));
  I = append_row(initial_state[3], to_vector(state_estimate[, 3]));
  T_ = append_row(initial_state[4], to_vector(state_estimate[, 4]));
  D = append_row(initial_state[5], to_vector(state_estimate[, 5]));

  daily_infections = S[:T] - S[2:] + 1E-4;
  daily_deaths = D[2:] - D[:T] + 1E-4; // I had problems with negativity
  effective_reproduction_number = (daily_infections ./ I[:T])*dI;
}

model {
  initial_state_raw[1] ~ beta(5.0, 0.5);
  initial_state_raw[2] ~ beta(1.1, 1.1);
  beta ~ normal(0, 0.5);
  dL ~ normal(mu_dL, sigma_dL);
  dI ~ normal(mu_dI, sigma_dI);
  dT ~ normal(mu_dT, sigma_dT);
  omega ~ beta(100, 9803);
  reciprocal_phi_deaths ~ exponential(5);

  for (i in 1:deaths_length) {
    target += neg_binomial_2_lpmf(deaths[i] |
              sum(daily_deaths[deaths_starts[i]:deaths_stops[i]]), phi_deaths);
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
