// https://discourse.mc-stan.org/t/codatmo-liverpool-uninove-models-slow-ode-implementation-and-trapezoidal-solver/22500/45

data {
  int no_days;
  int population;
  //observed new deaths per day
  int new_deaths[no_days];
  real beta_regularization;
  real likelihood;
}
parameters {
  /*
    The elements of a simplex sum to one.
    The first no_days entries represent daily new infections as a fraction
    of the total population while the last entry represents the proportion that
    remains susceptible at the end.
  */
  simplex[no_days+1] unit_dS;
  real<lower=0> dL;
  real<lower=0> dI;
  real<lower=0> dT;
  real<lower=0, upper=1> omega;
  real<lower=0> reciprocal_phi_deaths;
}
transformed parameters {
  vector[no_days] daily_infections = population * unit_dS[:no_days];
  vector[no_days] daily_deaths;
  vector[no_days] beta;
  vector[no_days] effective_reproduction_number;
  if(likelihood){
    vector[7] state = [
        0, 0,
        0, 0,
        0, 0,
        0
    ]';
    matrix[7, 7] transition_matrix = matrix_exp([
    //[E1   ,E2   ,I1   ,I2         ,T1   ,T2   ,D]
      [-2/dL,0    ,0    ,0          ,0    ,0    ,0],//E1
      [+2/dL,-2/dL,0    ,0          ,0    ,0    ,0],//E2
      [0    ,+2/dL,-2/dI,0          ,0    ,0    ,0],//I1
      [0    ,0    ,+2/dI,-2/dI      ,0    ,0    ,0],//I2
      [0    ,0    ,0    ,+2/dI*omega,-2/dT,0    ,0],//T1
      [0    ,0    ,0    ,0          ,+2/dT,-2/dT,0],//T2
      [0    ,0    ,0    ,0          ,0    ,+2/dT,0]//D
    ]);
    real S = population;
    real last_D;
    for(i in 1:no_days){
      last_D = state[7];
      S -= daily_infections[i];
      state[1] += daily_infections[i];
      state = transition_matrix * state;
      daily_deaths[i] = state[7] - last_D;
      beta[i] = daily_infections[i] * population / (S*(state[2]+state[3]));
      effective_reproduction_number[i] = daily_infections[i] / (state[3] + state[4]) * dI;
    }
  }
}
model {
  //One possible regularization
  if(beta_regularization){
    unit_dS[2:no_days] ~ lognormal(log(unit_dS[:no_days-1]), beta_regularization);
  }
  //This imposes a very wide prior on the proportion of still susceptible people!
  unit_dS[no_days+1] ~ uniform(0,1);
  dL ~ normal(4.0, 0.2);
  dI ~ normal(3.06, 0.21);
  dT ~ normal(16, 0.71);
  omega ~ beta(100, 9803);
  reciprocal_phi_deaths ~ exponential(5);
  if(likelihood){
     target += neg_binomial_2_lpmf(new_deaths | daily_deaths, 1/reciprocal_phi_deaths);
  }
}

generated quantities {
  vector[no_days-1] growth_rate = (log(daily_infections[2:]) - log(daily_infections[:no_days-1]))*100;
  int pred_deaths[no_days] = neg_binomial_2_rng(daily_deaths, 1 / reciprocal_phi_deaths);
}