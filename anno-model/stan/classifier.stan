data {
  int N_data_count;
  int M_feature_count;
  //row_vector[M_feature_count] x_features[N_data_count];
  matrix x_features[N_data_count,M_feature_count];
  int y_categories[N_data_count];
}

parameters {
  row_vector[M_feature_count] b_coefs;
  real intercept; 
}

model {
  for (n in 1:N) {
    y_categories[n] ~ binomial_logit(intercept + b_coefs[i,] * x_features[i])
  }
}