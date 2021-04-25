data {
  int N_data_count;
  int M_feature_count;
  //row_vector[M_feature_count] x_features[N_data_count];
  matrix [N_data_count,M_feature_count] x_features;
  int y_categories[N_data_count];
}

parameters {
  vector[M_feature_count] b_coefs;
  real intercept; 
}

model {
  intercept ~ normal(0,2);
  b_coefs ~ normal(0,2);
  y_categories ~ bernoulli_logit(x_features * b_coefs + intercept);
}