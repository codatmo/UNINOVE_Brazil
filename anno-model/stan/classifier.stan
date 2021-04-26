data {
  int N_data_count;
  int M_feature_count;
  int H_held_out;
  // row_vector[N_data_count] x_features[M_feature_count];
  matrix [N_data_count,M_feature_count] x_features;
  matrix [H_held_out,M_feature_count] x_held_out;
  int y_categories[N_data_count];
  int y_held_out[H_held_out];
}

parameters {
  vector[M_feature_count] b_coefs;
  real intercept; 
}

model {
  intercept ~ normal(0,2);
  b_coefs ~ normal(0,2);
  // y_categories ~ bernoulli_logit(x_features * b_coefs + intercept);
  target += bernoulli_logit_glm_lpmf(y_categories | x_features, intercept, b_coefs);
}

generated quantities {
  real precision;
  real recall;
  real tp = 0;
  real fp = 0;
  real tn = 0;
  real fn = 0;

  for (i in 1:H_held_out) {
    real prob = inv_logit(x_held_out[i] * b_coefs + intercept);
    print("prob=", prob);
    if (prob > .5) { 
      if (y_held_out[i] == 1) {
        tp += 1;
      }
      else {
        fp += 1;
      }
    }
    else {
      if (y_held_out[i] == 1) {
        fn += 1;
      }
      else {
        tn += 1;
      }
    }
  }
  precision = tp/(tp + fp);
  recall = tp/(tp + fn);
  print("precision=", precision);
  print("recall=", recall);
}
