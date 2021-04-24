data {
  int N_data_count;
  int M_feature_count;
  row_vector[M_feature_count] x_features[N_data_count];
  int y_categories[N_data_count];
}

parameters {
  row_vector[M_feature_count] b_coefs;
  real intercept; 
}

model {
  
  
  
}