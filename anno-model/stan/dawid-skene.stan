/**
 * Dawid-Skene model with long-form data input and fixed priors.
 * The true labels are marginalized out in the model block and then
 * reconstituted in the generated quantities block on the log scale.
 * The Stan User's Guide's chapter on Latent Discrete Models has an
 * explanation of the marginalization for full-panel data.
 *
 * Parameters and genrated quantities in output:
 *
 * pi[k]: prevalence of category k
 *
 * theta[j, k, k']: probability that annotator j provides label k' if
 * the true category is k
 *
 * log_p_cat[i, k]:  log probability that true category of item i is k
 */
data {
  int<lower = 0> N;

  int<lower = 0> I;
  int<lower = 1, upper = I> item[N];
  
  int<lower = 0> J;
  int<lower = 1, upper = J> rater[N];

  int<lower = 0> K;
  int<lower = 1, upper = K> label[N];
}
transformed data {
  vector<lower = 0>[K] alpha[K];
  for (k in 1:K) {
    alpha[k] = rep_vector(2, K);
    alpha[k, k] = 2 * K;
  }
}
parameters {
  simplex[K] pi; 
  simplex[K] theta[J, K];
}
model {
  vector[K] log_pi = log(pi);
  vector[K] log_theta[J, K] = log(theta);

  // likelihood
  vector[K] lp[I];
  for (i in 1:I)
    lp[i] = log_pi;
  for (n in 1:N)
    for (k in 1:K)
      lp[item[n], k] += log_theta[rater[n], k, label[n]];
  for (i in 1:I)
    target += log_sum_exp(lp[i]);

  // prior
  for (j in 1:J)
    for (k in 1:K)
      theta[j, k] ~ dirichlet(alpha[k]);
}
generated quantities {
  vector[K] log_p_cat[I];
  {
    vector[K] log_pi = log(pi);
    vector[K] log_theta[J, K] = log(theta);
    vector[K] lp[I];
    for (i in 1:I)
      lp[i] = log_pi;
    for (n in 1:N)
      for (k in 1:K)
        lp[item[n], k] += log_theta[rater[n], k, label[n]];
    for (i in 1:I)
      log_p_cat[i] = lp[i] - log_sum_exp(lp[i]);
  }
}
