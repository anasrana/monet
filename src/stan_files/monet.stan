data {
  // dimensions
  int<lower=1> no_gns; // # genes
  int<lower=1> no_tf; // # TFs
  int<lower=1> no_tpt; // # time points
  real<lower=0> de_b; //Laplace parameter
  real<lower=0> sigma;
  real<lower=0> sigma_x;
  // response variable
  vector[no_gns] y0;
  vector[no_gns] y[no_tpt]; // response variable
  vector[no_gns] a[no_tpt]; // masking matrix (e.g. ATACseq)
  real<lower=0> alpha;
}

parameters {
  // regression coefficient vector
  // real<lower=0> sigma;
  // real<lower=0> sigma_x;
  // restrict the model to take the choice of
  vector<lower=0>[no_tf] x[no_tpt];
  matrix[no_gns, no_tf] b;
  simplex[no_tf] w;
}

transformed parameters {
  //
  vector[no_tpt] tmp[no_gns];
  vector[no_tpt] mu[no_gns];
  for (i_tpt in 1:no_tpt) {
    tmp[i_tpt] = b * (w .* x[i_tpt]);
    mu[i_tpt] = y0 + a[i_tpt] .* tmp[i_tpt];
  }

}

model {
  // priors
  for (i_tpt in 1:no_tpt) {
    for (i_tf in 1:no_tf) {
      x[i_tpt, i_tf] ~ normal(0, sigma_x)T[0, ];
    }
  }
  // to_matrix(x) ~ normal(0, sigma_x)T[0, ];
  w ~ dirichlet(rep_vector(alpha, no_tf));

  for (i_tf in 1:no_tf) {
    b[, i_tf] ~ double_exponential(0, de_b);
  }
  // likelihood
  // y ~ multi_normal_cholesky(mu, diag_pre_multiply(sigma))
  for (i_tpt in 1:no_tpt) {
    y[i_tpt] ~ normal(mu[i_tpt], sigma);
  }
}

