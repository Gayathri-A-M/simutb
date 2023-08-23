data {
  int<lower=0>       N;
  int<lower=0>       n_dose;
  vector<lower=0>[N] y;
  int<lower=0>       DL[N];
  real<lower=0>      L_m;
  real<lower=0>      U_m;
  real<lower=0>      L_cv;
  real<lower=0>      U_cv;
}

parameters {
  vector<lower=L_m, upper = U_m>[n_dose] m;
  real<lower=L_cv,  upper = U_cv> cv;
}

transformed parameters {
  vector[n_dose] mu;
  real<lower=0>  sigma;

  for (j in 1:n_dose){
    mu[j] = log(m[j] / sqrt(1 + cv^2));
  }

  sigma = sqrt(log(1 + cv^2));
}

model {
  cv ~ uniform(L_cv, U_cv);
   m ~ uniform(L_m, U_m);

  for (i in 1:N){
    log(y[i]) ~ normal(mu[DL[i]], sigma);
  }
}

generated quantities {
  vector[n_dose]          x_tilde;
  vector<lower=0>[n_dose] y_tilde;

  for (j in 1:n_dose){
    x_tilde[j] = normal_rng(mu[j], sigma);
    y_tilde[j] = exp(x_tilde[j]);
  }
}
