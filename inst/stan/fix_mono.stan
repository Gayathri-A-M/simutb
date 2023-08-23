data {
  int<lower=0>       N;
  int<lower=0>       n_dose;
  vector<lower=0>[N] y;
  int<lower=0>       DL[N];
  real<lower=0>      L_m;
  real<lower=0>      U_m;
  real<lower=0>      L_m_inc;
  real<lower=0>      U_m_inc;
  real<lower=0>      L_cv;
  real<lower=0>      U_cv;
}

parameters {
  real<lower=L_m,  upper=U_m> alpha1;
  vector<lower=L_m_inc, upper=U_m_inc>[n_dose-1]  alpha;
  real<lower=L_cv, upper=U_cv> cv;
}

transformed parameters {
  vector<lower=L_m, upper=U_m>[n_dose] m;
  vector[n_dose] mu;
  real<lower=0>  sigma;

  m[1] = alpha1;
  for (j in 2:n_dose){
    m[j] = m[j-1] + alpha[j-1];
  }

  for (j in 1:n_dose){
    mu[j] = log(m[j] / sqrt(1 + cv^2));
  }
  sigma = sqrt(log(1 + cv^2));
}

model {
  cv     ~ uniform(L_cv, U_cv);
  alpha1 ~ uniform(L_m,U_m);
  alpha   ~ uniform(L_m_inc, U_m_inc);

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
