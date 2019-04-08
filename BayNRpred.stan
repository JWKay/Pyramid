//  Bayesian nonlinear regression for threshold model

data {                                  // comments in BayNR.stan are relevant here
  int<lower=0> N; 
  real apic[N]; 
  real thresh[N]; 
  real wts[N] ;
  real<lower=0> scale;
} 
transformed data{
 real a[N]; 
  real t[N]; 
  real w[N] ;
  for(i in 1:N){
  a[i] = apic[i]/scale;
  t[i] = thresh[i]/scale;
  w[i] =  wts[i]*scale;
  }
}
parameters {
  real b1; 
  real<lower=0>  b2;  
 real b3; 
  real<upper=0> b4;  
  real sigma;  
} 
transformed parameters{
real theta1;
real<lower=0> theta2;
real theta3;
real<upper=0> theta4;
theta1 = scale* b1;
theta2 = scale* b2;
theta3= scale* b3;
theta4 = scale* b4;
}
model {
  real m[N];
  for (i in 1:N) {
    m[i] = b1 + b2/( 1 +exp(-(a[i] -b3)/b4)) ;
  t[i] ~ normal(m[i], sigma/w[i]); 
  }
  b1 ~ normal(0.0, 10); 
  b2 ~ normal(0.0, 10); 
  b3 ~ normal(0.0, 10); 
  b4 ~ normal(0.0, 10); 
  sigma ~ uniform(0, 20); 
}
generated quantities{         // predicted values of the threshold are generated in this block, and the scaling undone
   real thrNew[N];
   for (i in 1:N){
   thrNew[i] = scale*normal_rng(b1 + b2/( 1 +exp(-(a[i] -b3)/b4)), sigma/w[i]);
   }
}
