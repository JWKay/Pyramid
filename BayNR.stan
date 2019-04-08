//  Bayesian nonlinear regression for the threshold model

data {
  int<lower=0> N;                 //  the number of combinations of basal and apical inputs
  vector[N] apic;                       //  the numbers of apical tuft inputs
  vector[N] thresh;                    //  the numbers of basal inputs
  vector[N] wts;                        // the weights
  real<lower=0> scale;         //   a scale factor to scale the inputs and the weights
} 
transformed data{                //  the scaling is performed in this block
 vector[N] a; 
  vector[N]  t; 
  vector[N] w;
   a = apic/scale;
  t = thresh/scale;
  w =  wts*scale;
}
parameters {                          // the parameters used to fit the model are defined in this block
  real b1; 
  real<lower=0>  b2;              // constraints are placed on b2 and b4 to ensure the the logistic function is increasing
 real b3; 
  real<upper=0> b4;  
  real sigma;  
} 
transformed parameters{       // the parameters used in the mathematical model are given here
real theta1;                             // scaling the input data results in a corresponding scaling of the b_i parameters
real<lower=0> theta2;            
real theta3;
real<upper=0> theta4;
theta1 = scale* b1;                 // the theta parameters are the ones used with the data on its original scale.
theta2 = scale* b2;
theta3= scale* b3;
theta4 = scale* b4;
}
model {
  real m[N];                                                       // the vector of mean thesholds
  for (i in 1:N) {
    m[i] = b1 + b2/( 1 +exp(-(a[i] -b3)/b4)) ;       // definition of the weighted nonlinear regression model
  t[i] ~ normal(m[i], sigma/w[i]);                        // Note Stan's notation here: N(m, s) denotes normal with mean m and variance s^2.
  }
  b1 ~ normal(0.0, 10);                                    // independent priors for the parameters b1 - b4 and sigma
  b2 ~ normal(0.0, 10);                                    // the normal(0, 10) priors are weakly informative, as is the 
  b3 ~ normal(0.0, 10);                                    // uniform prior on sigma
  b4 ~ normal(0.0, 10); 
  sigma ~ uniform(0, 20); 
}
