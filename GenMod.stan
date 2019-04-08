//  Stan code to fit the general model

data {                          
int<lower=0> N;                                      //   number of observations
int<lower=0,upper=1> ap[N];                 //   binary response variable
vector[N] basal;                                      //   number of basal inputs
vector[N] apical;                                     //   number of apical tuft  inputs
real<lower=0> scale;                             //    a scale factor for scaling the input explanatory variables
}
transformed data{                                //  the scaling of the input data is performed in this block
  vector[N] bas ; 
  vector[N] apic;

  bas = basal/scale;
  apic =  apical/scale;
}
parameters {                                      // the parameters used to fit the model are described in this block
real  b2;                   
real  b3;                
real  b4; 
}
transformed parameters{                    // the parameters used in the mathematical model are given here
real beta2;
real beta3;
real beta4;

beta2 = b2/100;                               // scaling the input data results in a corresponding scaling of the
beta3 = b3/10000;                           // b_i parameters
beta4 = b4/100;
}

model {
b2 ~ uniform(0, 10);                // mutually independent uniform priors on b2, b3, b4
b3 ~ uniform(0, 10);
b4 ~ uniform(0, 10);
ap ~ bernoulli_logit( -5.2933 + b2*bas .* (1 +exp(b3*bas .* apic .* (1+exp(b4*apic)))));  // definition of the nonlinear predictor
}
