# R code to simulate basal and apical data and compute the probabilities  of a 2nd AP
# given (a) basal and (b) basal and apical input.

# The MASS library is required for simulating from a bivariate Gaussian distribution.

library(MASS)

# Set the seed to make the results reporducible.

set.seed(1357531)

# The following commands are those required to generate data from a bivariate
# Gaussian distribution for the joint distributions of basal and apical inputs,
# separately for the two groups, Z_2 =0 and Z_2 =1.
# Since th emarginal and conditional distributions of a bivariate Gaussian are
# univariate Gaussian, we simulte the basal data first, and then the 
# apical data given the basal data.

# ngroup is the number of observations assumed in each group.

ngroup =500; ma1 = 15; ma0 =10;  sda1 =3; sda0 =2
mb1 = 12;  mb0 = 7; sdb1=2; sdb0 =2
corr1 = -0.7; corr0 = -0.7

# Simulate the basal data from two different univariate Gaussian distributions
# with means mb1 and mb0, respectively, and standard deviations, sdb1, sdb0.

b1 <- rnorm(ngroup, mb1, sdb1)
b0 <- rnorm(ngroup,mb0 , sdb0)


# Set up the parameters for the conditional distributions of apical input given basal input.


csd1 = sda1*sqrt(1-corr1^2) ; csd0 = sda0*sqrt(1-corr0^2)

beta1 = corr1*sda1/sdb1
alpha1 = ma1 - beta1*mean(b1)

beta0 = corr0*sda0/sdb0
alpha0 = ma0 - beta0*mean(b0)

# Simulate from the conditional distirbutions of apical input given basal input.

a1 <- rnorm(ngroup, alpha1+beta1*b1, csd1)
a0 <- rnorm(ngroup,  alpha0 +beta0*b0, csd0)


edge=2; xmax = max(a0, b0)  +edge; ymax = max(a1, b1)+edge
xmin=min(a0, b0) ; ymin=min(a1, b1)-edge

# Plot the data: 0s in blue, 1s in red

plot(b0, a0, col="blue", pch=16,xlim=c(0,xmax), ylim=c(0,ymax), pty="s")

points(b1, a1, col="red", pch=16)


# Sample from posteriors: Normal unknown variance

# First, b0. Simulate the unknown mean and variance parameters for group 0.
# sig2b0 is sigma^2_0 and meanb0 is mu_0.

bm0 = mean(b0); bv0 =var(b0); n0 = ngroup

nsim =10000

sig2b0 = (n0-1)*bv0 / rchisq(nsim, n0 -1)
meanb0 = rnorm(nsim, bm0, sqrt(sig2b0/n0))

# Second, b1. Simulate the unknown mean and variance parameters for group 1.
# sig2b1 is sigma^2_1 and meanb0 is mu_1.

bm1 = mean(b1); bv1 =var(b1); n1 = ngroup

nsim =10000

sig2b1 = (n1-1)*bv1 / rchisq(nsim, n1 -1)
meanb1 = rnorm(nsim, bm1, sqrt(sig2b1/n1))


# Third, simulate the parameters in the regression model for A_0 given B_0 (and Z_2 =0)
# sig2ab0 is sigma^2_3, and beta0 is psi_0 = (alpha_0, beta_0).


lm0 = lm(a0 ~ b0, x = TRUE)

coef0 = coef(lm0)
X0 = lm0$x
V0 = solve(t(X0)%*%X0)
df0 = ngroup-2
rss0 = sum(lm0$residuals^2)/df0

sig2ab0 = df0*rss0/rchisq(nsim, df0)
psi0 = matrix(0, nrow= nsim, ncol=2)
for (i in 1:nsim){SIG =  sig2ab0[i]*V0; 
psi0[i,] = mvrnorm(n = 1,mu = coef0, Sigma =SIG )}

# Fourth,  simulate the parameters in the regression model for A_0 given B_0 (and Z_2 =0)
# sig2ab1 is sigma^2_4, and beta1 is psi_1 = (alpha_1, beta_1).

lm1 = lm(a1 ~ b1, x = TRUE)

coef1 = coef(lm1)
X1 = lm1$x
V1 = solve(t(X1)%*%X1)
df1 = ngroup-2
rss1 = sum(lm1$residuals^2)/df1

sig2ab1 = df1*rss1/rchisq(nsim, df1)
psi1 = matrix(0, nrow= nsim, ncol=2)
for (i in 1:nsim){SIG =  sig2ab1[i]*V1; 
psi1[i,] = mvrnorm(n = 1,mu = coef1, Sigma =SIG )}

# Collect samples from posterior distributions

poster = cbind(meanb0, sig2b0, meanb1, sig2b1, psi0, sig2ab0, psi1, sig2ab1)



# Prediction of posterior  probability given the basal input
# newB  is a scalar here....the new value of b

predB = function(newB, x){
priorLO = -5.2933
temp = priorLO + 0.5*log(x[2]/x[4]) - (newB - x[3])^2/(2*x[4]) +  (newB - x[1])^2/(2*x[2]) 
1/(1 + exp(-temp))
 }
 
 # Make a prediction for a single new basal input
 
 



# Prediction of posterior  probability given the basal input and apical input
# newBA is a vector here....the new values of (b, a).

predBA = function(newBA, x){
priorLO = -5.2933
temp1 = priorLO + 0.5*log(x[2]/x[4]) - (newBA[1] - x[3])^2/(2*x[4]) +  (newBA[1] - x[1])^2/(2*x[2]) 
temp2 = 0.5*log(x[7]/x[10]) - (newBA[2] - x[8] - x[9]*newBA[1])^2/(2*x[10]) +  (newBA[2] - x[5] -x[6]*newBA[1])^2/(2*x[7]) 

1/(1 + exp(-(temp1 + temp2)))
 }



# Create output posterior probabilities.... NB.... basal first, apical second


# Define a grid of values of the basal and apical inputs

apic = seq(8, 17, .5)
bas = seq(8, 17, .5)
BAvals = expand.grid(bas, apic)
names(BAvals) = c("bas", "apic")

npts = nrow(BAvals)

# Initialise matrices to store the results.

PPest = matrix(0, nrow= npts, ncol=2)
PPstat = matrix(0, nrow=npts, ncol=4)

# For each new pair of basal and apical input, we apply the functions predB and predBA
# to compute Monte Carlo approximations for the posterior probabilities given (a) basal input
# and (b) both basal and apical input.

# PPest stores the two estimated posterior probabilities.
# PPstat stores the Monte Carlo standard error and three quantiles of the posterior
# distribution of P(S_2|b, a, theta).

for (i in 1:npts){
newvals =as.matrix(BAvals[i,])
ppvalsB = apply(poster[,1:4], 1, predB, newB = newvals[1])
ppvalsBA = apply(poster, 1, predBA, newBA = newvals)
PPest[i, 1] = mean(ppvalsB)
PPest[i, 2] = mean(ppvalsBA)
PPstat[i,] = c(sd(ppvalsBA)/sqrt(nsim), quantile(ppvalsBA, probs=c(0.05, 0.5, 0.95)))
}

# Set the working directory, in which the files will be placed.

setwd("/Users/jk/Desktop/KPAGL")

# Write to file the basal and apical inputs and the corresponding estimates of the 
# posterior probability of S_2 given both basal and apical inputs.

# This file will then be read into Mathematica to produce contour and surface plots
# for the supplementary file, S2 File.

PPdat = cbind(BAvals, PPest[,2])
write.table(PPdat, file ="PPdat.csv", row.names =FALSE, col.names =FALSE, sep =",")


PPdat2 = cbind(BAvals, PPest)
write.table(PPdat, file ="PPdat2.csv", row.names =FALSE, col.names =FALSE, sep =",")

# Write to file all the computed quantities for each pair of basal and apical inputs.

PPstats =cbind(PPdat2, PPstat)
write.table(PPstats, file ="PPstats.csv", row.names =FALSE, col.names =FALSE, sep =",")

# Make a plot of data for first figure in the supplementary, S2 File.		

plot(b0, a0, col="blue", pch=16, xlim=c(0,xmax), ylim=c(0,ymax), xlab="Basal input", ylab ="Apical input")
points(b1, a1, col="red", pch=16)
legend("bottomleft", legend =c("0", "1"), col=c("blue", "red"), pch =c(16, 16), inset=c(0.05, 0.05))

