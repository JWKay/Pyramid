#############################################################
#
#			R code for Figs 3(A, B) and Figs 4(A, C, E) & 10-fold Crossvalidation
#
##############################################################

#	Set the working directory from where files will read and to which files will be written

			setwd("/Users/jk/Desktop/KPAGL")
	
#	Read the 'spikes' data from Shai et al. (2015) which is in the current working directory

			spdat <- read.table("spikes_.dat")

#	Transpose the data, recode the data into binary format, and save to the current working directory
#  for plotting in Mathematica

			spdat1 <- t(spdat)
			
			spdat1[spdat1 == 1] <- 0
			spdat1[spdat1 == 2] <- 1
			spdat1[spdat1 == 3] <- 1
			spdat1[spdat1 == 4] <- 1
			
			write.table(spdat1, file = "spbin.csv", row.names =FALSE, col.names =FALSE, sep =",")
			
			
			
			
#	Penalised logistic regressions

#  Load the logistf and the MASS libraries. If not previously installed, the logist library can be installed using
#	the following R command.
#
#			install.packages("logistf", dependencies = TRUE)
#

			library(logistf)
			
#  Set the numbers of basal and apical tufts, and scale them.
			
			n_bas <- 31; bas <- seq(0, 300, length.out= n_bas)
			n_apic <- 21; apic <- seq(0, 200, length.out=n_apic)
			
# 	Initialise a matrix to store the threshold estimates and approximate estimated standard errors
# 	obtained using the delta method.

			thr <- matrix( rep(0, 2*n_apic), ncol =2)
			
# Run the penalised logistic regressions and store the results

			for(i in 1:n_apic){
			
				resp <-   as.numeric(spdat1[i,])
				mod <-   logistf(resp~bas)
				thr_est <-    - mod$coef[1]/mod$coef[2]
				m1 <-   mod$coef[1]; m2 = mod$coef[2]; v = vcov(mod)
				vv <-    (m1/m2)^2*(v[1,1]/m1^2 - 2*v[1,2]/(m1*m2) + v[2,2]/m2^2)
	
				thr[i,] <-  c( thr_est, sqrt(vv) )
				
				}
				
#	Save the threshold data to the current working directory

			thresh <-  cbind(thr[,1], apic)
			write.table(thresh, file = "threshdat.csv", row.names =FALSE, col.names =FALSE, sep =",")
			
			
			
			
			
			
#	Now load the rstan package, and set the options for multiple cores etc.

			library(rstan)
			options(mc.cores = parallel::detectCores())
			rstan_options(auto_write = TRUE)


#	Specify the input data to be pased to rstan.
			
			t_est <- thr[,1]
			wts = 1/thr[,2]
			
#  We feed the value of scale parameter into Stan so that the inputs and the weights will be scaled
#  internally within the Stan code in BayNR.stan. 

			data.list <- list(N = length(apic), thresh = t_est, apic = apic, wts = wts, scale =100)
			str(data.list)
			
# Fit the Stan model in file BayNR.stan, and examine the chain summaries, check Rhat and Neff

			NRfit <- stan(file = "BayNR.stan", data = data.list, seed = 357,  init=list( list(b1=0.5, b2=1, b3=1.5, b4=-1, sigma=1), 
		 list(b1=-1, b2=2, b3=2, b4=-2, sigma=1),  list(b1=1, b2=1.5, b3=3, b4=-0.5, sigma=1),  
		 list(b1=2, b2=0.5, b3=1.5, b4= -1.5, sigma=1)), 
		warmup =500, iter = 3000, chains = 4)
		
			NRout <- summary(NRfit)
			summ <- NRout$summary
			
#	Examine some graphical output: check the traces of chains for convergence and mixing, 
# and examine the density plots of the posterior distributions for the different chains.
# Check R-hat, N_eff and the Monte Carlo standard error (se_mean).

			

			stan_trace(NRfit, pars=c("theta1", "theta2", "theta3", "theta4", "sigma"), nrow=4, inc_warmup =FALSE)

			stan_dens(NRfit, pars=c("theta1", "theta2", "theta3", "theta4", "sigma"), separate_chains=TRUE)
			
#	Look at the autocorrelations, and check correlations using pairs(NRfit)

			stan_ac(NRfit, pars=c("theta1", "theta2", "theta3", "theta4", "sigma"), separate_chains=TRUE)
			
			
			
			
			
			
#	Create an array of size 651 by 3 to hold the binary AP repsonses, together with the numbers 
#	of basal and apical tuft  inputs.

			spdat <- read.table("spikes_.dat")

			spdat2 <- stack(spdat)
			
			spdat2[spdat2 == 1] <- 0
			spdat2[spdat2 == 2] <- 1
			spdat2[spdat2 == 3] <- 1
			spdat2[spdat2 == 4] <- 1
			
			apdat <- cbind(spdat2[,1], expand.grid(bas,apic))

			colnames(apdat) <- c("ap", "basal", "apical")
			

# Load a function required to compute the posterior predictive probabilities

			ProbEst<- function(a, b, x){
			lp = -5.2933 + b - x[1] - x[2]/(1 +exp(-(a -x[3])/x[4]))
			1/(1 +exp(-lp))
			}
			
#	Extract the 10000 by 4 matrix of simulations from the posterior distributions of
#	the model parameters, theta_1 etc. in columns 6 to 9.


			pdist <-  as.matrix(NRfit)[,6:9]
			
#	Compute the posterior predictive probabilities 
          
          
           	npts <- nrow(apdat)
			pest <- vector(mode ="numeric", length = npts)
			
		
			
			for (i in 1:npts){
			a = apdat$apical[i]
			b=apdat$basal[i]
			pest[i] = mean(apply(pdist, 1, ProbEst, a=a, b=b))
			}
			
#	Predict the class of each binary repsonse (0/1) and
#	compute the misclassification rate for the training data.

			class =pest >=0.5
			TMmisc_train = 100 * sum(class != apdat$ap)/npts
			
# Create a data matrix of the posterior predictive probabilities, and save it to the 
#	current working directory.

			TMdat = cbind(apdat$basal, apdat$apical, pest)
			
			write.table(TMdat, file ="TMdat.csv", row.names =FALSE, col.names =FALSE, sep =",")
			





# Now for 10 fold cross validation: The following code assumes ProbEst is currently loaded 
# and that the data matrix apdat is also available.
	
#  Randomize the indices
			
 #   		ind <- sample(1:21)
 
 # The random list of indices used:
    		
    		 ind <- c(3, 18, 21, 17, 15,  5,  7, 12, 10,  2,  6, 11, 16,  1, 19,  9,  8,  4, 13, 20, 14)
    		 
    		set1<- ind[1:2]
    		set2 = ind[3:4]
    		set3= ind[5:6]
    		set4= ind[7:8]
    		set5=ind[9:10]
    		set6<- ind[11:12]
    		set7 = ind[13:14]
    		set8= ind[15:16]
    		set9= ind[17:18]
    		set10=ind[19:21]
    		
    		
    		allset <- list(set1, set2, set3, set4, set5, set6, set7, set8, set9, set10)
    		
    		misc <-vector(mode ="numeric", length=10)
    		
    		
    		for(i in 1:10) {
    		
    		data.list <- list(N = nrow(thrdat[-allset[[i]], ] ), thresh = t_est[-allset[[i]] ], apic = apic[-allset[[i]] ], wts = wts[-allset[[i]] ], scale =100 )
			CVfit <- stan(file = "BayNR.stan", data = data.list, seed= 357,  init=list( list(b1=0.5, b2=1, b3=1.5, b4=-1, sigma=1)), warmup=1000, iter=11000, chains=1)
			s=summary(CVfit); summ =s$summary; print(summ)
			
			if (i < 10){
			dat = rbind(apdat[apdat$apical==10*allset[[i]][1], ], apdat[apdat$apical==10*allset[[i]][2], ])}
			if (i ==10){
			dat = rbind(apdat[apdat$apical==10*allset[[i]][1], ], apdat[apdat$apical==10*allset[[i]][2], ],  apdat[apdat$apical==10*allset[[i]][3], ])}
			pdist = as.matrix(CVfit)[,6:9]
			pest <- vector(mode ="numeric", length =nrow(dat))
			for (j in 1:nrow(dat)){
			a = dat$apical[j]
			b= dat$basal[j]
			pest[j] = mean(apply(pdist, 1, ProbEst, a=a, b=b))
			}

			class =pest >=0.5
			misc[i] =  sum(class != dat$ap)
			
			}
			
#  Misclassification percentage

    		100*sum(misc)/npts
    		
    		
			
			
			
#	Code for pointwise prediction intervals and for plotting them


			
			
			data.list <- list(N = length(apic), thresh = t_est, apic = apic, wts = wts, scale=100)
			str(data.list)

			NRpred <- stan(file = "BayNRpred.stan", data = data.list, seed =357, 
			init=list( list(b1=1, b2=1, b3=1, b4=-1, sigma=1)), warmup =1000, iter = 11000, chains = 1)
			
#	Examine the summary output

			PredOut = summary(NRpred)
			PredSumm = PredOut$summary
			
			PredSumm
			
			
#	The relevant output is in rows 10 to 30	of PredSumm. Columns 4, 6, 8 contain, respectively,
#	the .025th , .5th (median) and 0.975th quantiles of the posterior predictive distributions
			
#  Save the predictions and pointwise 95% prediction intervals for the 21 number of tufts used, 
# for plotting of Fig 3B in Mathematica.

	
			PredInt = cbind(	npts <- nrow(apdat)apic, low, fit, high)
			
            write.table(PredInt, file = "pred.csv", row.names =FALSE, col.names =FALSE, sep =",")
			
				

			

			

		
				
				
			
				
				 
	



			
			
			
			
			
			
			

		

			
