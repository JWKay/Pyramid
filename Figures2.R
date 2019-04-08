###################################################
#
#			R code for Figs. 4(B, D, F) &  10-fold Crossvalidation
#
##################################################

#	Set the working directory from where files will read and to which files will be written

			setwd("/Users/jk/Desktop/KPAGL")
	
#	Read the 'spikes' data from Shai et al. (2015) which is in the current working directory

			spdat <- read.table("spikes_.dat")

#	Transpose the data, recode the data into binary format, and stack it.

			spdat1 <- stack(spdat)
			
			spdat1[spdat1 == 1] <- 0
			spdat1[spdat1 == 2] <- 1
			spdat1[spdat1 == 3] <- 1
			spdat1[spdat1 == 4] <- 1
		

# Create lists with the numbers of basal and apical inputs

			n_bas <- 31; bas <- seq(0, 300, length.out= n_bas)
			n_apic <- 21; apic <- seq(0, 200, length.out=n_apic)
			
# Create a data matrix to hold the binary repsonses and corresponding numbers of basal
# and apical inputs.

			GMapdat <- cbind(spdat1[,1], expand.grid(bas,apic))

			colnames(GMapdat) <- c("ap", "basal", "apical")




#	Load the RStan library and set some options. If not previously installed, the RStan 
#	library can in installed using the R command
#
#				install.packages("rstan", dependencies = TRUE)
#

#	Now load the rstan package, and set the options for multiple cores etc.

			library(rstan)
			options(mc.cores = parallel::detectCores())
			rstan_options(auto_write = TRUE)
			
#	Specify the input data to be pased to rstan.

			data.list <- list(N = nrow(GMapdat), ap = GMapdat$ap, basal = GMapdat$basal, 
			apical = GMapdat$apical, scale =100)
			str(data.list)

# Fit the Stan model using file GenMod.stan.  The model code is in the working directory specified above.
#

			apfit <- stan(file = "GenMod.stan", data = data.list,  seed = 920389374,  init=list( list(b2=.5, b3=.5, b4=1.5), 
		 list(b2=1, b3=1, b4=1),  list(b2=1.5, b3=1.5, b4=1.5),   list(b2=2, b3=2, b4=2)), warmup =500, iter = 3000, chains = 4)
									
			summary(apfit)
			

			
			


#	Examine some graphical output: check the traces of chains for convergence and mixing, 
# and examine the density plots of the posterior distributions for the different chains.
# Check R-hat, N_eff and the Monte Carlo standard error (se_mean).



			stan_plot(apfit, show_density=TRUE, pars=c("beta2", "beta3", "beta4"), ci_level=0.95)

			stan_trace(apfit, pars=c("beta2", "beta3", "beta4"), nrow=3, inc_warmup =FALSE)


			stan_dens(apfit, pars=c("beta2", "beta3", "beta4"), separate_chains=TRUE)

			stan_ac(apfit, pars=c("beta2", "beta3", "beta4"), separate_chains=TRUE)




#  Load a function required to compute the posterior predictive probabilities

			GMProbEst<- function(a, b, x){
			lp = -5.2933 + x[1]*b*(1 + exp(x[2]*a*b*(1+exp(x[3]*a)))) 
			1/(1 +exp(-lp))
			}

#	Extract the 10000 by 6 matrix of simulations from the posterior distributions of
#	the model parameters, beta_1 etc. in columns 4 to 6.


			GMdist <- as.matrix(apfit)[,4:6]
			
			npts = nrow(GMapdat)	
			GMest <- vector(mode ="numeric", length = npts)
			for (i in 1:npts){
			a = GMapdat$apical[i]
			b= GMapdat$basal[i]
			GMest[i] = mean(apply( GMdist, 1, GMProbEst, a=a, b=b))
			}
			
#	Predict the class of each binary response (0/1) and
#	compute the misclassification rate for the training data.

			GMclass = GMest >=0.5
			GMmisc_train = 100 * sum(GMclass != GMapdat$ap)/npts

# Create a data matrix of the posterior predictive probabilities, and save it to the 
#	current working directory.

			GMdat = cbind(apdat$basal, apdat$apical, GMest)
			
			write.table(GMdat, file ="GMdat.csv", row.names =FALSE, col.names =FALSE, sep =",")
			


# Now for 10 fold cross validation, assuming that GMProbEst is loaded and the data matrix 
# GMapdat  is also loaded.
	
			
# Define the 10  index sets

    		ind <- sample(1:651)
    		
    		set1<- ind[1:65]
    		set2 = ind[66:130]
    		set3= ind[131:195]
    		set4= ind[196:260]
    		set5=ind[261:325]
    		set6<- ind[326:390]
    		set7 = ind[391:455]
    		set8= ind[456:520]
    		set9= ind[521:585]
    		set10=ind[586:651]
    		
    		
    		allset <- list(set1, set2, set3, set4, set5, set6, set7, set8, set9, set10)
    		
    		misc <-vector(mode ="numeric", length=10)

# Loop through the index sets, fit the model to 9 of the sets and predict the class
# of the observations in the set not used for model-fitting.
    		
    		for(i in 1:10) {
    		
    		data.list <- list(N = nrow(GMapdat[-allset[[i]], ] ),  ap = GMapdat[-allset[[i]],]$ap,  basal= GMapdat[-allset[[i]],]$basal,  
    		apical = GMapdat[-allset[[i]],]$apical, scale =100)
			GMCVfit <- stan(file = "GenMod.stan", data = data.list,  seed = 920389374,  init=list( list(b2=.5, b3=.5, b4=1.5)), warmup=1000, iter=11000, chains=1)
			s=summary(GMCVfit); summ =s$summary; print(summ)
			
			
			pdist = as.matrix(GMCVfit)[, 4:6]
			dat = GMapdat[allset[[i]], ]
			pest <- vector(mode ="numeric", length = nrow(dat))
			for (j in 1:nrow(dat)){
			a = dat$apical[j]
			b= dat$basal[j]
			pest[j] = mean(apply(pdist, 1, GMProbEst, a=a, b=b))
			}

			class =pest >=0.5
			misc[i] =  sum(class != dat$ap)
			
			}
			
#  Misclassification percentage

    		100*sum(misc)/651
    		
    		
    		
    		
    		
    		
    		
    		
    	


