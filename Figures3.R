#########################################
#
#			R code for Application of Information Theory
#
##########################################

#	Set the working directory from where files will read and to which files will be written

			setwd("/Users/jk/Desktop/KPAGL")

			spdat <- read.table("spikes_.dat")

#	 KPAL modelling the spikes data-- 4 by 4 by 3 system
 
			a <- seq(0, 200, length=21)
			b <- seq(0, 300, length=31)


			
			spdat1 <- stack( spdat)

			apdat <- cbind(  expand.grid(b,a), spdat1[,1])

#	Categorise the data

			ap <-  cut(apdat[,3], breaks=c(-0.5, 1.5, 2.5, 4.5), labels=c("0", "1", "2"))
			apical <-  cut(apdat[,2], breaks=c(-0.5, 50.5, 100.5, 150.5, 200.5), labels=c("0", "1", "2", "3"))
			basal <-  cut(apdat[,1], breaks=c(-0.5, 65, 145, 225, 305), labels=c("0", "1", "2", "3"))


			codes <- data.frame(ap, apical, basal)
			str(codes)

#	Tabulate the categorised data

			tab <- table(codes)
			
#	Print out the tabulated frequencies. They are in the Python script, cspid.py

			out <- data.frame(xtabs(~ ., codes))$Freq
			
			ind1 <- rep(c(0,1,2,3), each =12)
			ind2 <- rep(c(0,1,2,3), each =3, times=4)
			ind3 <- rep(c(0,1,2), times =16)
			
			cbind(ind1, ind2, ind3, out)
			
			


			
