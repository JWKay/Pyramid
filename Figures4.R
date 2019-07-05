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

#	Categorise the binary data

			ap <-  cut(apdat[,3], breaks=c(-0.5, 1.5, 2.5, 4.5), labels=c("0", "1", "2"))
			
# Define arrays to store the apical break points
   
             Abr1 <- c(-0.5, 55, 105, 155, 200.5)
             Abr2 <- c(-0.5, 45, 105, 155, 200.5)
            Abr3 <- c(-0.5, 45, 95, 155, 200.5)
            Abr4 <- c(-0.5, 45, 95, 145, 200.5)
            
            Abrks <- rbind(Abr1, Abr2, Abr3, Abr4)
            
            
# Define arrays to store the basal break points
            
            Bbr1 <- c(-0.5, 65, 145, 225, 300.5)
            Bbr2 <- c(-0.5, 75, 145, 225, 300.5)
            Bbr3 <- c(-0.5, 75, 155, 225, 300.5)
            Bbr4 <- c(-0.5, 75, 155, 235, 300.5)
            
            Bbrks <- rbind(Bbr1, Bbr2, Bbr3, Bbr4)
            
# Initialise an array for the output frequencies

            out <- matrix(0, nrow=16, ncol=48)
            
            for (ib in 1:4)
            {
            for(ia in 1:4){
            
        
			basal <-  cut(apdat[,1], breaks= Bbrks[ib,], labels=c("0", "1", "2", "3"))
            apical <-  cut(apdat[,2], breaks= Abrks[ia,], labels=c("0", "1", "2", "3"))

			codes <- data.frame(ap, apical, basal)
			

#	Tabulate the categorised data

			tab <- table(codes)
			
#	Save the tabulated frequencies. 

             ind = 4*(ib-1) + ia

			out[ind, ] <- data.frame(xtabs(~ ., codes))$Freq
            
            }
            }
            
# Write the frequencies to file in the current directory, for processing in Python.
            
          write.table(out, "freqs.txt", row.names =FALSE, col.names=FALSE)  
          
# Now read the results from Python

          	setwd("/Users/jk/Desktop/KPAGL")            # if required

          res = read.table("PIDout.txt")
        
# Examine summary statistics

         apply(res,2, summary)
       
       MEANS =  apply(res, 2, mean)
       
        SDS =  apply(res, 2, sd)
        
 # Produce the summary statistics for Table 4 and Figure 6.
 
 # Classical information measures
 
 MEANS[21:26]
 SDS[21:26]
 

# Partial information decompositions


# Ibroja

MEANS[1:4][c(4, 2, 3, 1)]
SDS[1:4][c(4, 2, 3, 1)]

# Imin

MEANS[5:8][c(4, 2, 3, 1)]
SDS[5:8][c(4, 2, 3, 1)]

# Iproj

MEANS[9:12][c(4, 2, 3, 1)]
SDS[9:12][c(4, 2, 3, 1)]

# Iccs

MEANS[13:16][c(4, 2, 3, 1)]
SDS[13:16][c(4, 2, 3, 1)]

# Idep

MEANS[17:20][c(4, 2, 3, 1)]
SDS[17:20][c(4, 2, 3, 1)]



       
       

