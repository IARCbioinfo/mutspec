#!/usr/bin/Rscript

#---------------------------------------#
# Author: Maude                         #
# Script: chi2test_MutSpecStat_Galaxy.r #
# Last update: 18/10/16                 #
#---------------------------------------#


#########################################################################################################################################
#                                  					 					Calculate the chi2 test for the strand bias            		          		  				#
#########################################################################################################################################

#-------------------------------------------------------------------------------
# Load library for recovering the arguments
#-------------------------------------------------------------------------------
suppressMessages(suppressWarnings(require("getopt")))



#-------------------------------------------------------------------------------
# Recover the arguments
#-------------------------------------------------------------------------------
spec = matrix(c(
                "folderChi2", "folderChi2", 1, "character",
                "help",       "h",          0, "logical"
               ),
               byrow=TRUE, ncol=4
            )

opt = getopt(spec)

# No argument is pass to the command line
if(length(opt) == 1)
{
  cat(paste("Usage:\n chi2test_MutSpecStat_Galaxy.r --folderChi2 <path_to_folder> \n",sep=""))
  q(status=1)
}

# Help was asked for.
if ( !is.null(opt$help) )
{
  # print a friendly message and exit with a non-zero error code
  cat(paste("Usage:\n  chi2test_MutSpecStat_Galaxy.r --folderChi2 <path_to_folder> \n",sep=""))
  q(status=1)
}





## Load the data. There is one column with the mutation type and the sample name but it's just for knowing what is corresponding to each line. The two columns with the number of variant per strand would be sufficient.
inputChi2 <- paste0(opt$folderChi2, "/Input_chi2_strandBias.txt")
strBias<-read.delim(inputChi2, dec=".")

# Chi2
pValChi2       <- c() # First I create an empty vector and then I apply a for on the data load
pValChi2_round <- c() # Empty vector with the rounded p-values
confInt        <- c() # Empty vector for the confident interval
proportion     <- c() # Empty vector for the proportion of NonTr compared to the (NonTr+Tr)
sampleSize     <- c() # Empty vector for the count of samples in NonTr and Tr
# For Pool_Data save the p-values in a different vector for not having them for the FDR
pValChi2_PoolData       <- c()
pValChi2_PoolData_Round <- c()

j = 1 # Timer for pValChi2_PoolData vector
k = 1 # Timer for pValChi2

for(i in 1:nrow(strBias))
{
	if(! sum(strBias[i,2:3]) == 0)
	{
		# For Pool_Data
		if(strBias[i,1] == "Pool_Data")
		{
			pValChi2_PoolData[j] <- prop.test(x=strBias[i,2],n=sum(strBias[i,2:3]),p=0.5)$p.value
			j <- j+1
		}
		# For the other sample(s)
		else
		{
			# Calculate the p-value
			pValChi2[k] <- prop.test(x=strBias[i,2],n=sum(strBias[i,2:3]),p=0.5)$p.value
			k <- k+1
		}

		# Calculate the confidence interval
		temp       <- prop.test(x=strBias[i,2],n=sum(strBias[i,2:3]),p=0.5)$conf.int
		confInt[i] <- paste0("[", round(temp[1],2), "-", round(temp[2],2), "]") # Same as paste(sep="")

		# Save the proportion
		proportion[i] <- strBias[i,2] / sum(strBias[i,2:3])

		# Save the sample size (count on NonTr and Tr)
		sampleSize[i]  <- paste(strBias[i,2], strBias[i,3], sep="-")
	} else
	{
		if(strBias[i,1] == "Pool_Data")
		{
			pValChi2_PoolData[j]       <- NA
			pValChi2_PoolData_Round[j] <- NA
			j <- j+1
		}
		else
		{
			# Not enough effective for the test
			pValChi2[k]       <- NA
			pValChi2_round[k] <- NA
			k <- k+1
		}

		confInt[i]        <- NA
		proportion[i]     <- NA
		sampleSize[i]     <- NA
	}
}

# Adjust with FDR
FDR<-p.adjust(pValChi2, method="BH")

# Rount the p-value
for(i in 1:nrow(strBias))
{
	pValChi2_round[i] <- format(pValChi2[i], scientific=T, digits=3)
}

# The option for the pool is specified
if(!is.null(pValChi2_PoolData))
{
	# Round the p-value for Pool_Data
	for(i in 1:6)
	{
		pValChi2_PoolData_Round[i] <- format(pValChi2_PoolData[i], scientific=T, digits=3)
	}
}


# I create a dataframe for add what I want
outputChi2 <- data.frame(round(strBias[,2]/strBias[,3], digits=2), sampleSize, round(proportion, 3), confInt)
outputChi2$Mut.type   <- strBias$Alteration
outputChi2$SampleName <- strBias$SampleName
colnames(outputChi2)[1:6]<-c("Strand_Bias", "NonTr-Tr", "Proportion", "Confidence Interval", "Mutation_Type", "SampleName")

# Transform the data frame into a matrix for adding the p-value for the samples and Pool_Data
matrix         <- as.matrix(outputChi2)
tempColPValFDR <- matrix(, nrow=length(sampleSize), ncol = 2) # Create an empty matrix with 2 columns for adding the p-value and the FDR
matrix         <- cbind(matrix, tempColPValFDR)
j = 1 # Timer for all the sample
k = 1 # Timer for Pool_Data
for(i in 1:nrow(matrix))
{
	if(matrix[i,6] == "Pool_Data")
	{
	  matrix[i,7] <- pValChi2_PoolData_Round[k]
	  matrix[i,8] <- "NA" # No FDR for Pool_Data
	  k = k+1
	}
	else
	{
	  matrix[i,7] <- pValChi2_round[j]
	  matrix[i,8] <- round(FDR[j], 3)
	  j = j+1
	}
}

# Reorder the columns
matrix <- cbind(matrix[,1:3], matrix[,7], matrix[,8], matrix[,4:6])
colnames(matrix)[4] <- "P-val-Chi2"
colnames(matrix)[5] <- "FDR"

# Export the file
# dec=".": Set the separator for the decimal by "."
outputFileChi2 <- paste0(opt$folderChi2, "/Output_chi2_strandBias.txt")
write.table(matrix,file=outputFileChi2,quote = FALSE,sep="\t",row.names = FALSE,dec=".")
