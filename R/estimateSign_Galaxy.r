#!/usr/bin/Rscript

#-----------------------------------#
# Author: Maude                     #
# Script: estimateSign_Galaxy.r     #
# Last update: 22/07/15             #
#-----------------------------------#

#########################################################################################################################################
#                                                   Estimate the number of signatures for NMF                                          #
#########################################################################################################################################

#-------------------------------------------------------------------------------
# Load library for recovering the arguments
#-------------------------------------------------------------------------------
suppressMessages(suppressWarnings(require("getopt")))


#-------------------------------------------------------------------------------
# Recover the arguments
#-------------------------------------------------------------------------------
spec = matrix(c(
                "input" ,      "i",      1, "character",
                "stop",        "stop",   1, "numeric",
                "cpu",         "cpu",    1, "integer",
                "output",      "o",      1, "character",
                "help",        "h",      0, "logical"
               ),
               byrow=TRUE, ncol=4
            )

opt = getopt(spec);

# No argument is pass to the command line
if(length(opt) == 1)
{
  cat(paste("Usage:\n estimateSign_Galaxy.r --input <matrix> --stop <maxNbSign> --cpu <cpu> --output <output_filename.png>\n",sep=""))
  q(status=1)
}

# Help was asked for.
if ( !is.null(opt$help) )
{
  # print a friendly message and exit with a non-zero error code
  cat(paste("Usage:\n estimateSign_Galaxy.r --input <matrix> --stop <maxNbSign> --cpu <cpu> --output <output_filename.png>\n",sep=""))
  q(status=1)
}



#-------------------------------------------------------------------------------
# Load library
#-------------------------------------------------------------------------------
suppressMessages(suppressWarnings(library(NMF)))


                                        ###############################################################################
                                        #                                 Load the functions                          #
                                        ###############################################################################

#-------------------------------------------------------------------------------
# Check the file doesn't have lines equal to zero
#-------------------------------------------------------------------------------
CheckFile <- function(rowsum, dataFrame, x)
{
  if(rowsum == 0)
  {
    write("\n\nERROR: There is not enough mutations for running NMF!!!", stderr())
    write(paste0("Input matrix contains at least one null row ", rownames(dataFrame)[x], "\n\n"), stderr())
    stop()
  }
}


                                        ###############################################################################
                                        #                                    Check file                               #
                                        ###############################################################################

# The input musn't contains lines equal to zero !!!
matrixNMF <- read.table(opt$input, header=T)
# suppresses the return of sapply function
invisible( sapply(1:nrow(matrixNMF), function(x) { CheckFile(rowSums(matrixNMF)[x], matrixNMF, x) } ) )



                                        ###############################################################################
                                        #                       Estimate the number of signatures                     #
                                        ###############################################################################
# Estimate the number of signatures with our data
nbCPU   <- paste0("vP", opt$cpu)
nbSign  <- 2:opt$stop # The minum number of signatures can't be lower than 2

estim_r <- nmf(matrixNMF, method="brunet", nbSign, nrun=50, .opt=nbCPU)

# Shuffle original data
v_random <- randomize(matrixNMF)
# Estimate quality measures from the shuffled data
estim_r_random <- nmf(v_random, method="brunet", nbSign, nrun=50, .opt=nbCPU)

# Plot the estimation for our data and the random ones
graphics.off()
options(bitmapType='cairo')
png(opt$output, width=3000, height=2000, res=300)
plot(estim_r, estim_r_random)
invisible( dev.off() )
