#!/usr/bin/Rscript

#-----------------------------------#
# Author: Maude                     #
# Script: compareSignature_Galaxy.r #
# Last update: 29/10/15             #
#-----------------------------------#


#########################################################################################################################################
#                                   Compare new sigantures with published one using the cosine similarity method                        #
#########################################################################################################################################


#-------------------------------------------------------------------------------
# Print a usage message if there is no argument pass to the command line
#-------------------------------------------------------------------------------
args <- commandArgs(TRUE)
usage <- function() 
{
  msg <- paste0('Usage:\n',
               ' compareSignature_Galaxy.r Published_Signature New_Signature Output_Folder\n'
              )
  cat(msg, '\n', file="/dev/stderr")
  quit(status=1)
}

input = args[length(args)]

if (length(args) == 0) { usage() }


#-------------------------------------------------------------------------------
# Load library
#-------------------------------------------------------------------------------
suppressMessages(suppressWarnings(library(lsa)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(reshape)))

#-------------------------------------------------------------------------------
# Recover the arguments
#-------------------------------------------------------------------------------
published_signature_file <- args[1] # The matrix with the published signatures
unknown_signature_file   <- args[2] # The matrix W from NMF from which we want to compare the signatures
dir								       <- args[3] # html directory


#-------------------------------------------------------------------------------
# Set the variables
#-------------------------------------------------------------------------------
# Create the outputs
output_cosineRes <- paste0(dir, "/Similarity_Matrix.txt")
output_png       <- paste0(dir, "/Similarity_Matrix.png")


#-------------------------------------------------------------------------------
# Calculate the cosine similarity and represent it with a heatmap
#-------------------------------------------------------------------------------
# Published signatures
dataFrame1 <- read.table(published_signature_file, header=T, sep="\t")
# Remove the first three colmumns (Substitution Type, Trinucleotide  Somatic, Mutation Type)
dataFrame1 <- dataFrame1[,4:ncol(dataFrame1)]
matrix1    <- as.matrix(dataFrame1)
  
# Unkown signatures
dataFrame2 <- read.table(unknown_signature_file, header=T, sep="\t")
# Remove the first two columns (alteration, context)
dataFrame2 <- dataFrame2[,3:ncol(dataFrame2)]
matrix2    <- as.matrix(dataFrame2)
# Recover the number of new signatures
NbNewSignature   <- ncol(dataFrame2) - 1

# Combined the two matrices (published and unknown signatures)
input_matrix_cos <- cbind(matrix1, matrix2)
# Calculate the cosine similarity
cosine_res <- cosine(input_matrix_cos)

# Keep only the comparison between the two matrices
nbSign            <- ncol(matrix1)+1 # +1 for havng the first signature of the matrix1
cosine_res_subset <- cosine_res[nbSign:nrow(cosine_res), 1:ncol(matrix1)]
  
# Save the matrix
write.table(cosine_res_subset, file=output_cosineRes, quote=F, sep="\t", col.names=T, row.names=T)

# Transform the matrix in a suitable format for ggplot2
cosineRes_subset_melt <- melt(cosine_res_subset)
# Rename the columns
colnames(cosineRes_subset_melt) <- c("Unknown_Signatures", "Published_Signatures", "Similarity")
# Reorder the Signature for having the same order as in the matrix. Turn your 'signature' column into a character vector
cosineRes_subset_melt$Published_Signatures <- as.character(cosineRes_subset_melt$Published_Signatures)
#Then turn it back into an ordered factor
cosineRes_subset_melt$Published_Signatures <- factor(cosineRes_subset_melt$Published_Signatures, levels=rev(unique(cosineRes_subset_melt$Published_Signature)))
  
# Base plot: heatmap
p1 <- ggplot(cosineRes_subset_melt, aes(x=Published_Signatures, y=Unknown_Signatures, fill=Similarity)) + geom_tile(colour="yellow") +scale_fill_gradientn(colours=c("yellow", "red")) + theme_classic()
  
# Rename the signatures
if(basename(published_signature_file) == "Frequency-COSMICv72-Hupki.txt")
{
  p1 <- p1 + scale_x_discrete(breaks = c("Signature.1", "Signature.2", "Signature.3", "Signature.4", "Signature.5", "Signature.6", "Signature.7", "Signature.8", "Signature.9",
                                         "Signature.10", "Signature.11", "Signature.12", "Signature.13", "Signature.14", "Signature.15", "Signature.16", "Signature.17",
                                         "Signature.18", "Signature.19", "Signature.20", "Signature.21", "Signature.22", "Signature.23", "Signature.24", "Signature.25",
                                         "Signature.26", "Signature.27", "Signature.28", "Signature.29", "Signature.30",
                                         "Signature.1.MEF", "Signature.2.MEF", "Signature.3.MEF", "Signature.5.MEF"),
                              labels = c("(Age) Sig 1", "(AID/APOBEC) Sig 2", "(BRCA1/2) Sig 3", "(Smoking) Sig 4", "Sig 5", "(DNA MMR deficiency) Sig 6", "(UV) Sig 7",
                                         "Sig 8", "(IgG) Sig 9", "(pol e) Sig 10", "(temozolomide) Sig 11", "Sig 12", "(AID/APOBEC) Sig 13", "Sig 14",
                                         "(DNA MMR deficiency) Sig 15", "Sig 16", "Sig 17", "Sig 18", "Sig 19", "(DNA MMR deficiency) Sig 20", "Sig 21", "(AA) Sig 22",
                                         "Sig 23", "(Aflatoxin) Sig 24", "Sig 25", "(DNA MMR deficiency) Sig 26", "Sig 27", "Sig 28", "(Tobacco chewing) Sig 29", "Sig 30",
                                         "(AA) Sig 1 MEF", "(AID) Sig 2 MEF", "(BaP) Sig 3 MEF", "(MNNG) Sig 5 MEF")
                              )
}

# Flipped cartesian coordinates so that horizontal becomes vertical, and vertical, horizontal
p1 <- p1 +  coord_flip()
# Remove the x axis line
p1 <- p1 + theme(axis.line.x=element_blank(), axis.line.y=element_blank())
# Add the cosine value only if >= 0.9
cosResLabel <- subset(cosineRes_subset_melt, round(cosineRes_subset_melt$Similarity, digits=2) >= 0.9) # Subset the data for keeping only the values greater thant 0.9
p1 <- p1 + geom_text(data = cosResLabel, aes(x = Published_Signatures, y = Unknown_Signatures, label = round(cosResLabel$Similarity, 2)))

graphics.off()
options(bitmapType='cairo')
png(output_png, width=3000, height=2000, res=300)
plot(p1)
invisible( dev.off() )
