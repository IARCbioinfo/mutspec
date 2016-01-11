#!/usr/bin/Rscript

#-----------------------------------#
# Author: Maude                     #
# Script: somaticSignature_Galaxy.r #
# Last update: 29/07/15             #
#-----------------------------------#


#########################################################################################################################################
#                  Run NMF algorithm and represent the composition of somatic signatures and the contribution in each samples           #
#########################################################################################################################################

#-------------------------------------------------------------------------------
# Load library for recovering the arguments
#-------------------------------------------------------------------------------
suppressMessages(suppressWarnings(require("getopt")))


#-------------------------------------------------------------------------------
# Recover the arguments
#-------------------------------------------------------------------------------
spec = matrix(c(
                "input" ,       "i",     1, "character",
                "nbSignature", "nbSign", 1, "integer",
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
  cat(paste("Usage:\n somaticSignature_Galaxy.r --input <matrix> --nbSignature <nbSign> --cpu <cpu> --output <outputdir>\n",sep=""))
  q(status=1)
}

# Help was asked for.
if ( !is.null(opt$help) )
{
  # print a friendly message and exit with a non-zero error code
  cat(paste("Usage:\n somaticSignature_Galaxy.r --input <matrix> --nbSignature <nbSign> --cpu <cpu> --output <outputdir>\n",sep=""))
  q(status=1)
}



#-------------------------------------------------------------------------------
# Load library
#-------------------------------------------------------------------------------
suppressMessages(suppressWarnings(library(NMF)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(reshape)))
suppressMessages(suppressWarnings(library(grid)))
suppressMessages(suppressWarnings(library(scales)))    # Set the maximum value to the y axis (graph composition somatic signature)
suppressMessages(suppressWarnings(library(gridExtra))) # function "unit"



                                        ###############################################################################
                                        #                                 Load the functions                          #
                                        ###############################################################################

#-------------------------------------------------------------------------------
# Set the font depending on X11 availability
#-------------------------------------------------------------------------------
font <- ""
# Check the device available
device <- capabilities()
# X11 is available
if(device[5]) { font <- "Helvetica" } else { font <- "Helvetica-Narrow" }

#-------------------------------------------------------------------------------
# My own theme
#-------------------------------------------------------------------------------
theme_custom <- function(base_size = 4, base_family = "")
{
  # Starts with theme_grey and then modify some parts
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.text        = element_text(size = rel(0.8), family=font),
      axis.ticks       = element_line(colour = "black", size=.2),
      axis.line        = element_line(colour = "black", size = .2),
      axis.ticks.length= unit(.05, "cm"),
      axis.ticks.margin= unit(.05, "cm"), # space between tick mark and tick label (‘unit’)
      legend.key.size  = unit(.2, "cm"),
      legend.margin    = unit(-.5, "cm"),
      panel.background = element_blank(),
      panel.border     = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text.y     = element_text(size = 3)
    )
}

#-------------------------------------------------------------------------------
# Customize the theme for adding a y axis
#-------------------------------------------------------------------------------
mytheme                    <- theme_custom()
mytheme$axis.line.x        <- mytheme$axis.line.y <- mytheme$axis.line
mytheme$axis.line.x$colour <- 'white'

#-------------------------------------------------------------------------------
# Replace the signature number by alphabet letter
#-------------------------------------------------------------------------------
ConvertNb2Aphabet <- function(c)
{
  if(c == "row1" || c == "col1") { c <- "A" } else
  if(c == "row2" || c == "col2") { c <- "B"}  else
  if(c == "row3" || c == "col3") { c <- "C"}  else
  if(c == "row4" || c == "col4") { c <- "D"}  else
  if(c == "row5" || c == "col5") { c <- "E"}  else
  if(c == "row6" || c == "col6") { c <- "F"}  else
  if(c == "row7" || c == "col7") { c <- "G"}  else
  if(c == "row8" || c == "col8") { c <- "H"}  else
  if(c == "row9" || c == "col9") { c <- "I"}  else
  if(c == "row10" || c == "col10") { c <- "J"}  else
  if(c == "row11" || c == "col11") { c <- "K"}  else
  if(c == "row12" || c == "col12") { c <- "L"}  else
  if(c == "row13" || c == "col13") { c <- "M"}  else
  if(c == "row14" || c == "col14") { c <- "N"}  else
  if(c == "row15" || c == "col15") { c <- "O"}  else
  if(c == "row16" || c == "col16") { c <- "P"}  else
  if(c == "row17" || c == "col17") { c <- "Q"}  else
  if(c == "row18" || c == "col18") { c <- "R"}  else
  if(c == "row19" || c == "col19") { c <- "S"}  else
  if(c == "row20" || c == "col20") { c <- "T"}  else
  if(c == "row21" || c == "col21") { c <- "U"}  else
  if(c == "row22" || c == "col22") { c <- "V"}  else
  if(c == "row23" || c == "col23") { c <- "W"}  else
  if(c == "row24" || c == "col24") { c <- "X"}  else
  if(c == "row25" || c == "col25") { c <- "Y"}  else
  if(c == "row26" || c == "col26") { c <- "Z"}  else { c <- c }
}

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

#-------------------------------------------------------------------------------
# Contribution to Signature as the number of SBS per sample
#-------------------------------------------------------------------------------
Contri2SignSBS <- function(Total_SBS, Percent)
{
  Total_SBS*Percent/100
}

#-------------------------------------------------------------------------------
# Combined two plots and share the legend
#-------------------------------------------------------------------------------
grid_arrange_shared_legend <- function(...)
{
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}

#-------------------------------------------------------------------------------
# Calculate the mean of each signatures in each cluster
#-------------------------------------------------------------------------------
meanCluster <- function(df)
{
  max <- opt$nbSignature+1
  sapply(2:max, function(x) { round(mean(as.numeric(as.matrix(df[,x]))), 2) } )
}




                                        ###############################################################################
                                        #                                    Check file                               #
                                        ###############################################################################

# The input musn't contains lines equal to zero !!!
matrixNMF <- read.table(opt$input, header=T)
# suppresses the return of sapply function
invisible( sapply(1:nrow(matrixNMF), function(x) { CheckFile(rowSums(matrixNMF)[x], matrixNMF, x) } ) )



                                        ###############################################################################
                                        #                                     Run NMF                                 #
                                        ###############################################################################

# Create the output directories
output_NMF  <- paste0(opt$output, "/NMF")
dir.create(output_NMF)
output_Figures <- paste0(output_NMF, "/", "Figures")
dir.create(output_Figures)
output_Files   <- paste0(output_NMF, "/", "Files")
dir.create(output_Files)

# Define the output filenames
output_cluster         <- paste0(output_Files,   "/", "Cluster_MixtureCoeff.txt")
figure_cluster         <- paste0(output_Figures, "/", "Heatmap_MixtureCoeff.png")
output_matrixW         <- paste0(output_Files,   "/", "MatrixW-Normto100.txt")
output_matrixW_ggplot2 <- paste0(output_Files,   "/", "MatrixW-Inputggplot2.txt")
output_matrixH_ggplot2 <- paste0(output_Files,   "/", "MatrixH-Inputggplot2.txt")
output_matrixH_cluster <- paste0(output_Files,   "/", "Average_ContriByCluster.txt")
figure_matrixW_png     <- paste0(output_Figures, "/", "CompositionSomaticMutation.png")
figure_matrixH_png     <- paste0(output_Figures, "/", "ContributionMutationSignature.png")
figure_matrixH_cluster <- paste0(output_Figures, "/", "Average_ContriByCluster.png")


# Run NMF
# request a certain number of cores to use .opt="vP4"
nbCPU     <- paste0("vP", opt$cpu)
res       <- nmf(matrixNMF, opt$nbSignature, "brunet", nrun=200, .opt=nbCPU)

# If there is more than 300 samples the creation of the heatmap returns an error
if(ncol(matrixNMF) <= 300)
{
  # Save the clustered heatmap generated by NMF
  graphics.off() # close graphics windows
  options(bitmapType='cairo')
  png(figure_cluster)
  coefmap(res, Colv="consensus")
  dev.off()
}

# Recover the matrix W and H
matrixW <- basis(res)
matrixH <- coef(res)

# Recover the cluster of the samples
matrix_cluster           <- cbind(as.numeric(predict(res, what="samples")), colnames(matrixNMF))
colnames(matrix_cluster) <- c("Cluster", "Samples")

## Save the cluster matrix
write.table(matrix_cluster, file=output_cluster, quote=F, sep="\t", col.names=T, row.names=F)



                                        ###############################################################################
                                        #                       Composition of somatic signatures                     #
                                        ###############################################################################

# Normalize to 100%
matrixW_norm <- t((t(matrixW)/colSums(matrixW))*100)
# Add a column name
colnames(matrixW_norm) <- colnames(matrixW_norm, do.NULL = FALSE, prefix = "col")
# Replace the name of the columns by the signature name
colnames(matrixW_norm) <- sapply(1:length(colnames(matrixW_norm)), function(x) { ConvertNb2Aphabet(colnames(matrixW_norm)[x]) } )

# Split the sequence context from the mutation type
context    <- c() # Create an empty vector for the sequence context
alteration <- c() # Create an empty vector for the mutation type
for(i in 1:nrow(matrixW_norm))
{
  temp <- strsplit((strsplit(rownames(matrixW_norm)[i], ""))[[1]], "")

  context[i]    <- paste0(temp[1], "_", temp[7])
  alteration[i] <- paste0(temp[3], temp[4], temp[5])
}

# Melt the matrix using the signatures as variable
matrixW_melt <- melt(matrixW_norm)

# Add columns for the mutation type and the sequence context
matrixW_norm <- cbind(matrixW_norm, alteration, context)
# Reorder (alteration) for having the same order as in the matrice of published signatures
matrixW_norm <- matrixW_norm[order(matrixW_norm[,"alteration"], matrixW_norm[,"context"]), ]
# Reorder (columns) for having the same order as in the matrice of published signatures
matrixW_norm <- cbind(matrixW_norm[,c("alteration", "context")], matrixW_norm[,1:(ncol(matrixW_norm)-2)]) # Put the column alteration and context at the begining
# Save the matrix
write.table(matrixW_norm, file=output_matrixW, quote=F, sep="\t", col.names=T, row.names=F)

# Add columns for the mutation type and the sequence context
matrixW_melt <- cbind(matrixW_melt, alteration)
matrixW_melt <- cbind(matrixW_melt, context)
# Rename the columns
colnames(matrixW_melt) <- c("", "Signature", "value", "alteration", "context")

# Save the input for ggplot2
input_ggplot2 <- as.matrix(matrixW_melt)
input_ggplot2 <- input_ggplot2[,2:ncol(input_ggplot2)]
write.table(input_ggplot2, file=output_matrixW_ggplot2, quote=F, sep="\t", col.names=T, row.names=F)

# Maximum value of the y axis
max_matrixW <- as.numeric(max(matrixW_melt$value))


# Base plot
p <- ggplot(matrixW_melt, aes(x=context, y=value, fill=alteration)) + geom_bar(stat="identity", width=0.5) + facet_grid(Signature ~ alteration, scales="free_y")
# Color the mutation types
p <- p + scale_fill_manual(values=c("blue", "black", "red", "#828282", "#00CC33", "pink"))
# Remove the legend
p <- p + guides(fill=FALSE)
# Customized theme (no background, no facet grid and strip, y axis only)
p <- p + mytheme
# Remove the title of the x facet strip
p <- p + theme(strip.text.x=element_blank())
# Remove the x axis ticks and title
p <- p + theme(axis.title.x=element_blank(), axis.ticks.x = element_blank(), axis.title.y=element_text(size=5))
# Rename the y axis
p <- p + ylab("% contribution to signatures")
# Set the maximum value of the y axis to the maximum value of the matrix W
p <- p + scale_y_continuous(limits=c(0,max_matrixW), oob=squish, breaks=c(0,round(max_matrixW)))
# Save some space for adding the sequence context at the bottom
p <- p + theme(plot.margin=unit(c(.3, 0, 0, 0), "cm"))
p <- p + scale_x_discrete(breaks = c("A_A","A_C","A_G","A_T", "C_A","C_C","C_G","C_T", "G_A","G_C","G_G","G_T", "T_A","T_C","T_G","T_T"),
                          labels =c('A\nA',"\nC","\nG","\nT", 'C\nA',"\nC","\nG","\nT",
                                    'G\nA',"\nC","\nG","\nT", 'T\nA',"\nC","\nG","\nT")
                         )


#------------------------------------------------------------------------------------------------------------------------------
# Change the color of the facets for the mutation type
#------------------------------------------------------------------------------------------------------------------------------
cols <- rep( c("blue", "black", "red", "#828282", "#00CC33", "pink")) # Facet strip colours

# Make a grob object
Pg <- ggplotGrob(p)
# To keep track of strip.background grobs
idx <- 0
# Find each strip.background and alter its backround colour
for( g in 1:length(Pg$grobs) )
{
  if( grepl( "strip.absoluteGrob" , Pg$grobs[[g]]$name ) )
  {
    idx <- idx + 1
    sb <- which( grepl( "strip\\.background" , names( Pg$grobs[[g]]$children ) ) )
    Pg$grobs[[g]]$children[[sb]][]$gp$fill <- cols[idx]
  }
}

# Reduce the size of the facet strip
Pg$heights[[3]] = unit(.05,"cm")


#------------------------------------------------------------------------------------------------------------------------------
# Save the graph in a png file
#------------------------------------------------------------------------------------------------------------------------------
options(bitmapType='cairo')
png(figure_matrixW_png, width=1300, heigh=500, res=300, pointsize = 4)
plot(Pg)
## Add label for the mutation type above the strip facet
grid.text(0.12, unit(1,"npc") - unit(1.4,"line"), label="C>A")
grid.text(0.27, unit(1,"npc") - unit(1.4,"line"), label="C>G")
grid.text(0.42, unit(1,"npc") - unit(1.4,"line"), label="C>T")
grid.text(0.58, unit(1,"npc") - unit(1.4,"line"), label="T>A")
grid.text(0.74, unit(1,"npc") - unit(1.4,"line"), label="T>C")
grid.text(0.89, unit(1,"npc") - unit(1.4,"line"), label="T>G")
invisible( dev.off() )



                                        ###############################################################################
                                        #            Contribution of mutational signature in each samples             #
                                        ###############################################################################

# Recover the total number of SBS per samples
NbSBS <- colSums(matrixNMF)
# Normalized matrix H  to 100%
matrixH_norm <- t((t(matrixH)/colSums(matrixH))*100)
# Add a row name
rownames(matrixH_norm) <- rownames(matrixH_norm, do.NULL = FALSE, prefix = "row")
# Replace the signature number by letter
rownames(matrixH_norm) <- sapply(1:length(rownames(matrixH_norm)), function(x) { ConvertNb2Aphabet(rownames(matrixH_norm)[x]) } )

## Combined the contribution with the total number of SBS
matrixH_norm_melt <- melt(matrixH_norm)
matrixH_norm_melt <- cbind(matrixH_norm_melt, rep(NbSBS, each = opt$nbSignature))
colnames(matrixH_norm_melt) <- c("Signature", "Sample", "Value", "Total_SBS")

# Calculate the contribution in number of SBS
matrixH_norm_melt$ContriSBS <- sapply(1:nrow(matrixH_norm_melt), function(x) { Contri2SignSBS(matrixH_norm_melt$Total_SBS[x], matrixH_norm_melt$Value[x]) } )


# Save the matrix
write.table(matrixH_norm_melt, file=output_matrixH_ggplot2, quote=F, sep="\t", col.names=T, row.names=F)

# Base plot for the contribution of each samples according the count of mutations
p2 <- ggplot(matrixH_norm_melt, aes(x=reorder(Sample, -ContriSBS), y=ContriSBS, fill=Signature)) + geom_bar(stat="identity") + theme_classic()
# Remove the name of samples
p2 <- p2 + theme(axis.text.x = element_blank()) 
# Reverse the y axis
p2 <- p2 + scale_y_reverse()
# Rename the y and x axis
p2 <- p2 + ylab("Number of mutations") + xlab("Samples")
# Remove the x axis line
p2 <- p2 + theme(axis.line.x=element_blank())

# Base plot for the contribution of each samples in percentages
p3 <- ggplot(matrixH_norm_melt, aes(x=reorder(Sample, -ContriSBS), y=Value, fill=Signature)) + geom_bar(stat="identity") + theme_classic() + theme(axis.text.x = element_blank()) + xlab("") + ylab("% of mutations")
# Remove the x axis line
p3 <- p3 + theme(axis.line.x=element_blank(), axis.ticks.x=element_blank())


# Plot PNG
png(figure_matrixH_png, width=3000, heigh=2000, res=300)
# Combined the two plots for the contribution of the samples
suppressWarnings( grid_arrange_shared_legend(p3, p2) )
invisible( dev.off() )


                                        ###############################################################################
                                        #            Average contributions of each signature in each cluster          #
                                        ###############################################################################

matrixH_cluster           <- cbind(matrix_cluster[,1], t(matrixH_norm))
colnames(matrixH_cluster) <- c("Cluster", colnames(t(matrixH_norm)))

df <- as.data.frame(matrixH_cluster)

tmp_mat <- sapply(1:opt$nbSignature, function(x) { meanCluster(df[df[,1] == x,]) } )
# Add a name for the row and the col
rownames(tmp_mat) <- sapply(1:opt$nbSignature, function(x) { paste0("Sig. ", x) } )
colnames(tmp_mat) <- sapply(1:opt$nbSignature, function(x) { paste0("Cluster ", x) } )
tmp_mat           <- t(tmp_mat)
# Recover the number of samples in each cluster
nbSampleByCluster <- sapply(1:opt$nbSignature, function(x) { as.numeric( strsplit( as.character(dim(df[df[,1] == x,])), " " ) )  } )
# Combined the average contribution and the number of samples
tmp_mat <- cbind(tmp_mat, nbSampleByCluster[1,])
# Add a name for the row and the col
colnames(tmp_mat)[opt$nbSignature+1] <- "Number of samples"
# Save the matrix
write.table(tmp_mat, file=output_matrixH_cluster, quote=F, sep="\t", col.names=T, row.names=T)

## Create an image of the table with ggplot2
# Dummy plot
p4 <- qplot(1:10, 1:10, geom = "blank") + 
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_rect(fill=NA,color="white", size=0.5, linetype="solid"),
                  axis.line = element_blank(),
                  axis.ticks = element_blank(),
                  panel.background = element_rect(fill="white"),
                  plot.background = element_rect(fill="white"),
                  legend.position = "none", 
                  axis.text = element_blank(),
                  axis.title  = element_blank()
                 )
# Adding a table
p4 <- p4 + annotation_custom(grob = tableGrob(tmp_mat),
                             xmin = 4, xmax = 7,
                             ymin = 0, ymax = 10)

# Save the table
png(figure_matrixH_cluster, width=2500, heigh=1000, res=300)
# Combined the two plots for the contribution of the samples
plot(p4)
invisible( dev.off() )


# Delete the empty plot created by the script
if (file.exists("Rplots.pdf")) invisible( file.remove("Rplots.pdf") )
