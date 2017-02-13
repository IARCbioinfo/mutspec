#!/usr/bin/Rscript

#-----------------------------------#
# Author: Maude                     #
# Script: figs_MutSpecStat_Galaxy.r #
# Last update: 18/10/16             #
#-----------------------------------#


#########################################################################################################################################
#                                     Create the figures for the report and the HTML page                          #
#########################################################################################################################################


#-------------------------------------------------------------------------------
# Load library for recovering the arguments
#-------------------------------------------------------------------------------
suppressMessages(suppressWarnings(require("getopt")))



#-------------------------------------------------------------------------------
# Recover the arguments
#-------------------------------------------------------------------------------
spec = matrix(c(
  "folderFigure", "folderFigure", 1, "character",
  "folderTemp",   "folderTemp",   1, "character",
  "filename",     "filename",     1, "character",
  "help",         "h",            0, "logical"
),
byrow=TRUE, ncol=4
)

opt = getopt(spec)

# No argument is pass to the command line
if(length(opt) == 1)
{
  cat(paste("Usage:\n figs_MutSpecStat_Galaxy.r --folderFigure <path_to_folder> --folderTemp <path_to_tempFolder> --filename <filename> \n",sep=""))
  q(status=1)
}

# Help was asked for.
if ( !is.null(opt$help) )
{
  # print a friendly message and exit with a non-zero error code
  cat(paste("Usage:\n  figs_MutSpecStat_Galaxy.r --folderFigure <path_to_folder> --filename <filename> \n",sep=""))
  q(status=1)
}



#-------------------------------------------------------------------------------
# Load library
#-------------------------------------------------------------------------------
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(gplots)))
suppressMessages(suppressWarnings(library(gtable)))




#-------------------------------------------------------------------------------
# OVERALL MUTATION DISTRIBUTION
#-------------------------------------------------------------------------------
inputDistrMut  <- paste0(opt$folderFigure, "/Overall_mutation_distribution/", opt$filename, "/", opt$filename, "-OverallMutationDistribution.txt")
outputDistrMut <- paste0(opt$folderFigure, "/Overall_mutation_distribution/", opt$filename, "/", opt$filename, "-OverallMutationDistribution.png")

# Load the input file
distrMut      <- read.table(inputDistrMut, header=T)

# Add the count of each category in the legend
distrMut$Legend[[1]] <- paste0(distrMut$Variant_type[[1]], " (", distrMut$Count[[1]], ")")
distrMut$Legend[[2]] <- paste0(distrMut$Variant_type[[2]], " (", distrMut$Count[[2]], ")")
distrMut$Legend[[3]] <- paste0(distrMut$Variant_type[[3]], " (", distrMut$Count[[3]], ")")

# Base plot
pie <- ggplot(distrMut, aes(x=factor(""), fill=Legend, weight=Count)) + geom_bar(width=1) + coord_polar(theta="y") + scale_x_discrete("", breaks=NULL) + scale_y_continuous("", breaks=NULL) + labs(fill="")
# Background of the plot entire white
pie <- pie + theme(panel.grid.major = element_line(colour="white"), panel.grid.minor = element_line(colour="white"), panel.background = element_rect(fill="white"))
# Legend on right in 3 rows
pie <- pie + theme(legend.position="bottom") +  guides(fill=guide_legend(nrow=3))
# Change the color and the title of the legend
pie <- pie + scale_fill_brewer("Variant type", palette="Set1")
# Remove all the margins
pie <- pie + theme(plot.margin=unit(c(-1, 0, -1.5, 0), "cm"))
# Save the pie chart for the HTML page (higher resolution)
options(bitmapType='cairo') # Use cairo device as isn't possible to install X11 on the server...
png(outputDistrMut, width=700, height=1100, res=300)
print(pie)
dev.off()



#-------------------------------------------------------------------------------
# SBS MUTATION DISTRIBUTION
#-------------------------------------------------------------------------------
inputDistrSBS        <- paste0(opt$folderFigure, "/SBS_distribution/", opt$filename, "/", opt$filename, "-SBS_distribution.txt")
outputDistrSBS       <- paste0(opt$folderFigure, "/SBS_distribution/", opt$filename, "/", opt$filename, "-SBS_distribution.png")
outputDistrSBSReport <- paste0(opt$folderTemp, "/", opt$filename, "-SBS_distribution-Report.png")


# Load the input file
distrSBS <- read.delim(inputDistrSBS)
distrSBS <- data.frame(distrSBS)
bar <- ggplot(distrSBS, aes(x=Mutation_Type, y=Percentage, fill=Mutation_Type))
bar <- bar + geom_bar(stat="identity", width=0.5)
# Theme classic
bar <- bar + theme_classic()
# Remove the axis legend
bar <- bar + xlab("")
# Set the color of the bars and Changing the labels in the legend
bar <- bar + scale_fill_manual(values=c("blue", "black", "red", "gray", "#00CC33", "pink"),
                               labels=c("C:G>A:T", "C:G>G:C", "C:G>T:A", "T:A>A:T", "T:A>C:G", "T:A>G:C")
)
# Remove the label in x axis
bar <- bar + theme(axis.text.x = element_blank())
# Change the name of the y label
bar <- bar + ylab("Percent")
# Save the plot for the HTML page (higher resolution)
png(outputDistrSBS, width=1800, height=1500, res=300)
print(bar)
dev.off()
# Save the plot for the report
bar
ggsave(outputDistrSBSReport)



#-------------------------------------------------------------------------------
# IMPACT ON PROTEIN
#-------------------------------------------------------------------------------
inputImpactProt        <- paste0(opt$folderFigure, "/Impact_protein_sequence/", opt$filename, "/", opt$filename, "-DistributionExoFunc.txt")
outputImpactProt       <- paste0(opt$folderFigure, "/Impact_protein_sequence/", opt$filename, "/", opt$filename, "-DistributionExoFunc.png")
outputImpactProtReport <- paste0(opt$folderTemp, "/", opt$filename, "-DistributionExoFunc-Report.png")

# Load the input file
impactProt <- read.table(inputImpactProt, header=T)
# Custom palette: black, orange, dark green, yellow, light blue, dark blue, darkslateblue, red, purple, pink, light green, turquoise, gray
cb_palette <- c("#000000", "#E69F00", "#006600", "#660000", "#F0E442", "#56B4E9", "#3300FF", "#483D8B", "#FF0000", "#9900CC", "#FF66CC", "#00CC00", "#66FFFF", "#C0C0C0")
pie <- ggplot(impactProt, aes(x=factor(""), fill=AA_Change, weight=Count)) + geom_bar(width=1) + coord_polar(theta="y") + scale_x_discrete("", breaks=NULL)+ scale_y_continuous("", breaks=NULL) + scale_fill_manual(values=cb_palette)
# Background of the plot entire white
pie <- pie + theme(panel.grid.major = element_line(colour="white"), panel.grid.minor = element_line(colour="white"), panel.background = element_rect(fill="white"))
# Legend in two column
pie <- pie + guides(fill=guide_legend(ncol=2)) + theme(legend.position="bottom")
# Remove the legend title
pie <- pie + labs(fill="")
# Save the plot for the HTML page (higher resolution)
png(outputImpactProt, width=1600, height=1800, res=300)
print(pie)
dev.off()
# Save the plot for the report
pie
ggsave(outputImpactProtReport)




#-------------------------------------------------------------------------------
# STRAND BIAS
#-------------------------------------------------------------------------------
inputSB        <- paste0(opt$folderFigure, "/Stranded_Analysis/", opt$filename, "/", opt$filename, "-StrandBias.txt")
outputSB       <- paste0(opt$folderFigure, "/Stranded_Analysis/", opt$filename, "/", opt$filename, "-StrandBias.png")
outputSBReport <- paste0(opt$folderTemp, "/", opt$filename, "-StrandBias-Report.png")

# Load the input file
file_sb       <- read.table(inputSB, header=T)
# Custom palette (blue, red)
cb_palette_SB <- c("#0072B2", "#CC0000")
# Base plot
p_sb          <- ggplot(file_sb, aes(x=Alteration, y=Count, fill=Strand)) + theme_classic() + geom_bar(stat="identity", position="dodge") + scale_fill_manual(values=cb_palette_SB) + theme(axis.text.x = element_text(angle=60, hjust=1)) + xlab("") + theme(legend.position="bottom")
# Save the plot for the HTML page (higher resolution)
png(outputSB, width=1000, height=1200, res=300)
print(p_sb)
dev.off()
# Save the plot for the report
p_sb
ggsave(outputSBReport)




#-------------------------------------------------------------------------------
# HEATMAP SEQUENCE CONTEXT - GENOMIC STRAND
#-------------------------------------------------------------------------------
inputHeatmapGenomic        <- paste0(opt$folderFigure, "/Trinucleotide_Sequence_Context/", opt$filename, "/", opt$filename, "-HeatmapCount-Genomic.txt")
outputHeatmapGenomic        <- paste0(opt$folderFigure, "/Trinucleotide_Sequence_Context/", opt$filename, "/", opt$filename, "-HeatmapCount-Genomic.png")
outputHeatmapGenomicReport <- paste0(opt$folderTemp, "/", opt$filename, "-HeatmapCount-Genomic-Report.png")

inputHeatmapGenomicPercent        <- paste0(opt$folderFigure, "/Trinucleotide_Sequence_Context/", opt$filename, "/", opt$filename, "-HeatmapPercent-Genomic.txt")
outputHeatmapGenomicPercent       <- paste0(opt$folderFigure, "/Trinucleotide_Sequence_Context/", opt$filename, "/", opt$filename, "-HeatmapPercent-Genomic.png")
outputHeatmapGenomicPercentReport <- paste0(opt$folderTemp, "/", opt$filename, "-HeatmapPercent-Genomic-Report.png")


## COUNT
heatmap_C <- read.table(inputHeatmapGenomic, header=T)
# Save the plot for the report
png(filename=outputHeatmapGenomicReport, bg="transparent", width=240, height=360)
# Heatmap with an absolute scale
heatmap.2(as.matrix(heatmap_C),Rowv=F,Colv=F,col=colorpanel(384,low="yellow",high="red"),dendrogram="none",scale="none",trace="none",key=F,labRow=rownames(as.matrix(heatmap_C)),labCol=colnames(as.matrix(heatmap_C)),lmat=rbind(c(5,1,4),c(3,1,2)), lhei=c(0.75,0.75),lwid=c(0.5,1.5,0.5))
dev.off()
# Save the plot for the HTML page (higher resolution)
png(filename=outputHeatmapGenomic, width=1100, height=1600, res=300)
heatmap.2(as.matrix(heatmap_C),Rowv=F,Colv=F,col=colorpanel(384,low="yellow",high="red"),dendrogram="none",scale="none",trace="none",key=F,labRow=rownames(as.matrix(heatmap_C)),labCol=colnames(as.matrix(heatmap_C)),lmat=rbind(c(5,1,4),c(3,1,2)), lhei=c(0.75,0.75),lwid=c(0.5,1.5,0.5))
dev.off()

## PERCENT
heatmap_P <- read.table(inputHeatmapGenomicPercent, header=T)
# Save the plot for the report
png(filename=outputHeatmapGenomicPercentReport,bg="transparent", width=240, height=360)
# Heatmap with an absolute scale
heatmap.2(as.matrix(heatmap_P),Rowv=F,Colv=F,col=colorpanel(384,low="yellow",high="red"),dendrogram="none",scale="none",trace="none",key=F,labRow=rownames(as.matrix(heatmap_P)),labCol=colnames(as.matrix(heatmap_P)),lmat=rbind(c(5,1,4),c(3,1,2)), lhei=c(0.75,0.75),lwid=c(0.5,1.5,0.5))
dev.off()
# Save the plot for the HTML page (higher resolution)
png(filename=outputHeatmapGenomicPercent, width=1100, height=1600, res=300)
heatmap.2(as.matrix(heatmap_P),Rowv=F,Colv=F,col=colorpanel(384,low="yellow",high="red"),dendrogram="none",scale="none",trace="none",key=F,labRow=rownames(as.matrix(heatmap_P)),labCol=colnames(as.matrix(heatmap_P)),lmat=rbind(c(5,1,4),c(3,1,2)), lhei=c(0.75,0.75),lwid=c(0.5,1.5,0.5))
dev.off()
