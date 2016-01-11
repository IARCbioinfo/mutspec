#!/usr/bin/Rscript

#---------------------------------------------#
# Author: Maude                               #
# Script: transcriptionalStrandBias_Galaxy.r  #
# Last update: 03/07/15                       #
#---------------------------------------------#

#########################################################################################################################################
#                                                             Transcriptional strand bias                                               #
#########################################################################################################################################

#-------------------------------------------------------------------------------
# Print a usage message if there is no argument pass to the command line
#-------------------------------------------------------------------------------
args <- commandArgs(TRUE)
usage <- function() 
{
  msg <- paste0('Usage:\n',
                ' transcriptionalStrandBias_Galaxy.r input Output_Folder_High_Resolution Output_Folder_Low_Resolution Label_Y_axis\n',
                '\ninput should be tab-separated: MutationTypeContext Strand Value Sample\n',
                '\nOutput_Folder_High_Resolution: Folder for saving the high resolution image (display on the HTML page)\n',
                '\nOutput_Folder_Low_Resolution:  Folder for saving the low resolution image (display on the Excel report)\n',
                '\nLabel_Y_axis: can be Count or Percent'
  )
  cat(msg, '\n', file="/dev/stderr")
  quit(status=1)
}

input = args[length(args)]

if (length(args) == 0) { usage() }



#-------------------------------------------------------------------------------
# Load library
#-------------------------------------------------------------------------------
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(gridExtra)))


#-------------------------------------------------------------------------------
# Recover the argument pass in the command line
#-------------------------------------------------------------------------------
input         <- args[1]
output        <- args[2]
output_temp   <- args[3] # Temp folder for the plot present in the Excel report
legend_y_axis <- args[4]


#-------------------------------------------------------------------------------
# Create the plot
#-------------------------------------------------------------------------------
## Load the data
txnSB <- read.table(input, header=T)
## Define the color for the transcribed (blue) and non-transcribed strand(red)
cb_palette_SB <- c("#0072B2", "#CC0000")
## Reorder the mutation on the x axis (same order as NMF)
txnSB$MutationTypeContext <- factor(txnSB$MutationTypeContext,
                                    levels=c(
            "C>A:A_A","C>A:A_C","C>A:A_G","C>A:A_T","C>A:C_A","C>A:C_C","C>A:C_G","C>A:C_T","C>A:G_A","C>A:G_C","C>A:G_G","C>A:G_T","C>A:T_A","C>A:T_C","C>A:T_G","C>A:T_T",
            "C>G:A_A","C>G:A_C","C>G:A_G","C>G:A_T","C>G:C_A","C>G:C_C","C>G:C_G","C>G:C_T","C>G:G_A","C>G:G_C","C>G:G_G","C>G:G_T","C>G:T_A","C>G:T_C","C>G:T_G","C>G:T_T",
            "C>T:A_A","C>T:A_C","C>T:A_G","C>T:A_T","C>T:C_A","C>T:C_C","C>T:C_G","C>T:C_T","C>T:G_A","C>T:G_C","C>T:G_G","C>T:G_T","C>T:T_A","C>T:T_C","C>T:T_G","C>T:T_T",
            "T>A:A_A","T>A:A_C","T>A:A_G","T>A:A_T","T>A:C_A","T>A:C_C","T>A:C_G","T>A:C_T","T>A:G_A","T>A:G_C","T>A:G_G","T>A:G_T","T>A:T_A","T>A:T_C","T>A:T_G","T>A:T_T",
            "T>C:A_A","T>C:A_C","T>C:A_G","T>C:A_T","T>C:C_A","T>C:C_C","T>C:C_G","T>C:C_T","T>C:G_A","T>C:G_C","T>C:G_G","T>C:G_T","T>C:T_A","T>C:T_C","T>C:T_G","T>C:T_T",
            "T>G:A_A","T>G:A_C","T>G:A_G","T>G:A_T","T>G:C_A","T>G:C_C","T>G:C_G","T>G:C_T","T>G:G_A","T>G:G_C","T>G:G_G","T>G:G_T","T>G:T_A","T>G:T_C","T>G:T_G","T>G:T_T"
            																)
																	 )
## Create a bar plot with custom color and classic theme
p_txnSB <- ggplot(txnSB, aes(x=MutationTypeContext, y=Value, fill=Strand))
# Add a background for better differentiate the different mutation types
p_txnSB <- p_txnSB + geom_rect(data=NULL,aes(xmin=0.25,xmax=16.5,ymin=-Inf,ymax=Inf), fill="#E5E5E5") +
                     geom_rect(data=NULL,aes(xmin=16.5,xmax=32.5,ymin=-Inf,ymax=Inf), fill="#EDEDED") +
                     geom_rect(data=NULL,aes(xmin=32.5,xmax=48.5,ymin=-Inf,ymax=Inf), fill="#E5E5E5") +
                     geom_rect(data=NULL,aes(xmin=48.5,xmax=80.5,ymin=-Inf,ymax=Inf), fill="#EDEDED") +
                     geom_rect(data=NULL,aes(xmin=64.5,xmax=80.5,ymin=-Inf,ymax=Inf), fill="#E5E5E5") +
                     geom_rect(data=NULL,aes(xmin=80.5,xmax=96.5,ymin=-Inf,ymax=Inf), fill="#EDEDED")
# Add the bar
p_txnSB <- p_txnSB + geom_bar(stat="identity", width=0.5) + theme_classic() + scale_fill_manual(values=cb_palette_SB)


# Rename the y axis
p_txnSB <- p_txnSB + ylab(legend_y_axis)
## Set the legend position to the top of plot and remove the legend title
p_txnSB <- p_txnSB + theme(legend.position="top") + labs(fill="")
## Add margins for having place to add the mutation type labels bellow the bar graph
p_txnSB <- p_txnSB + theme(plot.margin=unit(c(1,1,-.1,1.5), "cm"))
## Rename the x labels
p_txnSB <- p_txnSB + scale_x_discrete(name="",
                                      breaks=c(
              "T>G:A_A","T>G:A_C","T>G:A_G","T>G:A_T","T>G:C_A","T>G:C_C","T>G:C_G","T>G:C_T","T>G:G_A","T>G:G_C","T>G:G_G","T>G:G_T","T>G:T_A","T>G:T_C","T>G:T_G","T>G:T_T",
              "T>C:A_A","T>C:A_C","T>C:A_G","T>C:A_T","T>C:C_A","T>C:C_C","T>C:C_G","T>C:C_T","T>C:G_A","T>C:G_C","T>C:G_G","T>C:G_T","T>C:T_A","T>C:T_C","T>C:T_G","T>C:T_T",
              "T>A:A_A","T>A:A_C","T>A:A_G","T>A:A_T","T>A:C_A","T>A:C_C","T>A:C_G","T>A:C_T","T>A:G_A","T>A:G_C","T>A:G_G","T>A:G_T","T>A:T_A","T>A:T_C","T>A:T_G","T>A:T_T",
              "C>A:A_A","C>A:A_C","C>A:A_G","C>A:A_T","C>A:C_A","C>A:C_C","C>A:C_G","C>A:C_T","C>A:G_A","C>A:G_C","C>A:G_G","C>A:G_T","C>A:T_A","C>A:T_C","C>A:T_G","C>A:T_T",
              "C>G:A_A","C>G:A_C","C>G:A_G","C>G:A_T","C>G:C_A","C>G:C_C","C>G:C_G","C>G:C_T","C>G:G_A","C>G:G_C","C>G:G_G","C>G:G_T","C>G:T_A","C>G:T_C","C>G:T_G","C>G:T_T",
              "C>T:A_A","C>T:A_C","C>T:A_G","C>T:A_T","C>T:C_A","C>T:C_C","C>T:C_G","C>T:C_T","C>T:G_A","C>T:G_C","C>T:G_G","C>T:G_T","C>T:T_A","C>T:T_C","C>T:T_G","C>T:T_T",
              "G>A:A_A","G>A:A_C","G>A:A_G","G>A:A_T","G>A:C_A","G>A:C_C","G>A:C_G","G>A:C_T","G>A:G_A","G>A:G_C","G>A:G_G","G>A:G_T","G>A:T_A","G>A:T_C","G>A:T_G","G>A:T_T",
              "G>C:A_A","G>C:A_C","G>C:A_G","G>C:A_T","G>C:C_A","G>C:C_C","G>C:C_G","G>C:C_T","G>C:G_A","G>C:G_C","G>C:G_G","G>C:G_T","G>C:T_A","G>C:T_C","G>C:T_G","G>C:T_T",
              "G>T:A_A","G>T:A_C","G>T:A_G","G>T:A_T","G>T:C_A","G>T:C_C","G>T:C_G","G>T:C_T","G>T:G_A","G>T:G_C","G>T:G_G","G>T:G_T","G>T:T_A","G>T:T_C","G>T:T_G","G>T:T_T",
              "T>A:A_A","T>A:A_C","T>A:A_G","T>A:A_T","T>A:C_A","T>A:C_C","T>A:C_G","T>A:C_T","T>A:G_A","T>A:G_C","T>A:G_G","T>A:G_T","T>A:T_A","T>A:T_C","T>A:T_G","T>A:T_T",
              "T>C:A_A","T>C:A_C","T>C:A_G","T>C:A_T","T>C:C_A","T>C:C_C","T>C:C_G","T>C:C_T","T>C:G_A","T>C:G_C","T>C:G_G","T>C:G_T","T>C:T_A","T>C:T_C","T>C:T_G","T>C:T_T",
              "T>G:A_A","T>G:A_C","T>G:A_G","T>G:A_T","T>G:C_A","T>G:C_C","T>G:C_G","T>G:C_T","T>G:G_A","T>G:G_C","T>G:G_G","T>G:G_T","T>G:T_A","T>G:T_C","T>G:T_G","T>G:T_T"
              																),
                   										labels=c(
              "A_A","A_C","A_G","A_T","C_A","C_C","C_G","C_T","G_A","G_C","G_G","G_T","T_A","T_C","T_G","T_T",
              "A_A","A_C","A_G","A_T","C_A","C_C","C_G","C_T","G_A","G_C","G_G","G_T","T_A","T_C","T_G","T_T",
              "A_A","A_C","A_G","A_T","C_A","C_C","C_G","C_T","G_A","G_C","G_G","G_T","T_A","T_C","T_G","T_T",
              "A_A","A_C","A_G","A_T","C_A","C_C","C_G","C_T","G_A","G_C","G_G","G_T","T_A","T_C","T_G","T_T",
              "A_A","A_C","A_G","A_T","C_A","C_C","C_G","C_T","G_A","G_C","G_G","G_T","T_A","T_C","T_G","T_T",
              "A_A","A_C","A_G","A_T","C_A","C_C","C_G","C_T","G_A","G_C","G_G","G_T","T_A","T_C","T_G","T_T",
              "A_A","A_C","A_G","A_T","C_A","C_C","C_G","C_T","G_A","G_C","G_G","G_T","T_A","T_C","T_G","T_T",
              "A_A","A_C","A_G","A_T","C_A","C_C","C_G","C_T","G_A","G_C","G_G","G_T","T_A","T_C","T_G","T_T",
              "A_A","A_C","A_G","A_T","C_A","C_C","C_G","C_T","G_A","G_C","G_G","G_T","T_A","T_C","T_G","T_T",
              "A_A","A_C","A_G","A_T","C_A","C_C","C_G","C_T","G_A","G_C","G_G","G_T","T_A","T_C","T_G","T_T",
              "A_A","A_C","A_G","A_T","C_A","C_C","C_G","C_T","G_A","G_C","G_G","G_T","T_A","T_C","T_G","T_T",
              "A_A","A_C","A_G","A_T","C_A","C_C","C_G","C_T","G_A","G_C","G_G","G_T","T_A","T_C","T_G","T_T"
              																)
                  									 )
## Changing the appearance of x axis thicks
p_txnSB <- p_txnSB + theme(axis.text.x = element_text(angle=60, hjust=1, vjust=1))
## Close graphics windows
graphics.off()
## Save the plot for the HTML page (higher resolution)
options(bitmapType='cairo') # # Use cairo device as isn't possible to install X11 on the server...
png(paste0(output, ".png"), width=4000, height=1000, res=300)
plot(p_txnSB)
# Add a label bellow the bar graph for indicating the mutation type
grid.text(paste("C>A", sep=""), x=unit(.14, "npc"), y=unit(.7, "npc"), just=c("left", "bottom"), gp=gpar(fontface="bold",fontsize=10))
grid.text(paste("C>G", sep=""), x=unit(.29, "npc"),  y=unit(.7, "npc"), just=c("left", "bottom"), gp=gpar(fontface="bold",fontsize=10))
grid.text(paste("C>T", sep=""), x=unit(.45, "npc"), y=unit(.7, "npc"), just=c("left", "bottom"), gp=gpar(fontface="bold",fontsize=10))
grid.text(paste("T>A", sep=""), x=unit(.58, "npc"),  y=unit(.7, "npc"), just=c("left", "bottom"), gp=gpar(fontface="bold",fontsize=10))
grid.text(paste("T>C", sep=""), x=unit(.74, "npc"), y=unit(.7, "npc"), just=c("left", "bottom"), gp=gpar(fontface="bold",fontsize=10))
grid.text(paste("T>G", sep=""), x=unit(.9, "npc"),  y=unit(.7, "npc"), just=c("left", "bottom"), gp=gpar(fontface="bold",fontsize=10))
invisible( dev.off() )



# Save the plot for the report
p_txnSB
ggsave(paste0(output_temp, "-Report.png"), width=18)

# Delete the empty plot created by the script
if (file.exists("Rplots.pdf")) invisible( file.remove("Rplots.pdf") )
