#!/usr/bin/Rscript

#-----------------------------------#
# Author: Maude                     #
# Script: mutationSpectra_Galaxy.r  #
# Last update: 23/07/15             #
#-----------------------------------#

#########################################################################################################################################
#                                                   Represent the mutation spectra with a bar graph                                     #
#########################################################################################################################################

#-------------------------------------------------------------------------------
# Print a usage message if there is no argument pass to the command line
#-------------------------------------------------------------------------------
args <- commandArgs(TRUE)
usage <- function() 
{
  msg <- paste0('Usage:\n',
                ' mutationSpectra_Galaxy.r input_Mutation_Spectra Sample_Name Output_Folder_High_Resolution Output_Folder_Low_Resolution Count_ca Count_cg Count_ta Count_tc Count_tg\n',
                '\ninput_Mutation_Spectra should be tab-separated: alteration context value\n',
                '\nOutput_Folder_High_Resolution: Folder for saving the high resolution image (display on the HTML page)\n',
                '\nOutput_Folder_Low_Resolution:  Folder for saving the low resolution image (display on the Excel report)\n'
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
suppressMessages(suppressWarnings(library(reshape)))
suppressMessages(suppressWarnings(library(grid)))
suppressMessages(suppressWarnings(library(scales)))
suppressMessages(suppressWarnings(library(gridExtra)))



#-------------------------------------------------------------------------------
# Recover the arguments
#-------------------------------------------------------------------------------
input                <- args[1]
sampleName           <- args[2]
output_html          <- args[3]
output_report        <- args[4]
count_ca             <- as.numeric(args[5])
count_cg             <- as.numeric(args[6])
count_ct             <- as.numeric(args[7])
count_ta             <- as.numeric(args[8])
count_tc             <- as.numeric(args[9])
count_tg             <- as.numeric(args[10])

count_ca             <- paste("C>A (", count_ca, ")")
count_cg             <- paste("C>G (", count_cg, ")")
count_ct             <- paste("C>T (", count_ct, ")")
count_ta             <- paste("T>A (", count_ta, ")")
count_tc             <- paste("T>C (", count_tc, ")")
count_tg             <- paste("T>G (", count_tg, ")")



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
if(device[5]) { font <- "Helvetica" } else { font <- "mono" }

#-------------------------------------------------------------------------------
# My own thme
#-------------------------------------------------------------------------------
theme_custom <- function(base_size = 12, base_family = "")
{
  # Starts with theme_grey and then modify some parts
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
             theme(
                   axis.text         = element_text(size = rel(1), family=font),
                   axis.ticks        = element_line(colour = "black"),
                   axis.line         = element_line(colour = "black", size = .5),
                   legend.key        = element_blank(),
                   panel.background  = element_blank(),
                   panel.border      = element_blank(),
                   panel.grid.major  = element_blank(),
                   panel.grid.minor  = element_blank()
                  )
}


#-------------------------------------------------------------------------------
# Customize the theme for adding a y axis
#-------------------------------------------------------------------------------
mytheme                    <- theme_custom()
mytheme$axis.line.x        <- mytheme$axis.line.y <- mytheme$axis.line
mytheme$axis.line.x$colour <- 'white'

#-------------------------------------------------------------------------------
# Set the decimal precision to 0.0
#------------------------------------------------------------------------------
fmt <- function()
{
	function(x) format(x,nsmall = 1,scientific = FALSE, digits=1)
}



                                        ###############################################################################
                                        #                                      MAIN                                   #
                                        ###############################################################################

matrixW_inputggplot2 <- read.table(input, header=T)
matrixW_melt         <- melt(matrixW_inputggplot2)
max_matrixW          <- max(matrixW_inputggplot2[,3:ncol(matrixW_inputggplot2)])


p <- ggplot(matrixW_melt, aes(x=context, y=value, fill=alteration)) + geom_bar(stat="identity", width=0.5) + facet_grid(variable ~ alteration, scales="free_y")
# Color the mutation like Alexandrov et al.
p <- p + scale_fill_manual(values=c("blue", "black", "red", "#828282", "#00CC33", "pink"))
# Remove the legend
p <- p + guides(fill=FALSE)
# customized theme (no background, no facet grid and strip, y axis only)
p <- p + mytheme
# Remove the title of the x facet strip
p <- p + theme(strip.text.x=element_blank(), strip.text.y=element_blank())
# Remove the x axis label, thicks and title
p <- p + theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_text(size=15))
# Scale the y axis to the maximum value
p <- p + scale_y_continuous(limits=c(0,max_matrixW), oob=squish, breaks=c(0,max_matrixW), labels=fmt())
# Rename the y axis
p <- p + ylab("percent")
# Add a title to the plot
p <- p + ggtitle(sampleName) + theme(plot.title = element_text(vjust = 3.4, family=font))
# Add a top margin for writing the title of the plot
p <- p + theme(plot.margin=unit(c(.7,0,0,0), "cm"))
p <- p + scale_x_discrete(breaks = c("A_A","A_C","A_G","A_T", "C_A","C_C","C_G","C_T", "G_A","G_C","G_G","G_T", "T_A","T_C","T_G","T_T"),
		                      labels =c('A\nA',"\nC","\nG","\nT", 'C\nA',"\nC","\nG","\nT",
			                     	        'G\nA',"\nC","\nG","\nT", 'T\nA',"\nC","\nG","\nT"
                                   )
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
Pg$heights[[3]] = unit(.1,"cm")


# Save the plot for the HTML page (higher resolution)
graphics.off() # close graphics windows
# Use cairo device as isn't possible to install X11 on the server...
png(paste0(output_html, "/", sampleName, "-MutationSpectraPercent-Genomic.png"), width=3500, heigh=500, res=300, type=c("cairo-png"))
plot(Pg)
## Add label for the mutation type above the strip facet
grid.text(0.13, unit(0.90,"npc") - unit(1,"line"), label=count_ca)
grid.text(0.29, unit(0.90,"npc") - unit(1,"line"), label=count_cg)
grid.text(0.45, unit(0.90,"npc") - unit(1,"line"), label=count_ct)
grid.text(0.6, unit(0.90,"npc") - unit(1,"line"), label=count_ta)
grid.text(0.76, unit(0.90,"npc") - unit(1,"line"), label=count_tc)
grid.text(0.92, unit(0.90,"npc") - unit(1,"line"), label=count_tg)
invisible( dev.off() )

# Save the plot for the report
png(paste0(output_report, "/", sampleName, "-MutationSpectraPercent-Genomic-Report.png"), width=1000, heigh=150, type=c("cairo-png"))
plot(Pg)
## Add label for the mutation type above the strip facet
grid.text(0.13, unit(0.90,"npc") - unit(1,"line"), label=count_ca)
grid.text(0.29, unit(0.90,"npc") - unit(1,"line"), label=count_cg)
grid.text(0.45, unit(0.90,"npc") - unit(1,"line"), label=count_ct)
grid.text(0.6, unit(0.90,"npc") - unit(1,"line"), label=count_ta)
grid.text(0.76, unit(0.90,"npc") - unit(1,"line"), label=count_tc)
grid.text(0.92, unit(0.90,"npc") - unit(1,"line"), label=count_tg)
invisible( dev.off() )

# Delete the empty plot created by the script
if (file.exists("Rplots.pdf")) invisible( file.remove("Rplots.pdf") )
