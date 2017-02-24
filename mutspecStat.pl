#!/usr/bin/env perl

#-----------------------------------#
# Author: Maude                     #
# Script: mutspecStat.pl            #
# Last update: 09/02/17             #
#-----------------------------------#

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Basename; # my ($filename, $directories, $suffix) = fileparse($file, qr/\.[^.]*/);
use File::Path;
use Spreadsheet::WriteExcel;

our ($verbose, $man, $help) = (0, 0, 0); # Parse options and print usage if there is a syntax error, or if usage was explicitly requested.
our ($refGenome, $output, $folder_temp, $path_R_Scripts, $path_SeqrefGenome) = ("empty", "empty", "empty", "empty", "empty", "empty");  # The reference genome to use; The path for saving the result; The path for saving the temporary files; The path to R scripts; The path to the fasta reference sequences
our ($poolData, $oneReportPerSample) = (2, 2); # If a folder is pass as input file pool all the data and generate the report on the pool and for each samples; # Generate one report for each samples


GetOptions('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'refGenome=s'=>\$refGenome, 'outfile|o=s' => \$output, 'temp=s' => \$folder_temp, 'pathRscript=s' => \$path_R_Scripts, 'pathSeqRefGenome=s' => \$path_SeqrefGenome, 'poolData' => \$poolData, 'reportSample' => \$oneReportPerSample) or pod2usage(2);

our ($input) = @ARGV;

pod2usage(-verbose=>1, -exitval=>1, -output=>\*STDERR) if ($help);
pod2usage(-verbose=>2, -exitval=>1, -output=>\*STDERR) if ($man);
pod2usage(-verbose=>0, -exitval=>1, -output=>\*STDERR) if(@ARGV == 0); # No argument is pass to the command line print the usage of the script
pod2usage(-verbose=>0, -exitval=>1, -output=>\*STDERR) if(@ARGV == 2); # Only one argument is expected to be pass to @ARGV (the input)




####### The input must be a folder with one or several annotated files
if(!-d $input)
{
	print STDERR "Error: The input must be a Dataset List\n";
	print STDERR "Even for 1 file, please create a Dataset List\n";
	exit;
}


######################################################################################################################################################
#																																			GLOBAL VARIABLES																															 #
######################################################################################################################################################
# Recover the current path
our $pwd = `pwd`;
chomp($pwd);


# Path to R scripts
our $pathRscriptChi2test    = "$path_R_Scripts/R/chi2test_MutSpecStat_Galaxy.r";
our $pathRScriptFigs        = "$path_R_Scripts/R/figs_MutSpecStat_Galaxy.r";
our $pathRScriptTxnSB       = "$path_R_Scripts/R/transciptionalStrandBias.r";
our $pathRScriptMutSpectrum = "$path_R_Scripts/R/mutationSpectra_Galaxy.r";


# The path for saving the files with enough mutations for calculating the statistics;
our $folderCheckedForStat = "$pwd/folder_checked_$filename";
if(!-e $folderCheckedForStat) { mkdir($folderCheckedForStat) or die "$!: $folderCheckedForStat\n"; }

# Output dir with all the results
our $folderMutAnalysis = "";

# Hash table with the length of each chromosomes
our %chromosomes;
# Define the name of the column containing the chromosome, start, ref and alt alleles (based on Annovar output)
our ($chr_name, $start_name, $ref_name, $alt_name) = qw(Chr Start Ref Alt);
# Annovar annotation used
our $func_name       = "Func.refGene";
our $exonicFunc_name = "ExonicFunc.refGene";
our $strand_name     = "Strand";
our $context_name    = "context";
# Font formats
our ($format_A10, $format_A10Boldleft, $format_A10ItalicRed) = ("", "", "");
our ($formatT_left, $formatT_right, $formatT_bottomRight, $formatT_bottomLeft, $formatT_bottom, $formatT_bottomHeader, $formatT_bottomRightHeader, $formatT_bottomHeader2, $formatT_rightHeader);
our ($formatT_graphTitle);
our ($table_topleft, $table_topRight, $table_bottomleft, $table_bottomRight, $table_top, $table_right, $table_bottom, $table_bottomItalicRed, $table_left, $table_bottomrightHeader, $table_left2, $table_middleHeader, $table_middleHeader2);
# Hash table with the result of chi2 test for the strand bias
our %h_chi2 = ();
# For NMF input
our %h_inputNMF = ();


######################################################################################################################################################
#																																								MAIN 																																 #
######################################################################################################################################################
# Check the presence of the flags and create the output and temp directories
CheckFlags();

# First check if the files are annotated or not.
# If the files are annotated check there is enough mutations for generating the statistics, otherwise remove the samples from the analysis
checkVariants();

# Retrieve chromosomes length
checkChrDir();

# Calculate the statistics and generate the report
ReportMutDist();

# Remove the temporary directory
rmtree($folder_temp);
rmtree($folderCheckedForStat);


######################################################################################################################################################
#																																					 	FUNCTIONS																																 #
######################################################################################################################################################

# Check the presence of the flags and create the output and temp directories
sub CheckFlags
{
	# Check the reference genome
	if($refGenome eq "empty")
	{
		print STDERR "Missing flag !\n";
		print STDERR "You forget to specify the name for the reference genome!!!\nPlease specify it with the flag --refGenome\n";
		exit;
	}

	# If no output is specified write the result as the same place as the input file
	if($output eq "empty")
	{
		# The input is a folder with one or more annotated files
		my $directory = dirname( $input );

		$folderMutAnalysis = "$directory/Mutational_Analysis";
		if(!-e $folderMutAnalysis) { mkdir($folderMutAnalysis) or die "$!: $folderMutAnalysis\n"; }
	}
	else
	{
		if(!-e $output) { mkdir($output) or die "$!: $output\n"; }

		$folderMutAnalysis      = "$output/Mutational_Analysis";
		if(!-e $folderMutAnalysis) { mkdir($folderMutAnalysis) or die "$!: $folderMutAnalysis\n"; }
	}

	# If no temp folder is specified write the result in the current path
	my ($filename, $directories, $suffix) = fileparse($input, qr/\.[^.]*/);
	if($folder_temp eq "empty") { $folder_temp = "$pwd/TEMP_MutationalAnalysis_$filename"; }
	if(!-e $folder_temp)        { mkdir($folder_temp) or die "$!: $folder_temp\n"; }

	# Check the path to the R scripts
	if($path_R_Scripts eq "empty")
	{
		print STDERR "Missing flag !\n";
		print STDERR "You forget to specify the path for the R scripts!!!\nPlease specify it with the flag --pathRscript\n";
		exit;
	}


	foreach my $file (`ls $input/*`)
	{
		chomp($file);

		## Verify the name of file, must be <= 31 chars for the sheet name
		my ($filename, $directories, $suffix) = fileparse($file, qr/\.[^.]*/);

		if(length($filename) > 31)
		{
			print STDERR "Error: The filename of: $file\nMust be <= 31 chars\nPlease modify it before running the script\n";
			exit;
		}
	}
}

# Check input file(s)
sub checkVariants
{
	# Count the number of file(s) with enough mutations (at least 1 with a strand orientation)
	my $timerFile       = 0;
	my @listRemovedFile = ();

	foreach my $file (`ls $input/*`)
	{
		chomp($file);

		### Check if the file is annotated
		my $testAnnotation = "toto";
		$testAnnotation    = `grep 'Func.refGene' $file`;

		if($testAnnotation eq "toto")
		{
			print STDERR "Error: The input file you specify is not annotated!\nThe file concerned is: $file !!!!\nPlease first annotate your file before trying to generate the report on mutation spectra\n";
			exit;
		}
		else
		{
			### check if there is at least 1 mutation with a strand info
			my $strand_value = recoverNumCol($file, "Strand");
			my $nbSBScoding  = 0;

			open(F1, $file) or die "$!: $file\n";
			my $header = <F1>;
			while(<F1>)
			{
				$_      =~ s/[\r\n]+$//;
				my @tab = split("\t", $_);

				if($tab[$strand_value] ne "NA")
				{
					$nbSBScoding++;
				}
			}
			close F1;

			if($nbSBScoding != 0)
			{
				$timerFile++;
				`cp $file $folderCheckedForStat/`;
			}
			else
			{
				print STDOUT "\n\nWarning: There is no variant to compute statistics for $file\n\n";
				push(@listRemovedFile, $file);
			}
		}
	}

	if($timerFile == 0)
	{
		print STDERR "\n\nError: No variants to compute statistics for:\n";

		foreach (@listRemovedFile)
		{
			print STDERR $_."\n";
		}
		exit;
	}
}

# Retrieve chromosomes length
sub checkChrDir
{
	my @files = `ls $path_SeqrefGenome/$refGenome"_seq"/*.fa`;
	foreach my $file (@files)
	{
		if ($file !~ /chr(\d+|x|y)\.fa/i){next;}
		open(FILE,$file);
		<FILE>;
		my $seq="";
		while (<FILE>){ chomp; $seq.=$_;}
		$file =~ /chr(.*)\.fa/;
		$chromosomes{"chr".$1}=length($seq);
	}
}

# Calculate the statistics and generate the report
sub ReportMutDist
{
	print STDOUT "-----------------------------------------------------------------\n";
	print STDOUT "-----------------Report Mutational Analysis----------------------\n";
	print STDOUT "-----------------------------------------------------------------\n";

	my $folderFigure = "$folderMutAnalysis/Figures";
	if(-e $folderFigure) { rmtree($folderFigure); mkdir($folderFigure) or die "Can't create the directory $folderFigure\n"; }
	else { mkdir($folderFigure) or die "Can't create the directory $folderFigure\n"; }
	my $folderChi2 = "$folderFigure/Chi2";
	if(!-e $folderChi2) { mkdir($folderChi2) or die "Can't create the directory $folderChi2\n"; }
	my $folderWebLogo = "$folderFigure/WebLogo";
	if(!-e $folderWebLogo) { mkdir($folderWebLogo) or die "Can't create the directory $folderWebLogo\n"; }
	my $folderNMF    = "$folderFigure/Input_NMF";
	if(!-e $folderNMF) { mkdir($folderNMF) or die "Can't create the directory $folderNMF\n"; }


	################################################################################################
	###																	Calculates all the statistics														 ###
	################################################################################################

	########### Recover the functional region for all the samples. Allows to thave the same annotations for the pie chart "Impact on protein sequence"
	my @tab_func = recoverAnnovarAnnotation($func_name);
	if(@tab_func == 0)
	{
		print STDERR "Error: the table for the functional region is empty!!!!! check $folderCheckedForStat\n$func_name\n";
		exit;
	}

	############ Calculate the different statistics present in the report
	my %h_file   = ();
	CalculateStatistics(\%h_file, \@tab_func);

	############ Calculate the chi2 for the strand bias
	CalculateChi2(\%h_file, $folderChi2);

	############ Write the different statistics present in the report
	WriteStatistics(\%h_file, $#tab_func, $folderFigure, $folderChi2, $folderNMF);

	############ Create logo for studying the 10 flanking bases of the mutation
	CreateLogo(\%h_file, $folderWebLogo);
}


# Calculate the different statistics present in the report
sub CalculateStatistics
{
	my ($refH_file, $refT_func) = @_;

	our ($chr_value, $start_value, $ref_value, $alt_value, $func_value, $exonicFunc_value, $strand_value, $contextSeq_value) = ("", "", "", "", "", "", "", "", "", "");

	# Generate the pool of all the data
	if($poolData == 1)
	{
		my @listFile = `ls $folderCheckedForStat`;

		# For keeping the header only one time
		chomp($listFile[0]);
		system("cp $folderCheckedForStat/$listFile[0] $folderCheckedForStat/Pool_Data.txt");

		open(OUT, ">>", "$folderCheckedForStat/Pool_Data.txt") or die "$!: $folderCheckedForStat/Pool_Data.txt\n";

		for(my $i=1; $i<=$#listFile; $i++)
		{
			chomp($listFile[$i]);
			open(F1, "$folderCheckedForStat/$listFile[$i]") or die "$!: $folderCheckedForStat/$listFile[$i]\n";
			my $header = <F1>;
			while(<F1>) { print OUT $_; }
			close F1;
		}
		close OUT;
	}

	foreach my $file (`ls $folderCheckedForStat/*`)
	{
		chomp($file);
		############ Recover the number of the columns of interest
		$chr_value        = recoverNumCol($file, $chr_name);
		$start_value      = recoverNumCol($file, $start_name);
		$ref_value        = recoverNumCol($file, $ref_name);
		$alt_value        = recoverNumCol($file, $alt_name);
		$func_value       = recoverNumCol($file, $func_name);
		$exonicFunc_value = recoverNumCol($file, $exonicFunc_name);
		$strand_value     = recoverNumCol($file, $strand_name);
		$contextSeq_value = recoverNumCol($file, $context_name);
		############ Recover the number of the columns of interest

		############ Calculate the statistics for each file
		File2Hash($file, $func_value, $exonicFunc_value, $chr_value, $ref_value, $alt_value, $strand_value, $contextSeq_value, $refH_file, $refT_func);
	}
}

# Calculate the chi2 for the strand bias
sub CalculateChi2
{
	my ($refH_file, $folderChi2) = @_;

	# No value for the chi2
	if(scalar (keys %{$refH_file}) == 0)
	{
		print STDERR "Error: No value for calculating the chi2 for the strand bias\n";
		exit;
	}

	# Strand bias for one mutation type for all the samples
	my %h_tempchi2 = ();
	my ($ca_NonTr, $ca_Tr, $cg_NonTr, $cg_Tr, $ct_NonTr, $ct_Tr, $ta_NonTr, $ta_Tr, $tc_NonTr, $tc_Tr, $tg_NonTr, $tg_Tr) = (0,0,0,0,0,0, 0,0,0,0,0,0);

	my $nb_file = 0;

	foreach my $k_file (sort keys %{$refH_file})
	{
		$nb_file++;
		foreach my $k_func (sort keys %{$refH_file->{$k_file}{'6mutType'}})
		{
			foreach my $k_mutation (sort keys %{$refH_file->{$k_file}{'6mutType'}{$k_func}})
			{
				if($k_mutation eq "C:G>A:T")
				{
					$h_tempchi2{'C>A'}{$k_file}{'NonTr'} += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'};
					$h_tempchi2{'C>A'}{$k_file}{'Tr'}    += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'};
				}
				if($k_mutation eq "C:G>G:C")
				{
					$h_tempchi2{'C>G'}{$k_file}{'NonTr'} += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'};
					$h_tempchi2{'C>G'}{$k_file}{'Tr'}    += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'};
				}
				if($k_mutation eq "C:G>T:A")
				{
					$h_tempchi2{'C>T'}{$k_file}{'NonTr'} += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'};
					$h_tempchi2{'C>T'}{$k_file}{'Tr'}    += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'};
				}
				if($k_mutation eq "T:A>A:T")
				{
					$h_tempchi2{'T>A'}{$k_file}{'NonTr'} += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'};
					$h_tempchi2{'T>A'}{$k_file}{'Tr'}    += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'};
				}
				if($k_mutation eq "T:A>C:G")
				{
					$h_tempchi2{'T>C'}{$k_file}{'NonTr'} += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'};
					$h_tempchi2{'T>C'}{$k_file}{'Tr'}    += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'};
				}
				if($k_mutation eq "T:A>G:C")
				{
					$h_tempchi2{'T>G'}{$k_file}{'NonTr'} += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'};
					$h_tempchi2{'T>G'}{$k_file}{'Tr'}    += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'};
				}
			}
		}
	}

	# Create the input file for NMF
	open(CHI2, ">", "$folderChi2/Input_chi2_strandBias.txt") or die "$!: $folderChi2/Input_chi2_strandBias.txt\n";
	print CHI2 "SampleName\tNonTr\tTr\tAlteration\n";

	foreach my $k_mutation (sort keys %h_tempchi2)
	{
		foreach my $k_file (sort keys %{$h_tempchi2{$k_mutation}})
		{
			print CHI2 "$k_file\t$h_tempchi2{$k_mutation}{$k_file}{'NonTr'}\t$h_tempchi2{$k_mutation}{$k_file}{'Tr'}\t$k_mutation\n";
		}
	}
	close CHI2;


	`Rscript $pathRscriptChi2test --folderChi2 $folderChi2 2>&1`;
	# `Rscript $pathRscriptChi2test $folderChi2 2>&1`;


	if(!-e "$folderChi2/Output_chi2_strandBias.txt")
	{
		print STDERR "Error: Chi2 test didn't work !!!\n";
		exit;
	}
}

# Write the different statistics in the report
sub WriteStatistics
{
	my ($refH_file, $nb_func, $folderFigure, $folderChi2, $folderNMF) = @_;

	# Save the different graphs in specific folde
	if(!-e "$folderFigure/Overall_mutation_distribution") { mkdir("$folderFigure/Overall_mutation_distribution") or die "Can't create the directory $folderFigure/Overall_mutation_distribution\n"; }
	if(!-e "$folderFigure/Impact_protein_sequence") { mkdir("$folderFigure/Impact_protein_sequence") or die "Can't create the directory $folderFigure/Impact_protein_sequence\n"; }
	if(!-e "$folderFigure/SBS_distribution") { mkdir("$folderFigure/SBS_distribution") or die "Can't create the directory $folderFigure/SBS_distribution\n"; }
	if(!-e "$folderFigure/Stranded_Analysis") { mkdir("$folderFigure/Stranded_Analysis") or die "Can't create the directory $folderFigure/Stranded_Analysis\n"; }
	if(!-e "$folderFigure/Trinucleotide_Sequence_Context") { mkdir("$folderFigure/Trinucleotide_Sequence_Context") or die "Can't create the directory $folderFigure/Trinucleotide_Sequence_Context\n"; }
	if(!-e "$folderFigure/Distribution_SBS_Per_Chromosomes") { mkdir("$folderFigure/Distribution_SBS_Per_Chromosomes") or die "Can't create the directory $folderFigure/Distribution_SBS_Per_Chromosomes\n"; }


	# Create a workbook with all the samples
	my $wb = ""; my $ws_sum = "";
	my ($ws_inputNMF_count, $ws_inputNMF_percent) = ("", "");

	# Create one Excel file with all the samples
	if($oneReportPerSample == 2)
	{
		$wb = Spreadsheet::WriteExcel->new("$folderMutAnalysis/Report_Mutation_Spectra.xls");

		############## Set the variables for font formats in the Excel report
		Format_A10($wb, \$format_A10); # Text center in Arial 10
		Format_A10BoldLeft($wb, \$format_A10Boldleft); # Text on the left in Arial 10 bold
		Format_TextSection($wb, \$formatT_left, \$formatT_right, \$formatT_bottomRight, \$formatT_bottomLeft, \$formatT_bottom, \$formatT_bottomHeader, \$formatT_bottomRightHeader, \$formatT_bottomHeader2, \$formatT_rightHeader);
		Format_GraphTitle($wb, \$formatT_graphTitle);
		Format_Table($wb, \$table_topleft, \$table_topRight, \$table_bottomleft, \$table_bottomRight, \$table_top, \$table_right, \$table_bottom, \$table_bottomItalicRed, \$table_left, \$table_bottomrightHeader, \$table_left2, \$table_middleHeader, \$table_middleHeader2);
		Format_A10ItalicRed($wb, \$format_A10ItalicRed);

		############### Worksheet with a summary of the samples
		$ws_sum = $wb->add_worksheet("Sample_List");
		$ws_sum->write(0, 0, "Samples", $format_A10); $ws_sum->write(0, 1, "Total number SBS", $format_A10); $ws_sum->write(0, 2, "Total number of Indel", $format_A10); $ws_sum->write(0, 3, "Total number of mutations", $format_A10);
		$ws_sum->set_column(0,0, 50); $ws_sum->set_column(1,1, 20); $ws_sum->set_column(2,2, 20); $ws_sum->set_column(3,3, 22);

		############### Write the input matrix for NMF for the count and the un-normalized frequency
		$ws_inputNMF_count   = $wb->add_worksheet("Input_NMF_Count");
		$ws_inputNMF_percent = $wb->add_worksheet("Input_NMF_Percent");
	}


	################################################ Set the Rows and columns of the different part of the report
	my $row_SumSheet              = 1; # First row for the summary sheet of the report
	my $rowStart_SBSdistrBySeg    = 48; my $colStart_SBSdistrBySeg  = 0; # For the table SBS distribution by segment
	my $colStart_matrixSeqContext = 19; # Sequence context
	my $col_inputNMF              = 0;  # Write the names of the samples with at least 33 SBS


	## For each file
	foreach my $k_file (sort keys %{$refH_file})
	{
		print "File in process: $k_file\n";

		# Count the total of mutations for 6 mutation types on genomic strand
		my ($c_ca6_g, $c_cg6_g, $c_ct6_g, $c_ta6_g, $c_tc6_g, $c_tg6_g) = (0,0,0, 0,0,0);

		if($k_file ne "Pool_Data") { $col_inputNMF++; }

		############### Save the chi2 values into a hash table
		if(-e "$folderChi2/Output_chi2_strandBias.txt")
		{
			chi2hash("$folderChi2/Output_chi2_strandBias.txt", $k_file);
		}

		# Create one workbook for each sample
		if($oneReportPerSample == 1)
		{
			$wb = Spreadsheet::WriteExcel->new("$folderMutAnalysis/Report_Mutation_Spectra-$k_file.xls");

			############## Set the variables for font formats in the Excel report
			Format_A10($wb, \$format_A10); # Text center in Arial 10
			Format_A10BoldLeft($wb, \$format_A10Boldleft); # Text on the left in Arial 10 bold
			Format_TextSection($wb, \$formatT_left, \$formatT_right, \$formatT_bottomRight, \$formatT_bottomLeft, \$formatT_bottom, \$formatT_bottomHeader, \$formatT_bottomRightHeader, \$formatT_bottomHeader2, \$formatT_rightHeader);
			Format_GraphTitle($wb, \$formatT_graphTitle);
			Format_Table($wb, \$table_topleft, \$table_topRight, \$table_bottomleft, \$table_bottomRight, \$table_top, \$table_right, \$table_bottom, \$table_bottomItalicRed, \$table_left, \$table_bottomrightHeader, \$table_left2, \$table_middleHeader, \$table_middleHeader2);
			Format_A10ItalicRed($wb, \$format_A10ItalicRed);

		  ############### Worksheet with a summary of the samples
			$ws_sum = $wb->add_worksheet("Sample_List");
			$ws_sum->write(0, 0, "Samples", $format_A10); $ws_sum->write(0, 1, "Total number SBS", $format_A10); $ws_sum->write(0, 2, "Total number of Indel", $format_A10); $ws_sum->write(0, 3, "Total number of mutations", $format_A10);
			$ws_sum->set_column(0,0, 50); $ws_sum->set_column(1,1, 20); $ws_sum->set_column(2,2, 20); $ws_sum->set_column(3,3, 22);
			# Write in the Samples sheet the name and the total number of SBS
			$ws_sum->write(1, 0, "$k_file", $format_A10);
			$ws_sum->write(1, 1, $refH_file->{$k_file}{'TotalSBSGenomic'}, $format_A10); $ws_sum->write(1, 2, $refH_file->{$k_file}{'TotalIndelGenomic'}, $format_A10); $ws_sum->write($row_SumSheet, 3, $refH_file->{$k_file}{'TotalMutGenomic'}, $format_A10);
		}
		# One workbook with all the samples
		else
		{
			# Write in the Samples sheet the name and the total number of SBS
			$ws_sum->write($row_SumSheet, 0, $k_file, $format_A10);
			$ws_sum->write($row_SumSheet, 1, $refH_file->{$k_file}{'TotalSBSGenomic'}, $format_A10); $ws_sum->write($row_SumSheet, 2, $refH_file->{$k_file}{'TotalIndelGenomic'}, $format_A10); $ws_sum->write($row_SumSheet, 3, $refH_file->{$k_file}{'TotalMutGenomic'}, $format_A10);

			# For NMF don't consider the pool of the samples
			if($k_file ne "Pool_Data")
			{
				# Write in the input NMF sheet the name of the samples
				$ws_inputNMF_count->write(0, $col_inputNMF, $k_file);
				$ws_inputNMF_percent->write(0, $col_inputNMF, $k_file);
			}
		}

		# Calculate the correlation between the number of SBS and the size of the chromosome
		PearsonCoefficient($refH_file, $k_file);

		# Add a worksheet to the workbook
  	my $ws = $wb->add_worksheet($k_file);

  	# Write the titles of the different sections of the report
  	WriteBorderSection($wb, $ws, $rowStart_SBSdistrBySeg, $colStart_SBSdistrBySeg, $nb_func, $colStart_matrixSeqContext);

  	# Write the mutation types (6 types)
  	WriteHeaderSection($wb, $ws, $rowStart_SBSdistrBySeg, $colStart_SBSdistrBySeg, $nb_func, $colStart_matrixSeqContext);


  	# Save the figures of each samples in a different folder
		if(!-e "$folderFigure/Overall_mutation_distribution/$k_file") { mkdir("$folderFigure/Overall_mutation_distribution/$k_file") or die "Can't create the directory $folderFigure/Overall_mutation_distribution/$k_file\n"; }
		if(!-e "$folderFigure/Impact_protein_sequence/$k_file") { mkdir("$folderFigure/Impact_protein_sequence/$k_file") or die "Can't create the directory $folderFigure/Impact_protein_sequence/$k_file\n"; }
		if(!-e "$folderFigure/SBS_distribution/$k_file") { mkdir("$folderFigure/SBS_distribution/$k_file") or die "Can't create the directory $folderFigure/SBS_distribution\n"; }
		if(!-e "$folderFigure/Stranded_Analysis/$k_file") { mkdir("$folderFigure/Stranded_Analysis/$k_file") or die "Can't create the directory $folderFigure/Stranded_Analysis/$k_file\n"; }
		if(!-e "$folderFigure/Trinucleotide_Sequence_Context/$k_file") { mkdir("$folderFigure/Trinucleotide_Sequence_Context/$k_file") or die "Can't create the directory $folderFigure/Trinucleotide_Sequence_Context/$k_file\n"; }



  	##################################################################################################################################################
  	#################################################################	Write the statistics	##########################################################
  	##################################################################################################################################################
		my $row_SBSDistrBySegAndFunc_CG = $rowStart_SBSdistrBySeg+($nb_func*2)+16;


		######## Count of SBS by functional impact on the protein (Table 2) + Create the input for ggplot2 (pie chart with functional impact) + Create the input for ggplot2 (pie chart of SBS vs. Indels)
		writeDistrFuncImpact($ws, $refH_file, $k_file, "$folderFigure/Impact_protein_sequence/$k_file/$k_file-DistributionExoFunc.txt", "$folderFigure/Overall_mutation_distribution/$k_file/$k_file-OverallMutationDistribution.txt");


		######## Result of the chi2 for the strand bias (Table 3) + Create the input for ggplot2 (Strand bias bar graph)
		writeChi2result($wb, $ws, "$folderFigure/Stranded_Analysis/$k_file/$k_file-StrandBias.txt", $refH_file, $k_file);


		######## SBS distribution by functional region (Table 4) + Strand bias by functional region (Table 5) + Create the input for ggplot2 (SBS distribution) + Overall count and percent of SBS (Table 1)
		writeStatbyFuncRegion($refH_file, $k_file, $ws, $rowStart_SBSdistrBySeg, $colStart_SBSdistrBySeg, $nb_func, \$row_SBSDistrBySegAndFunc_CG, "$folderFigure/SBS_distribution/$k_file/$k_file-SBS_distribution.txt");


		######## Distribution of SBS per chromosomes and the result of Pearson test (Table 6)
		writeDistrByChr($ws, $refH_file, $k_file, $row_SBSDistrBySegAndFunc_CG, $colStart_SBSdistrBySeg, "$folderFigure/Distribution_SBS_Per_Chromosomes/$k_file-DistributionSNVS_per_chromosome.txt");


		######## Trinucleotide sequence context on genomic strand (Panel 1)
		# Represent the trinucleotide with a heatmap with count of SBS
		my $heatmapCountggplot2   = "$folderFigure/Trinucleotide_Sequence_Context/$k_file/$k_file-HeatmapCount-Genomic.txt";
		my $heatmapPercentggplot2 = "$folderFigure/Trinucleotide_Sequence_Context/$k_file/$k_file-HeatmapPercent-Genomic.txt";
		my $triNtBarChartggplot2  = "$folderFigure/Trinucleotide_Sequence_Context/$k_file/$k_file-MutationSpectraPercent-Genomic.txt";

		writeTriNtGenomic($ws, $refH_file, $k_file, $colStart_matrixSeqContext, $heatmapCountggplot2, $heatmapPercentggplot2, $triNtBarChartggplot2, \$c_ca6_g, \$c_cg6_g, \$c_ct6_g, \$c_ta6_g, \$c_tc6_g, \$c_tg6_g);

		# For the input matrix for NMF
		if($k_file ne "Pool_Data") { push(@{$h_inputNMF{'Sample'}}, $k_file); }


		######## Trinucleotide sequence context on genomic strand (Panel 2)
		my $triNtBarChartCodingCountggplot2   = "$folderFigure/Stranded_Analysis/$k_file/$k_file-StrandedSignatureCount.txt";
		my $triNtBarChartCodingPercentggplot2 = "$folderFigure/Stranded_Analysis/$k_file/$k_file-StrandedSignaturePercent.txt";

		writeTriNtCoding($ws, $rowStart_SBSdistrBySeg, $colStart_matrixSeqContext, $refH_file, $k_file, $triNtBarChartCodingCountggplot2, $triNtBarChartCodingPercentggplot2);


		######## Generate the figures and include them in the report
		createWriteFigs($ws, $rowStart_SBSdistrBySeg, $colStart_matrixSeqContext, $folderFigure, $k_file, $c_ca6_g, $c_cg6_g, $c_ct6_g, $c_ta6_g, $c_tc6_g, $c_tg6_g);


		# Next sample
		$row_SumSheet++;
	} # End $k_file

	######## Write the input matrix for NMF
	# One workbook with all the samples
	writeInputNMF($ws_inputNMF_count, $ws_inputNMF_percent, "$folderNMF/Input_NMF_Count.txt", "$folderNMF/Input_NMF_Frequency.txt");


	# Close the workbook
	$wb->close();
}

# Create logo for representing the sequence context with n bases
sub CreateLogo
{
	my ($refH_file, $folderWebLogo) = @_;

	my $folderSample = "";

	foreach my $k_file (sort keys %{$refH_file})
	{
		$folderSample      = "$folderWebLogo/$k_file";
		if(!-e $folderSample) { mkdir($folderSample) or die "Can't create the directory $folderSample\n"; }

		my $test_lengthSeqContext = 0;


		foreach my $k_mutation (sort keys %{$refH_file->{$k_file}{'WebLogo3'}})
		{
			$k_mutation     =~ /(\w)>(\w)/;
			my ($ref, $alt) = ($1, $2);

			open(WEBLOGO, ">", "$folderSample/$k_file-$ref$alt.fa") or die "$!: $folderSample/$k_file-$ref$alt.fa\n";
			foreach (@{$refH_file->{$k_file}{'WebLogo3'}{$k_mutation}})
			{
				print WEBLOGO ">$k_file\n$_\n";

				if(length($_) < 10) { $test_lengthSeqContext = 0; }
				else { $test_lengthSeqContext = 1; }
			}
			close WEBLOGO;
		}

		## Generate the logo
		foreach my $fastaFile (`ls $folderSample/*.fa`)
		{
			chomp($fastaFile);
			my ($filename, $directories, $suffix) = fileparse("$folderSample/$fastaFile", qr/\.[^.]*/);

			$filename =~ /(.+)\-/;
			my $title = $1;

			## Test if there is fasta sequences for the mutation type
			my $nbLigne_temp = `wc -l $fastaFile`;
			my @nbLigne = split(" ", $nbLigne_temp);


			if($nbLigne[0] == 0) { print "WARNING: No sequence for $filename\n"; next; }

			# When length sequence context is lower than 10 the image is to small for adding a title
			if($test_lengthSeqContext == 1) { system("weblogo -c classic -F png_print -U probability --title $title < $fastaFile > $folderSample/$filename-Probability.png"); }
			else { system("weblogo -c classic -F png_print -U probability < $fastaFile > $folderSample/$filename-Probability.png"); }
		}
	}
}



### Save the count of SBS for each file into a hash table
sub File2Hash
{
	my ($inputFile, $func_value, $exonicFunc_value, $chr_value, $ref_value, $alt_value, $strand_value, $contextSeq_value, $refH_file, $refT_func) = @_;

	my ($filename, $directories, $suffix) = fileparse($inputFile, qr/\.[^.]*/);

	# Initialisation of the hash
	my @tab_mutation   = qw(C:G>A:T C:G>G:C C:G>T:A T:A>A:T T:A>C:G T:A>G:C);
	my @tab_aaChange   = ("NonTr", "Tr", "TotalMutG");
	my @tabExoFunc     = ("frameshift insertion", "frameshift deletion", "frameshift block substitution", "frameshift substitution", "stopgain", "stoploss", "nonframeshift insertion", "nonframeshift deletion", "nonframeshift substitution", "nonframeshift block substitution", "nonsynonymous SNV", "synonymous SNV", "unknown", "NA");

	# Total number of SBS on the genomic strand
	$refH_file->{$filename}{'TotalSBSGenomic'} = 0;
	# Total number of Indel on the genomic strand
	$refH_file->{$filename}{'TotalIndelGenomic'} = 0;
	# Total number of SBS on the coding strand
	$refH_file->{$filename}{'TotalSBSCoding'} = 0;
	# Total number of SBS and Indel on the genomic strand
	$refH_file->{$filename}{'TotalMutGenomic'} = 0;


	#####################################################
	#		Initialisation of the tables and hash tables		#
	#####################################################

	## SBS by segment (6 mutation types)
	foreach my $elt_tabFunc (@$refT_func)
	{
	  foreach my $elt_tabMutation (@tab_mutation)
	  {
	   foreach my $elt_aaChange (@tab_aaChange)
	   {
	      $refH_file->{$filename}{'6mutType'}{$elt_tabFunc}{$elt_tabMutation}{$elt_aaChange} = 0;
	    }
	  }
	}

	## Pearson correlation
	$refH_file->{$filename}{'SBSPerChr'}{'AllMutType'} = 0;
	# Count of SBS per chromosome foreach mutation types
	foreach my $elt_tabMutation (@tab_mutation)
	{
		foreach my $chromosome (sort keys %chromosomes){ $refH_file->{$filename}{'SBSPerChr'}{$elt_tabMutation}{'CHR'}{$chromosome}{'chr'} = 0;}
 	 	$refH_file->{$filename}{'SBSPerChr'}{$elt_tabMutation}{'Pearson'} = 0;
	}
	foreach my $chromosome (sort keys %chromosomes)
	{
		$refH_file->{$filename}{'SBSPerChr'}{'TotalPerChr'}{$chromosome}{'chr'}=0;
	}

	## Impact of SBS on protein
	foreach my $elt_exoFunc (@tabExoFunc)
	{
	 	$refH_file->{$filename}{'ImpactSBS'}{$elt_exoFunc} = 0;
	}

	## Sequence context (genomic strand)
	my @tab_mutation2 = qw(C>A C>G C>T T>A T>C T>G);
	my @tab_context   = qw(A_A A_C A_G A_T C_A C_C C_G C_T G_A G_C G_G G_T T_A T_C T_G T_T);
	foreach my $elt_context (@tab_context)
	{
		foreach my $elt_mutation3 (@tab_mutation2)
	 	{
	 		$refH_file->{$filename}{'SeqContextG'}{$elt_context}{$elt_mutation3} = 0;
	 	}
	}

	## Sequence context (coding strand)
	my @tab_TrNonTr = qw(NonTr Tr);
	foreach my $elt_context (@tab_context)
	{
		foreach my $elt_mutation2 (@tab_mutation2)
	 	{
	  	foreach my $trNonTr (@tab_TrNonTr)
	  	{
	  		$refH_file->{$filename}{'SeqContextC'}{$elt_context}{$elt_mutation2}{$trNonTr} = 0;
	  	}
	  }
	}


	#####################################################
	#								Parse the intput file								#
	#####################################################

	open(F1,$inputFile) or die "$!: $inputFile\n";
	my $header = <F1>;
 	while(<F1>)
	{
		$_      =~ s/[\r\n]+$//;
		my @tab = split("\t", $_);

		### Don't consider random chromosomes and chromosome M
		if( ($tab[$chr_value] =~ /random/i) || ($tab[$chr_value] =~ /M/i) ) { next; }


		### Recover the trinucleotide sequence context: Extract the base just before and after the mutation
		my $context                   = "";
		my $contextSequence           = $tab[$contextSeq_value]; $contextSequence =~ tr/a-z/A-Z/;
		my @tempContextSequence       = split("", $contextSequence);
		my $total_nbBaseContext       = $#tempContextSequence;
		my $midlle_totalNbBaseContext = $total_nbBaseContext/2; # For having the middle of the sequence
		my $before                    = $midlle_totalNbBaseContext - 1; my $after = $midlle_totalNbBaseContext + 1;
		$context                      = $tempContextSequence[$before]."_".$tempContextSequence[$after];



		### Recover the annotations on the impact on the protein for creating the pie chart
	  my $exoFunc = "";
	  # Sometimes the annotation is repeated frameshift deletion;frameshift deletion
		if($tab[$exonicFunc_value] =~ /\;/)
		{
			my @temp = split(";", $tab[$exonicFunc_value]);
			if($temp[0] eq $temp[1]) { $exoFunc = $temp[0]; }
		}
		# The annotations have changed after MAJ Annovar 2014Jul22 (stopgain SNV => stopgain)
		elsif($tab[$exonicFunc_value] eq "stopgain SNV")      { $exoFunc = "stopgain"; }
		elsif($tab[$exonicFunc_value] eq "stoploss SNV")      { $exoFunc = "stoploss"; }
		elsif($tab[$exonicFunc_value] eq "nonsynonymous_SNV") { $exoFunc = "nonsynonymous SNV"; }
		elsif($tab[$exonicFunc_value] eq "stopgain_SNV")      { $exoFunc = "stopgain SNV"; }
		elsif($tab[$exonicFunc_value] eq "synonymous_SNV")    { $exoFunc = "synonymous SNV"; }
		else { $exoFunc = $tab[$exonicFunc_value]; }

		if(exists $refH_file->{$filename}{'ImpactSBS'}{$exoFunc})
		{
			# If the sequence context if not recovered correctly don't considered the variants
			if( ($context =~ /N/) || (length($context) != 3) ) { next; }

			$refH_file->{$filename}{'ImpactSBS'}{$exoFunc}++;
			$refH_file->{$filename}{'TotalMutGenomic'}++;
		}
		else { print "WARNING: Exonic function not considered: $exoFunc\n"; }

		#### Only SBS are considered for the statistics
		if( ($tab[$ref_value] =~ /^[ACGT]$/i) && ($tab[$alt_value] =~ /^[ACGT]$/i) )
		{
			# If the sequence context if not recovered correctly don't considered the variants
			if( ($context =~ /N/) || (length($context) != 3) ) { next; }

			# Total number of SBS on the genomic strand
			$refH_file->{$filename}{'TotalSBSGenomic'}++;

			# Total number of SBS on the coding strand with a sequence context
			if( ($tab[$strand_value] eq "+") || ($tab[$strand_value] eq "-") )
			{
				if( ($context ne "NA") && (($context =~ /N/) || (length($context) != 3)) ) { next; }
				$refH_file->{$filename}{'TotalSBSCoding'}++;
			}
		}
		else { $refH_file->{$filename}{'TotalIndelGenomic'}++; }

		### Number of SBS per chromosome: remove the "chr"
		my $chrNameForH=$tab[$chr_value];
		if(exists $refH_file->{$filename}{'SBSPerChr'}{'TotalPerChr'}{$chrNameForH}{'chr'})
		{
			$refH_file->{$filename}{'SBSPerChr'}{'TotalPerChr'}{$chrNameForH}{'chr'}++;
		}


		#### Some func value are repeated and separated by ";"
		my $funcSegment = "";
		if($tab[$func_value] =~ /;/) { my @temp = split(";", $tab[$func_value]); $funcSegment = $temp[0]; }
		else { $funcSegment = $tab[$func_value]; }



		#####################################################
		#		Calculate the statistics for each mutation type #
		#####################################################
	  if( (($tab[$ref_value] eq "C") && ($tab[$alt_value] eq "A")) || ( ($tab[$ref_value] eq "G") && ($tab[$alt_value] eq "T") ) )
	  {
	    my $mutation   = "C:G>A:T";

	    statPerMutType($filename, \@tab, $ref_value, $alt_value, \@tempContextSequence, $before, $after, $context, $funcSegment, $mutation, $refH_file, $chrNameForH, $strand_value, $midlle_totalNbBaseContext);
	  }
	  if( (($tab[$ref_value] eq "C") && ($tab[$alt_value] eq "G")) || ( ($tab[$ref_value] eq "G") && ($tab[$alt_value] eq "C") ) )
	  {
	    my $mutation   = "C:G>G:C";

	    statPerMutType($filename, \@tab, $ref_value, $alt_value, \@tempContextSequence, $before, $after, $context, $funcSegment, $mutation, $refH_file, $chrNameForH, $strand_value, $midlle_totalNbBaseContext);
	  }
	  if( (($tab[$ref_value] eq "C") && ($tab[$alt_value] eq "T")) || ( ($tab[$ref_value] eq "G") && ($tab[$alt_value] eq "A") ) )
	  {
	  	my $mutation   = "C:G>T:A";
	  	statPerMutType($filename, \@tab, $ref_value, $alt_value, \@tempContextSequence, $before, $after, $context, $funcSegment, $mutation, $refH_file, $chrNameForH, $strand_value, $midlle_totalNbBaseContext);
	  }
	  if( (($tab[$ref_value] eq "T") && ($tab[$alt_value] eq "A")) || ( ($tab[$ref_value] eq "A") && ($tab[$alt_value] eq "T") ) )
	  {
	  	my $mutation   = "T:A>A:T";
	  	statPerMutType($filename, \@tab, $ref_value, $alt_value, \@tempContextSequence, $before, $after, $context, $funcSegment, $mutation, $refH_file, $chrNameForH, $strand_value, $midlle_totalNbBaseContext);
	  }
	  if( (($tab[$ref_value] eq "T") && ($tab[$alt_value] eq "C")) || ( ($tab[$ref_value] eq "A") && ($tab[$alt_value] eq "G")) )
	  {
	  	my $mutation   = "T:A>C:G";
	  	statPerMutType($filename, \@tab, $ref_value, $alt_value, \@tempContextSequence, $before, $after, $context, $funcSegment, $mutation, $refH_file, $chrNameForH, $strand_value, $midlle_totalNbBaseContext);
	  }
	  if( (($tab[$ref_value] eq "T") && ($tab[$alt_value] eq "G")) || ( ($tab[$ref_value] eq "A") && ($tab[$alt_value] eq "C")) )
	  {
	  	my $mutation   = "T:A>G:C";
	  	statPerMutType($filename, \@tab, $ref_value, $alt_value, \@tempContextSequence, $before, $after, $context, $funcSegment, $mutation, $refH_file, $chrNameForH, $strand_value, $midlle_totalNbBaseContext);
	  }
	}
	close F1;
}

### Count the number of SBS for 12 and 6 categories
sub statPerMutType
{
	my ($filename, $refTab, $ref_value, $alt_value, $refTab_tempSeqContext, $before, $after, $context, $funcSegment, $mutation, $refH_file, $chrNameForH, $strand_value, $midlle_totalNbBaseContext) = @_;

	my @tab = @$refTab;
	my @tempContextSequence = @$refTab_tempSeqContext;


	# Split the mutations
	$mutation =~ /(\w)\:(\w)\>(\w)\:(\w)/;
	my ($ref1, $ref2, $alt1, $alt2) = ($1, $2, $3, $4);

	# Count the total number of mutations
	$refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'TotalMutG'}++;

	# Pearson correlation
	if(exists $refH_file->{$filename}{'SBSPerChr'}{$mutation}{'CHR'}{$chrNameForH}{'chr'})
	{
		$refH_file->{$filename}{'SBSPerChr'}{$mutation}{'CHR'}{$chrNameForH}{'chr'}++;
	}

	#### Sequence context - 6 mutation types - genomic strand
	my $mutationSeqContext6mutType = "$ref1>$alt1";
	# We want to express the mutation in C> or T>
	if( ($tab[$ref_value] eq $ref2) && ($tab[$alt_value] eq $alt2) )
	{
		my $base3 = complement($tempContextSequence[$before]);
		my $base5 = complement($tempContextSequence[$after]);
		my $context_reverse = $base5."_".$base3;
		if(exists $refH_file->{$filename}{'SeqContextG'}{$context_reverse}{$mutationSeqContext6mutType})
		{
			$refH_file->{$filename}{'SeqContextG'}{$context_reverse}{$mutationSeqContext6mutType}++;
		}
	}
	elsif(exists $refH_file->{$filename}{'SeqContextG'}{$context}{$mutationSeqContext6mutType})
	{
		$refH_file->{$filename}{'SeqContextG'}{$context}{$mutationSeqContext6mutType}++;
	}


	#### Strand analysis C>N and T>N on NonTr strand
	if( (($tab[$strand_value] eq "+") && (($tab[$ref_value] eq $ref1)&&($tab[$alt_value] eq $alt1))) || (($tab[$strand_value] eq "-") && (($tab[$ref_value] eq $ref2)&&($tab[$alt_value] eq $alt2))) )
	{
		if(exists $refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'NonTr'})
		{
			$refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'NonTr'}++;
		}

		# C>A With the sequence context (C>N and T>N on strand +)
		if( ($tab[$strand_value] eq "+") && (($tab[$ref_value] eq $ref1)&&($tab[$alt_value] eq $alt1)) )
		{
			if(exists $refH_file->{$filename}{'SeqContextC'}{$context}{$mutationSeqContext6mutType}{'NonTr'})
		  {
		  	$refH_file->{$filename}{'SeqContextC'}{$context}{$mutationSeqContext6mutType}{'NonTr'}++;
		  }
		}
		# C>A With the sequence context (G>N and A>N on strand -)
		else
		{
			my $base3           = complement($tempContextSequence[$before]);
		 	my $base5           = complement($tempContextSequence[$after]);
		 	my $context_reverse = $base5."_".$base3;

			if(exists $refH_file->{$filename}{'SeqContextC'}{$context_reverse}{$mutationSeqContext6mutType}{'NonTr'})
			{
		    $refH_file->{$filename}{'SeqContextC'}{$context_reverse}{$mutationSeqContext6mutType}{'NonTr'}++;
		  }
		}
	}

	#### Strand analysis C>N and T>N on Tr strand
	if( (($tab[$strand_value] eq "-") && (($tab[$ref_value] eq $ref1)&&($tab[$alt_value] eq $alt1))) || (($tab[$strand_value] eq "+") && (($tab[$ref_value] eq $ref2)&&($tab[$alt_value] eq $alt2))) )
	{
		if(exists $refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'Tr'})
		{
			$refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'Tr'}++;
		}

		# C>N and T>N With the sequence context (strand -)
		if( ($tab[$strand_value] eq "-") && (($tab[$ref_value] eq $ref1)&&($tab[$alt_value] eq $alt1)) )
		{
			if(exists $refH_file->{$filename}{'SeqContextC'}{$context}{$mutationSeqContext6mutType}{'Tr'})
		  {
		    $refH_file->{$filename}{'SeqContextC'}{$context}{$mutationSeqContext6mutType}{'Tr'}++;
		 	}
		}
		# C>N and T>N with the sequence context (strand +)
		if( ($tab[$strand_value] eq "+") && (($tab[$ref_value] eq $ref2)&&($tab[$alt_value] eq $alt2)) )
		{
			my $base3           = complement($tempContextSequence[$before]);
			my $base5           = complement($tempContextSequence[$after]);
			my $context_reverse = $base5."_".$base3;
		  if(exists $refH_file->{$filename}{'SeqContextC'}{$context_reverse}{$mutationSeqContext6mutType}{'Tr'})
		  {
		    $refH_file->{$filename}{'SeqContextC'}{$context_reverse}{$mutationSeqContext6mutType}{'Tr'}++;
		  }
		}
	}

	#### WebLogo-3
	if(($tab[$ref_value] eq $ref1) && ($tab[$alt_value] eq $alt1))
	{
		# For the logo all the sequences must have the same length
		if(scalar(@tempContextSequence) == 2)  { next; }
		my ($contextTemp1, $contextTemp2) = ("", "");
		for(my $i=0; $i<$midlle_totalNbBaseContext; $i++) { $contextTemp1 .= $tempContextSequence[$i]; }
		for(my $i=$midlle_totalNbBaseContext+1; $i<=$#tempContextSequence; $i++) { $contextTemp2 .= $tempContextSequence[$i]; }
		my $context = $contextTemp1.$ref1.$contextTemp2;
		push(@{$refH_file->{$filename}{'WebLogo3'}{$mutationSeqContext6mutType}}, $context);
	}
	else
	{

		if(scalar(@tempContextSequence) == 2)  { next; }
		my ($contextTemp1, $contextTemp2) = ("", "");
		for(my $i=0; $i<$midlle_totalNbBaseContext; $i++) { $contextTemp1 .= complement($tempContextSequence[$i]); }
		for(my $i=$midlle_totalNbBaseContext+1; $i<=$#tempContextSequence; $i++) { $contextTemp2 .= complement($tempContextSequence[$i]); }
		my $context = $contextTemp1.$ref1.$contextTemp2; $context = reverse $context;
		push(@{$refH_file->{$filename}{'WebLogo3'}{$mutationSeqContext6mutType}}, $context);
	}
}

# Calculate the correlation between the number of SBS and the size of the chromosome
sub PearsonCoefficient
{
	our ($refH_file, $filename) = @_;

	#### Calculate the Pearson coefficient
	my @total_SBS = (); # Pearson for all mutation types

	# Create a 2D array
	foreach my $k_mutation (sort keys %{$refH_file->{$filename}{'SBSPerChr'}})
	{
		my $x           = [];
		my $correlation = 0;

		if($k_mutation eq "AllMutType") { next; }
		elsif($k_mutation eq "TotalPerChr") { next; }
		elsif($k_mutation eq "ChrSize") { next; }
		else
		{
			my $testZero = 0; # The correlation function doesn't works if all the variables are equal to zero
			# generate an anonymous 2D array where $x->[1] is the row
			# $x->[1][1] is the value in row 1 column 1 and $x->[1][2] is the value of row 1 column 2
			# once you build the entire array, pass it to the correlation subroutine
			my $i=1;
			while ( my ($chromosome, $lenght) = each (%chromosomes))
			{
				$x->[$i][1] = $lenght; # First column contains the chromosome size
				$x->[$i][2] = $refH_file->{$filename}{'SBSPerChr'}{$k_mutation}{'CHR'}{$chromosome}{'chr'}; # Second column contains the count of SBS
				if($refH_file->{$filename}{'SBSPerChr'}{$k_mutation}{'CHR'}{$chromosome}{'chr'}==0) { $testZero++; }
				$i++;
			}
			if( $testZero == keys %{$refH_file->{$filename}{'SBSPerChr'}{$k_mutation}{'CHR'}} )
			{
				$correlation = 0;
			}
			# Pass the 2D array to the correlation subroutine
			else
			{
				$correlation = correlation($x);
			}

			$refH_file->{$filename}{'SBSPerChr'}{$k_mutation}{'Pearson'} = $correlation; # Pearson per mutation type
		}
	}

	#generate an anonymous 2D array for all mutation type
	my $testZero    = 0;
	my $x           = [];
	my $correlation = 0;
	my $i=1;
	while ( my ($chromosome, $lenght) = each (%chromosomes))
	{
		$x->[$i][1] = $lenght;
		$x->[$i][2] = $refH_file->{$filename}{'SBSPerChr'}{'TotalPerChr'}{$chromosome}{'chr'};
		$i++;
	}
	if($testZero == keys %{$refH_file->{$filename}{'SBSPerChr'}{'TotalPerChr'}} ) { $correlation = 0; }
	else { $correlation = correlation($x); }
	# Pass the 2D array to the correlation subroutine
	$refH_file->{$filename}{'SBSPerChr'}{'AllMutType'} = $correlation;

	sub correlation
	{
		my ($x) = @_;
		my ($mean_x,$mean_y) = mean($x);
		my $ssxx=ss($x,$mean_x,$mean_y,1,1);
		my $ssyy=ss($x,$mean_x,$mean_y,2,2);
		my $ssxy=ss($x,$mean_x,$mean_y,1,2);
		my $correl=correl($ssxx,$ssyy,$ssxy);;
		my $xcorrel=sprintf("%.2f",$correl);
		return($xcorrel);

		sub mean
		{
			my ($x)=@_;
			my $num = scalar(@{$x}) - 2;
			my $sum_x = '0';
			my $sum_y = '0';
			for (my $i = 2; $i < scalar(@{$x}); ++$i)
			{
			  $sum_x += $x->[$i][1];
			  $sum_y += $x->[$i][2];
			}
			 my $mu_x = $sum_x / $num;
			 my $mu_y = $sum_y / $num;
			 return($mu_x,$mu_y);
		}

		### ss = sum of squared (deviations to the mean)
		sub ss
		{
			my ($x,$mean_x,$mean_y,$one,$two)=@_;
			my $sum = '0';
			for (my $i=2;$i<scalar(@{$x});++$i)
			{
			  $sum += ($x->[$i][$one]-$mean_x)*($x->[$i][$two]-$mean_y);
			}
		 	return $sum;
		}

		sub correl
		{
			my($ssxx,$ssyy,$ssxy)=@_;

			my ($sign, $correl) = (0,0);
			if(abs($ssxy) == 0) { $sign = 0 }
			else { $sign=$ssxy/abs($ssxy); }

			 if( ($ssxx==0) || ($ssyy==0) ) { $correl = 0 }
			 else { $correl=$sign*sqrt($ssxy*$ssxy/($ssxx*$ssyy)); }

			 return $correl;
		}
	}
}

# Save the output of the chi2 into a hash table for writing the results in the Excel file
sub chi2hash
{
	my ($outputChi2, $k_file) = @_;

	open(F1, $outputChi2) or die "$!: $outputChi2\n";
	my $header = <F1>;
	# Strand_Bias($tab[0])	NonTr-Tr($tab[1])	Proportion($tab[2])	P-val-Chi2($tab[3])	FDR($tab[4])	Confidence Interval($tab[5])	Mutation_Type($tab[6])	SampleName($tab[7])
	while(<F1>)
	{
		$_      =~ s/[\r\n]+$//;
		my @tab = split("\t", $_);

		if($tab[7] eq $k_file)
		{
			if($tab[1] eq "NA")
			{
				$h_chi2{$tab[7]}{$tab[6]}{'NonTr'} = "NA";
				$h_chi2{$tab[7]}{$tab[6]}{'Tr'}    = "NA";
			}
			else
			{
				my ($nonTr, $tr) = split("-", $tab[1]);

				$h_chi2{$tab[7]}{$tab[6]}{'NonTr'} = $nonTr;
				$h_chi2{$tab[7]}{$tab[6]}{'Tr'}    = $tr;
			}


			$h_chi2{$tab[7]}{$tab[6]}{'p-value'} = $tab[3];
			$h_chi2{$tab[7]}{$tab[6]}{'ConfInt'} = $tab[5];

			# For the pool data the FDR isn't calculated so replace the NA (=Missing values) by "-"
			if($tab[7] eq "Pool_Data")
			{
				$h_chi2{$tab[7]}{$tab[6]}{'FDR'} = "-";
			}
			else
			{
				$h_chi2{$tab[7]}{$tab[6]}{'FDR'} = $tab[4];
			}
		}
	}
	close F1;
}

### Complement bases (for the sequence context)
sub complement
{
	if($_[0] eq "A") { return "T"; }
	if($_[0] eq "C") { return "G"; }
	if($_[0] eq "G") { return "C"; }
	if($_[0] eq "T") { return "A"; }
}

### Recover the functional region for all the samples. Allows to thave the same annotations for the pie chart "Impact on protein sequence"
sub recoverAnnovarAnnotation
{
	my ($func_name) = @_;

	my %hash = ();

	# The input is a folder
	foreach my $file (`ls $folderCheckedForStat/*`)
	{
		chomp($file);
		my $AV_annotation_value = recoverNumCol($file, $func_name);

		open(F1, $file) or die "$!: $file\n";
		my $header = <F1>;
		while(<F1>)
		{
			$_      =~ s/[\r\n]+$//;
			my @tab = split("\t", $_);

			# Some files can have an empty line at the end and WE DON'T WANT to consider it
			if(! defined $tab[0]) { next; }
			# Some func value are repeated and separated by ";"
			my $funcSegment = "";
			if($tab[$AV_annotation_value] =~ /;/)
			{
				my @temp = split(";", $tab[$AV_annotation_value]);
				$funcSegment = $temp[0];
			}
			else { $funcSegment = $tab[$AV_annotation_value]; }

			$hash{$funcSegment} = "";
		}
		close F1;
	}
	my @tab_AVAnnotation = ();
	foreach my $k (sort keys %hash) { push(@tab_AVAnnotation, $k); }
	return @tab_AVAnnotation;
}

sub recoverNumCol
{
	my ($input, $name_of_column) = @_;

  open(F1,$input) or die "recoverNumCol: $!: $input\n";
  # For having the name of the columns
  my $search_header = <F1>; $search_header =~ s/[\r\n]+$//; my @tab_search_header = split("\t",$search_header);
  close F1;
  # The number of the column
  my $name_of_column_NB  = "toto";
  for(my $i=0; $i<=$#tab_search_header; $i++)
  {
    if($tab_search_header[$i] eq $name_of_column) { $name_of_column_NB = $i; last; }
  }
  if($name_of_column_NB eq "toto") { print STDERR "Error recoverNumCol(): the column named $name_of_column doesn't exits in the input file $input!!!!!\n"; exit 3; }
  else                             { return $name_of_column_NB; }
}



######################################################################################################################################################
#																												Functions for writing in the Excel report																										 #
######################################################################################################################################################
# Write the header for the six mutation types
sub WriteHeaderSection
{
	our ($wb, $ws, $rowStart_SBSdistrBySeg, $colStart_SBSdistrBySeg, $nb_func, $colStart_matrixSeqContext) = @_;

	our ($format_CA, $format_CG, $format_CT, $format_TA, $format_TC, $format_TG, $format_TG2, $format_LeftHeader, $format_RightHeader, $format_LeftHeader2);
	Format_Header($wb, \$format_CA, \$format_CG, \$format_CT, \$format_TA, \$format_TC, \$format_TG, \$format_TG2, \$format_LeftHeader, \$format_RightHeader, \$format_LeftHeader2);

	our ($format_LeftCA, $format_LeftCG, $format_LeftCT, $format_LeftTA, $format_LeftTC, $format_LeftTG, $format_RightCA, $format_RightCG, $format_RightCT, $format_RightTA, $format_RightTC, $format_RightTG);
	Format_HeaderSBSDistrBySegAndFunc($wb, \$format_LeftCA, \$format_LeftCG, \$format_LeftCT, \$format_LeftTA, \$format_LeftTC, \$format_LeftTG, \$format_RightCA, \$format_RightCG, \$format_RightCT, \$format_RightTA, \$format_RightTC, \$format_RightTG);

	our $format_A11Bold      = ""; Format_A11Bold($wb, \$format_A11Bold); # Arial 11 bold and center
	our $format_A11BoldLeft  = ""; Format_A11BoldLeft($wb, \$format_A11BoldLeft); # Arial 11 bold and left

	our ($format_header12CA, $format_header12CG, $format_header12CT, $format_header12TA, $format_header12TC, $format_header12TG);
  Format_Header12MutType($wb, \$format_header12CA, \$format_header12CG, \$format_header12CT, \$format_header12TA, \$format_header12TC, \$format_header12TG);

	## Header for SBS distribution by segment
	HeaderMutTypeSBSDistrBySeg();

	## Header for strand bias by function
	$ws->set_column($colStart_SBSdistrBySeg+5, $colStart_SBSdistrBySeg+5, 11);

	my $row = $rowStart_SBSdistrBySeg+$nb_func+10; my $col = $colStart_SBSdistrBySeg;
	$ws->write($row, $col+1, ' ', $format_CA); $ws->write($row, $col+2, "C>A", $format_CA); $ws->write($row, $col+3, ' ', $format_CA);
	$ws->write($row, $col+5, ' ', $format_CG); $ws->write($row, $col+6, "C>G", $format_CG); $ws->write($row, $col+7, ' ', $format_CG);
	$ws->write($row, $col+9, ' ', $format_CT); $ws->write($row, $col+10, "C>T", $format_CT); $ws->write($row, $col+11, ' ', $format_RightCT);

	$row = $rowStart_SBSdistrBySeg+($nb_func*2)+14;
	$ws->write($row, $col+1, ' ', $format_TA); $ws->write($row, $col+2, "T>A", $format_TA); $ws->write($row, $col+3, ' ', $format_TA);
	$ws->write($row, $col+5, ' ', $format_TC); $ws->write($row, $col+6, "T>C", $format_TC); $ws->write($row, $col+7, ' ', $format_TC);
	$ws->write($row, $col+9, ' ', $format_TG2); $ws->write($row, $col+10, "T>G", $format_TG2); $ws->write($row, $col+11, ' ', $format_RightTG);

	$ws->set_row($rowStart_SBSdistrBySeg+$nb_func+11, 18); $ws->set_row($rowStart_SBSdistrBySeg+($nb_func*2)+15, 18);
	$ws->set_column($colStart_SBSdistrBySeg+5, $colStart_SBSdistrBySeg+5, 13); $ws->set_column($colStart_SBSdistrBySeg+9, $colStart_SBSdistrBySeg+9, 13);

	for(my $i=$rowStart_SBSdistrBySeg+$nb_func+10; $i<=$rowStart_SBSdistrBySeg+($nb_func*2)+14; $i+=$nb_func+4)
	{
		$ws->write($i+1, $colStart_SBSdistrBySeg, 'Segment', $format_LeftHeader); $ws -> write($i+1, $colStart_SBSdistrBySeg+1, 'Non-Tr/Tr', $format_A11Bold); $ws -> write($i+1, $colStart_SBSdistrBySeg+2, 'Non-Tr', $format_A11Bold); $ws -> write($i+1, $colStart_SBSdistrBySeg+3, 'Tr', $format_A11Bold);
		$ws -> write($i+1, $colStart_SBSdistrBySeg+5, 'Non-Tr/Tr', $format_A11Bold); $ws -> write($i+1, $colStart_SBSdistrBySeg+6, 'Non-Tr', $format_A11Bold); $ws -> write($i+1, $colStart_SBSdistrBySeg+7, 'Tr', $format_A11Bold);
		$ws -> write($i+1, $colStart_SBSdistrBySeg+9, 'Non-Tr/Tr', $format_A11Bold); $ws -> write($i+1, $colStart_SBSdistrBySeg+10, 'Non-Tr', $format_A11Bold); $ws -> write($i+1, $colStart_SBSdistrBySeg+11, 'Tr', $format_RightHeader);
	}


	## Header for Counts of SBS per chromosome and mutation type
	HeaderCountSBSPerChr();

	## Header for the short sequence context
	HeaderShortTriNtContext();

	## Header for the 12 mutation types with the sequence context (coding strand)
	HeaderLongTriNtContext();

	sub HeaderMutTypeSBSDistrBySeg
	{
		$ws->set_row($rowStart_SBSdistrBySeg+2, 18);
		$ws->write($rowStart_SBSdistrBySeg+2, $colStart_SBSdistrBySeg+2, "C:G>A:T", $format_CA); $ws->write_blank($rowStart_SBSdistrBySeg+2, $colStart_SBSdistrBySeg+3, $format_CA);
		$ws->write($rowStart_SBSdistrBySeg+2, $colStart_SBSdistrBySeg+4, "C:G>G:C", $format_CG); $ws->write_blank($rowStart_SBSdistrBySeg+2, $colStart_SBSdistrBySeg+5, $format_CG);
		$ws->write($rowStart_SBSdistrBySeg+2, $colStart_SBSdistrBySeg+6, "C:G>T:A", $format_CT); $ws->write_blank($rowStart_SBSdistrBySeg+2, $colStart_SBSdistrBySeg+7, $format_CT);
			$ws->write($rowStart_SBSdistrBySeg+2, $colStart_SBSdistrBySeg+8, "T:A>A:T", $format_TA); $ws->write_blank($rowStart_SBSdistrBySeg+2, $colStart_SBSdistrBySeg+9, $format_TA);
		$ws->write($rowStart_SBSdistrBySeg+2, $colStart_SBSdistrBySeg+10, "T:A>C:G", $format_TC); $ws->write_blank($rowStart_SBSdistrBySeg+2, $colStart_SBSdistrBySeg+11, $format_TC);
		$ws->write($rowStart_SBSdistrBySeg+2, $colStart_SBSdistrBySeg+12, "T:A>G:C", $format_TG); $ws->write_blank($rowStart_SBSdistrBySeg+2, $colStart_SBSdistrBySeg+13, $format_TG);

		$ws->write($rowStart_SBSdistrBySeg+3, $colStart_SBSdistrBySeg, "Segment", $format_LeftHeader); $ws->set_column($colStart_SBSdistrBySeg, $colStart_SBSdistrBySeg, 13); $ws->set_row($rowStart_SBSdistrBySeg+3, 18);
		$ws->write($rowStart_SBSdistrBySeg+3, $colStart_SBSdistrBySeg+1, "Total SBS", $format_A11Bold); $ws->set_column($colStart_SBSdistrBySeg+1, $colStart_SBSdistrBySeg+1, 11);
		$ws->write($rowStart_SBSdistrBySeg+3, $colStart_SBSdistrBySeg+2, "%", $format_A11Bold); $ws->write($rowStart_SBSdistrBySeg+3, $colStart_SBSdistrBySeg+3, "#", $format_A11Bold);
		$ws->write($rowStart_SBSdistrBySeg+3, $colStart_SBSdistrBySeg+4, "%", $format_A11Bold); $ws->write($rowStart_SBSdistrBySeg+3, $colStart_SBSdistrBySeg+5, "#", $format_A11Bold);
		$ws->write($rowStart_SBSdistrBySeg+3, $colStart_SBSdistrBySeg+6, "%", $format_A11Bold); $ws->write($rowStart_SBSdistrBySeg+3, $colStart_SBSdistrBySeg+7, "#", $format_A11Bold);
		$ws->write($rowStart_SBSdistrBySeg+3, $colStart_SBSdistrBySeg+8, "%", $format_A11Bold); $ws->write($rowStart_SBSdistrBySeg+3, $colStart_SBSdistrBySeg+9, "#", $format_A11Bold);
		$ws->write($rowStart_SBSdistrBySeg+3, $colStart_SBSdistrBySeg+10, "%", $format_A11Bold); $ws->write($rowStart_SBSdistrBySeg+3, $colStart_SBSdistrBySeg+11, "#", $format_A11Bold);
		$ws->write($rowStart_SBSdistrBySeg+3, $colStart_SBSdistrBySeg+12, "%", $format_A11Bold); $ws->write($rowStart_SBSdistrBySeg+3, 13, "#", $format_RightHeader);
	}

	sub HeaderCountSBSPerChr
	{
		$ws->set_column(3,3, 10); $ws->set_column(4,4, 10);
		$ws->set_row($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+8, 18);
		$ws->write($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+7, $colStart_SBSdistrBySeg+1, "Pearson", $format_A11Bold);
		$ws->write($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+8, $colStart_SBSdistrBySeg, "Chr", $format_LeftHeader);
		$ws->write($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+8, $colStart_SBSdistrBySeg+1, "Size", $format_A11Bold);
		$ws->write($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+8, $colStart_SBSdistrBySeg+2, "All SBS", $format_A11Bold);

		$ws->write($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+8, $colStart_SBSdistrBySeg+3, "C:G>A:T", $format_CA);
		$ws->write($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+8, $colStart_SBSdistrBySeg+4, "C:G>G:C", $format_CG);
		$ws->write($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+8, $colStart_SBSdistrBySeg+5, "C:G>T:A", $format_CT);
		$ws->write($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+8, $colStart_SBSdistrBySeg+6, "T:A>A:T", $format_TA);
		$ws->write($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+8, $colStart_SBSdistrBySeg+7, "T:A>C:G", $format_TC);
		$ws->write($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+8, $colStart_SBSdistrBySeg+8, "T:A>G:C", $format_TG);
	}

	sub HeaderShortTriNtContext
	{
		### GENOMIC STRAND
		$ws->write(2, $colStart_matrixSeqContext, 'Count matrix', $format_LeftHeader2);
		$ws->write(3, $colStart_matrixSeqContext+4, 'C>A', $format_CA); $ws->write(3, $colStart_matrixSeqContext+5, 'C>G', $format_CG); $ws->write(3, $colStart_matrixSeqContext+6, 'C>T', $format_CT); $ws->write(3, $colStart_matrixSeqContext+7, 'T>A', $format_TA); $ws->write(3, $colStart_matrixSeqContext+8, 'T>C', $format_TC); $ws->write(3, $colStart_matrixSeqContext+9, 'T>G', $format_TG2);

		$ws->write(2, $colStart_matrixSeqContext+11, 'Frequency matrix', $format_A11BoldLeft);
		$ws->write(3, $colStart_matrixSeqContext+14, 'C>A', $format_CA); $ws->write(3, $colStart_matrixSeqContext+15, 'C>G', $format_CG); $ws->write(3, $colStart_matrixSeqContext+16, 'C>T', $format_CT); $ws->write(3, $colStart_matrixSeqContext+17, 'T>A', $format_TA); $ws->write(3, $colStart_matrixSeqContext+18, 'T>C', $format_TC); $ws->write(3, $colStart_matrixSeqContext+19, 'T>G', $format_TG2);

			### sequence context with a bar graph
			$ws->write(25, $colStart_matrixSeqContext+10, "Mutation spectra frequency", $format_A11Bold);
		}

	sub HeaderLongTriNtContext
	{
		$ws->set_row($rowStart_SBSdistrBySeg+3, 15); $ws->set_row($rowStart_SBSdistrBySeg+4, 15); $ws->set_row($rowStart_SBSdistrBySeg+5, 15);
		$ws->write($rowStart_SBSdistrBySeg+3, $colStart_matrixSeqContext, "Count matrix", $format_LeftHeader2);
		$ws->write($rowStart_SBSdistrBySeg+4, $colStart_matrixSeqContext+1, "C>A", $format_CA); $ws->write_blank($rowStart_SBSdistrBySeg+4, $colStart_matrixSeqContext+2, $format_CA); $ws->write($rowStart_SBSdistrBySeg+5, $colStart_matrixSeqContext+1, "NonTr", $format_A11Bold); $ws->write($rowStart_SBSdistrBySeg+5, $colStart_matrixSeqContext+2, "Tr", $format_A11Bold);
		$ws->write($rowStart_SBSdistrBySeg+4, $colStart_matrixSeqContext+3, "C>G", $format_CG); $ws->write_blank($rowStart_SBSdistrBySeg+4, $colStart_matrixSeqContext+4, $format_CG); $ws->write($rowStart_SBSdistrBySeg+5, $colStart_matrixSeqContext+3, "NonTr", $format_A11Bold); $ws->write($rowStart_SBSdistrBySeg+5, $colStart_matrixSeqContext+4, "Tr", $format_A11Bold);
		$ws->write($rowStart_SBSdistrBySeg+4, $colStart_matrixSeqContext+5, "C>T", $format_CT); $ws->write_blank($rowStart_SBSdistrBySeg+4, $colStart_matrixSeqContext+6, $format_CT); $ws->write($rowStart_SBSdistrBySeg+5, $colStart_matrixSeqContext+5, "NonTr", $format_A11Bold); $ws->write($rowStart_SBSdistrBySeg+5, $colStart_matrixSeqContext+6, "Tr", $format_A11Bold);
		$ws->write($rowStart_SBSdistrBySeg+4, $colStart_matrixSeqContext+7, "T>A", $format_TA); $ws->write_blank($rowStart_SBSdistrBySeg+4, $colStart_matrixSeqContext+8, $format_TA); $ws->write($rowStart_SBSdistrBySeg+5, $colStart_matrixSeqContext+7, "NonTr", $format_A11Bold); $ws->write($rowStart_SBSdistrBySeg+5, $colStart_matrixSeqContext+8, "Tr", $format_A11Bold);
		$ws->write($rowStart_SBSdistrBySeg+4, $colStart_matrixSeqContext+9, "T>C", $format_TC); $ws->write_blank($rowStart_SBSdistrBySeg+4, $colStart_matrixSeqContext+10, $format_TC); $ws->write($rowStart_SBSdistrBySeg+5, $colStart_matrixSeqContext+9, "NonTr", $format_A11Bold); $ws->write($rowStart_SBSdistrBySeg+5, $colStart_matrixSeqContext+10, "Tr", $format_A11Bold);
		$ws->write($rowStart_SBSdistrBySeg+4, $colStart_matrixSeqContext+11, "T>G", $format_TG2); $ws->write_blank($rowStart_SBSdistrBySeg+4, $colStart_matrixSeqContext+12, $format_TG2); $ws->write($rowStart_SBSdistrBySeg+5, $colStart_matrixSeqContext+11, "NonTr", $format_A11Bold); $ws->write($rowStart_SBSdistrBySeg+5, $colStart_matrixSeqContext+12, "Tr", $format_A11Bold);


		$ws->set_row($rowStart_SBSdistrBySeg+24, 15); $ws->set_row($rowStart_SBSdistrBySeg+25, 15); $ws->set_row($rowStart_SBSdistrBySeg+26, 15);
		$ws->write($rowStart_SBSdistrBySeg+24, $colStart_matrixSeqContext, "Frequency matrix", $format_LeftHeader2);
		$ws->write($rowStart_SBSdistrBySeg+25, $colStart_matrixSeqContext+1, "C>A", $format_CA); $ws->write_blank($rowStart_SBSdistrBySeg+25, $colStart_matrixSeqContext+2, $format_CA); $ws->write($rowStart_SBSdistrBySeg+26, $colStart_matrixSeqContext+1, "NonTr", $format_A11Bold); $ws->write($rowStart_SBSdistrBySeg+26, $colStart_matrixSeqContext+2, "Tr", $format_A11Bold);
		$ws->write($rowStart_SBSdistrBySeg+25, $colStart_matrixSeqContext+3, "C>G", $format_CG); $ws->write_blank($rowStart_SBSdistrBySeg+25, $colStart_matrixSeqContext+4, $format_CG); $ws->write($rowStart_SBSdistrBySeg+26, $colStart_matrixSeqContext+3, "NonTr", $format_A11Bold); $ws->write($rowStart_SBSdistrBySeg+26, $colStart_matrixSeqContext+4, "Tr", $format_A11Bold);
		$ws->write($rowStart_SBSdistrBySeg+25, $colStart_matrixSeqContext+5, "C>T", $format_CT); $ws->write_blank($rowStart_SBSdistrBySeg+25, $colStart_matrixSeqContext+6, $format_CT); $ws->write($rowStart_SBSdistrBySeg+26, $colStart_matrixSeqContext+5, "NonTr", $format_A11Bold); $ws->write($rowStart_SBSdistrBySeg+26, $colStart_matrixSeqContext+6, "Tr", $format_A11Bold);
		$ws->write($rowStart_SBSdistrBySeg+25, $colStart_matrixSeqContext+7, "T>A", $format_TA); $ws->write_blank($rowStart_SBSdistrBySeg+25, $colStart_matrixSeqContext+8, $format_TA); $ws->write($rowStart_SBSdistrBySeg+26, $colStart_matrixSeqContext+7, "NonTr", $format_A11Bold); $ws->write($rowStart_SBSdistrBySeg+26, $colStart_matrixSeqContext+8, "Tr", $format_A11Bold);
		$ws->write($rowStart_SBSdistrBySeg+25, $colStart_matrixSeqContext+9, "T>C", $format_TC); $ws->write_blank($rowStart_SBSdistrBySeg+25, $colStart_matrixSeqContext+10, $format_TC); $ws->write($rowStart_SBSdistrBySeg+26, $colStart_matrixSeqContext+9, "NonTr", $format_A11Bold); $ws->write($rowStart_SBSdistrBySeg+26, $colStart_matrixSeqContext+10, "Tr", $format_A11Bold);
		$ws->write($rowStart_SBSdistrBySeg+25, $colStart_matrixSeqContext+11, "T>G", $format_TG2); $ws->write_blank($rowStart_SBSdistrBySeg+25, $colStart_matrixSeqContext+12, $format_TG2); $ws->write($rowStart_SBSdistrBySeg+26, $colStart_matrixSeqContext+11, "NonTr", $format_A11Bold); $ws->write($rowStart_SBSdistrBySeg+26, $colStart_matrixSeqContext+12, "Tr", $format_A11Bold);
	}
}
# Write the titles of the different sections of the report
sub WriteBorderSection
{
	our ($wb, $ws, $rowStart_SBSdistrBySeg, $colStart_SBSdistrBySeg, $nb_func, $colStart_matrixSeqContext) = @_;

	our ($format_topLeft, $format_topRight, $format_bottomLeft, $format_bottomRight, $format_top, $format_right, $format_bottom, $format_left);
	Format_section($wb, \$format_topLeft, \$format_topRight, \$format_bottomLeft, \$format_bottomRight, \$format_top, \$format_right, \$format_bottom, \$format_left);

	TableSBSDistrBySeg();
	TableStrandBiasBySegment();
	CountSBSPerChr();
	ShortTriNtContext(); # 6 mut type
	LongTriNtContext();  # 12 mut type

	sub TableSBSDistrBySeg
	{
		# Top-Left
		$ws->write($rowStart_SBSdistrBySeg, $colStart_SBSdistrBySeg, "Table 4. SBS distribution by functional region", $format_topLeft); $ws->set_row($rowStart_SBSdistrBySeg, 18); # Set the height of the row to 0.25"
		# Top
		for(my $i=1; $i<=13; $i++) { $ws->write_blank($rowStart_SBSdistrBySeg, $colStart_SBSdistrBySeg+$i, $format_top); }
		# Top-Right
		$ws->write_blank($rowStart_SBSdistrBySeg, $colStart_SBSdistrBySeg+13, $format_topRight);
		# Right
		$ws->write_blank($rowStart_SBSdistrBySeg+1, $colStart_SBSdistrBySeg+13, $format_right);
		# Bottom-left
		$ws->write_blank($rowStart_SBSdistrBySeg+$nb_func+5, $colStart_SBSdistrBySeg, $format_bottomLeft);
		# Left
		$ws->write_blank($rowStart_SBSdistrBySeg+1, $colStart_SBSdistrBySeg, $format_left); $ws->write_blank($rowStart_SBSdistrBySeg+2, $colStart_SBSdistrBySeg, $format_left);
	}

	sub TableStrandBiasBySegment
	{
		# Top-Left
		$ws->write($rowStart_SBSdistrBySeg+$nb_func+8, $colStart_SBSdistrBySeg, "Table 5. Strand bias by functional region", $format_topLeft); $ws->set_row($rowStart_SBSdistrBySeg+$nb_func+8, 18); # Set the height of the row to 0.25"
		# Top
		for(my $i=1; $i<=10; $i++) { $ws->write_blank($rowStart_SBSdistrBySeg+$nb_func+8, $colStart_SBSdistrBySeg+$i, $format_top); }
		# Top-Right
		$ws->write_blank($rowStart_SBSdistrBySeg+$nb_func+8, $colStart_SBSdistrBySeg+11, $format_topRight);
		# Right
		$ws->write_blank($rowStart_SBSdistrBySeg+$nb_func+9, $colStart_SBSdistrBySeg+11, $format_right); $ws->write_blank($rowStart_SBSdistrBySeg+($nb_func*2)+13, $colStart_SBSdistrBySeg+11, $format_right);
		# Left
		$ws->write_blank($rowStart_SBSdistrBySeg+$nb_func+9, $colStart_SBSdistrBySeg, $format_left); $ws->write_blank($rowStart_SBSdistrBySeg+$nb_func+10, $colStart_SBSdistrBySeg, $format_left); $ws->write_blank($rowStart_SBSdistrBySeg+($nb_func*2)+13, $colStart_SBSdistrBySeg, $format_left); $ws->write_blank($rowStart_SBSdistrBySeg+($nb_func*2)+14, $colStart_SBSdistrBySeg, $format_left);
		# Bottom
		$ws->write_blank($rowStart_SBSdistrBySeg+($nb_func*3)+16, $colStart_SBSdistrBySeg+4, $format_bottom); $ws->write_blank($rowStart_SBSdistrBySeg+($nb_func*3)+16, $colStart_SBSdistrBySeg+8, $format_bottom);
	}

	sub CountSBSPerChr
	{
		#### Top-Left
		$ws->write($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+4, $colStart_SBSdistrBySeg, "Table 6. SBS distribution per chromosome", $format_topLeft); $ws->set_row($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+4, 18); # Set the height of the row to 0.25"
		#### Top
		for(my $i=1; $i<8; $i++) { $ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+4, $colStart_SBSdistrBySeg+$i, $format_top); }
		#### Top-Right
		$ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+4, $colStart_SBSdistrBySeg+8, $format_topRight);
		#### Right
		$ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+5, $colStart_SBSdistrBySeg+8, $format_right); $ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+6, $colStart_SBSdistrBySeg+8, $format_right);

		#### Bottom-Right
		# Human genome = 24 chromosomes
		if($refGenome =~ /hg/)    { $ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+33, $colStart_SBSdistrBySeg+8, $format_bottomRight); }
		# Mouse genome = 21 chromosomes
		if($refGenome =~ /mm/) { $ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+30, $colStart_SBSdistrBySeg+8, $format_bottomRight); }
		# Rat genome = 22 chromosomes
		if($refGenome =~ /rn/) { $ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+31, $colStart_SBSdistrBySeg+8, $format_bottomRight); }

		#### Bottom
		if($refGenome =~ /hg/)
		{
			$ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+33, $colStart_SBSdistrBySeg+1, $format_bottom);
			for(my $i=3; $i<=7; $i++) { $ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+33, $colStart_SBSdistrBySeg+$i, $format_bottom); }
		}
		if($refGenome =~ /mm/)
		{
			$ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+30, $colStart_SBSdistrBySeg+1, $format_bottom);
			for(my $i=3; $i<=7; $i++) { $ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+30, $colStart_SBSdistrBySeg+$i, $format_bottom); }
		}
		if($refGenome =~ /rn/)
		{
			$ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+31, $colStart_SBSdistrBySeg+1, $format_bottom);
			for(my $i=3; $i<=7; $i++) { $ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+31, $colStart_SBSdistrBySeg+$i, $format_bottom); }
		}

		#### Left
		$ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+5, $colStart_SBSdistrBySeg, $format_left); $ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+6, $colStart_SBSdistrBySeg, $format_left); $ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+7, $colStart_SBSdistrBySeg, $format_left);

		#### Bottom-left
		if($refGenome =~ /hg/) { $ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+33, $colStart_SBSdistrBySeg, $format_bottomLeft);  }
		if($refGenome =~ /mm/) { $ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+30, $colStart_SBSdistrBySeg, $format_bottomLeft); }
		if($refGenome =~ /rn/) { $ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+31, $colStart_SBSdistrBySeg, $format_bottomLeft); }
	}

	sub ShortTriNtContext
	{
		my $format_headerSection = $wb->add_format(valign  => 'left', bold => 1, font => 'Arial', size => 12);
		$format_headerSection->set_left(2); $format_headerSection->set_left_color('blue');

		# Top-left
		$ws->write(0, $colStart_matrixSeqContext, "Panel 1. Trinucleotide sequence context of SBS on the genomic sequence", $format_topLeft);
		# Top
		for(my $i=1; $i<=19; $i++) { $ws->write_blank(0, $colStart_matrixSeqContext+$i, $format_top); }
		# Top-right
		$ws->write_blank(0, $colStart_matrixSeqContext+20, $format_topRight);
		# Right
		for(my $i=1; $i<=37; $i++)  { $ws->write_blank($i, $colStart_matrixSeqContext+20, $format_right); }
		# Bottom-right
		$ws->write_blank(37, $colStart_matrixSeqContext+20, $format_bottomRight);
		# Bottom
		for(my $i=1; $i<=19; $i++) { $ws->write_blank(38, $colStart_matrixSeqContext+$i, $format_top); }
		# Bottom-left
		$ws->write_blank(37, $colStart_matrixSeqContext, $format_bottomLeft);
		# Left
		$ws->write(1, $colStart_matrixSeqContext, "", $format_left);
		for(my $i=3; $i<=36; $i++)  { $ws->write_blank($i, $colStart_matrixSeqContext, $format_left); }
	}

	sub LongTriNtContext
	{
		# Top-left
		$ws->write($rowStart_SBSdistrBySeg, $colStart_matrixSeqContext, "Panel 2. Stranded analysis of trinucleotide sequence context of SBS", $format_topLeft);
		# Top
		for(my $i=1; $i<=28; $i++) { $ws->write_blank($rowStart_SBSdistrBySeg, $colStart_matrixSeqContext+$i, $format_top); }
		# Top-right
		$ws->write_blank($rowStart_SBSdistrBySeg, $colStart_matrixSeqContext+29, $format_topRight);
		# Right
		for(my $i=1; $i<=42; $i++)  { $ws->write_blank($rowStart_SBSdistrBySeg+$i, $colStart_matrixSeqContext+29, $format_right); }
		# Bottom-right
		$ws->write_blank(91, $colStart_matrixSeqContext+29, $format_bottomRight);
		# Bottom
		for(my $i=13; $i<=28; $i++) { $ws->write_blank(92, $colStart_matrixSeqContext+$i, $format_top); }
		# Bottom-left
		$ws->write_blank(91, $colStart_matrixSeqContext, $format_bottomLeft);
		# Left
		$ws->write_blank($rowStart_SBSdistrBySeg+1, $colStart_matrixSeqContext, $format_left); $ws->write_blank($rowStart_SBSdistrBySeg+2, $colStart_matrixSeqContext, $format_left); $ws->write_blank($rowStart_SBSdistrBySeg+4, $colStart_matrixSeqContext, $format_left); $ws->write_blank($rowStart_SBSdistrBySeg+5, $colStart_matrixSeqContext, $format_left);
		$ws->write_blank($rowStart_SBSdistrBySeg+22, $colStart_matrixSeqContext, $format_left); $ws->write_blank($rowStart_SBSdistrBySeg+23, $colStart_matrixSeqContext, $format_left); $ws->write_blank($rowStart_SBSdistrBySeg+25, $colStart_matrixSeqContext, $format_left); $ws->write_blank($rowStart_SBSdistrBySeg+26, $colStart_matrixSeqContext, $format_left);
	}
}

############################################################################################################
# Write count of SBS by functional impact on the protein (Table 2) + Create the input for ggplot2 (pie chart with functional impact) + Create the input for ggplot2 (pie chart of SBS vs. Indels)
sub writeDistrFuncImpact
{
	my ($ws, $refH_file, $sample, $funcImpactggplot2, $overallDistrggplot2) = @_;

	# Set the row for the table 2
	my $lImpactSBS = 31;
	my ($deletion, $insertion) = (0, 0);


	$ws->write(29, 6, "Table 2. Frequency and counts of functional impact", $format_A10Boldleft);
  $ws->set_column(6, 6, 13); $ws->set_column(10, 10, 15);
  $ws->write(30, 6, "RefSeq gene", $table_topleft);
  $ws->write(30, 7, "", $table_top);
  $ws->write(30, 8, "Percent", $table_top);
  $ws->write(30, 9, "Count", $table_topRight);



 	open(IMPACTSBS, ">", $funcImpactggplot2) or die "$!: $funcImpactggplot2\n";
  print IMPACTSBS "AA_Change\tCount\tPercent\n";

  # Pie chart with the distribution of SBS vs Indel
  open(SBSINDEL, ">", $overallDistrggplot2) or die "$!: $overallDistrggplot2\n";
  print SBSINDEL "Variant_type\tCount\tPercent\n";


  foreach my $k_exoFunc(sort keys %{$refH_file->{$sample}{'ImpactSBS'}})
  {
  	my $percent = ($refH_file->{$sample}{'ImpactSBS'}{$k_exoFunc} / $refH_file->{$sample}{'TotalMutGenomic'})*100;
  	$percent    = sprintf("%.2f", $percent);

  	if($k_exoFunc eq "NA")
  	{
  		print IMPACTSBS "Not_Applicable\t$percent\t$refH_file->{$sample}{'ImpactSBS'}{$k_exoFunc}\n";
  	}
  	else
  	{
  		my $temp = $k_exoFunc;
  		$temp =~ s/ /_/g;
  		print IMPACTSBS "$temp\t$percent\t$refH_file->{$sample}{'ImpactSBS'}{$k_exoFunc}\n";
  	}

  	$ws->write($lImpactSBS, 6, $k_exoFunc, $table_left2);
  	$ws->write($lImpactSBS, 8, $percent, $format_A10);
  	$ws->write($lImpactSBS, 9, $refH_file->{$sample}{'ImpactSBS'}{$k_exoFunc}, $table_right);

  	$lImpactSBS++;

  	# Pie chart with the distribution of SBS vs Indel
  	if($k_exoFunc    =~ /deletion/i)  { $deletion  += $refH_file->{$sample}{'ImpactSBS'}{$k_exoFunc}; }
		elsif($k_exoFunc =~ /insertion/i) { $insertion += $refH_file->{$sample}{'ImpactSBS'}{$k_exoFunc}; }
  }
  close IMPACTSBS;

  $ws->write($lImpactSBS, 9, $refH_file->{$sample}{'TotalMutGenomic'}, $table_bottomrightHeader);
  $ws->write($lImpactSBS, 6, "", $table_bottomleft); $ws->write($lImpactSBS, 7, "", $table_bottom); $ws->write($lImpactSBS, 8, "", $table_bottom);

  # Pie chart with the distribution of SBS vs Indel
  my $percentSBSIndel = ($deletion/$refH_file->{$sample}{'TotalMutGenomic'})*100; $percentSBSIndel = sprintf("%.2f", $percentSBSIndel);
  print SBSINDEL "Deletion\t$deletion\t$percentSBSIndel\n";
  $percentSBSIndel = ($insertion/$refH_file->{$sample}{'TotalMutGenomic'})*100; $percentSBSIndel = sprintf("%.2f", $percentSBSIndel);
  print SBSINDEL "Insertion\t$insertion\t$percentSBSIndel\n";
  $percentSBSIndel = ($refH_file->{$sample}{TotalSBSGenomic}/$refH_file->{$sample}{'TotalMutGenomic'})*100; $percentSBSIndel = sprintf("%.2f", $percentSBSIndel);
  print SBSINDEL "SBS\t$refH_file->{$sample}{TotalSBSGenomic}\t$percentSBSIndel\n";
  close SBSINDEL;
}


############################################################################################################
# Write the result of the chi2 for the strand bias (Table 3)
sub writeChi2result
{
	my ($wb, $ws, $strandBiasggplot2, $refH_file, $k_file) = @_;

	# Define the header of the table 3
	$ws->write(28, 11, "Table 3. Significance of the strand biases", $format_A10Boldleft);
	$ws->set_column(11, 11, 13); $ws->set_column(16, 16, 15); $ws->set_column(17, 17, 10);
	$ws->write(29, 11, "Mutation Type", $table_topleft); $ws->write(29, 12, "Non-Tr/Tr", $table_top); $ws->write(29, 13, "Non-Tr", $table_top); $ws->write(29, 14, "Tr", $table_top); $ws->write(29, 15, "P-value", $table_top); $ws->write(29, 16, "FDR q value", $table_top); $ws->write(29, 17, "95% CI", $table_topRight);
	$ws->write(39, 11, "Table 3. Significance of the strand biases", $format_A10Boldleft);
	$ws->write(40, 11, "Mutation Type", $table_topleft); $ws->write(40, 12, "Non-Tr/Tr", $table_top); $ws->write(40, 13, "Non-Tr", $table_top); $ws->write(40, 14, "Tr", $table_top); $ws->write(40, 15, "P-value", $table_top); $ws->write(40, 16, "FDR q value", $table_top); $ws->write(40, 17, "95% CI", $table_topRight);


	# Define the count on non-transcribed and transcribed strand
	# C>A
	my ($ca_NonTr, $ca_Tr) = (0,0);
	if( ($h_chi2{$k_file}{"C>A"}{'NonTr'} eq "NA") || ($h_chi2{$k_file}{"C>A"}{'Tr'} eq "NA") )
	{
		$ca_NonTr = 0;
		$ca_Tr    = 0;
	}
	else
	{
		$ca_NonTr = $h_chi2{$k_file}{"C>A"}{'NonTr'};
		$ca_Tr    = $h_chi2{$k_file}{"C>A"}{'Tr'};
	}
	# C>G
	my ($cg_NonTr, $cg_Tr) = (0,0);
	if( ($h_chi2{$k_file}{"C>G"}{'NonTr'} eq "NA") || ($h_chi2{$k_file}{"C>G"}{'Tr'} eq "NA") )
	{
		$cg_NonTr = 0;
		$cg_Tr    = 0;
	}
	else
	{
		$cg_NonTr = $h_chi2{$k_file}{"C>G"}{'NonTr'};
		$cg_Tr    = $h_chi2{$k_file}{"C>G"}{'Tr'};
	}
	# C>T
	my ($ct_NonTr, $ct_Tr) = (0,0);
	if( ($h_chi2{$k_file}{"C>T"}{'NonTr'} eq "NA") || ($h_chi2{$k_file}{"C>T"}{'Tr'} eq "NA") )
	{
		$ct_NonTr = 0;
		$ct_Tr    = 0;
	}
	else
	{
		$ct_NonTr = $h_chi2{$k_file}{"C>T"}{'NonTr'};
		$ct_Tr    = $h_chi2{$k_file}{"C>T"}{'Tr'};
	}
	# T>A
	my ($ta_NonTr, $ta_Tr) = (0,0);
	if( ($h_chi2{$k_file}{"T>A"}{'NonTr'} eq "NA") || ($h_chi2{$k_file}{"T>A"}{'Tr'} eq "NA") )
	{
		$ta_NonTr = 0;
		$ta_Tr    = 0;
	}
	else
	{
		$ta_NonTr = $h_chi2{$k_file}{"T>A"}{'NonTr'};
		$ta_Tr    = $h_chi2{$k_file}{"T>A"}{'Tr'};
	}
	# T>C
	my ($tc_NonTr, $tc_Tr) = (0,0);
	if( ($h_chi2{$k_file}{"T>C"}{'NonTr'} eq "NA") || ($h_chi2{$k_file}{"T>C"}{'Tr'} eq "NA") )
	{
		$tc_NonTr = 0;
		$tc_Tr    = 0;
	}
	else
	{
		$tc_NonTr = $h_chi2{$k_file}{"T>C"}{'NonTr'};
		$tc_Tr    = $h_chi2{$k_file}{"T>C"}{'Tr'};
	}
	# T>G
	my ($tg_NonTr, $tg_Tr) = (0,0);
	if( ($h_chi2{$k_file}{"T>G"}{'NonTr'} eq "NA") || ($h_chi2{$k_file}{"T>G"}{'Tr'} eq "NA") )
	{
		$tg_NonTr = 0;
		$tg_Tr    = 0;
	}
	else
	{
		$tg_NonTr = $h_chi2{$k_file}{"T>G"}{'NonTr'};
		$tg_Tr    = $h_chi2{$k_file}{"T>G"}{'Tr'};
	}



	# Create an input for representing the strand bias with ggplot2
	open(SB, ">", $strandBiasggplot2) or die "$!: $strandBiasggplot2\n";
	print SB "Alteration\tStrand\tCount\n";

	print SB "C>A\tNonTranscribed\t".$ca_NonTr."\n"."C>A\tTranscribed\t".$ca_Tr."\n";
	print SB "C>G\tNonTranscribed\t".$cg_NonTr."\n"."C>G\tTranscribed\t".$cg_Tr."\n";
	print SB "C>T\tNonTranscribed\t".$ct_NonTr."\n"."C>T\tTranscribed\t".$ct_Tr."\n";
	print SB "T>A\tNonTranscribed\t".$ta_NonTr."\n"."T>A\tTranscribed\t".$ta_Tr."\n";
	print SB "T>C\tNonTranscribed\t".$tc_NonTr."\n"."T>C\tTranscribed\t".$tc_Tr."\n";
	print SB "T>G\tNonTranscribed\t".$tg_NonTr."\n"."T>G\tTranscribed\t".$tg_Tr."\n";
	close SB;


	### Calcul the ratio for the strand bias and write it on the Excel file (Table 3)
	writeRatioSB($ca_NonTr, $ca_Tr, $ws, 30, $refH_file, $k_file, "C>A", "G>T", $format_A10, $table_left, $table_middleHeader, $table_right);
	writeRatioSB($cg_NonTr, $cg_Tr, $ws, 31, $refH_file, $k_file, "C>G", "G>C", $format_A10, $table_left, $table_middleHeader, $table_right);
	writeRatioSB($ct_NonTr, $ct_Tr, $ws, 32, $refH_file, $k_file, "C>T", "G>A", $format_A10, $table_left, $table_middleHeader, $table_right);
	writeRatioSB($ta_NonTr, $ta_Tr, $ws, 33, $refH_file, $k_file, "T>A", "A>T", $format_A10, $table_left, $table_middleHeader, $table_right);
	writeRatioSB($tc_NonTr, $tc_Tr, $ws, 34, $refH_file, $k_file, "T>C", "A>G", $format_A10, $table_left, $table_middleHeader, $table_right);
	writeRatioSB($tg_NonTr, $tg_Tr, $ws, 35, $refH_file, $k_file, "T>G", "A>C", $table_bottom, $table_bottomleft, $table_middleHeader2, $table_bottomRight);


	### Write a warning message when NonTr+Tr < 10
	my $format_italic_red = $wb->add_format(font=>'Arial', size=>10, italic=>1, color => 'red');

	if( (($ca_NonTr+$ca_Tr)< 10) || (($cg_NonTr+$cg_Tr)< 10) || (($ct_NonTr+$ct_Tr)< 10) || (($ta_NonTr+$ta_Tr)< 10) || (($tc_NonTr+$tc_Tr)< 10) || (($tg_NonTr+$tg_Tr)< 10) )
	{
		$ws->write(36, 11, "Warning message: chi-squared approximation may be incorrect because the number of SBS", $format_italic_red);
		$ws->write(37, 11, "on Non-transcribed and transcribed strand is lower than 10", $format_italic_red);
	}
}
# Write values in Table 3 (Sub function of writeChi2result)
sub writeRatioSB
{
	my ($count_NonTr, $count_Tr, $ws, $row, $refH_file, $k_file, $mut1, $mut2, $formatText, $formatTextMut1, $formatTextRatio, $formatTable) = @_;

	my ($ratio_mut1, $ratio_mut2, $percent_NonTr, $percent_Tr) = (0, 0, 0, 0, 0);
	if( ($count_NonTr==0) || ($count_Tr==0) ) { $ratio_mut1 = 0; $ratio_mut2 = 0; $percent_NonTr = 0; $percent_Tr = 0; }
	else
	{
		$ratio_mut1    = $count_NonTr/$count_Tr; $ratio_mut1 = sprintf("%.2f", $ratio_mut1);
		$ratio_mut2    = $count_Tr/$count_NonTr; $ratio_mut2 = sprintf("%.2f", $ratio_mut2);
		$percent_NonTr = ($count_NonTr/$refH_file->{$k_file}{'TotalSBSGenomic'})*100;
		$percent_Tr    = ($count_Tr/$refH_file->{$k_file}{'TotalSBSGenomic'})*100;
	}

	# C>N and T>N
	$ws->write($row, 11, $mut1, $formatTextMut1);
	$ws->write($row, 12, $ratio_mut1, $formatTextRatio);
	$ws->write($row, 13, $count_NonTr, $formatText);
	$ws->write($row, 14, $count_Tr, $formatText);
	# Write in italic and red (= warning message) when the count of NonTr + Tr is lower than 10
	if( ($count_NonTr+$count_Tr) < 10 )
	{
		if($mut1 eq "T>G")
		{
			if(! exists $h_chi2{$k_file}{$mut1}{'p-value'})   { $ws->write_string($row, 15, ""); }
			elsif($h_chi2{$k_file}{$mut1}{'p-value'} eq "NA") { $ws->write_string($row, 15, $h_chi2{$k_file}{$mut1}{'p-value'}, $table_bottom); }
			else { $ws->write_string($row, 15, $h_chi2{$k_file}{$mut1}{'p-value'}, $table_bottomItalicRed); }
		}
		else
		{
			if(! exists $h_chi2{$k_file}{$mut1}{'p-value'})   { $ws->write_string($row, 15, ""); }
			elsif($h_chi2{$k_file}{$mut1}{'p-value'} eq "NA") { $ws->write_string($row, 15, $h_chi2{$k_file}{$mut1}{'p-value'}, $formatText); }
			else { $ws->write_string($row, 15, $h_chi2{$k_file}{$mut1}{'p-value'}, $format_A10ItalicRed); }
		}
	}
	else
	{
		$ws->write_string($row, 15, $h_chi2{$k_file}{$mut1}{'p-value'}, $formatText);
	}

	$ws->write($row, 16, $h_chi2{$k_file}{$mut1}{'FDR'}, $formatText);
	$ws->write($row, 17, $h_chi2{$k_file}{$mut1}{'ConfInt'}, $formatTable);


	# G>N and A>N
	$ws->write($row+11, 11, $mut2, $formatTextMut1);
	$ws->write($row+11, 12, $ratio_mut2, $formatTextRatio);
	$ws->write($row+11, 13, $count_Tr, $formatText);
	$ws->write($row+11, 14, $count_NonTr, $formatText);
	if( ($count_NonTr+$count_Tr) < 10 )
	{
		if($mut1 eq "T>G")
		{
			if(! exists $h_chi2{$k_file}{$mut1}{'p-value'})   { $ws->write_string($row+11, 15, ""); }
			elsif($h_chi2{$k_file}{$mut1}{'p-value'} eq "NA") { $ws->write_string($row+11, 15, $h_chi2{$k_file}{$mut1}{'p-value'}, $table_bottom); }
			else { $ws->write_string($row+11, 15, $h_chi2{$k_file}{$mut1}{'p-value'}, $table_bottomItalicRed); }
		}
		else
		{
			if(! exists $h_chi2{$k_file}{$mut1}{'p-value'})   { $ws->write_string($row+11, 15, ""); }
			elsif($h_chi2{$k_file}{$mut1}{'p-value'} eq "NA") { $ws->write_string($row+11, 15, $h_chi2{$k_file}{$mut1}{'p-value'}, $formatText); }
			else { $ws->write_string($row+11, 15, $h_chi2{$k_file}{$mut1}{'p-value'}, $format_A10ItalicRed); }
		}
	}
	else
	{
		$ws->write_string($row+11, 15, $h_chi2{$k_file}{$mut1}{'p-value'}, $formatText);
	}

	$ws->write($row+11, 16, $h_chi2{$k_file}{$mut1}{'FDR'}, $formatText);
	$ws->write($row+11, 17, $h_chi2{$k_file}{$mut1}{'ConfInt'}, $formatTable);
}


############################################################################################################
# SBS distribution by functional region (Table 4) & Strand bias by functional region (Table 5) & Overall count and percent of SBS (Table 1)
sub writeStatbyFuncRegion
{
	my ($refH_file, $sample, $ws, $rowStart_SBSdistrBySeg, $colStart_SBSdistrBySeg, $nb_func, $ref_RowSBSDistrBySegAndFuncCG, $mutDistrggplot2) = @_;

	my $row_SBSdistrBySeg                = $rowStart_SBSdistrBySeg+4;
	my $row_SBSDistrBySegAndFunc_CA      = $rowStart_SBSdistrBySeg + $nb_func + 12;
	my $rowEndCG_SBSDistrBySegAndFunc_CG = $$ref_RowSBSDistrBySegAndFuncCG + $nb_func;
	my $row_SBSDistrBySegAndFunc_CT      = $rowStart_SBSdistrBySeg + ($nb_func*3) + 20;
	my $colTable4                        = 0;

	my ($count_ca, $count_cg, $count_ct, $count_ta, $count_tc, $count_tg) = (0,0,0,0,0,0);


	## 6 mutation types by segment
	foreach my $k_func (sort keys %{$refH_file->{$sample}{'6mutType'}})
	{
		my $totalSBS_bySegment = 0;

		# Write the functional region for the section SBS distribution by segment (Table 4)
		$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg, $k_func, $formatT_left);
		# Write the exonic func for the section strand bias by segment (Table 5)
		$ws->write($row_SBSDistrBySegAndFunc_CA, $colStart_SBSdistrBySeg, $k_func, $formatT_left);
		# Write the last functional element in the table
		if($$ref_RowSBSDistrBySegAndFuncCG == $rowEndCG_SBSDistrBySegAndFunc_CG)
		{
			$ws->write($$ref_RowSBSDistrBySegAndFuncCG, $colStart_SBSdistrBySeg, $k_func, $formatT_bottomLeft);
		}
		else
		{
			$ws->write($$ref_RowSBSDistrBySegAndFuncCG, $colStart_SBSdistrBySeg, $k_func, $formatT_left);
		}

		foreach my $k_mutation (sort keys %{$refH_file->{$sample}{'6mutType'}{$k_func}})
		{
			### Write the count of SBS per mutation on genomic (Table 4) and coding strand (Table 5)
			if($k_mutation eq "C:G>A:T")
			{
				writeCountSBS($refH_file, $sample, $k_func, $k_mutation, $ws, $row_SBSDistrBySegAndFunc_CA, $colStart_SBSdistrBySeg, $colStart_SBSdistrBySeg+3, $row_SBSdistrBySeg, $rowEndCG_SBSDistrBySegAndFunc_CG, \$count_ca);
			}
			if($k_mutation eq "C:G>G:C")
			{
				writeCountSBS($refH_file, $sample, $k_func, $k_mutation, $ws, $row_SBSDistrBySegAndFunc_CA, $colStart_SBSdistrBySeg+4, $colStart_SBSdistrBySeg+5, $row_SBSdistrBySeg, $rowEndCG_SBSDistrBySegAndFunc_CG, \$count_cg);
			}
			if($k_mutation eq "C:G>T:A")
			{
				writeCountSBS($refH_file, $sample, $k_func, $k_mutation, $ws, $row_SBSDistrBySegAndFunc_CA, $colStart_SBSdistrBySeg+8, $colStart_SBSdistrBySeg+7, $row_SBSdistrBySeg, $rowEndCG_SBSDistrBySegAndFunc_CG, \$count_ct);
			}
			if($k_mutation eq "T:A>A:T")
			{
				writeCountSBS($refH_file, $sample, $k_func, $k_mutation, $ws, $$ref_RowSBSDistrBySegAndFuncCG, $colStart_SBSdistrBySeg, $colStart_SBSdistrBySeg+9, $row_SBSdistrBySeg, $rowEndCG_SBSDistrBySegAndFunc_CG, \$count_ta);
			}
			if($k_mutation eq "T:A>C:G")
			{
				writeCountSBS($refH_file, $sample, $k_func, $k_mutation, $ws, $$ref_RowSBSDistrBySegAndFuncCG, $colStart_SBSdistrBySeg+4, $colStart_SBSdistrBySeg+11, $row_SBSdistrBySeg, $rowEndCG_SBSDistrBySegAndFunc_CG, \$count_tc);
			}
			if($k_mutation eq "T:A>G:C")
			{
				writeCountSBS($refH_file, $sample, $k_func, $k_mutation, $ws, $$ref_RowSBSDistrBySegAndFuncCG, $colStart_SBSdistrBySeg+8, $colStart_SBSdistrBySeg+13, $row_SBSdistrBySeg, $rowEndCG_SBSDistrBySegAndFunc_CG, \$count_tg);
			}

			# Calculate the total number of SBS on the genomic strand for each mutation types by exonic region
			$totalSBS_bySegment += $refH_file->{$sample}{'6mutType'}{$k_func}{$k_mutation}{'TotalMutG'};
		} # End $k_mutation

		$row_SBSDistrBySegAndFunc_CA++; $$ref_RowSBSDistrBySegAndFuncCG++; #$row_SBSDistrBySegAndFunc_CT++;

		# Write the percent by exonic region
		writePercentSBS($refH_file, $sample, $k_func, "C:G>A:T", $row_SBSdistrBySeg, $colStart_SBSdistrBySeg+2, $ws, $totalSBS_bySegment);
		writePercentSBS($refH_file, $sample, $k_func, "C:G>G:C", $row_SBSdistrBySeg, $colStart_SBSdistrBySeg+4, $ws, $totalSBS_bySegment);
		writePercentSBS($refH_file, $sample, $k_func, "C:G>T:A", $row_SBSdistrBySeg, $colStart_SBSdistrBySeg+6, $ws, $totalSBS_bySegment);
		writePercentSBS($refH_file, $sample, $k_func, "T:A>A:T", $row_SBSdistrBySeg, $colStart_SBSdistrBySeg+8, $ws, $totalSBS_bySegment);
		writePercentSBS($refH_file, $sample, $k_func, "T:A>C:G", $row_SBSdistrBySeg, $colStart_SBSdistrBySeg+10, $ws, $totalSBS_bySegment);
		writePercentSBS($refH_file, $sample, $k_func, "T:A>G:C", $row_SBSdistrBySeg, $colStart_SBSdistrBySeg+12, $ws, $totalSBS_bySegment);

		# Write the count of SBS by segment
		$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+1, $totalSBS_bySegment, $format_A10);

		$row_SBSdistrBySeg++;
	} # End $k_func


	# Write the total number of SBS on the genomic strand (Table 4)
	$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+1, $refH_file->{$sample}{'TotalSBSGenomic'}, $formatT_bottomHeader);


	##### Calculate the total percentages by mutation type
	my $percent_ca = ($count_ca / $refH_file->{$sample}{'TotalSBSGenomic'}) * 100; $percent_ca = sprintf("%.2f", $percent_ca);
	my $percent_cg = ($count_cg / $refH_file->{$sample}{'TotalSBSGenomic'}) * 100; $percent_cg = sprintf("%.2f", $percent_cg);
	my $percent_ct = ($count_ct / $refH_file->{$sample}{'TotalSBSGenomic'}) * 100; $percent_ct = sprintf("%.2f", $percent_ct);
	my $percent_ta = ($count_ta / $refH_file->{$sample}{'TotalSBSGenomic'}) * 100; $percent_ta = sprintf("%.2f", $percent_ta);
	my $percent_tc = ($count_tc / $refH_file->{$sample}{'TotalSBSGenomic'}) * 100; $percent_tc = sprintf("%.2f", $percent_tc);
	my $percent_tg = ($count_tg / $refH_file->{$sample}{'TotalSBSGenomic'}) * 100; $percent_tg = sprintf("%.2f", $percent_tg);


	# Write the total percentage (Table 4)
	$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+2, $percent_ca, $formatT_bottom);
	$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+3, $count_ca, $formatT_bottomHeader);
	$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+4, $percent_cg, $formatT_bottom);
	$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+5, $count_cg, $formatT_bottomHeader);
	$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+6, $percent_ct, $formatT_bottom);
	$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+7, $count_ct, $formatT_bottomHeader);
	$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+8, $percent_ta, $formatT_bottom);
	$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+9, $count_ta, $formatT_bottomHeader);
	$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+10, $percent_tc, $formatT_bottom);
	$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+11, $count_tc, $formatT_bottomHeader);
	$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+12, $percent_tg, $formatT_bottom);
	$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+13, $count_tg, $formatT_bottomRightHeader);



	# Overall distribution of SBS (Table 1)
	$ws->write(0, 0, "Graph 1. SBS distribution", $formatT_graphTitle); $ws->set_row(0, 18);
	$ws->write(29, 0, "Table 1. Frequency and counts of all SBS", $format_A10Boldleft);
	$ws->write(30, 0, "Mutation type", $table_topleft);
	$ws->write(30, 1, "Percentage", $table_top);
	$ws->write(30, 2, "Count", $table_topRight);
	$ws->write(31, 0, "C:G>A:T", $table_left); $ws->write(31, 1, $percent_ca, $format_A10); $ws->write(31, 2, $count_ca, $table_right);
	$ws->write(32, 0, "C:G>G:C", $table_left); $ws->write(32, 1, $percent_cg, $format_A10); $ws->write(32, 2, $count_cg, $table_right);
	$ws->write(33, 0, "C:G>T:A", $table_left); $ws->write(33, 1, $percent_ct, $format_A10); $ws->write(33, 2, $count_ct, $table_right);
	$ws->write(34, 0, "T:A>A:T", $table_left); $ws->write(34, 1, $percent_ta, $format_A10); $ws->write(34, 2, $count_ta, $table_right);
	$ws->write(35, 0, "T:A>C:G", $table_left); $ws->write(35, 1, $percent_tc, $format_A10); $ws->write(35, 2, $count_tc, $table_right);
	$ws->write(36, 0, "T:A>G:C", $table_left); $ws->write(36, 1, $percent_tg, $format_A10); $ws->write(36, 2, $count_tg, $table_right);
	$ws->write(37, 0, "", $table_bottomleft); $ws->write(37, 1, "", $table_bottom); $ws->write(37, 2, $refH_file->{$sample}{'TotalSBSGenomic'}, $table_bottomrightHeader);

	# Create an input for ggplot2 for representing the distribution of SBS for each mutation types (Figure 1)
	open(DISTRSBS, ">", $mutDistrggplot2) or die "$!: $mutDistrggplot2\n";
	print DISTRSBS "Mutation_Type\tCount\tPercentage\tSample\n";
	print DISTRSBS "C:G>A:T\t$count_ca\t$percent_ca\t$sample\n";
	print DISTRSBS "C:G>G:C\t$count_cg\t$percent_cg\t$sample\n";
	print DISTRSBS "C:G>T:A\t$count_ct\t$percent_ct\t$sample\n";
	print DISTRSBS "T:A>A:T\t$count_ta\t$percent_ta\t$sample\n";
	print DISTRSBS "T:A>C:G\t$count_tc\t$percent_tc\t$sample\n";
	print DISTRSBS "T:A>G:C\t$count_tg\t$percent_tg\t$sample\n";
	close DISTRSBS;
}
# Write the percentage in table 4 of the Excel report (Sub function of writeStatbyFuncRegion)
sub writePercentSBS
{
	my ($refH_file, $k_file, $k_func, $mutation, $row, $col, $ws, $totalSBS) = @_;

	my $percent = 0;
	if($refH_file->{$k_file}{'6mutType'}{$k_func}{$mutation}{'TotalMutG'} == 0)
	{
		$percent = 0;
	}
	else
	{
		$percent = ($refH_file->{$k_file}{'6mutType'}{$k_func}{$mutation}{'TotalMutG'} / $totalSBS ) * 100;
		$percent = sprintf("%.2f", $percent);
	}
	$ws->write($row, $col, $percent, $format_A10);
}
# Write the count in table 4 and table 5 of the Excel report (Sub function of writeStatbyFuncRegion)
sub writeCountSBS
{
	my ($refH_file, $k_file, $k_func, $k_mutation, $ws, $row_SBSDistrBySegAndFunc_CA, $colStart_SBSdistrBySeg, $colTable4, $row_SBSdistrBySeg, $rowEndCG_SBSDistrBySegAndFunc_CG, $refCount) = @_;

	my $ratioSB = 0;
	if( ($refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'} == 0) || ($refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'} == 0) )
	{
		$ratioSB = 0;
	}
	else
	{
		$ratioSB = $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'} / $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'};
	}
	$ratioSB = sprintf("%.2f", $ratioSB);


	if($row_SBSDistrBySegAndFunc_CA == $rowEndCG_SBSDistrBySegAndFunc_CG)
	{
		# Write the ratio of NonTr / Tr (Table 5)
		$ws->write($row_SBSDistrBySegAndFunc_CA, $colStart_SBSdistrBySeg+1, $ratioSB, $formatT_bottom);
		# Write the count of SBS in the NonTr and Tr strand
		$ws->write($row_SBSDistrBySegAndFunc_CA, $colStart_SBSdistrBySeg+2, $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'}, $formatT_bottom);
		$ws->write($row_SBSDistrBySegAndFunc_CA, $colStart_SBSdistrBySeg+3, $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'}, $formatT_bottom);

		if( ($k_mutation eq "C:G>T:A") || ($k_mutation eq "T:A>G:C") )
		{
			$ws->write($row_SBSDistrBySegAndFunc_CA, $colStart_SBSdistrBySeg+3, $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'}, $formatT_bottomRight);
		}
	}
	else
	{
		# Write the ratio of NonTr / Tr (Table 5)
		$ws->write($row_SBSDistrBySegAndFunc_CA, $colStart_SBSdistrBySeg+1, $ratioSB, $format_A10);
		# Write the count of SBS in the NonTr and Tr strand
		$ws->write($row_SBSDistrBySegAndFunc_CA, $colStart_SBSdistrBySeg+2, $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'}, $format_A10);
		$ws->write($row_SBSDistrBySegAndFunc_CA, $colStart_SBSdistrBySeg+3, $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'}, $format_A10);

		if( ($k_mutation eq "C:G>T:A") || ($k_mutation eq "T:A>G:C") )
		{
			$ws->write($row_SBSDistrBySegAndFunc_CA, $colStart_SBSdistrBySeg+3, $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'}, $formatT_right);
		}
	}

	if($k_mutation eq "C:G>A:T")
	{
		# Calculate the total number of SBS per mut type (genomic strand)
		$$refCount += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'TotalMutG'};
	}
	elsif($k_mutation eq "C:G>G:C")
	{
		$$refCount += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'TotalMutG'};
	}
	elsif($k_mutation eq "C:G>T:A")
	{
		$$refCount += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'TotalMutG'};
	}
	elsif($k_mutation eq "T:A>A:T")
	{
		$$refCount += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'TotalMutG'};
	}
	elsif($k_mutation eq "T:A>C:G")
	{
		$$refCount += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'TotalMutG'};
	}
	elsif($k_mutation eq "T:A>G:C")
	{
		$$refCount += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'TotalMutG'};
	}


	# Write the count by exonic region (Table 4)
	$ws->write($row_SBSdistrBySeg, $colTable4, $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'TotalMutG'}, $format_A10);

	if($k_mutation eq "T:A>G:C")
	{
		$ws->write($row_SBSdistrBySeg, $colTable4, $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'TotalMutG'}, $formatT_right);
	}
}


############################################################################################################
# SBS distribution by chromosomes (Table 6)
sub writeDistrByChr
{
	my ($ws, $refH_file, $sample, $row, $col, $distrByChrggplot2) = @_;


	# For the HTML report
  open(SBSPerChr, ">", $distrByChrggplot2) or die "$!: $distrByChrggplot2\n";
  print SBSPerChr "\tPearson\t$refH_file->{$sample}{'SBSPerChr'}{'AllMutType'}\t", $refH_file->{$sample}{'SBSPerChr'}{"C:G>A:T"}{'Pearson'},"\t", $refH_file->{$sample}{'SBSPerChr'}{"C:G>G:C"}{'Pearson'},"\t", $refH_file->{$sample}{'SBSPerChr'}{"C:G>T:A"}{'Pearson'},"\t", $refH_file->{$sample}{'SBSPerChr'}{"T:A>A:T"}{'Pearson'},"\t", $refH_file->{$sample}{'SBSPerChr'}{"T:A>C:G"}{'Pearson'},"\t", $refH_file->{$sample}{'SBSPerChr'}{"T:A>G:C"}{'Pearson'},"\n";
  print SBSPerChr "Chr\tSize\tAll SBS\tC:G>A:T\tC:G>G:C\tC:G>T:A\tT:A>A:T\tT:A>C:G\tT:A>G:C\n";


  my $row_SBSPerChr = $row + 8;


  # Write the Pearson coefficient
  $ws->write($row+6, $col+3, $refH_file->{$sample}{'SBSPerChr'}{"C:G>A:T"}{'Pearson'}, $format_A10);
  $ws->write($row+6, $col+4, $refH_file->{$sample}{'SBSPerChr'}{"C:G>G:C"}{'Pearson'}, $format_A10);
  $ws->write($row+6, $col+5, $refH_file->{$sample}{'SBSPerChr'}{"C:G>T:A"}{'Pearson'}, $format_A10);
  $ws->write($row+6, $col+6, $refH_file->{$sample}{'SBSPerChr'}{"T:A>A:T"}{'Pearson'}, $format_A10);
  $ws->write($row+6, $col+7, $refH_file->{$sample}{'SBSPerChr'}{"T:A>C:G"}{'Pearson'}, $format_A10);
  $ws->write($row+6, $col+8, $refH_file->{$sample}{'SBSPerChr'}{"T:A>G:C"}{'Pearson'}, $formatT_right);

  # Write the chromosome number and their sizes / Write count SBS per chromosomes
  my $line=0;

  foreach my $chromosome (sort keys %chromosomes)
	{
		$ws->write($row_SBSPerChr+($line), $col, $chromosome, $formatT_left);
		$ws->write($row_SBSPerChr+($line), $col+1, $chromosomes{$chromosome}, $format_A10);
		$ws->write($row_SBSPerChr+($line), $col+2, $refH_file->{$sample}{'SBSPerChr'}{'TotalPerChr'}{$chromosome}{'chr'}, $format_A10);

		# Write the count per mutation
		$ws->write($row_SBSPerChr+($line), $col+3, $refH_file->{$sample}{'SBSPerChr'}{"C:G>A:T"}{'CHR'}{$chromosome}{'chr'}, $format_A10);
		$ws->write($row_SBSPerChr+($line), $col+4, $refH_file->{$sample}{'SBSPerChr'}{"C:G>G:C"}{'CHR'}{$chromosome}{'chr'}, $format_A10);
		$ws->write($row_SBSPerChr+($line), $col+5, $refH_file->{$sample}{'SBSPerChr'}{"C:G>T:A"}{'CHR'}{$chromosome}{'chr'}, $format_A10);
		$ws->write($row_SBSPerChr+($line), $col+6, $refH_file->{$sample}{'SBSPerChr'}{"T:A>A:T"}{'CHR'}{$chromosome}{'chr'}, $format_A10);
		$ws->write($row_SBSPerChr+($line), $col+7, $refH_file->{$sample}{'SBSPerChr'}{"T:A>C:G"}{'CHR'}{$chromosome}{'chr'}, $format_A10);
		$ws->write($row_SBSPerChr+($line), $col+8, $refH_file->{$sample}{'SBSPerChr'}{"T:A>G:C"}{'CHR'}{$chromosome}{'chr'}, $formatT_right);


		# For the HTML report
		print SBSPerChr "$chromosome\t", $chromosomes{$chromosome},"\t", $refH_file->{$sample}{'SBSPerChr'}{'TotalPerChr'}{$chromosome}{'chr'},"\t", $refH_file->{$sample}{'SBSPerChr'}{"C:G>A:T"}{'CHR'}{$chromosome}{'chr'},"\t", $refH_file->{$sample}{'SBSPerChr'}{"C:G>G:C"}{'CHR'}{$chromosome}{'chr'},"\t", $refH_file->{$sample}{'SBSPerChr'}{"C:G>T:A"}{'CHR'}{$chromosome}{'chr'},"\t", $refH_file->{$sample}{'SBSPerChr'}{"T:A>A:T"}{'CHR'}{$chromosome}{'chr'},"\t", $refH_file->{$sample}{'SBSPerChr'}{"T:A>C:G"}{'CHR'}{$chromosome}{'chr'},"\t", $refH_file->{$sample}{'SBSPerChr'}{"T:A>G:C"}{'CHR'}{$chromosome}{'chr'},"\n";

		$line++;
	}

	# Write the Pearson coefficient for the total number of SBS
	$ws->write($row+6, $col+2, $refH_file->{$sample}{'SBSPerChr'}{'AllMutType'}, $format_A10);
	$ws->write($row_SBSPerChr+(keys %chromosomes), $col+2, $refH_file->{$sample}{'TotalSBSGenomic'}, $formatT_bottomHeader);

	print SBSPerChr "\t\t$refH_file->{$sample}{'TotalSBSGenomic'}\n";
	close SBSPerChr;
}


############################################################################################################
# Trinucleotide sequence context on genomic strand (Panel 1)
sub writeTriNtGenomic
{
	my ($ws, $refH_file, $sample, $col, $heatmapCountggplot2, $heatmapPercentggplot2, $triNtBarChartggplot2, $ref_c_ca6_g, $ref_c_cg6_g, $ref_c_ct6_g, $ref_c_ta6_g, $ref_c_tc6_g, $ref_c_tg6_g) = @_;

	# Initialise the row of the panel 1
	my $row_SeqContext6 = 4;
  # Percent total of mutations for 6 mutation types on genomic strand
  my ($p_ca6_g, $p_cg6_g, $p_ct6_g, $p_ta6_g, $p_tc6_g, $p_tg6_g) = (0,0,0, 0,0,0);
  my $maxValue = 0; # For the heatmap

  # For checking if the total number of SBS is correct
  my $total_SBS_genomic = 0;


  open(HEATMAPCGENOMIC, ">", $heatmapCountggplot2) or die "$!: $heatmapCountggplot2\n";
	print HEATMAPCGENOMIC "\tC>A\tC>G\tC>T\tT>A\tT>C\tT>G\n";
	open(HEATMAPPGENOMIC, ">", $heatmapPercentggplot2) or die "$!: $heatmapPercentggplot2\n";
	print HEATMAPPGENOMIC "\tC>A\tC>G\tC>T\tT>A\tT>C\tT>G\n";

	## Bar plot NMF like
	open(BARPLOTNMFLIKE, ">", $triNtBarChartggplot2) or die "$!: $triNtBarChartggplot2\n";
	print BARPLOTNMFLIKE "alteration\tcontext\tvalue\n";


	foreach my $k_context (sort keys %{$refH_file->{$sample}{'SeqContextG'}})
	{
		if( ($k_context =~ /N/) || (length($k_context) != 3) ) { next; }

		# Write the context: 6 mut type on genomic strand
		$ws->write($row_SeqContext6 , $col+3, $k_context, $format_A10);
		$ws->write($row_SeqContext6 , $col+13, $k_context, $format_A10);

  	# Count for the heatmap
  	print HEATMAPCGENOMIC $k_context."\t";
  	print HEATMAPPGENOMIC $k_context."\t";

  	foreach my $k_mutation (sort keys %{$refH_file->{$sample}{'SeqContextG'}{$k_context}})
  	{
  		# For checking the total number of SBS
  		$total_SBS_genomic += $refH_file->{$sample}{'SeqContextG'}{$k_context}{$k_mutation};

  		# Calculate the percentages
  		my $percent = 0;
  		if($refH_file->{$sample}{'SeqContextG'}{$k_context}{$k_mutation} == 0) { $percent = 0; }
  		else
  		{
  			$percent = ($refH_file->{$sample}{'SeqContextG'}{$k_context}{$k_mutation} / $refH_file->{$sample}{'TotalSBSGenomic'}) * 100;
  			$percent = sprintf("%.2f", $percent);
  		}

  		# For representing the sequence context with a bar plot (NMF like style)
  		print BARPLOTNMFLIKE $k_mutation,"\t", $k_context,"\t", $percent,"\n";

  		# Write the count for the heatmap
  		print HEATMAPCGENOMIC $refH_file->{$sample}{'SeqContextG'}{$k_context}{$k_mutation}."\t";
  		print HEATMAPPGENOMIC "$percent\t";


  		# For NMF input
		  my $count = $refH_file->{$sample}{'SeqContextG'}{$k_context}{$k_mutation};
			if($sample ne "Pool_Data") { push(@{$h_inputNMF{'Count'}{$k_context}{$k_mutation}}, $count); }
			if($sample ne "Pool_Data") { push(@{$h_inputNMF{'Percent'}{$k_context}{$k_mutation}}, $percent); }


  		if($k_mutation eq "C>A")
  		{
  			triNtByMut($ws, $row_SeqContext6, $col+4, $col+14, $refH_file, $sample, $k_context, $k_mutation, $percent, $maxValue, $ref_c_ca6_g, \$p_ca6_g);
  		}
  		elsif($k_mutation eq "C>G")
  		{
  			triNtByMut($ws, $row_SeqContext6, $col+5, $col+15, $refH_file, $sample, $k_context, $k_mutation, $percent, $maxValue, $ref_c_cg6_g, \$p_cg6_g);
  		}
  		elsif($k_mutation eq "C>T")
  		{
  			triNtByMut($ws, $row_SeqContext6, $col+6, $col+16, $refH_file, $sample, $k_context, $k_mutation, $percent, $maxValue, $ref_c_ct6_g, \$p_ct6_g);
  		}
  		elsif($k_mutation eq "T>A")
  		{
  			triNtByMut($ws, $row_SeqContext6, $col+7, $col+17, $refH_file, $sample, $k_context, $k_mutation, $percent, $maxValue, $ref_c_ta6_g, \$p_ta6_g);
  		}
  		elsif($k_mutation eq "T>C")
  		{
  			triNtByMut($ws, $row_SeqContext6, $col+8, $col+18, $refH_file, $sample, $k_context, $k_mutation, $percent, $maxValue, $ref_c_tc6_g, \$p_tc6_g);
  		}
  		elsif($k_mutation eq "T>G")
  		{
  			triNtByMut($ws, $row_SeqContext6, $col+9, $col+19, $refH_file, $sample, $k_context, $k_mutation, $percent, $maxValue, $ref_c_tg6_g, \$p_tg6_g);
  		}
  		else
  		{
  			print STDERR "Error: Mutation type not considered for: $k_mutation\n";
  			exit;
  		}
  	}
  	$row_SeqContext6++;

  	print HEATMAPCGENOMIC "\n";
  	print HEATMAPPGENOMIC "\n";
  }
	close HEATMAPCGENOMIC; close HEATMAPPGENOMIC;
	close BARPLOTNMFLIKE;


	# Write the total number of SBS per mutation type: COUNT
	$ws->write($row_SeqContext6, $col+4, $$ref_c_ca6_g, $formatT_bottomHeader2);
	$ws->write($row_SeqContext6, $col+5, $$ref_c_cg6_g, $formatT_bottomHeader2);
	$ws->write($row_SeqContext6, $col+6, $$ref_c_ct6_g, $formatT_bottomHeader2);
	$ws->write($row_SeqContext6, $col+7, $$ref_c_ta6_g, $formatT_bottomHeader2);
	$ws->write($row_SeqContext6, $col+8, $$ref_c_tc6_g, $formatT_bottomHeader2);
	$ws->write($row_SeqContext6, $col+9, $$ref_c_tg6_g, $formatT_bottomHeader2);
  if($total_SBS_genomic != $refH_file->{$sample}{'TotalSBSGenomic'})
  {
  	print STDERR "Error in the calculation of the total number of SBS on the genomic strand!!!!\n";
  	print STDERR "From hash table $refH_file->{$sample}{'TotalSBSGenomic'}\tVS\t$total_SBS_genomic\n";
  	exit;
  }

	# Write the total number of SBS per mutation type: PERCENT
	$ws->write($row_SeqContext6, $col+14, $p_ca6_g, $formatT_bottomHeader2);
	$ws->write($row_SeqContext6, $col+15, $p_cg6_g, $formatT_bottomHeader2);
	$ws->write($row_SeqContext6, $col+16, $p_ct6_g, $formatT_bottomHeader2);
	$ws->write($row_SeqContext6, $col+17, $p_ta6_g, $formatT_bottomHeader2);
	$ws->write($row_SeqContext6, $col+18, $p_tc6_g, $formatT_bottomHeader2);
	$ws->write($row_SeqContext6, $col+19, $p_tg6_g, $formatT_bottomHeader2);

	my $totalPercent_genomic = $p_ca6_g + $p_cg6_g + $p_ct6_g + $p_ta6_g + $p_tc6_g + $p_tg6_g;
	$totalPercent_genomic    = sprintf("%.0f", $totalPercent_genomic);

	if($totalPercent_genomic != 100)
	{
		print STDERR "Error in the calculation of the total percentages on the genomic strand!!!\n";
		print STDERR "The total is equal to=\t$totalPercent_genomic\n";
		exit;
	}
}
# Trinucleotide count and percentage by mutation type (Sub function of writeTriNtGenomic)
sub triNtByMut
{
	my ($ws, $row, $colC, $colP, $refH_file, $sample, $context, $mutation, $percent, $maxValue, $refCountG, $refPercentG) = @_;

	### COUNT
	$ws->write($row, $colC, $refH_file->{$sample}{'SeqContextG'}{$context}{$mutation}, $format_A10);

	### PERCENTAGE
	$ws->write($row, $colP, $percent, $format_A10);

	# For the heatmap
	if($percent >= $maxValue) { $maxValue = $percent; }

	# For the total amount per mutation types
	$$refCountG += $refH_file->{$sample}{'SeqContextG'}{$context}{$mutation};
	$$refPercentG += $percent;
}


############################################################################################################
# Trinucleotide sequence context on coding strand (Panel 2)
sub writeTriNtCoding
{
	my ($ws, $row, $col, $refH_file, $sample, $triNtBarChartCodingCountggplot2, $triNtBarChartCodingPercentggplot2) = @_;

	# Initialise the row
	my $row_SeqContext12        = $row+6;
	my $row_SeqContext12Percent = $row+27;

  # Total count and percent calculated for the strand bias
  my ($ca_NonTr, $ca_Tr, $cg_NonTr, $cg_Tr, $ct_NonTr, $ct_Tr, $ta_NonTr, $ta_Tr, $tc_NonTr, $tc_Tr, $tg_NonTr, $tg_Tr) = (0,0,0, 0,0,0, 0,0,0, 0,0,0);
  my ($percent_ca_NonTr, $percent_ca_Tr, $percent_cg_NonTr, $percent_cg_Tr, $percent_ct_NonTr, $percent_ct_Tr, $percent_ta_NonTr, $percent_ta_Tr, $percent_tc_NonTr, $percent_tc_Tr, $percent_tg_NonTr, $percent_tg_Tr) = (0,0,0, 0,0,0, 0,0,0, 0,0,0);

  # For checking if the total number of SBS is correct
  my $total_SBS_coding  = 0;


  open(COUNT, ">", $triNtBarChartCodingCountggplot2) or die "$!: $triNtBarChartCodingCountggplot2\n";
	print COUNT "MutationTypeContext\tStrand\tValue\tSample\n";
	open(PERCENT, ">", $triNtBarChartCodingPercentggplot2) or die "$!: $triNtBarChartCodingPercentggplot2\n";
	print PERCENT "MutationTypeContext\tStrand\tValue\tSample\n";

	foreach my $k_context (sort keys %{$refH_file->{$sample}{'SeqContextC'}})
	{
		if( ($k_context =~ /N/) || (length($k_context) != 3) ) { next; }

		# Write the context: 12 mut type on coding strand
		$ws->write($row_SeqContext12 , $col, $k_context, $formatT_left);
		$ws->write($row_SeqContext12Percent , $col, $k_context, $formatT_left);

  	foreach my $k_mutation (sort keys %{$refH_file->{$sample}{'SeqContextC'}{$k_context}})
  	{
  		# Percent: 12 mut type on coding strand
  		my ($percent_NonTr, $percent_Tr) = (0, 0);

  		if($refH_file->{$sample}{'SeqContextC'}{$k_context}{$k_mutation}{'NonTr'} != 0)
  		{
  			$percent_NonTr = ( $refH_file->{$sample}{'SeqContextC'}{$k_context}{$k_mutation}{'NonTr'} / $refH_file->{$sample}{'TotalSBSCoding'} ) * 100;
  		}

  		if($refH_file->{$sample}{'SeqContextC'}{$k_context}{$k_mutation}{'Tr'} != 0)
  		{
  			$percent_Tr = ( $refH_file->{$sample}{'SeqContextC'}{$k_context}{$k_mutation}{'Tr'} / $refH_file->{$sample}{'TotalSBSCoding'} ) * 100;
  		}


  		# Counts
  		print COUNT "$k_mutation:$k_context\tNonTranscribed\t$refH_file->{$sample}{'SeqContextC'}{$k_context}{$k_mutation}{'NonTr'}\t$sample\n";
  		print COUNT "$k_mutation:$k_context\tTranscribed\t$refH_file->{$sample}{'SeqContextC'}{$k_context}{$k_mutation}{'Tr'}\t$sample\n";

  		# Percentages
			$percent_NonTr = sprintf("%.2f", $percent_NonTr);
			$percent_Tr    = sprintf("%.2f", $percent_Tr);
			print PERCENT "$k_mutation:$k_context\tNonTranscribed\t$percent_NonTr\t$sample\n";
			print PERCENT "$k_mutation:$k_context\tTranscribed\t$percent_Tr\t$sample\n";

  		# Calculate the total number for each mutation types
  		if($k_mutation eq "C>A")
  		{
  			triNtByMutCoding($refH_file, $sample, $k_context, $k_mutation, $ws, $row_SeqContext12, $col+1, $row_SeqContext12Percent, \$ca_NonTr, \$ca_Tr, $percent_NonTr, $percent_Tr);

  			$percent_ca_NonTr += $percent_NonTr;
  			$percent_ca_Tr    += $percent_Tr;
  		}
  		if($k_mutation eq "C>G")
  		{
  			triNtByMutCoding($refH_file, $sample, $k_context, $k_mutation, $ws, $row_SeqContext12, $col+3, $row_SeqContext12Percent, \$cg_NonTr, \$cg_Tr, $percent_NonTr, $percent_Tr);

  			$percent_cg_NonTr += $percent_NonTr;
  			$percent_cg_Tr    += $percent_Tr;
  		}
  		if($k_mutation eq "C>T")
  		{
  			triNtByMutCoding($refH_file, $sample, $k_context, $k_mutation, $ws, $row_SeqContext12, $col+5, $row_SeqContext12Percent, \$ct_NonTr, \$ct_Tr, $percent_NonTr, $percent_Tr);

  			$percent_ct_NonTr += $percent_NonTr;
  			$percent_ct_Tr    += $percent_Tr;
  		}
  		if($k_mutation eq "T>A")
  		{
  			triNtByMutCoding($refH_file, $sample, $k_context, $k_mutation, $ws, $row_SeqContext12, $col+7, $row_SeqContext12Percent, \$ta_NonTr, \$ta_Tr, $percent_NonTr, $percent_Tr);

  			$percent_ta_NonTr += $percent_NonTr;
  			$percent_ta_Tr    += $percent_Tr;
  		}
  		if($k_mutation eq "T>C")
  		{
  			triNtByMutCoding($refH_file, $sample, $k_context, $k_mutation, $ws, $row_SeqContext12, $col+9, $row_SeqContext12Percent, \$tc_NonTr, \$tc_Tr, $percent_NonTr, $percent_Tr);

  			$percent_tc_NonTr += $percent_NonTr;
  			$percent_tc_Tr    += $percent_Tr;
  		}
  		if($k_mutation eq "T>G")
  		{
  			triNtByMutCoding($refH_file, $sample, $k_context, $k_mutation, $ws, $row_SeqContext12, $col+11, $row_SeqContext12Percent, \$tg_NonTr, \$tg_Tr, $percent_NonTr, $percent_Tr);

  			$percent_tg_NonTr += $percent_NonTr;
  			$percent_tg_Tr    += $percent_Tr;
  		}

  		# For checking if the total number of SBS is correct
  		$total_SBS_coding  += $refH_file->{$sample}{'SeqContextC'}{$k_context}{$k_mutation}{'NonTr'} + $refH_file->{$sample}{'SeqContextC'}{$k_context}{$k_mutation}{'Tr'};
  	}
  	$row_SeqContext12++; $row_SeqContext12Percent++;
  }
	close COUNT; close PERCENT;


	## Write the total of each mutation types: 12 mut type on coding strand
	$ws->write($row_SeqContext12, $col+1, $ca_NonTr, $formatT_bottomHeader2); $ws->write($row_SeqContext12, $col+2, $ca_Tr, $formatT_bottomHeader2);
	$ws->write($row_SeqContext12, $col+3, $cg_NonTr, $formatT_bottomHeader2); $ws->write($row_SeqContext12, $col+4, $cg_Tr, $formatT_bottomHeader2);
	$ws->write($row_SeqContext12, $col+5, $ct_NonTr, $formatT_bottomHeader2); $ws->write($row_SeqContext12, $col+6, $ct_Tr, $formatT_bottomHeader2);
	$ws->write($row_SeqContext12, $col+7, $ta_NonTr, $formatT_bottomHeader2); $ws->write($row_SeqContext12, $col+8, $ta_Tr, $formatT_bottomHeader2);
	$ws->write($row_SeqContext12, $col+9, $tc_NonTr, $formatT_bottomHeader2); $ws->write($row_SeqContext12, $col+10, $tc_Tr, $formatT_bottomHeader2);
	$ws->write($row_SeqContext12, $col+11, $tg_NonTr, $formatT_bottomHeader2); $ws->write($row_SeqContext12, $col+12, $tg_Tr, $formatT_bottomHeader2);
	# Write the total percentages of each mutation types: 12 mut type on coding strand
	$ws->write($row_SeqContext12Percent, $col+1, $percent_ca_NonTr, $formatT_bottomHeader); $ws->write($row_SeqContext12Percent, $col+2, $percent_ca_Tr, $formatT_bottomHeader);
	$ws->write($row_SeqContext12Percent, $col+3, $percent_cg_NonTr, $formatT_bottomHeader); $ws->write($row_SeqContext12Percent, $col+4, $percent_cg_Tr, $formatT_bottomHeader);
	$ws->write($row_SeqContext12Percent, $col+5, $percent_ct_NonTr, $formatT_bottomHeader); $ws->write($row_SeqContext12Percent, $col+6, $percent_ct_Tr, $formatT_bottomHeader);
	$ws->write($row_SeqContext12Percent, $col+7, $percent_ta_NonTr, $formatT_bottomHeader); $ws->write($row_SeqContext12Percent, $col+8, $percent_ta_Tr, $formatT_bottomHeader);
	$ws->write($row_SeqContext12Percent, $col+9, $percent_tc_NonTr, $formatT_bottomHeader); $ws->write($row_SeqContext12Percent, $col+10, $percent_tc_Tr, $formatT_bottomHeader);
	$ws->write($row_SeqContext12Percent, $col+11, $percent_tg_NonTr, $formatT_bottomHeader); $ws->write($row_SeqContext12Percent, $col+12, $percent_tg_Tr, $formatT_bottomHeader);

	if($total_SBS_coding == $refH_file->{$sample}{'TotalSBSCoding'})
	{
		$ws->write($row_SeqContext12, $col+13, $refH_file->{$sample}{'TotalSBSCoding'}, $formatT_bottomHeader2)
	}
	else
	{
		print STDERR "Error: in the calculation of the total number of SBS on the coding strand!!!!\n";
		print STDERR "From hash table $refH_file->{$sample}{'TotalSBSCoding'}\tVS\t$total_SBS_coding\n";
		exit;
	}

	my $totalP_SBS_coding = $percent_ca_NonTr + $percent_ca_Tr + $percent_cg_NonTr + $percent_cg_Tr + $percent_ct_NonTr + $percent_ct_Tr + $percent_ta_NonTr + $percent_ta_Tr + $percent_tc_NonTr + $percent_tc_Tr + $percent_tg_NonTr + $percent_tg_Tr;
  		$totalP_SBS_coding = sprintf("%.0f", $totalP_SBS_coding);

	if($totalP_SBS_coding != 100)
	{
		print STDERR "Error: The percentages for the trinucleotide sequence context on the coding strand for 12 mutation types is not equal to 100!!!\n$totalP_SBS_coding\n";
		exit;
	}
}
# Trinucleotide count and percentage by mutation type on Coding strand (Sub function of writeTriNtCoding)
sub triNtByMutCoding
{
	my ($refH_file, $sample, $context, $mutation, $ws, $row, $col, $rowP, $refNonTr, $refTr, $percent_NonTr, $percent_Tr) = @_;

	$$refNonTr += $refH_file->{$sample}{'SeqContextC'}{$context}{$mutation}{'NonTr'};
	$$refTr    += $refH_file->{$sample}{'SeqContextC'}{$context}{$mutation}{'Tr'};

	# COUNT : 12 mutation type (stranded bar graph)
  $ws->write($row, $col, $refH_file->{$sample}{'SeqContextC'}{$context}{$mutation}{'NonTr'}, $format_A10);
  $ws->write($row, $col+1, $refH_file->{$sample}{'SeqContextC'}{$context}{$mutation}{'Tr'}, $format_A10);


	## PERCENT : 12 mutation type (stranded bar graph)
	$ws->write($rowP, $col, $percent_NonTr, $format_A10);
	$ws->write($rowP, $col+1, $percent_Tr, $format_A10);
}


############################################################################################################
# Create and write the figures on the Excel report
sub createWriteFigs
{
	my ($ws, $row, $col, $folderFigure, $sample, $c_ca6_g, $c_cg6_g, $c_ct6_g, $c_ta6_g, $c_tc6_g, $c_tg6_g) = @_;

	######## Create figures
	# Bar char for SBS distribution (Figure 1)
	# Pie char for Impact on protein sequence (Figure 2)
	# Stranded distribution of SBS (Figure 3)
	# Heatmaps for trinucleotide context
	`Rscript $pathRScriptFigs --folderFigure $folderFigure --folderTemp $folder_temp --filename $sample`;

	# Bar chart for trinucleotide context on coding strand
	`Rscript $pathRScriptTxnSB $folderFigure/Stranded_Analysis/$sample/$sample-StrandedSignatureCount.txt $folderFigure/Stranded_Analysis/$sample/$sample-StrandedSignatureCount $folder_temp/$sample-StrandedSignatureCount Count`;

	`Rscript $pathRScriptTxnSB $folderFigure/Stranded_Analysis/$sample/$sample-StrandedSignaturePercent.txt $folderFigure/Stranded_Analysis/$sample/$sample-StrandedSignaturePercent $folder_temp/$sample-StrandedSignaturePercent Percent`;

	# Bar plot for representing the sequence context (NMF like style)
	`Rscript $pathRScriptMutSpectrum $folderFigure/Trinucleotide_Sequence_Context/$sample/$sample-MutationSpectraPercent-Genomic.txt $sample $folderFigure/Trinucleotide_Sequence_Context/$sample $folder_temp $c_ca6_g $c_cg6_g $c_ct6_g $c_ta6_g $c_tc6_g $c_tg6_g`;


	######## Write the figures in the Excel report
	# Bar char for SBS distribution (Figure 1)
	$ws->insert_image(1, 0, "$folder_temp/$sample-SBS_distribution-Report.png", 0, 0, .2, .2);

	# Impact of the SBS on the protein (Figure 2)
	$ws->write(0, 6, "Graph 2. Impact on protein sequence", $formatT_graphTitle);
	$ws->insert_image(1, 6, "$folder_temp/$sample-DistributionExoFunc-Report.png", 0, 0, .2, .2);

	# Stranded distribution of SBS (Figure 3)
	$ws->write(0, 11, "Graph 3. Stranded distribution of SBS", $formatT_graphTitle);
	$ws->insert_image(1, 11, "$folder_temp/$sample-StrandBias-Report.png", 0, 0, .2, .2);

	## Trinucleotide context on coding strand (Scale the inserted image: width x 0.7, height x 0.8)
	$ws->insert_image($row+3, $col+15, "$folder_temp/$sample-StrandedSignatureCount-Report.png", 0, 0, .16, .16);
	$ws->insert_image($row+24, $col+15, "$folder_temp/$sample-StrandedSignaturePercent-Report.png", 0, 0, .16, .16);

	# Heatamp for the sequence context on the genomic strand (6 mutation types)
	$ws->insert_image(4, $col, "$folder_temp/$sample-HeatmapCount-Genomic-Report.png");
	$ws->insert_image(4, $col+10, "$folder_temp/$sample-HeatmapPercent-Genomic-Report.png");

	# Bar plot for the sequence context on the genomic strand (6 mutation types)
	$ws->insert_image(27, $col+3, "$folder_temp/$sample-MutationSpectraPercent-Genomic-Report.png");
}


############################################################################################################
# Write NMF input for count and percentages in the Excel report
sub writeInputNMF
{
	my ($ws_inputNMF_count, $ws_inputNMF_percent, $outCount, $outPercent) = @_;


	open(OUTINPUTNMFC, ">", $outCount) or die "$!: $outCount\n"; # with the count
	open(OUTINPUTNMFP, ">", $outPercent) or die "$!: $outPercent\n"; # With the frequency un-normalized

	foreach my $k_sample (@{$h_inputNMF{'Sample'}})
	{
		print OUTINPUTNMFC "\t$k_sample";
		print OUTINPUTNMFP "\t$k_sample";
	}
	print OUTINPUTNMFC "\n"; print OUTINPUTNMFP "\n";

	my $row_inputNMF = 1;
	foreach my $k_context (sort keys %{$h_inputNMF{'Count'}})
	{
		$k_context =~ /(\w)_(\w)/; my ($base5, $base3) = ($1, $2);
		foreach my $k_mutation (sort keys %{$h_inputNMF{'Count'}{$k_context}})
		{
			my ($col_inputNMF_Count, $col_inputNMF_Percent) = (1, 1);
			my $contextNMF = $base5."[$k_mutation]".$base3;

			# Write the input in the Excel report, only when all the samples are in the same workbook
			if($oneReportPerSample == 2)
			{
				$ws_inputNMF_count->write($row_inputNMF, 0, $contextNMF);
				$ws_inputNMF_percent->write($row_inputNMF, 0, $contextNMF);
			}

			print OUTINPUTNMFC $contextNMF,"\t"; print OUTINPUTNMFP $contextNMF,"\t";

			foreach (@{$h_inputNMF{'Count'}{$k_context}{$k_mutation}})   { print OUTINPUTNMFC "$_\t"; } print OUTINPUTNMFC "\n";
			foreach (@{$h_inputNMF{'Percent'}{$k_context}{$k_mutation}}) { print OUTINPUTNMFP "$_\t"; } print OUTINPUTNMFP "\n";

			foreach (@{$h_inputNMF{'Count'}{$k_context}{$k_mutation}})
			{
				if($oneReportPerSample == 2)
				{
					$ws_inputNMF_count->write($row_inputNMF, $col_inputNMF_Count, $_);
				}
			 	$col_inputNMF_Count++;
			}
			foreach (@{$h_inputNMF{'Percent'}{$k_context}{$k_mutation}})
			{
				if($oneReportPerSample == 2)
				{
					$ws_inputNMF_percent->write($row_inputNMF, $col_inputNMF_Percent, $_);
				}
			  $col_inputNMF_Percent++;
			 }
			 $row_inputNMF++;
		}
	}
	close OUTINPUTNMFP; close OUTINPUTNMFC;
}


######################################################################################################################################################
#																								Define format and background colors for the Excel report																						 #
######################################################################################################################################################
# Font: Arial size 10
sub Format_A10
{
	my ($wb, $format) = @_;
	$$format = $wb->add_format(font=>'Arial', size=>10); $$format->set_align('center');
}
# Font: Arial size 11 bold and center
sub Format_A11Bold
{
	my ($wb, $format) = @_;
	$$format = $wb->add_format(font=>'Arial', size=>11, bold=>1); $$format->set_align('center');
}
# Font: Arial size 10 italic red and center
sub Format_A10ItalicRed
{
	my ($wb, $format) = @_;
	$$format = $wb->add_format(font=>'Arial', size=>10, italic=>1, color => 'red'); $$format->set_align('center');
}
# Format: Arialt size 11 bold and left
sub Format_A11BoldLeft
{
	my ($wb, $format) = @_;
	$$format = $wb->add_format(valign =>'left', font=>'Arial', size=>11, bold=>1);
}
# Font: Arialt size 10 bold and left
sub Format_A10BoldLeft
{
	my ($wb, $format) = @_;
	$$format = $wb->add_format(valign =>'left', font=>'Arial', size=>10, bold=>1);
}
# Define the format of the border of the section (for delimiting the different section of the report)
sub Format_section
{
	my ($wb, $format_topLeft, $format_topRight, $format_bottomLeft, $format_bottomRight, $format_top, $format_right, $format_bottom, $format_left) = @_;

	$$format_topLeft = $wb->add_format(valign  => 'left', bold => 1, font => 'Arial', size => 12);
	$$format_topLeft->set_top(2); $$format_topLeft->set_top_color('blue');
	$$format_topLeft->set_left(2); $$format_topLeft->set_left_color('blue');

	$$format_topRight = $wb->add_format(valign  => 'left', bold => 1, font => 'Arial', size => 12);
	$$format_topRight->set_top(2); $$format_topRight->set_top_color('blue');
	$$format_topRight->set_right(2); $$format_topRight->set_right_color('blue');

	$$format_bottomLeft = $wb->add_format(valign  => 'left', bold => 1, font => 'Arial', size => 12);
	$$format_bottomLeft->set_bottom(2); $$format_bottomLeft->set_bottom_color('blue');
	$$format_bottomLeft->set_left(2); $$format_bottomLeft->set_left_color('blue');

	$$format_bottomRight = $wb->add_format(valign  => 'left', bold => 1, font => 'Arial', size => 12);
	$$format_bottomRight->set_bottom(2); $$format_bottomRight->set_bottom_color('blue');
	$$format_bottomRight->set_right(2); $$format_bottomRight->set_right_color('blue');

	$$format_top    = $wb->add_format(); $$format_top->set_top(2);       $$format_top->set_top_color('blue');
	$$format_right  = $wb->add_format(); $$format_right->set_right(2);   $$format_right->set_right_color('blue');
	$$format_bottom = $wb->add_format(); $$format_bottom->set_bottom(2); $$format_bottom->set_bottom_color('blue');
	$$format_left   = $wb->add_format(); $$format_left->set_left(2);     $$format_left->set_left_color('blue');
}
# Define the header
sub Format_Header
{
	my ($wb, $format_CA, $format_CG, $format_CT, $format_TA, $format_TC, $format_TG, $format_TG2, $format_LeftHeader, $format_RightHeader, $format_LeftHeader2) = @_;

	my ($blue, $black, $red, $gray, $green, $pink);
	Color($wb, \$blue, \$black, \$red, \$gray, \$green, \$pink);

	my ($bgColor_blue, $bgColor_black, $bgColor_red, $bgColor_gray, $bgColor_green, $bgColor_pink);
	BackgroundColor($wb, \$bgColor_blue, \$bgColor_black, \$bgColor_red, \$bgColor_gray, \$bgColor_green, \$bgColor_pink);


	$$format_CA = $wb->add_format(bg_color => $blue, font=>'Arial', bold=>1, size=>11, color=>'white'); $$format_CA->set_align('center'); $$format_CA->set_center_across();
	$$format_CG = $wb->add_format(bg_color => $black, font=>'Arial', bold=>1, size=>11, color=>'white'); $$format_CG->set_align('center'); $$format_CG->set_center_across();
	$$format_CT = $wb->add_format(bg_color => $red, font=>'Arial', bold=>1, size=>11, color=>'white'); $$format_CT->set_align('center'); $$format_CT->set_center_across();
	$$format_TA = $wb->add_format(bg_color => $gray, font=>'Arial', bold=>1, size=>11, color=>'white'); $$format_TA->set_align('center'); $$format_TA->set_center_across();
	$$format_TC = $wb->add_format(bg_color => $green, font=>'Arial', bold=>1, size=>11, color=>'white'); $$format_TC->set_align('center'); $$format_TC->set_center_across();
	$$format_TG = $wb->add_format(bg_color=>$bgColor_pink, font=>'Arial', bold=>1, size=>11, color=>'white'); $$format_TG->set_align('center'); $$format_TG->set_center_across();
	$$format_TG->set_right(2); $$format_TG->set_right_color('blue');

	$$format_TG2 = $wb->add_format(bg_color => $pink, font=>'Arial', bold=>1, size=>11, color=>'white'); $$format_TG2->set_align('center'); $$format_TG2->set_center_across();

	$$format_LeftHeader = $wb->add_format(bold=>1, font=>'Arial', size=>11); $$format_LeftHeader->set_align('center'); $$format_LeftHeader->set_left(2); $$format_LeftHeader->set_left_color('blue');
	$$format_LeftHeader2 = $wb->add_format(bold=>1, font=>'Arial', size=>11); $$format_LeftHeader2->set_left(2); $$format_LeftHeader2->set_left_color('blue');
	$$format_RightHeader = $wb->add_format(bold=>1, font=>'Arial', size=>11); $$format_RightHeader->set_align('center'); $$format_RightHeader->set_right(2); $$format_RightHeader->set_right_color('blue');
}
# Define the header for the part "Strand bias by segment"
sub Format_HeaderSBSDistrBySegAndFunc
{
	my ($wb, $format_LeftCA, $format_LeftCG, $format_LeftCT, $format_LeftTA, $format_LeftTC, $format_LeftTG, $format_RightCA, $format_RightCG, $format_RightCT, $format_RightTA, $format_RightTC, $format_RightTG) = @_;

	my ($bgColor_blue, $bgColor_black, $bgColor_red, $bgColor_gray, $bgColor_green, $bgColor_pink);
	BackgroundColor($wb, \$bgColor_blue, \$bgColor_black, \$bgColor_red, \$bgColor_gray, \$bgColor_green, \$bgColor_pink);

	$$format_LeftCA = $wb->add_format(bg_color=>$bgColor_blue, font=>'Arial', bold=>1, size=>11, color=>'white'); $$format_LeftCA->set_align('center'); $$format_LeftCA->set_left(2); $$format_LeftCA->set_left_color('blue');
	$$format_LeftCG = $wb->add_format(bg_color=>$bgColor_black, font=>'Arial', bold=>1, size=>11, color=>'white'); $$format_LeftCG->set_align('center'); $$format_LeftCG->set_left(2); $$format_LeftCG->set_left_color('blue');
	$$format_LeftCT = $wb->add_format(bg_color=>$bgColor_red, font=>'Arial', bold=>1, size=>11, color=>'white'); $$format_LeftCT->set_align('center'); $$format_LeftCT->set_left(2); $$format_LeftCT->set_left_color('blue');
	$$format_LeftTA = $wb->add_format(bg_color=>$bgColor_gray, font=>'Arial', bold=>1, size=>11, color=>'white'); $$format_LeftTA->set_align('center'); $$format_LeftTA->set_left(2); $$format_LeftTA->set_left_color('blue');
	$$format_LeftTC = $wb->add_format(bg_color=>$bgColor_green, font=>'Arial', bold=>1, size=>11, color=>'white'); $$format_LeftTC->set_align('center'); $$format_LeftTC->set_left(2); $$format_LeftTC->set_left_color('blue');
	$$format_LeftTG = $wb->add_format(bg_color=>$bgColor_pink, font=>'Arial', bold=>1, size=>11, color=>'white'); $$format_LeftTG->set_align('center'); $$format_LeftTG->set_left(2); $$format_LeftTG->set_left_color('blue');


	$$format_RightCA = $wb->add_format(bg_color=>$bgColor_blue, font=>'Arial', bold=>1, size=>11, color=>'white'); $$format_RightCA->set_align('center'); $$format_RightCA->set_right(2); $$format_RightCA->set_right_color('blue');
	$$format_RightCG = $wb->add_format(bg_color=>$bgColor_black, font=>'Arial', bold=>1, size=>11, color=>'white'); $$format_RightCG->set_align('center'); $$format_RightCG->set_right(2); $$format_RightCG->set_right_color('blue');
	$$format_RightCT = $wb->add_format(bg_color=>$bgColor_red, font=>'Arial', bold=>1, size=>11, color=>'white'); $$format_RightCT->set_align('center'); $$format_RightCT->set_right(2); $$format_RightCT->set_right_color('blue');
	$$format_RightTA = $wb->add_format(bg_color=>$bgColor_gray, font=>'Arial', bold=>1, size=>11, color=>'white'); $$format_RightTA->set_align('center'); $$format_RightTA->set_right(2); $$format_RightTA->set_right_color('blue');
	$$format_RightTC = $wb->add_format(bg_color=>$bgColor_green, font=>'Arial', bold=>1, size=>11, color=>'white'); $$format_RightTC->set_align('center'); $$format_RightTC->set_right(2); $$format_RightTC->set_right_color('blue');
	$$format_RightTG = $wb->add_format(bg_color=>$bgColor_pink, font=>'Arial', bold=>1, size=>11, color=>'white'); $$format_RightTG->set_align('center'); $$format_RightTG->set_right(2); $$format_RightTG->set_right_color('blue');
}
# Define the header for the part "Trinucleotide sequence context on the coding strand"
sub Format_Header12MutType
{
	my ($wb, $format_CA, $format_CG, $format_CT, $format_TA, $format_TC, $format_TG) = @_;

	my ($bgColor_blue, $bgColor_black, $bgColor_red, $bgColor_gray, $bgColor_green, $bgColor_pink);
	BackgroundColor($wb, \$bgColor_blue, \$bgColor_black, \$bgColor_red, \$bgColor_gray, \$bgColor_green, \$bgColor_pink);

	$$format_CA = $wb->add_format(bg_color=>$bgColor_blue, font=>'Arial', bold=>1, size=>11, color=>'white'); $$format_CA->set_align('center');
	$$format_CG = $wb->add_format(bg_color=>$bgColor_black, font=>'Arial', bold=>1, size=>11, color=>'white'); $$format_CG->set_align('center');
	$$format_CT = $wb->add_format(bg_color=>$bgColor_red, font=>'Arial', bold=>1, size=>11, color=>'white'); $$format_CT->set_align('center');
	$$format_TA = $wb->add_format(bg_color=>$bgColor_gray, font=>'Arial', bold=>1, size=>11, color=>'white'); $$format_TA->set_align('center');
	$$format_TC = $wb->add_format(bg_color=>$bgColor_green, font=>'Arial', bold=>1, size=>11, color=>'white'); $$format_TC->set_align('center');
	$$format_TG = $wb->add_format(bg_color=>$bgColor_pink, font=>'Arial', bold=>1, size=>11, color=>'white'); $$format_TG->set_align('center');
}
# Define the format for the text that needs a section border
sub Format_TextSection
{
	my ($wb, $formatT_left, $formatT_right, $formatT_bottomRight, $formatT_bottomLeft, $formatT_bottom, $formatT_bottomHeader, $formatT_bottomRightHeader, $formatT_bottomHeader2, $formatT_rightHeader) = @_;

	$$formatT_left = $wb->add_format(valign=>'center', font=>'Arial', size=>10);
	$$formatT_left->set_left(2); $$formatT_left->set_left_color('blue');

	$$formatT_right = $wb->add_format(valign=>'center', font=>'Arial', size=>10);
	$$formatT_right->set_right(2); $$formatT_right->set_right_color('blue');

	$$formatT_bottomRight = $wb->add_format(valign=>'center', font=>'Arial', size=>10);
	$$formatT_bottomRight->set_bottom(2); $$formatT_bottomRight->set_bottom_color('blue');
	$$formatT_bottomRight->set_right(2); $$formatT_bottomRight->set_right_color('blue');

	$$formatT_bottomLeft = $wb->add_format(valign=>'center', font=>'Arial', size=>10);
	$$formatT_bottomLeft->set_bottom(2); $$formatT_bottomLeft->set_bottom_color('blue');
	$$formatT_bottomLeft->set_left(2); $$formatT_bottomLeft->set_left_color('blue');

	$$formatT_bottom = $wb->add_format(valign=>'center', font=>'Arial', size=>10);
	$$formatT_bottom->set_bottom(2); $$formatT_bottom->set_bottom_color('blue');

	my $bgColor_totallighGray = $wb->set_custom_color(54, 230, 230, 230);
	$$formatT_bottomHeader = $wb->add_format(bg_color=>$bgColor_totallighGray, font=>'Arial', bold=>1, size=>11); $$formatT_bottomHeader->set_align('center');
	$$formatT_bottomHeader->set_bottom(2); $$formatT_bottomHeader->set_bottom_color('blue');

	$$formatT_bottomRightHeader = $wb->add_format(bg_color=>$bgColor_totallighGray, font=>'Arial', bold=>1, size=>11); $$formatT_bottomRightHeader->set_align('center');
	$$formatT_bottomRightHeader->set_bottom(2); $$formatT_bottomRightHeader->set_bottom_color('blue');
	$$formatT_bottomRightHeader->set_right(2); $$formatT_bottomRightHeader->set_right_color('blue');

	$$formatT_bottomHeader2 = $wb->add_format(bg_color=>$bgColor_totallighGray, font=>'Arial', bold=>1, size=>11); $$formatT_bottomHeader2->set_align('center');

	$$formatT_rightHeader = $wb->add_format(bg_color=>$bgColor_totallighGray, font=>'Arial', bold=>1, size=>11); $$formatT_rightHeader->set_align('center');
	$$formatT_rightHeader->set_right(2); $$formatT_rightHeader->set_right_color('blue');
}
# Define the format for the graphs titles
sub Format_GraphTitle
{
	my ($wb, $formatT_graphTitle) = @_;

	$$formatT_graphTitle = $wb->add_format(font=>'Arial', size=>12, bold=>1);
}
# Define the format of the border of the tables
sub Format_Table
{
	my ($wb, $table_topleft, $table_topRight, $table_bottomleft, $table_bottomRight, $table_top, $table_right, $table_bottom, $table_bottomItalicRed, $table_left, $table_bottomrightHeader, $table_left2, $table_middleHeader, $table_middleHeader2) = @_;

	$$table_topleft = $wb->add_format(valign=>'center', bold=>1, font=>'Arial', size=>10); $$table_topleft->set_top(1); $$table_topleft->set_left(1);
	$$table_topRight = $wb->add_format(valign=>'center', bold=>1, font=>'Arial', size=>10); $$table_topRight->set_top(1); $$table_topRight->set_right(1);
	$$table_bottomleft = $wb->add_format(valign=>'center', bold=>1, font=>'Arial', size=>10); $$table_bottomleft->set_bottom(1); $$table_bottomleft->set_left(1);
	$$table_bottomRight = $wb->add_format(valign=>'center', font=>'Arial', size=>10); $$table_bottomRight->set_bottom(1); $$table_bottomRight->set_right(1);

	$$table_top          = $wb->add_format(valign=>'center', bold=>1, font=>'Arial', size=>10);   $$table_top->set_top(1);
	$$table_right        = $wb->add_format(valign=>'center', font=>'Arial', size=>10);            $$table_right->set_right(1);
	$$table_bottom       = $wb->add_format(valign=>'center', font=>'Arial', size=>10);            $$table_bottom->set_bottom(1);
	$$table_bottomItalicRed = $wb->add_format(valign=>'center', font=>'Arial', size=>10, italic=>1, color => 'red'); $$table_bottomItalicRed->set_bottom(1);
		$$table_left         = $wb->add_format(valign=>'center', bold=>1, font=>'Arial', size=>10);   $$table_left->set_left(1);

	my $bgColor_totallighGray = $wb->set_custom_color(54, 230, 230, 230);
	$$table_bottomrightHeader = $wb->add_format(bg_color=>$bgColor_totallighGray, font=>'Arial', bold=>1, size=>10); $$table_bottomrightHeader->set_bottom(1); $$table_bottomrightHeader->set_right(1);

	$$table_left2 =  $wb->add_format(valign=>'left', font=>'Arial', size=>10); $$table_left2->set_left(1);

	$$table_middleHeader  = $wb->add_format(valign=>'center', bg_color=>$bgColor_totallighGray, font=>'Arial', bold=>1, size=>10);
	$$table_middleHeader2 = $wb->add_format(valign=>'center', bg_color=>$bgColor_totallighGray, font=>'Arial', bold=>1, size=>10); $$table_middleHeader2->set_bottom(1);
}
# Define the color
sub Color
{
	my ($wb, $blue, $black, $red, $gray, $green, $pink) = @_;

	$$blue     = $wb->set_custom_color(40, 0, 0, 204);# C:G>A:T in blue
	$$black    = $wb->set_custom_color(41, 0, 0, 0);# C:G>G:C in black
	$$red      = $wb->set_custom_color(42, 255, 0, 0);# C:G>T:A in red
	$$gray     = $wb->set_custom_color(43, 205, 205, 205); # T:A>A:T in light gray
	$$green    = $wb->set_custom_color(44, 0, 204, 51);# T:A>C:G in green
	$$pink     = $wb->set_custom_color(45, 255, 192, 203);# T:A>G:C in pink
}
sub BackgroundColor
{
	my ($wb, $bgColor_blue, $bgColor_black, $bgColor_red, $bgColor_gray, $bgColor_green, $bgColor_pink) = @_;

	$$bgColor_blue  = $wb->set_custom_color(48, 0, 0, 204);
	$$bgColor_black = $wb->set_custom_color(49, 0, 0, 0);
	$$bgColor_red   = $wb->set_custom_color(50, 255, 0, 0);
	$$bgColor_gray  = $wb->set_custom_color(51, 205, 205, 205);
	$$bgColor_green = $wb->set_custom_color(52, 0, 204, 51);
	$$bgColor_pink  = $wb->set_custom_color(53, 255, 192, 203);
}




=head1 NAME

mutSpec-Stat

=head1 SYNOPSIS

	mutSpecstat.pl [arguments] <query-file>

  <query-file>                                   a folder with one or multiple VCFs

  Arguments:
        -h,        --help                        print help message
        -m,        --man                         print complete documentation
        -v,        --verbose                     use verbose output
                   --refGenome                   the reference genome to use (human, mouse or rat genomes)
        -o,        --outfile <string>            output directory for the result. If none is specify the result will be write in the same directory as the input file
                   --temp <string>               the path for saving the temporary files
                   --pathSeqRefGenome            the path to the fasta reference sequences
                   --poolData                    generate the pool of all the samples (optional)
                   --reportSample                generate a report for each sample (optional)


Function: automatically run a pipeline and calculate various statistics on mutations

 Example: mutSpecstat.pl --refGenome hg19 --outfile output_directory --temp path_to_temporary_directory --pathRscript path_to_R_scripts --pathSeqRefGenome path_fasta_ref_seq --poolData --reportSample inputFolder

 Version: 02-2017 (February 2016)


=head1 OPTIONS

=over 8

=item B<--help>

print a brief usage message and detailed explanation of options.

=item B<--man>

print the complete manual of the program.

=item B<--verbose>

use verbose output.

=item B<--refGenome>

the reference genome to use, could be human, mouse or rat genomes.

=item B<--outfile>

the directory of output file names. If it is nor specify the same directory as the input file is used.

=item B<--temp>

the path for saving temporary files generated by the script.
If any is specify a temporary folder is created in the same directory where the script is running.
Deleted when the script is finish

=item B<--pathSeqRefGenome>

The path to the fasta reference sequences

=item B<--poolData only for the report>

calculate the statistics on the pool of all the data pass in input

=item B<--reportSample only for the report>

generate a report for each samples

=head1 DESCRIPTION

mutSpecstat is a perl script for calculated various statistics on mutations
An Excel report containing the mutation type distribution per functional region, the strand bias and the sequence context on genomic and coding sequence is created.
The different statistics are illustrated using ggplot2.

=cut
