#!/usr/bin/env perl

#-----------------------------------#
# Author: Maude                     #
# Script: mutspecStat.pl            #
# Last update: 01/12/15             #
#-----------------------------------#

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Basename; # my ($filename, $directories, $suffix) = fileparse($file, qr/\.[^.]*/);
use File::Path;
use Statistics::R;
use Spreadsheet::WriteExcel;

our ($verbose, $man, $help) = (0, 0, 0); # Parse options and print usage if there is a syntax error, or if usage was explicitly requested.
our ($refGenome, $output, $folder_temp, $path_R_Scripts, $path_SeqrefGenome) = ("empty", "empty", "empty", "empty", "empty");  # The reference genome to use; The path for saving the result; The path for saving the temporary files; The path to R scripts; The path to the fasta reference sequences
our ($poolData, $oneReportPerSample) = (2, 2); # If a folder is pass as input file pool all the data and generate the report on the pool and for each samples; # Generate one report for each samples


GetOptions('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'refGenome=s'=>\$refGenome, 'outfile|o=s' => \$output, 'pathTemporary|temp=s' => \$folder_temp, 'pathRscript=s' => \$path_R_Scripts, 'pathSeqRefGenome=s' => \$path_SeqrefGenome, 'poolData' => \$poolData, 'reportSample' => \$oneReportPerSample) or pod2usage(2);

our ($input) = @ARGV;

pod2usage(-verbose=>1, -exitval=>1, -output=>\*STDERR) if ($help);
pod2usage(-verbose=>2, -exitval=>1, -output=>\*STDERR) if ($man);
pod2usage(-verbose=>0, -exitval=>1, -output=>\*STDERR) if(@ARGV == 0); # No argument is pass to the command line print the usage of the script
pod2usage(-verbose=>0, -exitval=>1, -output=>\*STDERR) if(@ARGV == 2); # Only one argument is expected to be pass to @ARGV (the input)



######################################################################################################################################################
#																																			GLOBAL VARIABLES																															 #
######################################################################################################################################################
# Recover the current path
our $pwd = `pwd`;
chomp($pwd);

# Path to R scripts
our $pathRScriptTxnSB       = "$path_R_Scripts/R/transciptionalStrandBias.r";
our $pathRScriptMutSpectrum = "$path_R_Scripts/R/mutationSpectra_Galaxy.r";

our $folderMutAnalysis = "";
our @pathInput         = split("/", $input);

# Hash table with the length of each chromosomes
our %chromosomes;

######################################################################################################################################################
#																																								MAIN 																																 #
######################################################################################################################################################
# Check the presence of the flags and create the output and temp directories
CheckFlags();

# Retrieve chromosomes length
checkChrDir();


print "-----------------------------------------------------------------\n";
print "-----------------Report Mutational Analysis----------------------\n";
print"-----------------------------------------------------------------\n";

# First check if the file is annotated or not
CheckAnnotationFile($input);

# Calculate the statistics and generate the report
my @colInfoAV = qw(Chr Start Ref Alt);
ReportMutDist($input, $folderMutAnalysis, $folder_temp, \@colInfoAV, $refGenome);

# Remove the temporary directory
rmtree($folder_temp);


######################################################################################################################################################
#																																							FUNCTIONS																															 #
######################################################################################################################################################

# Check the presence of the flags and create the output and temp directories
sub CheckFlags
{
	# Check the reference genome
	if($refGenome eq "empty") { print STDERR "You forget to specify the name for the reference genome!!!\nPlease specify it with the flag --refGenome\n"; exit; }

	# If no output is specified write the result as the same place as the input file
	if($output eq "empty")
	{
		my $folderRes         = "";
		for(my $i=0; $i<$#pathInput; $i++) { $folderRes .= "$pathInput[$i]/"; }

		$folderMutAnalysis = "$folderRes/Mutational_Analysis";
		if(!-e $folderMutAnalysis) { mkdir($folderMutAnalysis) or die "$!: $folderMutAnalysis\n"; }
	}
	else
	{
		if(!-e $output) { mkdir($output) or die "$!: $output\n"; }

		$folderMutAnalysis      = "$output/Mutational_Analysis";
		if(!-e $folderMutAnalysis) { mkdir($folderMutAnalysis) or die "$!: $folderMutAnalysis\n"; }
	}

	# If no temp folder is specified write the result in the current path
	if($folder_temp eq "empty") { $folder_temp   = "$pwd/TEMP_MutationalAnalysis_$pathInput[$#pathInput]"; }
	if(!-e $folder_temp)        { mkdir($folder_temp) or die "$!: $folder_temp\n"; }

	# Check the path to the R scripts
	if($path_R_Scripts eq "empty") { print STDERR "You forget to specify the path for the R scripts!!!\nPlease specify it with the flag --pathRscript\n"; exit; }


	# The input is a folder
	if(-d $input) { foreach my $file (`ls $input`) { CheckLengthFilename("$input/$file"); } }
	# The input is one file
	else { CheckLengthFilename($input); }
}
# Check the length of the file, must be < 32 characters for the Excel sheet
sub CheckLengthFilename
{
	my ($inputFile) = @_;

	## Verify the name of file, must be <= 31 chars for the sheet name
	my ($filename, $directories, $suffix) = fileparse($inputFile, qr/\.[^.]*/);

	if(length($filename) > 31) { print STDERR "The file: $inputFile must be <= 31 chars\nPlease modify it before running the script\n"; exit; }
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

# Check if the file is annotated or not
sub CheckAnnotationFile
{
	my ($inputFile)  = @_;

	# A folder is pass in argument
	if(-d $inputFile)
	{
		foreach my $file (`ls $inputFile`)
		{
			chomp($file);

			open(F1, "$inputFile/$file") or die "$!: $inputFile/$file\n";
			my $search_header = <F1>; $search_header =~ s/[\r\n]+$//; my @tab_search_header = split("\t",$search_header);
			close F1;
			# The number of the column
		  my $value_of_column_NB  = "toto";
		  for(my $i=0; $i<=$#tab_search_header; $i++)
		  {
		    if($tab_search_header[$i] eq "Func.refGene") { $value_of_column_NB = $i; }
		  }
		  if($value_of_column_NB eq "toto") { print STDERR "Error the input file you specify is not annotated! $inputFile/$file !!!!\nPlease first annotate your file before trying to generate the report on mutation patterns\n"; exit; }
		}
	}
	else
	{
		open(F1, $inputFile) or die "$!: $inputFile\n";
		my $search_header = <F1>; $search_header =~ s/[\r\n]+$//; my @tab_search_header = split("\t",$search_header);
		close F1;
		# The number of the column
		my $value_of_column_NB  = "toto";
		for(my $i=0; $i<=$#tab_search_header; $i++)
		{
			if($tab_search_header[$i] eq "Func.refGene") { $value_of_column_NB = $i; }
		}
		if($value_of_column_NB eq "toto") { print STDERR "Error the input file you specify is not annotated! $inputFile !!!!\nPlease first annotate your file before trying to generate the report on mutation patterns\n"; exit; }
	}
}

# Calculate the statistics and generate the report
sub ReportMutDist
{
	our ($input, $output, $folder_temp, $refTab_colInfo, $refGenome) = @_;

	my @column = @$refTab_colInfo;

	our ($chr_name, $start_name, $ref_name, $alt_name) = split(/,/, join(',', @column)); # Separe each element

	our $func_name       = "Func.refGene";
	our $exonicFunc_name = "ExonicFunc.refGene";
	our $strand_name     = "Strand";
	our $context_name    = "context";

	my $folderFigure = "$output/Figures";
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
	############ Recover Annovar annotations (for having the save number of functional regions for each samples)
	my @tab_func = recoverAnnovarAnnotation($input, $func_name);
	if(@tab_func == 0) { print STDERR "Error the table for the functional region is empty!!!!! check $input $func_name\n"; exit; }

	############ Calculate the different statistics present in the report
	my %h_file   = ();
	CalculateStatistics(\%h_file, \@tab_func);

	############ Calculate the chi2 for the strand bias
	CalculateChi2(\%h_file, $folderChi2);

	############ Write the different statistics present in the report
	WriteStatistics(\%h_file, $#tab_func, $folderFigure, $folderChi2, $folderNMF);

	############ Create logo for studying the 10 flanking bases of the mutation
	CreateLogo(\%h_file, $folderWebLogo);


	################### Subroutines for generating the report for the mutational analysis
	sub recoverAnnovarAnnotation
	{
		my ($input, $AV_annotation) = @_;

		my %hash = ();

		# The input is a folder
		if(-d $input)
		{
			foreach my $file (`ls $input`)
			{
				$file =~ s/[\r\n]+$//;
				my $AV_annotation_value = recoverNumCol("$input/$file", $AV_annotation);

				open(F1, "$input/$file") or die "$!: $input/$file\n";
				my $header = <F1>;
				while(<F1>)
				{
					$_ =~ s/[\r\n]+$//;
					my @tab = split("\t", $_);

					# Some files can have an empty line at the end and WE DON'T WANT to consider it
					if(! defined $tab[0]) { next; }
					# Some func value are repeated and separated by ";"
					my $funcSegment = "";
					if($tab[$AV_annotation_value] =~ /;/) { my @temp = split(";", $tab[$AV_annotation_value]); $funcSegment = $temp[0]; }
					else { $funcSegment = $tab[$AV_annotation_value]; }

					$hash{$funcSegment} = "";
				}
				close F1;
			}
			my @tab_AVAnnotation = ();
			foreach my $k (sort keys %hash) { push(@tab_AVAnnotation, $k); }
			return @tab_AVAnnotation;
		}
		# The input is a file
		else
		{
			my $AV_annotation_value = recoverNumCol($input, $AV_annotation);

			open(F1, $input) or die "$!: $input\n";
			my $header = <F1>;
			while(<F1>)
			{
				$_ =~ s/[\r\n]+$//;
				my @tab = split("\t", $_);

				# Some func value are repeated and separated by ";"
				my $funcSegment = "";
				if($tab[$AV_annotation_value] =~ /;/) { my @temp = split(";", $tab[$AV_annotation_value]); $funcSegment = $temp[0]; }
				else { $funcSegment = $tab[$AV_annotation_value]; }

				$hash{$funcSegment} = "";
			}
			close F1;
			my @tab_AVAnnotation = ();
			foreach my$k (sort keys %hash) { push(@tab_AVAnnotation, $k); }
			return @tab_AVAnnotation;
		}
	}
	# Calculate the different statistics present in the report
	sub CalculateStatistics
	{
		my ($refH_file, $refT_func) = @_;

		my ($chr_value, $start_value, $ref_value, $alt_value, $func_value, $exonicFunc_value, $strand_value, $contextSeq_value) = ("", "", "", "", "", "", "", "", "", "");

		# If the input is a folder
		if(-d $input)
		{
			my $folderPool = "$folder_temp/Pool";
			if(!-e $folderPool) { mkdir($folderPool) or die "Can't create the directory $folderPool\n"; }

			# Copy each sample
			foreach my $file (`ls $input`) { chomp($file); system("cp $input/$file $folderPool"); }

			# Generate the pool of all the data
			if($poolData == 1)
			{
				my @listFile = `ls $input`;

				# For keeping the header only one time
				chomp($listFile[0]);
				system("cp $input/$listFile[0] $folderPool/Pool_Data.txt");

				open(OUT, ">>", "$folderPool/Pool_Data.txt") or die "$!: $folderPool/Pool_Data.txt\n";

				for(my $i=1; $i<=$#listFile; $i++)
				{
					chomp($listFile[$i]);
					open(F1, "$input/$listFile[$i]") or die "$!: $input/$listFile[$i]\n";
					my $header = <F1>;
					while(<F1>) { print OUT $_; }
					close F1;
				}
				close OUT;
			}

			foreach my $file (`ls $folderPool`)
			{
				chomp($file);
				############ Recover the number of the columns of interest
				$chr_value        = recoverNumCol("$folderPool/$file", $chr_name);
				$start_value      = recoverNumCol("$folderPool/$file", $start_name);
				$ref_value        = recoverNumCol("$folderPool/$file", $ref_name);
				$alt_value        = recoverNumCol("$folderPool/$file", $alt_name);
				$func_value       = recoverNumCol("$folderPool/$file", $func_name);
				$exonicFunc_value = recoverNumCol("$folderPool/$file", $exonicFunc_name);
				$strand_value     = recoverNumCol("$folderPool/$file", $strand_name);
				$contextSeq_value = recoverNumCol("$folderPool/$file", $context_name);
				############ Recover the number of the columns of interest

				############ Control the annotated file pass in argument
				## Check if the files have variants
				my $nbLines_originalFile = `wc -l $folderPool/$file`; $nbLines_originalFile =~ /(\d+) /;
				if($1==1) { print STDERR "\n\nNo line in the file $folderPool/$file\n\n"; exit; }
				## Check if there is variant with strand information. If not the rest of the script generates errors
				my $testFile = 0;
				CheckVariantReport("$folderPool/$file", $strand_value, \$testFile);
				if($testFile==0) { print STDERR "\n\nNo strand information for the file $folderPool/$file\n\n"; exit; }
				############ Control the annotated file pass in argument

				############ Calculate the statistics
				File2Hash("$folderPool/$file", $func_value, $exonicFunc_value, $chr_value, $ref_value, $alt_value, $strand_value, $contextSeq_value, $refH_file, $refT_func);
			}
		}
		# If the input is a file
		else
		{
			############ Recover the number of the columns of interest
			$chr_value        = recoverNumCol($input, $chr_name);
			$start_value      = recoverNumCol($input, $start_name);
			$ref_value        = recoverNumCol($input, $ref_name);
			$alt_value        = recoverNumCol($input, $alt_name);
			$func_value       = recoverNumCol($input, $func_name);
			$exonicFunc_value = recoverNumCol($input, $exonicFunc_name);
			$strand_value     = recoverNumCol($input, $strand_name);
			$contextSeq_value = recoverNumCol($input, $context_name);
			############ Recover the number of the columns of interest

			############ Control the annotated file pass in argument
			## Check if the files have variants
			my $nbLines_originalFile = `wc -l $input`; $nbLines_originalFile =~ /(\d+) /;
			if($1==1) { print STDERR "\n\nNo line in the file $input\n\n"; exit; }
			## Check if there is variant with strand information. If not the rest of the script generates errors
			my $testFile = 0;
			CheckVariantReport($input, $strand_value, \$testFile);
			if($testFile==0) { print STDERR "\n\nNo strand information for the file $input\n\n"; exit; }
			############ Control the annotated file pass in argument

			############ Calculate the statistics
			File2Hash($input, $func_value, $exonicFunc_value, $chr_value, $ref_value, $alt_value, $strand_value, $contextSeq_value, $refH_file, $refT_func);
		}
	}
	# Check if there is at least one variant with a strand information
	sub CheckVariantReport
	{
		my ($file, $strand_value, $refS_testFile) = @_;

		open(F1, $file) or die "$!: $file\n";
		my $header = <F1>;
		while(<F1>)
		{
			$_      =~ s/[\r\n]+$//;
			my @tab = split("\t", $_);

			if( ($tab[$strand_value] eq "+") || ($tab[$strand_value] eq "-") ) { $$refS_testFile++; }
		}
		close F1;
	}
	# Convert the annotated VCF into a hash table
	sub File2Hash
	{
		my ($inputFile, $func_value, $exonicFunc_value, $chr_value, $ref_value, $alt_value, $strand_value, $contextSeq_value, $refH_file, $refT_func) = @_;
		my ($filename, $directories, $suffix) = fileparse($inputFile, qr/\.[^.]*/);

		# Initialisation of the hash
		my @tab_mutation   = qw(C:G>A:T C:G>G:C C:G>T:A T:A>A:T T:A>C:G T:A>G:C);
		my @tab_aaChange   = ("NonTr", "Tr", "TotalMutG");
		my @tab_humanChrom = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y);
		my @tab_mouseChrom = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y);
		my @tabExoFunc     = ("frameshift insertion", "frameshift deletion", "frameshift block substitution", "frameshift substitution", "stopgain", "stoploss", "nonframeshift insertion", "nonframeshift deletion", "nonframeshift substitution", "nonframeshift block substitution", "nonsynonymous SNV", "synonymous SNV", "unknown", "NA");

		# Total number of SBS on the genomic strand
		$refH_file->{$filename}{'TotalSBSGenomic'} = 0;
		# Total number of Indel on the genomic strand
		$refH_file->{$filename}{'TotalIndelGenomic'} = 0;
		# Total number of SBS on the coding strand
		$refH_file->{$filename}{'TotalSBSCoding'} = 0;
		# Total number of SBS and Indel on the genomic strand
		$refH_file->{$filename}{'TotalMutGenomic'} = 0;

		#####################################
		# SBS by segment (6 mutation types) #
		#####################################
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

	  #######################
	  # Pearson correlation #
	  #######################
	  $refH_file->{$filename}{'SBSPerChr'}{'AllMutType'} = 0;
	  # Count of SBS per chromosome foreach mutation types
	  foreach my $elt_tabMutation (@tab_mutation)
	  {
	  	  foreach my $chromosome (sort keys %chromosomes){ $refH_file->{$filename}{'SBSPerChr'}{$elt_tabMutation}{'CHR'}{$chromosome}{'chr'} = 0;}
		  $refH_file->{$filename}{'SBSPerChr'}{$elt_tabMutation}{'Pearson'} = 0;
	  }

	  foreach my $chromosome (sort keys %chromosomes){
	  	$refH_file->{$filename}{'SBSPerChr'}{'TotalPerChr'}{$chromosome}{'chr'}=0;
	  }

	  ############################
	  # Impact of SBS on protein #
	  ############################
	  foreach my $elt_exoFunc (@tabExoFunc)
	  {
	  	$refH_file->{$filename}{'ImpactSBS'}{$elt_exoFunc} = 0;
	  }

	  #####################################
	  # Sequence context (genomic strand) #
	  #####################################
	  my @tab_mutation2 = qw(C>A C>G C>T T>A T>C T>G);
	  my @tab_context   = qw(A_A A_C A_G A_T C_A C_C C_G C_T G_A G_C G_G G_T T_A T_C T_G T_T);
	  foreach my $elt_context (@tab_context)
	  {
	  	foreach my $elt_mutation3 (@tab_mutation2)
	  	{
	  		$refH_file->{$filename}{'SeqContextG'}{$elt_context}{$elt_mutation3} = 0;
	  	}
	  }

	  ####################################
	  # Sequence context (coding strand) #
	  ####################################
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

	  open(F1,$inputFile) or die "$!: $inputFile\n";
	  my $header = <F1>;
	  while(<F1>)
	  {
			$_      =~ s/[\r\n]+$//;
			my @tab = split("\t", $_);

			# Random chromosome and chromosome M
			if( ($tab[$chr_value] =~ /random/i) || ($tab[$chr_value] =~ /M/i) ) { next; }

			############################################## Extract the base just before and after the mutation ##############################################
			my $context                   = "";
			my $contextSequence           = $tab[$contextSeq_value]; $contextSequence =~ tr/a-z/A-Z/;
			my @tempContextSequence       = split("", $contextSequence);
			my $total_nbBaseContext       = $#tempContextSequence;
			my $midlle_totalNbBaseContext = $total_nbBaseContext/2; # For having the middle of the sequence
			my $before                    = $midlle_totalNbBaseContext - 1; my $after = $midlle_totalNbBaseContext + 1;
			$context                      = $tempContextSequence[$before]."_".$tempContextSequence[$after];
			############################################## Extract the base just before and after the mutation ##############################################


			############################################################### Impact on protein ###############################################################
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
			############################################################### Impact on protein ###############################################################

			################################################### Only SBS are considered for the statistics ##################################################
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
			################################################### Only SBS are considered for the statistics ##################################################

			# Number of SBS per chromosome: remove the "chr"
			my $chrNameForH=$tab[$chr_value];
			if(exists $refH_file->{$filename}{'SBSPerChr'}{'TotalPerChr'}{$chrNameForH}{'chr'}) { $refH_file->{$filename}{'SBSPerChr'}{'TotalPerChr'}{$chrNameForH}{'chr'}++; }


			################################################### Some func value are repeated and separated by ";" ##################################################
			my $funcSegment = "";
			if($tab[$func_value] =~ /;/) { my @temp = split(";", $tab[$func_value]); $funcSegment = $temp[0]; }
			else { $funcSegment = $tab[$func_value]; }


			############################################################### MUTATION C> #############################################
	    ###################################### C:G>A:T
	    if( (($tab[$ref_value] eq "C") && ($tab[$alt_value] eq "A")) || ( ($tab[$ref_value] eq "G") && ($tab[$alt_value] eq "T") ) )
	    {
	    	my $mutation   = "C:G>A:T";
	    	$refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'TotalMutG'}++; # Count the total number of mutations

	    	# Pearson correlation
	    	if(exists $refH_file->{$filename}{'SBSPerChr'}{$mutation}{'CHR'}{$chrNameForH}{'chr'}) { $refH_file->{$filename}{'SBSPerChr'}{$mutation}{'CHR'}{$chrNameForH}{'chr'}++; }

	    	# Sequence context - 6 mutation types - genomic strand
	    	my $mutationSeqContext6mutType = "C>A";
	    	# We want to express the mutation in C>
	    	if( ($tab[$ref_value] eq "G") && ($tab[$alt_value] eq "T") )
	    	{
	    		my $base3 = complement($tempContextSequence[$before]); my $base5 = complement($tempContextSequence[$after]);
	    		my $context_reverse = $base5."_".$base3;
	    		if(exists $refH_file->{$filename}{'SeqContextG'}{$context_reverse}{$mutationSeqContext6mutType}) { $refH_file->{$filename}{'SeqContextG'}{$context_reverse}{$mutationSeqContext6mutType}++; }
	    	}
	    	elsif(exists $refH_file->{$filename}{'SeqContextG'}{$context}{$mutationSeqContext6mutType}) { $refH_file->{$filename}{'SeqContextG'}{$context}{$mutationSeqContext6mutType}++; }

	    	# Strand analysis C>A on NonTr strand
	    	if( (($tab[$strand_value] eq "+") && (($tab[$ref_value] eq "C")&&($tab[$alt_value] eq "A"))) || (($tab[$strand_value] eq "-") && (($tab[$ref_value] eq "G")&&($tab[$alt_value] eq "T"))) )
	    	{
	    		if(exists $refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'NonTr'}) { $refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'NonTr'}++; }

	    		# C>A With the sequence context (C>A strand  = +)
	    		if( ($tab[$strand_value] eq "+") && (($tab[$ref_value] eq "C")&&($tab[$alt_value] eq "A")) )
	    		{
	    			if(exists $refH_file->{$filename}{'SeqContextC'}{$context}{'C>A'}{'NonTr'}) { $refH_file->{$filename}{'SeqContextC'}{$context}{'C>A'}{'NonTr'}++; }
	    		}
	    		# C>A With the sequence context (G>T strand = -)
	    		else
	    		{
	    			my $base3 = complement($tempContextSequence[$before]); my $base5 = complement($tempContextSequence[$after]);
	    			my $context_reverse = $base5."_".$base3;
	    			if(exists $refH_file->{$filename}{'SeqContextC'}{$context_reverse}{'C>A'}{'NonTr'}) { $refH_file->{$filename}{'SeqContextC'}{$context_reverse}{'C>A'}{'NonTr'}++; }
	    		}
	    	}
	    	# Strand analysis C>A on Tr strand
	    	if( (($tab[$strand_value] eq "-") && (($tab[$ref_value] eq "C")&&($tab[$alt_value] eq "A"))) || (($tab[$strand_value] eq "+") && (($tab[$ref_value] eq "G")&&($tab[$alt_value] eq "T"))) )
	    	{
	    		if(exists $refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'Tr'}) { $refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'Tr'}++; }

	    		# C>A With the sequence context (C>A strand = -)
	    		if( ($tab[$strand_value] eq "-") && (($tab[$ref_value] eq "C")&&($tab[$alt_value] eq "A")) )
	    		{
	    			if(exists $refH_file->{$filename}{'SeqContextC'}{$context}{'C>A'}{'Tr'}) { { $refH_file->{$filename}{'SeqContextC'}{$context}{'C>A'}{'Tr'}++; } }
	    		}
	    		# C>A with the sequence context (G>T strand = +)
	    		if( ($tab[$strand_value] eq "+") && (($tab[$ref_value] eq "G")&&($tab[$alt_value] eq "T")) )
	    		{
	    			my $base3 = complement($tempContextSequence[$before]); my $base5 = complement($tempContextSequence[$after]);
	    			my $context_reverse = $base5."_".$base3;
	    			if(exists $refH_file->{$filename}{'SeqContextC'}{$context_reverse}{'C>A'}{'Tr'}) { $refH_file->{$filename}{'SeqContextC'}{$context_reverse}{'C>A'}{'Tr'}++; }
	    		}
	    	}
	    	# WebLogo-3
	    	if(($tab[$ref_value] eq "C") && ($tab[$alt_value] eq "A"))
	    	{
	    		# For the logo all the sequences must have the same length
	    		if(scalar(@tempContextSequence) == 2)  { next; }
	    		my ($contextTemp1, $contextTemp2) = ("", "");
					for(my $i=0; $i<$midlle_totalNbBaseContext; $i++) { $contextTemp1 .= $tempContextSequence[$i]; }
					for(my $i=$midlle_totalNbBaseContext+1; $i<=$#tempContextSequence; $i++) { $contextTemp2 .= $tempContextSequence[$i]; }
					my $context = $contextTemp1."C".$contextTemp2;
					push(@{$refH_file->{$filename}{'WebLogo3'}{'CA'}}, $context);
	    	}
	    	else
	    	{
	    		if(scalar(@tempContextSequence) == 2)  { next; }
	    		my ($contextTemp1, $contextTemp2) = ("", "");
					for(my $i=0; $i<$midlle_totalNbBaseContext; $i++) { $contextTemp1 .= complement($tempContextSequence[$i]); }
					for(my $i=$midlle_totalNbBaseContext+1; $i<=$#tempContextSequence; $i++) { $contextTemp2 .= complement($tempContextSequence[$i]); }
					my $context = $contextTemp1."C".$contextTemp2; $context = reverse $context;
					push(@{$refH_file->{$filename}{'WebLogo3'}{'CA'}}, $context);
	    	}
	    }
	    ###################################### C:G>G:C
	   	if( (($tab[$ref_value] eq "C") && ($tab[$alt_value] eq "G")) || ( ($tab[$ref_value] eq "G") && ($tab[$alt_value] eq "C") ) )
	    {
	    	my $mutation   = "C:G>G:C";
	    	$refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'TotalMutG'}++; # Count the total number of mutations

		    # Pearson correlation
	    	if(exists $refH_file->{$filename}{'SBSPerChr'}{$mutation}{'CHR'}{$chrNameForH}{'chr'}) { $refH_file->{$filename}{'SBSPerChr'}{$mutation}{'CHR'}{$chrNameForH}{'chr'}++; }

	    	# Sequence context - 6 mutation types - genomic strand
	    	my $mutationSeqContext6mutType = "C>G";
	    	# We want to express the mutation in C>
	    	if( ($tab[$ref_value] eq "G") && ($tab[$alt_value] eq "C") )
	    	{
	    		my $base3 = complement($tempContextSequence[$before]); my $base5 = complement($tempContextSequence[$after]);
	    		my $context_reverse = $base5."_".$base3;
	    		if(exists $refH_file->{$filename}{'SeqContextG'}{$context_reverse}{$mutationSeqContext6mutType}) { $refH_file->{$filename}{'SeqContextG'}{$context_reverse}{$mutationSeqContext6mutType}++; }
	    	}
	    	elsif(exists $refH_file->{$filename}{'SeqContextG'}{$context}{$mutationSeqContext6mutType}) { $refH_file->{$filename}{'SeqContextG'}{$context}{$mutationSeqContext6mutType}++; }

	    	# Strand analysis C>G on NonTr strand
	    	if( (($tab[$strand_value] eq "+") && (($tab[$ref_value] eq "C")&&($tab[$alt_value] eq "G"))) || (($tab[$strand_value] eq "-") && (($tab[$ref_value] eq "G")&&($tab[$alt_value] eq "C"))) )
	    	{
	    		if(exists $refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'NonTr'}) { $refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'NonTr'}++; }

	    		# C>G with the sequence context (C>G strand = +)
	    		if( ($tab[$strand_value] eq "+") && (($tab[$ref_value] eq "C")&&($tab[$alt_value] eq "G")) )
	    		{
	    			if(exists $refH_file->{$filename}{'SeqContextC'}{$context}{'C>G'}{'NonTr'}) { $refH_file->{$filename}{'SeqContextC'}{$context}{'C>G'}{'NonTr'}++; }
	    		}
	    		# C>G with the sequence context (G>C strand = -)
	    		if( ($tab[$strand_value] eq "-") && (($tab[$ref_value] eq "G")&&($tab[$alt_value] eq "C")) )
	    		{
	    			my $base3           = complement($tempContextSequence[$before]); my $base5 = complement($tempContextSequence[$after]);
	    			my $context_reverse = $base5."_".$base3;
	    			if(exists $refH_file->{$filename}{'SeqContextC'}{$context_reverse}{'C>G'}{'NonTr'}) { $refH_file->{$filename}{'SeqContextC'}{$context_reverse}{'C>G'}{'NonTr'}++; }
	    		}
	    	}
	    	# Strand analysis C>G on Tr strand
	    	if( (($tab[$strand_value] eq "-") && (($tab[$ref_value] eq "C")&&($tab[$alt_value] eq "G"))) || (($tab[$strand_value] eq "+") && (($tab[$ref_value] eq "G")&&($tab[$alt_value] eq "C"))) )
	    	{
	    		if(exists $refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'Tr'}) { $refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'Tr'}++; }

	    		# C>G with the sequence context (C>G strand = -)
	    		if( ($tab[$strand_value] eq "-") && (($tab[$ref_value] eq "C")&&($tab[$alt_value] eq "G")) )
	    		{
	    			if(exists $refH_file->{$filename}{'SeqContextC'}{$context}{'C>G'}{'Tr'}) { $refH_file->{$filename}{'SeqContextC'}{$context}{'C>G'}{'Tr'}++; }
	    		}
	    		# C>G with the sequence context (G>C strand = +)
	    		if( ($tab[$strand_value] eq "+") && (($tab[$ref_value] eq "G")&&($tab[$alt_value] eq "C")) )
	    		{
	    			my $base3           = complement($tempContextSequence[$before]); my $base5 = complement($tempContextSequence[$after]);
	    			my $context_reverse = $base5."_".$base3;
	    			if(exists $refH_file->{$filename}{'SeqContextC'}{$context_reverse}{'C>G'}{'Tr'}) { $refH_file->{$filename}{'SeqContextC'}{$context_reverse}{'C>G'}{'Tr'}++; }
	    		}
	    	}
	    	# WebLogo-3
	    	if(($tab[$ref_value] eq "C") && ($tab[$alt_value] eq "G"))
	    	{
	    		if(scalar(@tempContextSequence) == 2)  { next; }
	    		my ($contextTemp1, $contextTemp2) = ("", "");
					for(my $i=0; $i<$midlle_totalNbBaseContext; $i++) { $contextTemp1 .= $tempContextSequence[$i]; }
					for(my $i=$midlle_totalNbBaseContext+1; $i<=$#tempContextSequence; $i++) { $contextTemp2 .= $tempContextSequence[$i]; }
					my $context = $contextTemp1."C".$contextTemp2;
					push(@{$refH_file->{$filename}{'WebLogo3'}{'CG'}}, $context);
	    	}
	    	else
	    	{
	    		if(scalar(@tempContextSequence) == 2)  { next; }
	    		my ($contextTemp1, $contextTemp2) = ("", "");
					for(my $i=0; $i<$midlle_totalNbBaseContext; $i++) { $contextTemp1 .= complement($tempContextSequence[$i]); }
					for(my $i=$midlle_totalNbBaseContext+1; $i<=$#tempContextSequence; $i++) { $contextTemp2 .= complement($tempContextSequence[$i]); }
					my $context = $contextTemp1."C".$contextTemp2; $context = reverse $context;
					push(@{$refH_file->{$filename}{'WebLogo3'}{'CG'}}, $context);
	    	}
	    }
	    ###################################### C:G>T:A
	   	if( (($tab[$ref_value] eq "C") && ($tab[$alt_value] eq "T")) || ( ($tab[$ref_value] eq "G") && ($tab[$alt_value] eq "A") ) )
	    {
	    	my $mutation   = "C:G>T:A";
	    	$refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'TotalMutG'}++; # Count the total number of mutations

	    	# Pearson correlation
	    	if(exists $refH_file->{$filename}{'SBSPerChr'}{$mutation}{'CHR'}{$chrNameForH}{'chr'}) { $refH_file->{$filename}{'SBSPerChr'}{$mutation}{'CHR'}{$chrNameForH}{'chr'}++; }

	    	# Sequence context - 6 mutation types - genomic strand
	    	my $mutationSeqContext6mutType = "C>T";
	    	# We want to express the mutation in C>
	    	if( ($tab[$ref_value] eq "G") && ($tab[$alt_value] eq "A") )
	    	{
	    		my $base3 = complement($tempContextSequence[$before]); my $base5 = complement($tempContextSequence[$after]);
	    		my $context_reverse = $base5."_".$base3;
	    		if(exists $refH_file->{$filename}{'SeqContextG'}{$context_reverse}{$mutationSeqContext6mutType}) { $refH_file->{$filename}{'SeqContextG'}{$context_reverse}{$mutationSeqContext6mutType}++; }
	    	}
	    	elsif(exists $refH_file->{$filename}{'SeqContextG'}{$context}{$mutationSeqContext6mutType}) { $refH_file->{$filename}{'SeqContextG'}{$context}{$mutationSeqContext6mutType}++; }

	    	# Strand analysis C>T on NonTr strand
	    	if( (($tab[$strand_value] eq "+") && (($tab[$ref_value] eq "C")&&($tab[$alt_value] eq "T"))) || (($tab[$strand_value] eq "-") && (($tab[$ref_value] eq "G")&&($tab[$alt_value] eq "A"))) )
	    	{
	    		if(exists $refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'NonTr'}) { $refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'NonTr'}++; }

	    		# C>T with the sequence context (C>T strand = +)
	    		if( ($tab[$strand_value] eq "+") && (($tab[$ref_value] eq "C")&&($tab[$alt_value] eq "T")) )
	    		{
	    			if(exists $refH_file->{$filename}{'SeqContextC'}{$context}{'C>T'}{'NonTr'}) { $refH_file->{$filename}{'SeqContextC'}{$context}{'C>T'}{'NonTr'}++; }
	    		}
	    		# C>T with the sequence context (G>A strand = -)
	    		if( ($tab[$strand_value] eq "-") && (($tab[$ref_value] eq "G")&&($tab[$alt_value] eq "A")) )
	    		{
	    			my $base3           = complement($tempContextSequence[$before]); my $base5 = complement($tempContextSequence[$after]);
	    			my $context_reverse = $base5."_".$base3;
	    			if(exists $refH_file->{$filename}{'SeqContextC'}{$context_reverse}{'C>T'}{'NonTr'}) { $refH_file->{$filename}{'SeqContextC'}{$context_reverse}{'C>T'}{'NonTr'}++; }
	    		}
	    	}
	    	# Strand analysis C>T on Tr strand
	    	if( (($tab[$strand_value] eq "-") && (($tab[$ref_value] eq "C")&&($tab[$alt_value] eq "T"))) || (($tab[$strand_value] eq "+") && (($tab[$ref_value] eq "G")&&($tab[$alt_value] eq "A"))) )
	    	{
	    		if(exists $refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'Tr'}) { $refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'Tr'}++; }

	    		# C>T with the sequence context (C>T strand = -)
	    		if( ($tab[$strand_value] eq "-") && (($tab[$ref_value] eq "C")&&($tab[$alt_value] eq "T")) )
	    		{
	    			if(exists $refH_file->{$filename}{'SeqContextC'}{$context}{'C>T'}{'Tr'}) { $refH_file->{$filename}{'SeqContextC'}{$context}{'C>T'}{'Tr'}++; }
	    		}
	    		# C>T with the sequence context (G>A strand = +)
	    		if( ($tab[$strand_value] eq "+") && (($tab[$ref_value] eq "G")&&($tab[$alt_value] eq "A")) )
	    		{
	    			my $base3           = complement($tempContextSequence[$before]); my $base5 = complement($tempContextSequence[$after]);
	    			my $context_reverse = $base5."_".$base3;
	    			if(exists $refH_file->{$filename}{'SeqContextC'}{$context_reverse}{'C>T'}{'Tr'}) { $refH_file->{$filename}{'SeqContextC'}{$context_reverse}{'C>T'}{'Tr'}++; }
	    		}
	    	}
	    	# WebLogo-3
	    	if(($tab[$ref_value] eq "C") && ($tab[$alt_value] eq "T"))
	    	{
	    		if(scalar(@tempContextSequence) == 2)  { next; }
	    		my ($contextTemp1, $contextTemp2) = ("", "");
					for(my $i=0; $i<$midlle_totalNbBaseContext; $i++) { $contextTemp1 .= $tempContextSequence[$i]; }
					for(my $i=$midlle_totalNbBaseContext+1; $i<=$#tempContextSequence; $i++) { $contextTemp2 .= $tempContextSequence[$i]; }
					my $context = $contextTemp1."C".$contextTemp2;
					push(@{$refH_file->{$filename}{'WebLogo3'}{'CT'}}, $context);
	    	}
	    	else
	    	{
	    		if(scalar(@tempContextSequence) == 2)  { next; }
	    		my ($contextTemp1, $contextTemp2) = ("", "");
					for(my $i=0; $i<$midlle_totalNbBaseContext; $i++) { $contextTemp1 .= complement($tempContextSequence[$i]); }
					for(my $i=$midlle_totalNbBaseContext+1; $i<=$#tempContextSequence; $i++) { $contextTemp2 .= complement($tempContextSequence[$i]); }
					my $context = $contextTemp1."C".$contextTemp2; $context = reverse $context;
					push(@{$refH_file->{$filename}{'WebLogo3'}{'CT'}}, $context);
	    	}
	    }

	    ############################################################### MUTATION T> #############################################
	    ###################################### T:A>A:T
	    if( (($tab[$ref_value] eq "T") && ($tab[$alt_value] eq "A")) || ( ($tab[$ref_value] eq "A") && ($tab[$alt_value] eq "T") ) )
	    {
	    	my $mutation   = "T:A>A:T";
	    	$refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'TotalMutG'}++; # Count the total number of mutations

	    	# Pearson correlation
	    	if(exists $refH_file->{$filename}{'SBSPerChr'}{$mutation}{'CHR'}{$chrNameForH}{'chr'}) { $refH_file->{$filename}{'SBSPerChr'}{$mutation}{'CHR'}{$chrNameForH}{'chr'}++; }

	    	# Sequence context - 6 mutation types - genomic strand
	    	my $mutationSeqContext6mutType = "T>A";
	    	# We want to express the mutation in T>
	    	if( ($tab[$ref_value] eq "A") && ($tab[$alt_value] eq "T") )
	    	{
	    		my $base3 = complement($tempContextSequence[$before]); my $base5 = complement($tempContextSequence[$after]);
	    		my $context_reverse = $base5."_".$base3;
	    		if(exists $refH_file->{$filename}{'SeqContextG'}{$context_reverse}{$mutationSeqContext6mutType}) { $refH_file->{$filename}{'SeqContextG'}{$context_reverse}{$mutationSeqContext6mutType}++; }
	    	}
	    	elsif(exists $refH_file->{$filename}{'SeqContextG'}{$context}{$mutationSeqContext6mutType}) { $refH_file->{$filename}{'SeqContextG'}{$context}{$mutationSeqContext6mutType}++; }

	    	# Strand analysis T>A on NonTr stand
	    	if( (($tab[$strand_value] eq "+") && (($tab[$ref_value] eq "T")&&($tab[$alt_value] eq "A"))) || (($tab[$strand_value] eq "-") && (($tab[$ref_value] eq "A")&&($tab[$alt_value] eq "T"))) )
	    	{
	    		if(exists $refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'NonTr'}) { $refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'NonTr'}++; }

	    		# T>A with the sequence context (T>A strand = +)
	    		if( ($tab[$strand_value] eq "+") && (($tab[$ref_value] eq "T")&&($tab[$alt_value] eq "A")) )
	    		{
	    			if(exists $refH_file->{$filename}{'SeqContextC'}{$context}{'T>A'}{'NonTr'}) { $refH_file->{$filename}{'SeqContextC'}{$context}{'T>A'}{'NonTr'}++; }
	    		}
	    		# T>A with the sequence context (A>T strand = -)
	    		else
	    		{
	    			my $base3           = complement($tempContextSequence[$before]); my $base5 = complement($tempContextSequence[$after]);
	    			my $context_reverse = $base5."_".$base3;
	    			if(exists $refH_file->{$filename}{'SeqContextC'}{$context_reverse}{'T>A'}{'NonTr'}) { $refH_file->{$filename}{'SeqContextC'}{$context_reverse}{'T>A'}{'NonTr'}++; }
	    		}
	    	}
	    	# Strand analysis T>A on Tr strand
	    	if( (($tab[$strand_value] eq "-") && (($tab[$ref_value] eq "T")&&($tab[$alt_value] eq "A"))) || (($tab[$strand_value] eq "+") && (($tab[$ref_value] eq "A")&&($tab[$alt_value] eq "T"))) )
	    	{
	    		if(exists $refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'Tr'}) { $refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'Tr'}++; }

	    		# T>A <ith the sequence context (T>A strand = -)
	    		if( ($tab[$strand_value] eq "-") && (($tab[$ref_value] eq "T")&&($tab[$alt_value] eq "A")) )
	    		{
	    			if(exists $refH_file->{$filename}{'SeqContextC'}{$context}{'T>A'}{'Tr'}) { $refH_file->{$filename}{'SeqContextC'}{$context}{'T>A'}{'Tr'}++; }
	    		}
	 				# T>A with the sequence context (A>T strand = +)
	    		else
	    		{
	    			my $base3           = complement($tempContextSequence[$before]); my $base5 = complement($tempContextSequence[$after]);
	    			my $context_reverse = $base5."_".$base3;
	    			if(exists $refH_file->{$filename}{'SeqContextC'}{$context_reverse}{'T>A'}{'Tr'}) { $refH_file->{$filename}{'SeqContextC'}{$context_reverse}{'T>A'}{'Tr'}++; }
	    		}
	    	}
	    	# WebLogo-3
	    	if(($tab[$ref_value] eq "T") && ($tab[$alt_value] eq "A"))
	    	{
	    		if(scalar(@tempContextSequence) == 2)  { next; }
	    		my ($contextTemp1, $contextTemp2) = ("", "");
					for(my $i=0; $i<$midlle_totalNbBaseContext; $i++)                        { $contextTemp1 .= $tempContextSequence[$i]; }
					for(my $i=$midlle_totalNbBaseContext+1; $i<=$#tempContextSequence; $i++) { $contextTemp2 .= $tempContextSequence[$i]; }
					my $context = $contextTemp1."T".$contextTemp2;
					push(@{$refH_file->{$filename}{'WebLogo3'}{'TA'}}, $context);
	    	}
	    	else
	    	{
	    		if(scalar(@tempContextSequence) == 2)  { next; }
	    		my ($contextTemp1, $contextTemp2) = ("", "");
					for(my $i=0; $i<$midlle_totalNbBaseContext; $i++) { $contextTemp1 .= complement($tempContextSequence[$i]); }
					for(my $i=$midlle_totalNbBaseContext+1; $i<=$#tempContextSequence; $i++) { $contextTemp2 .= complement($tempContextSequence[$i]); }
					my $context = $contextTemp1."T".$contextTemp2; $context = reverse $context;
					push(@{$refH_file->{$filename}{'WebLogo3'}{'TA'}}, $context);
	    	}
	    }
	    ###################################### T:A>C:G
	    if( (($tab[$ref_value] eq "T") && ($tab[$alt_value] eq "C")) || ( ($tab[$ref_value] eq "A") && ($tab[$alt_value] eq "G")) )
	    {
	    	my $mutation   = "T:A>C:G";
	    	$refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'TotalMutG'}++; # Count the total number of mutations

	    	# Pearson correlation
	    	if(exists $refH_file->{$filename}{'SBSPerChr'}{$mutation}{'CHR'}{$chrNameForH}{'chr'}) { $refH_file->{$filename}{'SBSPerChr'}{$mutation}{'CHR'}{$chrNameForH}{'chr'}++; }

	    	# Sequence context - 6 mutation types - genomic strand
	    	my $mutationSeqContext6mutType = "T>C";
	    	# We want to express the mutation in T>
	    	if( ($tab[$ref_value] eq "A") && ($tab[$alt_value] eq "T") )
	    	{
	    		my $base3 = complement($tempContextSequence[$before]); my $base5 = complement($tempContextSequence[$after]);
	    		my $context_reverse = $base5."_".$base3;
	    		if(exists $refH_file->{$filename}{'SeqContextG'}{$context_reverse}{$mutationSeqContext6mutType}) { $refH_file->{$filename}{'SeqContextG'}{$context_reverse}{$mutationSeqContext6mutType}++; }
	    	}
	    	elsif(exists $refH_file->{$filename}{'SeqContextG'}{$context}{$mutationSeqContext6mutType}) { $refH_file->{$filename}{'SeqContextG'}{$context}{$mutationSeqContext6mutType}++; }

	    	# Strand analysis T>C on NonTr strand
	    	if( (($tab[$strand_value] eq "+") && (($tab[$ref_value] eq "T")&&($tab[$alt_value] eq "C"))) || (($tab[$strand_value] eq "-") && (($tab[$ref_value] eq "A")&&($tab[$alt_value] eq "G"))) )
	    	{
	    		if(exists $refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'NonTr'}) { $refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'NonTr'}++; }

	    		# T>C (T>C strand = +)
	    		if( ($tab[$strand_value] eq "+") && (($tab[$ref_value] eq "T")&&($tab[$alt_value] eq "C")) )
	    		{
	    			if(exists $refH_file->{$filename}{'SeqContextC'}{$context}{'T>C'}{'NonTr'}) { $refH_file->{$filename}{'SeqContextC'}{$context}{'T>C'}{'NonTr'}++; }
	    		}
	    		# T>C (A>G strand = -)
	    		else
	    		{
	    			my $base3           = complement($tempContextSequence[$before]); my $base5 = complement($tempContextSequence[$after]);
	    			my $context_reverse = $base5."_".$base3;
	    			if(exists $refH_file->{$filename}{'SeqContextC'}{$context_reverse}{'T>C'}{'NonTr'}) { $refH_file->{$filename}{'SeqContextC'}{$context_reverse}{'T>C'}{'NonTr'}++; }
	    		}
	    	}
	    	# Strand analysis T>C on Tr strand
	    	if( (($tab[$strand_value] eq "-") && (($tab[$ref_value] eq "T")&&($tab[$alt_value] eq "C"))) || (($tab[$strand_value] eq "+") && (($tab[$ref_value] eq "A")&&($tab[$alt_value] eq "G"))) )
	    	{
	    		if(exists $refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'Tr'}) { $refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'Tr'}++; }

	    		# T>C (T>C strand = -)
	    		if( ($tab[$strand_value] eq "-") && (($tab[$ref_value] eq "T")&&($tab[$alt_value] eq "C")) )
	    		{
	    			if(exists $refH_file->{$filename}{'SeqContextC'}{$context}{'T>C'}{'Tr'}) { $refH_file->{$filename}{'SeqContextC'}{$context}{'T>C'}{'Tr'}++; }
	    		}
	 				# T>C (A>G strand = +)
	    		else
	    		{
	    			my $base3           = complement($tempContextSequence[$before]); my $base5 = complement($tempContextSequence[$after]);
	    			my $context_reverse = $base5."_".$base3;
	    			if(exists $refH_file->{$filename}{'SeqContextC'}{$context_reverse}{'T>C'}{'Tr'}) { $refH_file->{$filename}{'SeqContextC'}{$context_reverse}{'T>C'}{'Tr'}++; }
	    		}
	    	}
	    	# WebLogo-3
	    	if(($tab[$ref_value] eq "T") && ($tab[$alt_value] eq "C"))
	    	{
	    		if(scalar(@tempContextSequence) == 2)  { next; }
	    		my ($contextTemp1, $contextTemp2) = ("", "");
					for(my $i=0; $i<$midlle_totalNbBaseContext; $i++) { $contextTemp1 .= $tempContextSequence[$i]; }
					for(my $i=$midlle_totalNbBaseContext+1; $i<=$#tempContextSequence; $i++) { $contextTemp2 .= $tempContextSequence[$i]; }
					my $context = $contextTemp1."T".$contextTemp2; $context = reverse $context;
					push(@{$refH_file->{$filename}{'WebLogo3'}{'TC'}}, $context);
	    	}
	    	else
	    	{
	    		if(scalar(@tempContextSequence) == 2)  { next; }
	    		my ($contextTemp1, $contextTemp2) = ("", "");
	    		for(my $i=0; $i<$midlle_totalNbBaseContext; $i++) { $contextTemp1 .= complement($tempContextSequence[$i]); }
					for(my $i=$midlle_totalNbBaseContext+1; $i<=$#tempContextSequence; $i++) { $contextTemp2 .= complement($tempContextSequence[$i]); }
					my $context = $contextTemp1."T".$contextTemp2;
					push(@{$refH_file->{$filename}{'WebLogo3'}{'TC'}}, $context);
	    	}
	    }
	    ###################################### T:A>G:C
	    if( (($tab[$ref_value] eq "T") && ($tab[$alt_value] eq "G")) || ( ($tab[$ref_value] eq "A") && ($tab[$alt_value] eq "C")) )
	    {
	    	my $mutation   = "T:A>G:C";
	    	$refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'TotalMutG'}++; # Count the total number of mutations

	    	# Pearson correlation
	    	if(exists $refH_file->{$filename}{'SBSPerChr'}{$mutation}{'CHR'}{$chrNameForH}{'chr'}) { $refH_file->{$filename}{'SBSPerChr'}{$mutation}{'CHR'}{$chrNameForH}{'chr'}++; }

	    	# Sequence context - 6 mutation types - genomic strand
	    	my $mutationSeqContext6mutType = "T>G";
	    	# We want to express the mutation in T>
	    	if( ($tab[$ref_value] eq "A") && ($tab[$alt_value] eq "T") )
	    	{
	    		my $base3 = complement($tempContextSequence[$before]); my $base5 = complement($tempContextSequence[$after]);
	    		my $context_reverse = $base5."_".$base3;
	    		if(exists $refH_file->{$filename}{'SeqContextG'}{$context_reverse}{$mutationSeqContext6mutType}) { $refH_file->{$filename}{'SeqContextG'}{$context_reverse}{$mutationSeqContext6mutType}++; }
	    	}
	    	elsif(exists $refH_file->{$filename}{'SeqContextG'}{$context}{$mutationSeqContext6mutType}) { $refH_file->{$filename}{'SeqContextG'}{$context}{$mutationSeqContext6mutType}++; }

	    	# Strand analysis T>G on NonTr strand
	    	if( (($tab[$strand_value] eq "+") && (($tab[$ref_value] eq "T")&&($tab[$alt_value] eq "G"))) || (($tab[$strand_value] eq "-") && (($tab[$ref_value] eq "A")&&($tab[$alt_value] eq "C"))) )
	    	{
	    		if(exists $refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'NonTr'}) { $refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'NonTr'}++; }

	    		# T>G (T>G strand = +)
	    		if( ($tab[$strand_value] eq "+") && (($tab[$ref_value] eq "T")&&($tab[$alt_value] eq "G")) )
	    		{
	    			if(exists $refH_file->{$filename}{'SeqContextC'}{$context}{'T>G'}{'NonTr'}) { $refH_file->{$filename}{'SeqContextC'}{$context}{'T>G'}{'NonTr'}++; }
	    		}
	    		# T>G (A>C strand = -)
	    		else
	    		{
	    			my $base3           = complement($tempContextSequence[$before]); my $base5 = complement($tempContextSequence[$after]);
	    			my $context_reverse = $base5."_".$base3;
	    			if(exists $refH_file->{$filename}{'SeqContextC'}{$context_reverse}{'T>G'}{'NonTr'}) { $refH_file->{$filename}{'SeqContextC'}{$context_reverse}{'T>G'}{'NonTr'}++; }
	    		}
	    	}
	    	# Strand analysis T>G on Tr strand
	    	if( (($tab[$strand_value] eq "-") && (($tab[$ref_value] eq "T")&&($tab[$alt_value] eq "G"))) || (($tab[$strand_value] eq "+") && (($tab[$ref_value] eq "A")&&($tab[$alt_value] eq "C"))) )
	    	{
	    		if(exists $refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'Tr'}) { $refH_file->{$filename}{'6mutType'}{$funcSegment}{$mutation}{'Tr'}++; }

	    		# T>G (T>G strand = -)
	    		if( ($tab[$strand_value] eq "-") && (($tab[$ref_value] eq "T")&&($tab[$alt_value] eq "G")) )
	    		{
	    			if(exists $refH_file->{$filename}{'SeqContextC'}{$context}{'T>G'}{'Tr'}) { $refH_file->{$filename}{'SeqContextC'}{$context}{'T>G'}{'Tr'}++; }
	    		}
	 				# T>G (A>C strand = +)
	    		else
	    		{
	    			my $base3           = complement($tempContextSequence[$before]); my $base5 = complement($tempContextSequence[$after]);
	    			my $context_reverse = $base5."_".$base3;
	    			if(exists $refH_file->{$filename}{'SeqContextC'}{$context_reverse}{'T>G'}{'Tr'}) { $refH_file->{$filename}{'SeqContextC'}{$context_reverse}{'T>G'}{'Tr'}++; }
	    		}
	    	}
	    	# WebLogo-3
	    	if(($tab[$ref_value] eq "T") && ($tab[$alt_value] eq "G"))
	    	{
	    		if(scalar(@tempContextSequence) == 2)  { next; }
	    		my ($contextTemp1, $contextTemp2) = ("", "");
	    		for(my $i=0; $i<$midlle_totalNbBaseContext; $i++) { $contextTemp1 .= $tempContextSequence[$i]; }
					for(my $i=$midlle_totalNbBaseContext+1; $i<=$#tempContextSequence; $i++) { $contextTemp2 .= $tempContextSequence[$i]; }
					my $context = $contextTemp1."T".$contextTemp2; $context = reverse $context;
					push(@{$refH_file->{$filename}{'WebLogo3'}{'TG'}}, $context);
	    	}
	    	else
	    	{
	    		if(scalar(@tempContextSequence) == 2)  { next; }
	    		my ($contextTemp1, $contextTemp2) = ("", "");
					for(my $i=0; $i<$midlle_totalNbBaseContext; $i++) { $contextTemp1 .= complement($tempContextSequence[$i]); }
					for(my $i=$midlle_totalNbBaseContext+1; $i<=$#tempContextSequence; $i++) { $contextTemp2 .= complement($tempContextSequence[$i]); }
					my $context = $contextTemp1."T".$contextTemp2;
					push(@{$refH_file->{$filename}{'WebLogo3'}{'TG'}}, $context);
	    	}
	    }
	  }
	  close F1;
	}
	# Write the different statistics in the report
	sub WriteStatistics
	{
		my ($refH_file, $nb_func, $folderFigure, $folderChi2, $folderNMF) = @_;

		# Save the different graphs in specific folder instead of in a general one.
		if(!-e "$folderFigure/Overall_mutation_distribution") { mkdir("$folderFigure/Overall_mutation_distribution") or die "Can't create the directory $folderFigure/Overall_mutation_distribution\n"; }
		if(!-e "$folderFigure/Impact_protein_sequence") { mkdir("$folderFigure/Impact_protein_sequence") or die "Can't create the directory $folderFigure/Impact_protein_sequence\n"; }
		if(!-e "$folderFigure/SBS_distribution") { mkdir("$folderFigure/SBS_distribution") or die "Can't create the directory $folderFigure/SBS_distribution\n"; }
		if(!-e "$folderFigure/Stranded_Analysis") { mkdir("$folderFigure/Stranded_Analysis") or die "Can't create the directory $folderFigure/Stranded_Analysis\n"; }
		if(!-e "$folderFigure/Trinucleotide_Sequence_Context") { mkdir("$folderFigure/Trinucleotide_Sequence_Context") or die "Can't create the directory $folderFigure/Trinucleotide_Sequence_Context\n"; }
		if(!-e "$folderFigure/Distribution_SBS_Per_Chromosomes") { mkdir("$folderFigure/Distribution_SBS_Per_Chromosomes") or die "Can't create the directory $folderFigure/Distribution_SBS_Per_Chromosomes\n"; }


		# Create a workbook with all the samples
		my $wb = ""; my $ws_sum = ""; my %h_chi2 = ();
		my ($ws_inputNMF_count, $ws_inputNMF_percent) = ("", "");
		############### Define the format
		my ($format_A10, $format_A10Boldleft, $format_A10Italic) = ("", "", "");
		my ($formatT_left, $formatT_right, $formatT_bottomRight, $formatT_bottomLeft, $formatT_bottom, $formatT_bottomHeader, $formatT_bottomRightHeader, $formatT_bottomHeader2, $formatT_rightHeader);
		my ($formatT_graphTitle);
		my ($table_topleft, $table_topRight, $table_bottomleft, $table_bottomRight, $table_top, $table_right, $table_bottom, $table_bottomItalic, $table_left, $table_bottomrightHeader, $table_left2, $table_middleHeader, $table_middleHeader2);

		if($oneReportPerSample == 2)
		{
			$wb = Spreadsheet::WriteExcel->new("$output/Report_Mutation_Spectra.xls");

			############### Define the format
			Format_A10($wb, \$format_A10); # Text center in Arial 10
			Format_A10BoldLeft($wb, \$format_A10Boldleft); # Text on the left in Arial 10 bold
			Format_TextSection($wb, \$formatT_left, \$formatT_right, \$formatT_bottomRight, \$formatT_bottomLeft, \$formatT_bottom, \$formatT_bottomHeader, \$formatT_bottomRightHeader, \$formatT_bottomHeader2, \$formatT_rightHeader);
			Format_GraphTitle($wb, \$formatT_graphTitle);
			Format_Table($wb, \$table_topleft, \$table_topRight, \$table_bottomleft, \$table_bottomRight, \$table_top, \$table_right, \$table_bottom, \$table_bottomItalic, \$table_left, \$table_bottomrightHeader, \$table_left2, \$table_middleHeader, \$table_middleHeader2);
			Format_A10Italic($wb, \$format_A10Italic);


			############### Worksheet with a summary of the samples
			$ws_sum = $wb->add_worksheet("Sample_List");
			$ws_sum->write(0, 0, "Samples", $format_A10); $ws_sum->write(0, 1, "Total number SBS", $format_A10); $ws_sum->write(0, 2, "Total number of Indel", $format_A10); $ws_sum->write(0, 3, "Total number of mutations", $format_A10);
			$ws_sum->set_column(0,0, 50); $ws_sum->set_column(1,1, 20); $ws_sum->set_column(2,2, 20); $ws_sum->set_column(3,3, 22);

			############### Save the chi2 values into a hash table
			if(-e "$folderChi2/Output_chi2_strandBias.txt")
			{
				open(F1, "$folderChi2/Output_chi2_strandBias.txt") or die "$!: $folderChi2/Output_chi2_strandBias.txt\n";
				my $header = <F1>; # Strand_Bias($tab[0])	NonTr-Tr($tab[1])	Proportion($tab[2])	P-val-Chi2($tab[3])	FDR($tab[4])	Confidence Interval($tab[5])	Mutation_Type($tab[6])	SampleName($tab[7])
				while(<F1>)
				{
					$_      =~ s/[\r\n]+$//;
					my @tab = split("\t", $_);

					$h_chi2{$tab[7]}{$tab[6]}{'p-value'} = $tab[3]; $h_chi2{$tab[7]}{$tab[6]}{'ConfInt'} = $tab[5];

					# For the pool data the FDR isn't calculated so replace the NA (=Missing values) by "-"
					if($tab[7] eq "Pool_Data") { $h_chi2{$tab[7]}{$tab[6]}{'FDR'} = "-"; }
					else { $h_chi2{$tab[7]}{$tab[6]}{'FDR'} = $tab[4]; }
				}
				close F1;
			}
			############### Write the input matrix for NMF for the count and the un-normalized frequency
			$ws_inputNMF_count   = $wb->add_worksheet("Input_NMF_Count");
			$ws_inputNMF_percent = $wb->add_worksheet("Input_NMF_Percent");
		}


		################################################ Set the Rows and columns of the different part of the report ################################################
		my $row_SumSheet              = 1; # First row for the summary sheet of the report
		my $rowStart_SBSdistrBySeg    = 48; my $colStart_SBSdistrBySeg  = 0; # For the table SBS distribution by segment
		my $colStart_matrixSeqContext = 19; # Sequence context
		my $col_inputNMF              = 0;  # Write the names of the samples with at least 33 SBS

		# For NMF input
  	my %h_inputNMF = ();

		## For each file
		foreach my $k_file (sort keys $refH_file)
		{
			print "File in process: $k_file\n";
			if($k_file ne "Pool_Data") { $col_inputNMF++; }

			# Create one workbook for each sample
			if($oneReportPerSample == 1)
			{
				$wb = Spreadsheet::WriteExcel->new("$output/Report_Mutation_Spectra-$k_file.xls");

				############### Define the format
				Format_A10($wb, \$format_A10); # Text center in Arial 10
				Format_A10BoldLeft($wb, \$format_A10Boldleft); # Text on the left in Arial 10 bold
				Format_TextSection($wb, \$formatT_left, \$formatT_right, \$formatT_bottomRight, \$formatT_bottomLeft, \$formatT_bottom, \$formatT_bottomHeader, \$formatT_bottomRightHeader, \$formatT_bottomHeader2, \$formatT_rightHeader);
				Format_GraphTitle($wb, \$formatT_graphTitle);
				Format_Table($wb, \$table_topleft, \$table_topRight, \$table_bottomleft, \$table_bottomRight, \$table_top, \$table_right, \$table_bottom, \$table_bottomItalic, \$table_left, \$table_bottomrightHeader, \$table_left2, \$table_middleHeader, \$table_middleHeader2);
				Format_A10Italic($wb, \$format_A10Italic);


		  	############### Worksheet with a summary of the samples
				$ws_sum = $wb->add_worksheet("Sample_List");
				$ws_sum->write(0, 0, "Samples", $format_A10); $ws_sum->write(0, 1, "Total number SBS", $format_A10); $ws_sum->write(0, 2, "Total number of Indel", $format_A10); $ws_sum->write(0, 3, "Total number of mutations", $format_A10);
				$ws_sum->set_column(0,0, 50); $ws_sum->set_column(1,1, 20); $ws_sum->set_column(2,2, 20); $ws_sum->set_column(3,3, 22);
				# Write in the Samples sheet the name and the total number of SBS
				$ws_sum->write(1, 0, "$k_file", $format_A10);
				$ws_sum->write(1, 1, $refH_file->{$k_file}{'TotalSBSGenomic'}, $format_A10); $ws_sum->write(1, 2, $refH_file->{$k_file}{'TotalIndelGenomic'}, $format_A10); $ws_sum->write($row_SumSheet, 3, $refH_file->{$k_file}{'TotalMutGenomic'}, $format_A10);

				############### Save the chi2 values into a hash table
				if(-e "$folderChi2/Output_chi2_strandBias.txt")
				{
					open(F1, "$folderChi2/Output_chi2_strandBias.txt") or die "$!: $folderChi2/Output_chi2_strandBias.txt\n";
					my $header = <F1>; # Strand_Bias($tab[0])	NonTr-Tr($tab[1])	Proportion($tab[2])	P-val-Chi2($tab[3])	FDR($tab[4])	Confidence Interval($tab[5])	Mutation_Type($tab[6])	SampleName($tab[7])
					while(<F1>)
					{
						$_      =~ s/[\r\n]+$//;
						my @tab = split("\t", $_);

						if($tab[7] eq $k_file)
						{
							$h_chi2{$tab[7]}{$tab[6]}{'p-value'} = $tab[3]; $h_chi2{$tab[7]}{$tab[6]}{'ConfInt'} = $tab[5];

							# For the pool data the FDR isn't calculated so replace the NA (=Missing values) by "-"
							if($tab[7] eq "Pool_Data") { $h_chi2{$tab[7]}{$tab[6]}{'FDR'} = "-"; }
							else { $h_chi2{$tab[7]}{$tab[6]}{'FDR'} = $tab[4]; }
						}
					}
					close F1;
				}

				############### Write the input matrix for NMF
				if($k_file ne "Pool_Data")
				{
					# For NMF don't consider the pool of the samples
					$ws_inputNMF_count   = $wb->add_worksheet("Input_NMF_Count");
					$ws_inputNMF_percent = $wb->add_worksheet("Input_NMF_Percent");
					# Write in the input NMF sheet the name of the samples
					$ws_inputNMF_count->write(0, 1, $k_file);
					$ws_inputNMF_percent->write(0, 1, $k_file);
				}
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

			# Count of SBS per chromosome
			PearsonCoefficient($refH_file, $k_file);

			# Add a worksheet to the workbook
  		my $ws = $wb->add_worksheet($k_file);

  		# Write the titles of the different sections of the report
  		WriteBoderSection($wb, $ws, $rowStart_SBSdistrBySeg, $colStart_SBSdistrBySeg, $nb_func, $colStart_matrixSeqContext);

  		# Write the mutation types (6 types)
  		WriteHeaderSection($wb, $ws, $rowStart_SBSdistrBySeg, $colStart_SBSdistrBySeg, $nb_func, $colStart_matrixSeqContext);


  		# Save the figures of each samples in a different folder
			if(!-e "$folderFigure/Overall_mutation_distribution/$k_file") { mkdir("$folderFigure/Overall_mutation_distribution/$k_file") or die "Can't create the directory $folderFigure/Overall_mutation_distribution/$k_file\n"; }
			if(!-e "$folderFigure/Impact_protein_sequence/$k_file") { mkdir("$folderFigure/Impact_protein_sequence/$k_file") or die "Can't create the directory $folderFigure/Impact_protein_sequence/$k_file\n"; }
			if(!-e "$folderFigure/SBS_distribution/$k_file") { mkdir("$folderFigure/SBS_distribution/$k_file") or die "Can't create the directory $folderFigure/SBS_distribution\n"; }
			if(!-e "$folderFigure/Stranded_Analysis/$k_file") { mkdir("$folderFigure/Stranded_Analysis/$k_file") or die "Can't create the directory $folderFigure/Stranded_Analysis/$k_file\n"; }
			if(!-e "$folderFigure/Trinucleotide_Sequence_Context/$k_file") { mkdir("$folderFigure/Trinucleotide_Sequence_Context/$k_file") or die "Can't create the directory $folderFigure/Trinucleotide_Sequence_Context/$k_file\n"; }



  		###########################################################################################################################################################
  		#################################################################	Write the statistics	###################################################################
  		###########################################################################################################################################################
  		my ($ca_genomique, $cg_genomique, $ct_genomique, $ta_genomique, $tc_genomique, $tg_genomique) = (0,0,0,0,0,0);
  		my ($ca_NonTr, $ca_Tr, $cg_NonTr, $cg_Tr, $ct_NonTr, $ct_Tr, $ta_NonTr, $ta_Tr, $tc_NonTr, $tc_Tr, $tg_NonTr, $tg_Tr) = (0,0,0,0,0,0, 0,0,0,0,0,0);

			my $row_SBSdistrBySeg           = $rowStart_SBSdistrBySeg+4;
			my $row_SBSDistrBySegAndFunc_CA = $rowStart_SBSdistrBySeg+$nb_func+12;
			my $row_SBSDistrBySegAndFunc_CG = $rowStart_SBSdistrBySeg+($nb_func*2)+16; my $rowEndCG_SBSDistrBySegAndFunc_CG = $row_SBSDistrBySegAndFunc_CG+$nb_func;
			my $row_SBSDistrBySegAndFunc_CT = $rowStart_SBSdistrBySeg+($nb_func*3)+20;

			## 6 mutation types by segment
			foreach my $k_func (sort keys $refH_file->{$k_file}{'6mutType'})
			{
				my $totalSBS_bySegment = 0;

				# Write the functional region for the section SBS distribution by segment
				$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg, $k_func, $formatT_left);

				# Write the exonic func for the section strand bias by segment
				$ws->write($row_SBSDistrBySegAndFunc_CA, $colStart_SBSdistrBySeg, $k_func, $formatT_left);

				if($row_SBSDistrBySegAndFunc_CG == $rowEndCG_SBSDistrBySegAndFunc_CG)
				{
					$ws->write($row_SBSDistrBySegAndFunc_CG, $colStart_SBSdistrBySeg, $k_func, $formatT_bottomLeft);
				}
				else
				{
					$ws->write($row_SBSDistrBySegAndFunc_CG, $colStart_SBSdistrBySeg, $k_func, $formatT_left);
				}

				foreach my $k_mutation (sort keys $refH_file->{$k_file}{'6mutType'}{$k_func})
				{
					if($k_mutation eq "C:G>A:T")
					{
						# Write the ratio NonTr(CA)/Tr(GT)
						my $ratioSB = 0;
						if( ($refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'} == 0) || ($refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'} == 0) ) { $ratioSB = 0;  }
						else { $ratioSB = $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'} / $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'}; }
						$ratioSB = sprintf("%.2f", $ratioSB);
						$ws->write($row_SBSDistrBySegAndFunc_CA, $colStart_SBSdistrBySeg+1, $ratioSB, $format_A10);

						# Write the count of SBS in the NonTr and Tr strand
						$ws->write($row_SBSDistrBySegAndFunc_CA, $colStart_SBSdistrBySeg+2, $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'}, $format_A10);
						$ws->write($row_SBSDistrBySegAndFunc_CA, $colStart_SBSdistrBySeg+3, $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'}, $format_A10);

						# Calculate the total number of SBS per mut type (genomic strand)
						$ca_genomique += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'TotalMutG'};
						# Calculate the total number of SBS by NonTr / Tr strand
						$ca_NonTr += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'}; $ca_Tr += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'};

						# Write the count by exonic region
						$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+3, $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'TotalMutG'}, $format_A10);
					}
					if($k_mutation eq "C:G>G:C")
					{
						# Write the ratio NonTr(CG)/Tr(GC)
						my $ratioSB = 0;
						if( ($refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'} == 0) || ($refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'} == 0) ) { $ratioSB = 0;  }
						else { $ratioSB = $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'} / $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'}; }
						$ratioSB = sprintf("%.2f", $ratioSB);
						$ws->write($row_SBSDistrBySegAndFunc_CA, $colStart_SBSdistrBySeg+5, $ratioSB, $format_A10);

						# Write the count of SBS in the NonTr and Tr strand
						$ws->write($row_SBSDistrBySegAndFunc_CA, $colStart_SBSdistrBySeg+6, $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'}, $format_A10);
						$ws->write($row_SBSDistrBySegAndFunc_CA, $colStart_SBSdistrBySeg+7, $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'}, $format_A10);

						# Calculate the total number of SBS per mut type (genomic strand)
						$cg_genomique += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'TotalMutG'};
						# Calculate the total number of SBS by NonTr / Tr strand
						$cg_NonTr += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'}; $cg_Tr += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'};

						# Write the count by exonic region
						$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+5, $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'TotalMutG'}, $format_A10);
					}
					if($k_mutation eq "C:G>T:A")
					{
						# Write the ratio NonTr(CT)/Tr(GA)
						my $ratioSB = 0;
						if( ($refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'} == 0) || ($refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'} == 0) ) { $ratioSB = 0;  }
						else { $ratioSB = $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'} / $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'}; }
						$ratioSB = sprintf("%.2f", $ratioSB);
						$ws->write($row_SBSDistrBySegAndFunc_CA, $colStart_SBSdistrBySeg+9, $ratioSB, $format_A10);

						# Write the count of SBS in the NonTr and Tr strand
						$ws->write($row_SBSDistrBySegAndFunc_CA, $colStart_SBSdistrBySeg+10, $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'}, $format_A10);
						$ws->write($row_SBSDistrBySegAndFunc_CA, $colStart_SBSdistrBySeg+11, $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'}, $formatT_right);

						# Calculate the total number of SBS per mut type (genomic strand)
						$ct_genomique += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'TotalMutG'};
						# Calculate the total number of SBS by NonTr / Tr strand
						$ct_NonTr += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'}; $ct_Tr += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'};

						# Write the count by exonic region
						$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+7, $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'TotalMutG'}, $format_A10);
					}
					if($k_mutation eq "T:A>A:T")
					{
						# Write the ratio NonTr(AT)/Tr(TA)
						my $ratioSB = 0;
						if( ($refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'} == 0) || ($refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'} == 0) ) { $ratioSB = 0;  }
						else { $ratioSB = $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'} / $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'}; }
						$ratioSB = sprintf("%.2f", $ratioSB);


						if($row_SBSDistrBySegAndFunc_CG == $rowEndCG_SBSDistrBySegAndFunc_CG)
						{
							# Write the ratio NonTr(AC)/Tr(TG)
							$ws->write($row_SBSDistrBySegAndFunc_CG, $colStart_SBSdistrBySeg+1, $ratioSB, $formatT_bottom);
							# Write the count of SBS in the NonTr and Tr strand
							$ws->write($row_SBSDistrBySegAndFunc_CG, $colStart_SBSdistrBySeg+2, $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'}, $formatT_bottom);
							$ws->write($row_SBSDistrBySegAndFunc_CG, $colStart_SBSdistrBySeg+3, $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'}, $formatT_bottom);
						}
						else
						{
							# Write the ratio NonTr(AC)/Tr(TG)
							$ws->write($row_SBSDistrBySegAndFunc_CG, $colStart_SBSdistrBySeg+1, $ratioSB, $format_A10);
							# Write the count of SBS in the NonTr and Tr strand
							$ws->write($row_SBSDistrBySegAndFunc_CG, $colStart_SBSdistrBySeg+2, $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'}, $format_A10);
							$ws->write($row_SBSDistrBySegAndFunc_CG, $colStart_SBSdistrBySeg+3, $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'}, $format_A10);
						}


						# Calculate the total number of SBS per mut type (genomic strand)
						$ta_genomique += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'TotalMutG'};
						# Calculate the total number of SBS by NonTr / Tr strand
						$ta_NonTr += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'}; $ta_Tr += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'};

						# Write the count by exonic region
						$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+9, $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'TotalMutG'}, $format_A10);
					}
					if($k_mutation eq "T:A>C:G")
					{
						# Write the ratio NonTr(AG)/Tr(TC)
						my $ratioSB = 0;
						if( ($refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'} == 0) || ($refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'} == 0) ) { $ratioSB = 0;  }
						else { $ratioSB = $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'} / $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'}; }
						$ratioSB = sprintf("%.2f", $ratioSB);

						if($row_SBSDistrBySegAndFunc_CG == $rowEndCG_SBSDistrBySegAndFunc_CG)
						{
							# Write the ratio NonTr(AC)/Tr(TG)
							$ws->write($row_SBSDistrBySegAndFunc_CG, $colStart_SBSdistrBySeg+5, $ratioSB, $formatT_bottom);
							# Write the count of SBS in the NonTr and Tr strand
							$ws->write($row_SBSDistrBySegAndFunc_CG, $colStart_SBSdistrBySeg+6, $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'}, $formatT_bottom);
							$ws->write($row_SBSDistrBySegAndFunc_CG, $colStart_SBSdistrBySeg+7, $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'}, $formatT_bottom);
						}
						else
						{
							# Write the ratio NonTr(AC)/Tr(TG)
							$ws->write($row_SBSDistrBySegAndFunc_CG, $colStart_SBSdistrBySeg+5, $ratioSB, $format_A10);
							# Write the count of SBS in the NonTr and Tr strand
							$ws->write($row_SBSDistrBySegAndFunc_CG, $colStart_SBSdistrBySeg+6, $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'}, $format_A10);
							$ws->write($row_SBSDistrBySegAndFunc_CG, $colStart_SBSdistrBySeg+7, $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'}, $format_A10);
						}

						# Calculate the total number of SBS per mut type (genomic strand)
						$tc_genomique += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'TotalMutG'};
						# Calculate the total number of SBS by NonTr / Tr strand
						$tc_NonTr += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'}; $tc_Tr += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'};

						# Write the count by exonic region
						$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+11, $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'TotalMutG'}, $format_A10);
					}
					if($k_mutation eq "T:A>G:C")
					{
						# Calculate the ratio for the strand bias
						my $ratioSB = 0;
						if( ($refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'} == 0) || ($refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'} == 0) ) { $ratioSB = 0;  }
						else { $ratioSB = $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'} / $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'}; }
						$ratioSB = sprintf("%.2f", $ratioSB);

						if($row_SBSDistrBySegAndFunc_CG == $rowEndCG_SBSDistrBySegAndFunc_CG)
						{
							# Write the ratio NonTr(AC)/Tr(TG)
							$ws->write($row_SBSDistrBySegAndFunc_CG, $colStart_SBSdistrBySeg+9, $ratioSB, $formatT_bottom);
							# Write the count of SBS in the NonTr and Tr strand
							$ws->write($row_SBSDistrBySegAndFunc_CG, $colStart_SBSdistrBySeg+10, $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'}, $formatT_bottom);
							$ws->write($row_SBSDistrBySegAndFunc_CG, $colStart_SBSdistrBySeg+11, $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'}, $formatT_bottomRight);
						}
						else
						{
							# Write the ratio NonTr(AC)/Tr(TG)
							$ws->write($row_SBSDistrBySegAndFunc_CG, $colStart_SBSdistrBySeg+9, $ratioSB, $format_A10);
							# Write the count of SBS in the NonTr and Tr strand
							$ws->write($row_SBSDistrBySegAndFunc_CG, $colStart_SBSdistrBySeg+10, $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'}, $format_A10);
							$ws->write($row_SBSDistrBySegAndFunc_CG, $colStart_SBSdistrBySeg+11, $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'}, $formatT_right);
						}

						# Calculate the total number of SBS per mut type (genomic strand)
						$tg_genomique += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'TotalMutG'};
						# Calculate the total number of SBS by NonTr / Tr strand
						$tg_NonTr += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'NonTr'}; $tg_Tr += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'Tr'};

						# Write the count by exonic region
						$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+13, $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'TotalMutG'}, $formatT_right);
					}

					# Calculate the total number of SBS on the genomic strand for each mutation types by exonic region
					$totalSBS_bySegment += $refH_file->{$k_file}{'6mutType'}{$k_func}{$k_mutation}{'TotalMutG'};
				} # End $k_mutation

				$row_SBSDistrBySegAndFunc_CA++; $row_SBSDistrBySegAndFunc_CG++; $row_SBSDistrBySegAndFunc_CT++;

				# Write the percent by exonic region
				my $percent_ca = 0;
				if($refH_file->{$k_file}{'6mutType'}{$k_func}{'C:G>A:T'}{'TotalMutG'} == 0) { $percent_ca = 0; }
				else { $percent_ca = ($refH_file->{$k_file}{'6mutType'}{$k_func}{'C:G>A:T'}{'TotalMutG'} / $totalSBS_bySegment ) * 100; $percent_ca = sprintf("%.2f", $percent_ca); }
				$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+2, $percent_ca, $format_A10);
				my $percent_cg = 0;
				if($refH_file->{$k_file}{'6mutType'}{$k_func}{'C:G>A:T'}{'TotalMutG'} == 0) { $percent_cg = 0; }
				else { $percent_cg = ($refH_file->{$k_file}{'6mutType'}{$k_func}{'C:G>G:C'}{'TotalMutG'} / $totalSBS_bySegment ) * 100; $percent_cg = sprintf("%.2f", $percent_cg); }
				$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+4, $percent_cg, $format_A10);
				my $percent_ct = 0;
				if($refH_file->{$k_file}{'6mutType'}{$k_func}{'C:G>A:T'}{'TotalMutG'} == 0) { $percent_ct = 0; }
				else { $percent_ct = ($refH_file->{$k_file}{'6mutType'}{$k_func}{'C:G>T:A'}{'TotalMutG'} / $totalSBS_bySegment ) * 100; $percent_ct = sprintf("%.2f", $percent_ct); }
				$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+6, $percent_ct, $format_A10);
				my $percent_ta = 0;
				if($refH_file->{$k_file}{'6mutType'}{$k_func}{'C:G>A:T'}{'TotalMutG'} == 0) { $percent_ta = 0; }
				else { $percent_ta = ($refH_file->{$k_file}{'6mutType'}{$k_func}{'T:A>A:T'}{'TotalMutG'} / $totalSBS_bySegment ) * 100; $percent_ta = sprintf("%.2f", $percent_ta); }
				$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+8, $percent_ta, $format_A10);
				my $percent_tc = 0;
				if($refH_file->{$k_file}{'6mutType'}{$k_func}{'C:G>A:T'}{'TotalMutG'} == 0) { $percent_tc = 0; }
				else { $percent_tc = ($refH_file->{$k_file}{'6mutType'}{$k_func}{'T:A>C:G'}{'TotalMutG'} / $totalSBS_bySegment ) * 100; $percent_tc = sprintf("%.2f", $percent_tc); }
				$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+10, $percent_tc, $format_A10);
				my $percent_tg = 0;
				if($refH_file->{$k_file}{'6mutType'}{$k_func}{'C:G>A:T'}{'TotalMutG'} == 0) { $percent_tg = 0; }
				else { $percent_tg = ($refH_file->{$k_file}{'6mutType'}{$k_func}{'T:A>G:C'}{'TotalMutG'} / $totalSBS_bySegment ) * 100; $percent_tg = sprintf("%.2f", $percent_tg); }
				$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+12, $percent_tg, $format_A10);

				# Write the count of SBS by segment
				$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+1, $totalSBS_bySegment, $format_A10);

				$row_SBSdistrBySeg++;
			} # End $k_func

			# Write the total number of SBS on the genomic strand
			$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+1, $refH_file->{$k_file}{'TotalSBSGenomic'}, $formatT_bottomHeader);

			# Write the total and the percentage of SBS for each mutation types and save it to a text file
			open(DISTRSBS, ">", "$folderFigure/SBS_distribution/$k_file/$k_file-SBS_distribution.txt") or die "$!: $folderFigure/SBS_distribution/$k_file/$k_file-SBS_distribution.txt\n";
			print DISTRSBS "Mutation_Type\tCount\tPercentage\tSample\n";
			my $percent_ca = 0;
			if($ca_genomique == 0) { $percent_ca = 0; }
			else { $percent_ca = ($ca_genomique/$refH_file->{$k_file}{'TotalSBSGenomic'})*100; $percent_ca = sprintf("%.2f", $percent_ca); }
			$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+2, $percent_ca, $formatT_bottom); print DISTRSBS "C:G>A:T\t$ca_genomique\t$percent_ca\t$k_file\n";
			$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+3, $ca_genomique, $formatT_bottomHeader);
			my $percent_cg = 0;
			if($cg_genomique == 0) { $percent_cg = 0; }
			else { $percent_cg = ($cg_genomique/$refH_file->{$k_file}{'TotalSBSGenomic'})*100; $percent_cg = sprintf("%.2f", $percent_cg); }
			$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+4, $percent_cg, $formatT_bottom); print DISTRSBS "C:G>G:C\t$cg_genomique\t$percent_cg\t$k_file\n";
			$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+5, $cg_genomique, $formatT_bottomHeader);
			my $percent_ct = 0;
			if($ct_genomique == 0) { $percent_ct = 0; }
			else { $percent_ct = ($ct_genomique/$refH_file->{$k_file}{'TotalSBSGenomic'})*100; $percent_ct = sprintf("%.2f", $percent_ct); }
			$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+6, $percent_ct, $formatT_bottom); print DISTRSBS "C:G>T:A\t$ct_genomique\t$percent_ct\t$k_file\n";
			$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+7, $ct_genomique, $formatT_bottomHeader);
			my $percent_ta = 0;
			if($ta_genomique == 0) { $percent_ta = 0; }
			else { $percent_ta = ($ta_genomique/$refH_file->{$k_file}{'TotalSBSGenomic'})*100; $percent_ta = sprintf("%.2f", $percent_ta); }
			$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+8, $percent_ta, $formatT_bottom); print DISTRSBS "T:A>A:T\t$ta_genomique\t$percent_ta\t$k_file\n";
			$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+9, $ta_genomique, $formatT_bottomHeader);
			my $percent_tc = 0;
			if($tc_genomique == 0) { $percent_tc = 0; }
			else { $percent_tc = ($tc_genomique/$refH_file->{$k_file}{'TotalSBSGenomic'})*100; $percent_tc = sprintf("%.2f", $percent_tc); }
			$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+10, $percent_tc, $formatT_bottom); print DISTRSBS "T:A>C:G\t$tc_genomique\t$percent_tc\t$k_file\n";
			$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+11, $tc_genomique, $formatT_bottomHeader);
			my $percent_tg = 0;
			if($tg_genomique == 0) { $percent_tg = 0; }
			else { $percent_tg = ($tg_genomique/$refH_file->{$k_file}{'TotalSBSGenomic'})*100; $percent_tg = sprintf("%.2f", $percent_tg); }
			$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+12, $percent_tg, $formatT_bottom); print DISTRSBS "T:A>G:C\t$tg_genomique\t$percent_tg\t$k_file\n";
			$ws->write($row_SBSdistrBySeg, $colStart_SBSdistrBySeg+13, $tg_genomique, $formatT_bottomRightHeader);
			close DISTRSBS;

			###########################################################################################################################################################
  		###################################################################	Write Strand BIAS	#####################################################################
  		###########################################################################################################################################################
			# Write the SB for each mutation type (table 3)
			$ws->write(28, 11, "Table 3. Significance of the strand biases", $format_A10Boldleft);
			$ws->set_column(11, 11, 13); $ws->set_column(16, 16, 15); $ws->set_column(17, 17, 10);
			$ws->write(29, 11, "Mutation Type", $table_topleft); $ws->write(29, 12, "Non-Tr/Tr", $table_top); $ws->write(29, 13, "Non-Tr", $table_top); $ws->write(29, 14, "Tr", $table_top); $ws->write(29, 15, "P-value", $table_top); $ws->write(29, 16, "FDR q value", $table_top); $ws->write(29, 17, "95% CI", $table_topRight);

			$ws->write(39, 11, "Table 3. Significance of the strand biases", $format_A10Boldleft);
			$ws->write(40, 11, "Mutation Type", $table_topleft); $ws->write(40, 12, "Non-Tr/Tr", $table_top); $ws->write(40, 13, "Non-Tr", $table_top); $ws->write(40, 14, "Tr", $table_top); $ws->write(40, 15, "P-value", $table_top); $ws->write(40, 16, "FDR q value", $table_top); $ws->write(40, 17, "95% CI", $table_topRight);

			# For ggplot2
			open(SB, ">", "$folderFigure/Stranded_Analysis/$k_file/$k_file-StrandBias.txt") or die "$!: $folderFigure/Stranded_Analysis/$k_file/$k_file-StrandBias.txt\n";
			print SB "Alteration\tStrand\tCount\n";


			#-----------------------------------------------------------------------------------------------------#
			my ($ratio_ca, $ratio_gt, $percent_ca_NonTr, $percent_ca_Tr) = (0, 0, 0, 0, 0);
			if( ($ca_NonTr==0) || ($ca_Tr==0) ) { $ratio_ca = 0; $ratio_gt = 0; $percent_ca_NonTr = 0; $percent_ca_Tr = 0; }
			else
			{
				$ratio_ca = $ca_NonTr/$ca_Tr; $ratio_ca = sprintf("%.2f", $ratio_ca);
				$ratio_gt = $ca_Tr/$ca_NonTr; $ratio_gt = sprintf("%.2f", $ratio_gt);
				$percent_ca_NonTr = ($ca_NonTr/$refH_file->{$k_file}{'TotalSBSGenomic'})*100; $percent_ca_Tr = ($ca_Tr/$refH_file->{$k_file}{'TotalSBSGenomic'})*100;
			}
			print SB "C>A\tNonTranscribed\t$ca_NonTr\n", "C>A\tTranscribed\t$ca_Tr\n";
			# C>A
			$ws->write(30, 11, "C>A", $table_left); $ws->write(30, 12, $ratio_ca, $table_middleHeader); $ws->write(30, 13, $ca_NonTr, $format_A10); $ws->write(30, 14, $ca_Tr, $format_A10);
			# Write in italic (= warning message) when the count of NonTr + Tr is lower than 10
			if(($ca_NonTr+$ca_Tr)< 10) { $ws->write_string(30, 15, $h_chi2{$k_file}{'C>A'}{'p-value'}, $format_A10Italic); }
			else { $ws->write_string(30, 15, $h_chi2{$k_file}{'C>A'}{'p-value'}, $format_A10); }
			$ws->write(30, 16, $h_chi2{$k_file}{'C>A'}{'FDR'}, $format_A10); $ws->write(30, 17, $h_chi2{$k_file}{'C>A'}{'ConfInt'}, $table_right);
			# G>T
			$ws->write(41, 11, "G>T", $table_left); $ws->write(41, 12, $ratio_gt, $table_middleHeader); $ws->write(41, 13, $ca_Tr, $format_A10); $ws->write(41, 14, $ca_NonTr, $format_A10);
			if(($ca_NonTr+$ca_Tr)< 10) { $ws->write_string(41, 15, $h_chi2{$k_file}{'C>A'}{'p-value'}, $format_A10Italic); }
			else { $ws->write_string(41, 15, $h_chi2{$k_file}{'C>A'}{'p-value'}, $format_A10); }
			$ws->write(41, 16, $h_chi2{$k_file}{'C>A'}{'FDR'}, $format_A10); $ws->write(41, 17, $h_chi2{$k_file}{'C>A'}{'ConfInt'}, $table_right);

			#-----------------------------------------------------------------------------------------------------#
			my ($ratio_cg, $ratio_gc, $percent_cg_NonTr, $percent_cg_Tr) = (0, 0, 0, 0, 0);
			if( ($cg_NonTr==0) || ($cg_Tr==0) ) { $ratio_cg = 0; $ratio_gc = 0; $percent_cg_NonTr = 0; $percent_cg_Tr = 0; }
			else
			{
				$ratio_cg = $cg_NonTr/$cg_Tr; $ratio_cg = sprintf("%.2f", $ratio_cg);
				$ratio_gc = $cg_Tr/$cg_NonTr; $ratio_gc = sprintf("%.2f", $ratio_gc);
				$percent_cg_NonTr = ($cg_NonTr/$refH_file->{$k_file}{'TotalSBSGenomic'})*100; $percent_cg_Tr = ($cg_Tr/$refH_file->{$k_file}{'TotalSBSGenomic'})*100;
			}
			print SB "C>G\tNonTranscribed\t$cg_NonTr\n", "C>G\tTranscribed\t$cg_Tr\n";
			# C>G
			$ws->write(31, 11, "C>G", $table_left); $ws->write(31, 12, $ratio_cg, $table_middleHeader); $ws->write(31, 13, $cg_NonTr, $format_A10); $ws->write(31, 14, $cg_Tr, $format_A10);
			# Write in italic (= warning message) when the count of NonTr + Tr is lower than 10
			if(($cg_NonTr+$cg_Tr)< 10) { $ws->write_string(31, 15, $h_chi2{$k_file}{'C>G'}{'p-value'}, $format_A10Italic); }
			else { $ws->write_string(31, 15, $h_chi2{$k_file}{'C>G'}{'p-value'}, $format_A10); }
			$ws->write(31, 16, $h_chi2{$k_file}{'C>G'}{'FDR'}, $format_A10); $ws->write(31, 17, $h_chi2{$k_file}{'C>G'}{'ConfInt'}, $table_right);
			# G>C
			$ws->write(42, 11, "G>C", $table_left); $ws->write(42, 12, $ratio_gc, $table_middleHeader); $ws->write(42, 13, $cg_Tr, $format_A10); $ws->write(42, 14, $cg_NonTr, $format_A10);
			if(($cg_NonTr+$cg_Tr)< 10) { $ws->write_string(42, 15, $h_chi2{$k_file}{'C>G'}{'p-value'}, $format_A10Italic); }
			else { $ws->write_string(42, 15, $h_chi2{$k_file}{'C>G'}{'p-value'}, $format_A10); }
			$ws->write(42, 16, $h_chi2{$k_file}{'C>G'}{'FDR'}, $format_A10); $ws->write(42, 17, $h_chi2{$k_file}{'C>G'}{'ConfInt'}, $table_right);

			#-----------------------------------------------------------------------------------------------------#
			my ($ratio_ct, $ratio_ga, $percent_ct_NonTr, $percent_ct_Tr) = (0, 0, 0, 0, 0);
			if( ($ct_NonTr==0) || ($ct_Tr==0) ) { $ratio_ct = 0; $ratio_ga = 0; $percent_ct_NonTr = 0; $percent_ct_Tr = 0; }
			else
			{
				$ratio_ct = $ct_NonTr/$ct_Tr; $ratio_ct = sprintf("%.2f", $ratio_ct);
				$ratio_ga = $ct_Tr/$ct_NonTr; $ratio_ga = sprintf("%.2f", $ratio_ga);
				$percent_ct_NonTr = ($ct_NonTr/$refH_file->{$k_file}{'TotalSBSGenomic'})*100; $percent_ct_Tr = ($ct_Tr/$refH_file->{$k_file}{'TotalSBSGenomic'})*100;
			}
			print SB "C>T\tNonTranscribed\t$ct_NonTr\n", "C>T\tTranscribed\t$ct_Tr\n";
			# C>T
			$ws->write(32, 11, "C>T", $table_left); $ws->write(32, 12, $ratio_ct, $table_middleHeader); $ws->write(32, 13, $ct_NonTr, $format_A10); $ws->write(32, 14, $ct_Tr, $format_A10);
			# Write in italic (= warning message) when the count of NonTr + Tr is lower than 10
			if(($ct_NonTr+$ct_Tr)< 10) { $ws->write_string(32, 15, $h_chi2{$k_file}{'C>T'}{'p-value'}, $format_A10Italic); }
			else { $ws->write_string(32, 15, $h_chi2{$k_file}{'C>T'}{'p-value'}, $format_A10); }
			$ws->write(32, 16, $h_chi2{$k_file}{'C>T'}{'FDR'}, $format_A10); $ws->write(32, 17, $h_chi2{$k_file}{'C>T'}{'ConfInt'}, $table_right);
			# G>A
			$ws->write(43, 11, "G>A", $table_left); $ws->write(43, 12, $ratio_ga, $table_middleHeader); $ws->write(43, 13, $ct_Tr, $format_A10); $ws->write(43, 14, $ct_NonTr, $format_A10);
			if(($ct_NonTr+$ct_Tr)< 10) { $ws->write_string(43, 15, $h_chi2{$k_file}{'C>T'}{'p-value'}, $format_A10Italic); }
			else { $ws->write_string(43, 15, $h_chi2{$k_file}{'C>T'}{'p-value'}, $format_A10); }
			$ws->write(43, 16, $h_chi2{$k_file}{'C>T'}{'FDR'}, $format_A10); $ws->write(43, 17, $h_chi2{$k_file}{'C>T'}{'ConfInt'}, $table_right);

			#-----------------------------------------------------------------------------------------------------#
			my ($ratio_ta, $ratio_at, $percent_ta_NonTr, $percent_ta_Tr) = (0, 0, 0, 0, 0);
			if( ($ta_NonTr==0) || ($ta_Tr==0) ) { $ratio_ta = 0; $ratio_at = 0; $percent_ta_NonTr = 0; $percent_ta_Tr = 0; }
			else
			{
				$ratio_ta = $ta_NonTr/$ta_Tr; $ratio_ta = sprintf("%.2f", $ratio_ta);
				$ratio_at = $ta_Tr/$ta_NonTr; $ratio_at = sprintf("%.2f", $ratio_at);
				$percent_ta_NonTr = ($ta_NonTr/$refH_file->{$k_file}{'TotalSBSGenomic'})*100; $percent_ta_Tr = ($ta_Tr/$refH_file->{$k_file}{'TotalSBSGenomic'})*100;
			}
			print SB "T>A\tNonTranscribed\t$ta_NonTr\n", "T>A\tTranscribed\t$ta_Tr\n";
			# T>A
			$ws->write(33, 11, "T>A", $table_left); $ws->write(33, 12, $ratio_ta, $table_middleHeader); $ws->write(33, 13, $ta_NonTr, $format_A10); $ws->write(33, 14, $ta_Tr, $format_A10);
			# Write in italic (= warning message) when the count of NonTr + Tr is lower than 10
			if(($ta_NonTr+$ta_Tr)< 10) { $ws->write_string(33, 15, $h_chi2{$k_file}{'T>A'}{'p-value'}, $format_A10Italic); }
			else { $ws->write_string(33, 15, $h_chi2{$k_file}{'T>A'}{'p-value'}, $format_A10); }
			$ws->write(33, 16, $h_chi2{$k_file}{'T>A'}{'FDR'}, $format_A10); $ws->write(33, 17, $h_chi2{$k_file}{'T>A'}{'ConfInt'}, $table_right);
			# A>T
			$ws->write(44, 11, "A>T", $table_left); $ws->write(44, 12, $ratio_at, $table_middleHeader); $ws->write(44, 13, $ta_Tr, $format_A10); $ws->write(44, 14, $ta_NonTr, $format_A10);
			if(($ta_NonTr+$ta_Tr)< 10) { $ws->write_string(44, 15, $h_chi2{$k_file}{'T>A'}{'p-value'}, $format_A10Italic); }
			else { $ws->write_string(44, 15, $h_chi2{$k_file}{'T>A'}{'p-value'}, $format_A10); }
			$ws->write(44, 16, $h_chi2{$k_file}{'T>A'}{'FDR'}, $format_A10); $ws->write(44, 17, $h_chi2{$k_file}{'T>A'}{'ConfInt'}, $table_right);

			#-----------------------------------------------------------------------------------------------------#
			my ($ratio_tc, $ratio_ag, $percent_tc_NonTr, $percent_tc_Tr) = (0, 0, 0, 0, 0);
			if( ($tc_NonTr==0) || ($tc_Tr==0) ) { $ratio_tc = 0; $ratio_ag = 0; $percent_tc_NonTr = 0; $percent_tc_Tr = 0; }
			else
			{
				$ratio_tc = $tc_NonTr/$tc_Tr; $ratio_tc = sprintf("%.2f", $ratio_tc);
				$ratio_ag = $tc_Tr/$tc_NonTr; $ratio_ag = sprintf("%.2f", $ratio_ag);
				$percent_tc_NonTr = ($tc_NonTr/$refH_file->{$k_file}{'TotalSBSGenomic'})*100; $percent_tc_Tr = ($tc_Tr/$refH_file->{$k_file}{'TotalSBSGenomic'})*100;
			}
			print SB "T>C\tNonTranscribed\t$tc_NonTr\n", "T>C\tTranscribed\t$tc_Tr\n";
			# T>C
			$ws->write(34, 11, "T>C", $table_left); $ws->write(34, 12, $ratio_tc, $table_middleHeader); $ws->write(34, 13, $tc_NonTr, $format_A10); $ws->write(34, 14, $tc_Tr, $format_A10);
			# Write in italic (= warning message) when the count of NonTr + Tr is lower than 10
			if(($tc_NonTr+$tc_Tr)< 10) { $ws->write_string(34, 15, $h_chi2{$k_file}{'T>C'}{'p-value'}, $format_A10Italic); }
			else { $ws->write_string(34, 15, $h_chi2{$k_file}{'T>C'}{'p-value'}, $format_A10); }
			$ws->write(34, 16, $h_chi2{$k_file}{'T>C'}{'FDR'}, $format_A10); $ws->write(34, 17, $h_chi2{$k_file}{'T>C'}{'ConfInt'}, $table_right);
			# A>G
			$ws->write(45, 11, "A>G", $table_left); $ws->write(45, 12, $ratio_ag, $table_middleHeader); $ws->write(45, 13, $tc_Tr, $format_A10); $ws->write(45, 14, $tc_NonTr, $format_A10);
			if(($tc_NonTr+$tc_Tr)< 10) { $ws->write_string(45, 15, $h_chi2{$k_file}{'T>C'}{'p-value'}, $format_A10Italic); }
			else { $ws->write_string(45, 15, $h_chi2{$k_file}{'T>C'}{'p-value'}, $format_A10); }
			$ws->write(45, 16, $h_chi2{$k_file}{'T>C'}{'FDR'}, $format_A10); $ws->write(45, 17, $h_chi2{$k_file}{'T>C'}{'ConfInt'}, $table_right);

			#-----------------------------------------------------------------------------------------------------#
			my ($ratio_tg, $ratio_ac, $percent_tg_NonTr, $percent_tg_Tr) = (0, 0, 0, 0, 0);
			if( ($tg_NonTr==0) || ($tg_Tr==0) ) { $ratio_tg = 0; $ratio_ac = 0; $percent_tg_NonTr = 0; $percent_tg_Tr = 0; }
			else
			{
				$ratio_tg = $tg_NonTr/$tg_Tr; $ratio_tg = sprintf("%.2f", $ratio_tg);
				$ratio_ac = $tg_Tr/$tg_NonTr; $ratio_ac = sprintf("%.2f", $ratio_ac);
				$percent_tg_NonTr = ($tg_NonTr/$refH_file->{$k_file}{'TotalSBSGenomic'})*100; $percent_tg_Tr = ($tg_Tr/$refH_file->{$k_file}{'TotalSBSGenomic'})*100;
			}
			print SB "T>G\tNonTranscribed\t$tg_NonTr\n", "T>G\tTranscribed\t$tg_Tr\n";
			# T>G
			$ws->write(35, 11, "T>G", $table_bottomleft); $ws->write(35, 12, $ratio_tg, $table_middleHeader2); $ws->write(35, 13, $tg_NonTr, $table_bottom); $ws->write(35, 14, $tg_Tr, $table_bottom);
			# Write in italic (= warning message) when the count of NonTr + Tr is lower than 10
			if(($tg_NonTr+$tg_Tr)< 10) { $ws->write_string(35, 15, $h_chi2{$k_file}{'T>G'}{'p-value'}, $table_bottomItalic); }
			else { $ws->write_string(35, 15, $h_chi2{$k_file}{'T>G'}{'p-value'}, $table_bottom); }
			$ws->write(35, 16, $h_chi2{$k_file}{'T>G'}{'FDR'}, $table_bottom); $ws->write(35, 17, $h_chi2{$k_file}{'T>G'}{'ConfInt'}, $table_bottomRight);
			# A>C
			$ws->write(46, 11, "A>C", $table_bottomleft); $ws->write(46, 12, $ratio_ac, $table_middleHeader2); $ws->write(46, 13, $tg_Tr, $table_bottom); $ws->write(46, 14, $tg_NonTr, $table_bottom);
			if(($tg_NonTr+$tg_Tr)< 10) { $ws->write_string(46, 15, $h_chi2{$k_file}{'T>G'}{'p-value'}, $table_bottomItalic); }
			else { $ws->write_string(46, 15, $h_chi2{$k_file}{'T>G'}{'p-value'}, $table_bottom); }
			$ws->write(46, 16, $h_chi2{$k_file}{'T>G'}{'FDR'}, $table_bottom); $ws->write(46, 17, $h_chi2{$k_file}{'T>G'}{'ConfInt'}, $table_bottomRight);

			### Write a warning message when NonTr+Tr < 10
			my $format_italic  = $wb->add_format(font=>'Arial', size=>10, italic=>1);
			$ws->write(36, 11, "Warning message: chi-squarred approximation may be incorrect because the number of SBS", $format_italic);
			$ws->write(37, 11, "on Non-transcribed and transcribed strand is lower than 10", $format_italic);

			close SB;


			###########################################################################################################################################################
  		###################################################################	Write SBS Per Chr	#####################################################################
  		###########################################################################################################################################################
  		# For the HTML report
  		open(SBSPerChr, ">", "$folderFigure/Distribution_SBS_Per_Chromosomes/$k_file-DistributionSNVS_per_chromosome.txt") or die "$!: $folderFigure/Distribution_SBS_Per_Chromosomes/$k_file-DistributionSNVS_per_chromosome.txt\n";
  		print SBSPerChr "\tPearson\t$refH_file->{$k_file}{'SBSPerChr'}{'AllMutType'}\t", $refH_file->{$k_file}{'SBSPerChr'}{"C:G>A:T"}{'Pearson'},"\t", $refH_file->{$k_file}{'SBSPerChr'}{"C:G>G:C"}{'Pearson'},"\t", $refH_file->{$k_file}{'SBSPerChr'}{"C:G>T:A"}{'Pearson'},"\t", $refH_file->{$k_file}{'SBSPerChr'}{"T:A>A:T"}{'Pearson'},"\t", $refH_file->{$k_file}{'SBSPerChr'}{"T:A>C:G"}{'Pearson'},"\t", $refH_file->{$k_file}{'SBSPerChr'}{"T:A>G:C"}{'Pearson'},"\n";
  		print SBSPerChr "Chr\tSize\tAll SBS\tC:G>A:T\tC:G>G:C\tC:G>T:A\tT:A>A:T\tT:A>C:G\tT:A>G:C\n";

  		my $row_SBSPerChr = $row_SBSDistrBySegAndFunc_CG + 8; # Line 158

  		# Write the Pearson coefficient
  		$ws->write($row_SBSDistrBySegAndFunc_CG+6, $colStart_SBSdistrBySeg+3, $refH_file->{$k_file}{'SBSPerChr'}{"C:G>A:T"}{'Pearson'}, $format_A10);
  		$ws->write($row_SBSDistrBySegAndFunc_CG+6, $colStart_SBSdistrBySeg+4, $refH_file->{$k_file}{'SBSPerChr'}{"C:G>G:C"}{'Pearson'}, $format_A10);
  		$ws->write($row_SBSDistrBySegAndFunc_CG+6, $colStart_SBSdistrBySeg+5, $refH_file->{$k_file}{'SBSPerChr'}{"C:G>T:A"}{'Pearson'}, $format_A10);
  		$ws->write($row_SBSDistrBySegAndFunc_CG+6, $colStart_SBSdistrBySeg+6, $refH_file->{$k_file}{'SBSPerChr'}{"T:A>A:T"}{'Pearson'}, $format_A10);
  		$ws->write($row_SBSDistrBySegAndFunc_CG+6, $colStart_SBSdistrBySeg+7, $refH_file->{$k_file}{'SBSPerChr'}{"T:A>C:G"}{'Pearson'}, $format_A10);
  		$ws->write($row_SBSDistrBySegAndFunc_CG+6, $colStart_SBSdistrBySeg+8, $refH_file->{$k_file}{'SBSPerChr'}{"T:A>G:C"}{'Pearson'}, $formatT_right);

  		# Write the chromosome number and their sizes / Write the total of SBS per chromosome
  		my $line=0;

  		foreach my $chromosome (sort keys %chromosomes)
			{
				$ws->write($row_SBSPerChr+($line), $colStart_SBSdistrBySeg, $chromosome, $formatT_left);
				$ws->write($row_SBSPerChr+($line), $colStart_SBSdistrBySeg+1, $chromosomes{$chromosome}, $format_A10);
				$ws->write($row_SBSPerChr+($line), $colStart_SBSdistrBySeg+2, $refH_file->{$k_file}{'SBSPerChr'}{'TotalPerChr'}{$chromosome}{'chr'}, $format_A10);

				# Write the count per mutation
				$ws->write($row_SBSPerChr+($line), $colStart_SBSdistrBySeg+3, $refH_file->{$k_file}{'SBSPerChr'}{"C:G>A:T"}{'CHR'}{$chromosome}{'chr'}, $format_A10);
				$ws->write($row_SBSPerChr+($line), $colStart_SBSdistrBySeg+4, $refH_file->{$k_file}{'SBSPerChr'}{"C:G>G:C"}{'CHR'}{$chromosome}{'chr'}, $format_A10);
				$ws->write($row_SBSPerChr+($line), $colStart_SBSdistrBySeg+5, $refH_file->{$k_file}{'SBSPerChr'}{"C:G>T:A"}{'CHR'}{$chromosome}{'chr'}, $format_A10);
				$ws->write($row_SBSPerChr+($line), $colStart_SBSdistrBySeg+6, $refH_file->{$k_file}{'SBSPerChr'}{"T:A>A:T"}{'CHR'}{$chromosome}{'chr'}, $format_A10);
				$ws->write($row_SBSPerChr+($line), $colStart_SBSdistrBySeg+7, $refH_file->{$k_file}{'SBSPerChr'}{"T:A>C:G"}{'CHR'}{$chromosome}{'chr'}, $format_A10);
				$ws->write($row_SBSPerChr+($line), $colStart_SBSdistrBySeg+8, $refH_file->{$k_file}{'SBSPerChr'}{"T:A>G:C"}{'CHR'}{$chromosome}{'chr'}, $formatT_right);


				# For the HTML report
				print SBSPerChr "$chromosome\t", $chromosomes{$chromosome},"\t", $refH_file->{$k_file}{'SBSPerChr'}{'TotalPerChr'}{$chromosome}{'chr'},"\t", $refH_file->{$k_file}{'SBSPerChr'}{"C:G>A:T"}{'CHR'}{$chromosome}{'chr'},"\t", $refH_file->{$k_file}{'SBSPerChr'}{"C:G>G:C"}{'CHR'}{$chromosome}{'chr'},"\t", $refH_file->{$k_file}{'SBSPerChr'}{"C:G>T:A"}{'CHR'}{$chromosome}{'chr'},"\t", $refH_file->{$k_file}{'SBSPerChr'}{"T:A>A:T"}{'CHR'}{$chromosome}{'chr'},"\t", $refH_file->{$k_file}{'SBSPerChr'}{"T:A>C:G"}{'CHR'}{$chromosome}{'chr'},"\t", $refH_file->{$k_file}{'SBSPerChr'}{"T:A>G:C"}{'CHR'}{$chromosome}{'chr'},"\n";
				$line++;
			}

			# Write the Pearson coefficient for the total number of SBS
			$ws->write($row_SBSDistrBySegAndFunc_CG+6, $colStart_SBSdistrBySeg+2, $refH_file->{$k_file}{'SBSPerChr'}{'AllMutType'}, $format_A10);
			$ws->write($row_SBSPerChr+(keys %chromosomes), $colStart_SBSdistrBySeg+2, $refH_file->{$k_file}{'TotalSBSGenomic'}, $formatT_bottomHeader);

			print SBSPerChr "\t\t$refH_file->{$k_file}{'TotalSBSGenomic'}\n";
			close SBSPerChr;



			###########################################################################################################################################################
  		####################################################################### Impact on protein #################################################################
  		###########################################################################################################################################################
  		$ws->write(29, 6, "Table 2. Frequency and counts of functional impact", $format_A10Boldleft);
  		$ws->set_column(6, 6, 13); $ws->set_column(10, 10, 15);
  		$ws->write(30, 6, "RefSeq gene", $table_topleft); $ws->write(30, 7, "", $table_top); $ws->write(30, 8, "Percent", $table_top); $ws->write(30, 9, "Count", $table_topRight);
  		my $lImpactSBS = 31;
  		open(IMPACTSBS, ">", "$folderFigure/Impact_protein_sequence/$k_file/$k_file-DistributionExoFunc.txt") or die "$!: $folderFigure/Impact_protein_sequence/$k_file/$k_file-DistributionExoFunc.txt\n";
  		print IMPACTSBS "AA_Change\tCount\tPercent\n";

  		# Pie chart with the distribution of SBS vs Indel
  		open(SBSINDEL, ">", "$folderFigure/Overall_mutation_distribution/$k_file/$k_file-OverallMutationDistribution.txt") or die "$!: $folderFigure/Overall_mutation_distribution/$k_file/$k_file-OverallMutationDistribution.txt\n";
  		print SBSINDEL "Variant_type\tCount\tPercent\n";
  		my ($deletion, $insertion) = (0, 0);


  		foreach my $k_exoFunc(sort keys $refH_file->{$k_file}{'ImpactSBS'})
  		{
  			my $percent = ($refH_file->{$k_file}{'ImpactSBS'}{$k_exoFunc} / $refH_file->{$k_file}{'TotalMutGenomic'})*100;
  			$percent = sprintf("%.2f", $percent);

  			if($k_exoFunc eq "NA") { print IMPACTSBS "Not_Applicable\t$percent\t$refH_file->{$k_file}{'ImpactSBS'}{$k_exoFunc}\n"; }
  			else { my $temp = $k_exoFunc; $temp =~ s/ /_/g;  print IMPACTSBS "$temp\t$percent\t$refH_file->{$k_file}{'ImpactSBS'}{$k_exoFunc}\n"; }

  			$ws->write($lImpactSBS, 6, $k_exoFunc, $table_left2); $ws->write($lImpactSBS, 8, $percent, $format_A10); $ws->write($lImpactSBS, 9, $refH_file->{$k_file}{'ImpactSBS'}{$k_exoFunc}, $table_right);
  			$lImpactSBS++;

  			# Pie chart with the distribution of SBS vs Indel
  			if($k_exoFunc    =~ /deletion/i)  { $deletion  += $refH_file->{$k_file}{'ImpactSBS'}{$k_exoFunc}; }
				elsif($k_exoFunc =~ /insertion/i) { $insertion += $refH_file->{$k_file}{'ImpactSBS'}{$k_exoFunc}; }
  		}
  		close IMPACTSBS;
  		$ws->write($lImpactSBS, 9, $refH_file->{$k_file}{'TotalMutGenomic'}, $table_bottomrightHeader);
  		$ws->write($lImpactSBS, 6, "", $table_bottomleft); $ws->write($lImpactSBS, 7, "", $table_bottom); $ws->write($lImpactSBS, 8, "", $table_bottom);

  		# Pie chart with the distribution of SBS vs Indel
  		my $percentSBSIndel = ($deletion/$refH_file->{$k_file}{'TotalMutGenomic'})*100; $percentSBSIndel = sprintf("%.2f", $percentSBSIndel);
  		print SBSINDEL "Deletion\t$deletion\t$percentSBSIndel\n";
  		$percentSBSIndel = ($insertion/$refH_file->{$k_file}{'TotalMutGenomic'})*100; $percentSBSIndel = sprintf("%.2f", $percentSBSIndel);
  		print SBSINDEL "Insertion\t$insertion\t$percentSBSIndel\n";
  		$percentSBSIndel = ($refH_file->{$k_file}{TotalSBSGenomic}/$refH_file->{$k_file}{'TotalMutGenomic'})*100; $percentSBSIndel = sprintf("%.2f", $percentSBSIndel);
  		print SBSINDEL "SBS\t$refH_file->{$k_file}{TotalSBSGenomic}\t$percentSBSIndel\n";
  		close SBSINDEL;

  		###########################################################################################################################################################
  		######################################################## SEQUENCE CONTEXT ON GENOMIC STRAND	###############################################################
  		###########################################################################################################################################################
  		my $row_SeqContext6 = 4;
  		# Count the total of mutations for 6 mutation types on genomic strand
  		my ($c_ca6_g, $c_cg6_g, $c_ct6_g, $c_ta6_g, $c_tc6_g, $c_tg6_g) = (0,0,0, 0,0,0);
  		my ($p_ca6_g, $p_cg6_g, $p_ct6_g, $p_ta6_g, $p_tc6_g, $p_tg6_g) = (0,0,0, 0,0,0);
  		my $maxValue = 0; # For the heatmap

  		# For checking if the total number of SBS is correct
  		my $total_SBS_genomic = 0;


  		open(HEATMAPCGENOMIC, ">", "$folderFigure/Trinucleotide_Sequence_Context/$k_file/$k_file-HeatmapCount-Genomic.txt") or die "$!: $folderFigure/Trinucleotide_Sequence_Context/$k_file/$k_file-HeatmapCount-Genomic.txt\n";
		  print HEATMAPCGENOMIC "\tC>A\tC>G\tC>T\tT>A\tT>C\tT>G\n";
		  open(HEATMAPPGENOMIC, ">", "$folderFigure/Trinucleotide_Sequence_Context/$k_file/$k_file-HeatmapPercent-Genomic.txt") or die "$!: $folderFigure/Trinucleotide_Sequence_Context/$k_file/$k_file-HeatmapPercent-Genomic.txt\n";
		  print HEATMAPPGENOMIC "\tC>A\tC>G\tC>T\tT>A\tT>C\tT>G\n";

		  ## Bar plot NMF like
		  open(BARPLOTNMFLIKE, ">", "$folderFigure/Trinucleotide_Sequence_Context/$k_file/$k_file-MutationSpectraPercent-Genomic.txt") or die "$!: $folderFigure/Trinucleotide_Sequence_Context/$k_file/$k_file-MutationSpectraPercent-Genomic.txt\n";
		  print BARPLOTNMFLIKE "alteration\tcontext\tvalue\n";

		  foreach my $k_context (sort keys $refH_file->{$k_file}{'SeqContextG'})
		  {
		  	if( ($k_context =~ /N/) || (length($k_context) != 3) ) { next; }

		  	# Write the context: 6 mut type on genomic strand
  			$ws->write($row_SeqContext6 , $colStart_matrixSeqContext+3, $k_context, $format_A10); $ws->write($row_SeqContext6 , $colStart_matrixSeqContext+13, $k_context, $format_A10);

  			foreach my $k_mutation (sort keys $refH_file->{$k_file}{'SeqContextG'}{$k_context})
  			{
  				# For checking the total number of SBS
  				$total_SBS_genomic += $refH_file->{$k_file}{'SeqContextG'}{$k_context}{$k_mutation};

  				# Calculate the percentages
  				my $percent = 0;
  				if($refH_file->{$k_file}{'SeqContextG'}{$k_context}{$k_mutation} == 0) { $percent = 0; }
  				else
  				{
  					$percent = ($refH_file->{$k_file}{'SeqContextG'}{$k_context}{$k_mutation} / $refH_file->{$k_file}{'TotalSBSGenomic'}) * 100;
  					$percent = sprintf("%.2f", $percent);
  				}

  				# For representing the sequence context with a bar plot (NMF like style)
  				print BARPLOTNMFLIKE $k_mutation,"\t", $k_context,"\t", $percent,"\n";

  				if($k_mutation eq "C>A")
  				{
  					### COUNT
  					$ws->write($row_SeqContext6, $colStart_matrixSeqContext+4, $refH_file->{$k_file}{'SeqContextG'}{$k_context}{$k_mutation}, $format_A10);
  					# Write the count for the heatmap
  					print HEATMAPCGENOMIC "$k_context\t$refH_file->{$k_file}{'SeqContextG'}{$k_context}{$k_mutation}\t";

  					### PERCENTAGE
  					$ws->write($row_SeqContext6, $colStart_matrixSeqContext+14, $percent, $format_A10);
  					print HEATMAPPGENOMIC "$k_context\t$percent\t";

  					# For NMF input
  					my $count = $refH_file->{$k_file}{'SeqContextG'}{$k_context}{$k_mutation};
						if($k_file ne "Pool_Data") { push(@{$h_inputNMF{'Count'}{$k_context}{'C>A'}}, $count); }
						if($k_file ne "Pool_Data") { push(@{$h_inputNMF{'Percent'}{$k_context}{'C>A'}}, $percent); }

  					# For the heatmap
  					if($percent >= $maxValue) { $maxValue = $percent; }

  					# For the total amount per mutation types
  					$c_ca6_g += $refH_file->{$k_file}{'SeqContextG'}{$k_context}{$k_mutation};
  					$p_ca6_g += $percent;
  				}
  				if($k_mutation eq "C>G")
  				{
  					### COUNT
  					$ws->write($row_SeqContext6, $colStart_matrixSeqContext+5, $refH_file->{$k_file}{'SeqContextG'}{$k_context}{$k_mutation}, $format_A10);
  					# Write the count for the heatmap
  					print HEATMAPCGENOMIC "$refH_file->{$k_file}{'SeqContextG'}{$k_context}{$k_mutation}\t";

  					### PERCENTAGE
  					$ws->write($row_SeqContext6, $colStart_matrixSeqContext+15, $percent, $format_A10);
  					print HEATMAPPGENOMIC "$percent\t";

  					# For NMF input
  					my $count = $refH_file->{$k_file}{'SeqContextG'}{$k_context}{$k_mutation};
						if($k_file ne "Pool_Data") { push(@{$h_inputNMF{'Count'}{$k_context}{'C>G'}}, $count); }
						if($k_file ne "Pool_Data") { push(@{$h_inputNMF{'Percent'}{$k_context}{'C>G'}}, $percent); }

  					# For the heatmap
  					if($percent >= $maxValue) { $maxValue = $percent; }

  					# For the total amount per mutation types
  					$c_cg6_g += $refH_file->{$k_file}{'SeqContextG'}{$k_context}{$k_mutation};
  					$p_cg6_g += $percent;
  				}
  				if($k_mutation eq "C>T")
  				{
  					### COUNT
  					$ws->write($row_SeqContext6, $colStart_matrixSeqContext+6, $refH_file->{$k_file}{'SeqContextG'}{$k_context}{$k_mutation}, $format_A10);
  					# Write the count for the heatmap
  					print HEATMAPCGENOMIC "$refH_file->{$k_file}{'SeqContextG'}{$k_context}{$k_mutation}\t";

  					### PERCENTAGE
  					$ws->write($row_SeqContext6, $colStart_matrixSeqContext+16, $percent, $format_A10);
  					print HEATMAPPGENOMIC "$percent\t";

  					# For NMF input
  					my $count = $refH_file->{$k_file}{'SeqContextG'}{$k_context}{$k_mutation};
						if($k_file ne "Pool_Data") { push(@{$h_inputNMF{'Count'}{$k_context}{'C>T'}}, $count); }
						if($k_file ne "Pool_Data") { push(@{$h_inputNMF{'Percent'}{$k_context}{'C>T'}}, $percent); }

  					# For the heatmap
  					if($percent >= $maxValue) { $maxValue = $percent; }

  					# For the total amount per mutation types
  					$c_ct6_g += $refH_file->{$k_file}{'SeqContextG'}{$k_context}{$k_mutation};
  					$p_ct6_g += $percent;
  				}
  				if($k_mutation eq "T>A")
  				{
  					### COUNT
  					$ws->write($row_SeqContext6, $colStart_matrixSeqContext+7, $refH_file->{$k_file}{'SeqContextG'}{$k_context}{$k_mutation}, $format_A10);
  					# Write the count for the heatmap
  					print HEATMAPCGENOMIC "$refH_file->{$k_file}{'SeqContextG'}{$k_context}{$k_mutation}\t";

  					### PERCENTAGE
  					$ws->write($row_SeqContext6, $colStart_matrixSeqContext+17, $percent, $format_A10);
  					print HEATMAPPGENOMIC "$percent\t";

  					# For NMF input
  					my $count = $refH_file->{$k_file}{'SeqContextG'}{$k_context}{$k_mutation};
						if($k_file ne "Pool_Data") { push(@{$h_inputNMF{'Count'}{$k_context}{'T>A'}}, $count); }
						if($k_file ne "Pool_Data") { push(@{$h_inputNMF{'Percent'}{$k_context}{'T>A'}}, $percent); }

  					# For the heatmap
  					if($percent >= $maxValue) { $maxValue = $percent; }

  					# For the total amount per mutation types
  					$c_ta6_g += $refH_file->{$k_file}{'SeqContextG'}{$k_context}{$k_mutation};
  					$p_ta6_g += $percent;
  				}
  				if($k_mutation eq "T>C")
  				{
  					### COUNT
  					$ws->write($row_SeqContext6, $colStart_matrixSeqContext+8, $refH_file->{$k_file}{'SeqContextG'}{$k_context}{$k_mutation}, $format_A10);
  					# Write the count for the heatmap
  					print HEATMAPCGENOMIC "$refH_file->{$k_file}{'SeqContextG'}{$k_context}{$k_mutation}\t";

  					### PERCENTAGE
  					$ws->write($row_SeqContext6, $colStart_matrixSeqContext+18, $percent, $format_A10);
  					print HEATMAPPGENOMIC "$percent\t";

  					# For NMF input
  					my $count = $refH_file->{$k_file}{'SeqContextG'}{$k_context}{$k_mutation};
						if($k_file ne "Pool_Data") { push(@{$h_inputNMF{'Count'}{$k_context}{'T>C'}}, $count); }
						if($k_file ne "Pool_Data") { push(@{$h_inputNMF{'Percent'}{$k_context}{'T>C'}}, $percent); }

  					# For the heatmap
  					if($percent >= $maxValue) { $maxValue = $percent; }

  					# For the total amount per mutation types
  					$c_tc6_g += $refH_file->{$k_file}{'SeqContextG'}{$k_context}{$k_mutation};
  					$p_tc6_g += $percent;
  				}
  				if($k_mutation eq "T>G")
  				{
  					### COUNT
  					$ws->write($row_SeqContext6, $colStart_matrixSeqContext+9, $refH_file->{$k_file}{'SeqContextG'}{$k_context}{$k_mutation}, $format_A10);
  					# Write the count for the heatmap
  					print HEATMAPCGENOMIC "$refH_file->{$k_file}{'SeqContextG'}{$k_context}{$k_mutation}\n";

  					### PERCENTAGE
  					$ws->write($row_SeqContext6, $colStart_matrixSeqContext+19, $percent, $format_A10);
  					print HEATMAPPGENOMIC "$percent\n";

  					# For NMF input
  					my $count = $refH_file->{$k_file}{'SeqContextG'}{$k_context}{$k_mutation};
						if($k_file ne "Pool_Data") { push(@{$h_inputNMF{'Count'}{$k_context}{'T>G'}}, $count); }
						if($k_file ne "Pool_Data") { push(@{$h_inputNMF{'Percent'}{$k_context}{'T>G'}}, $percent); }

  					# For the heatmap
  					if($percent >= $maxValue) { $maxValue = $percent; }

  					# For the total amount per mutation types
  					$c_tg6_g += $refH_file->{$k_file}{'SeqContextG'}{$k_context}{$k_mutation};
  					$p_tg6_g += $percent;
  				}
  			}
  			$row_SeqContext6++;
		  }
		  close HEATMAPCGENOMIC; close HEATMAPPGENOMIC;
		  close BARPLOTNMFLIKE;


		  # Write the total number of SBS per mutation type: COUNT
		  $ws->write($row_SeqContext6, $colStart_matrixSeqContext+4, $c_ca6_g, $formatT_bottomHeader2); $ws->write($row_SeqContext6, $colStart_matrixSeqContext+5, $c_cg6_g, $formatT_bottomHeader2); $ws->write($row_SeqContext6, $colStart_matrixSeqContext+6, $c_ct6_g, $formatT_bottomHeader2);
		  $ws->write($row_SeqContext6, $colStart_matrixSeqContext+7, $c_ta6_g, $formatT_bottomHeader2); $ws->write($row_SeqContext6, $colStart_matrixSeqContext+8, $c_tc6_g, $formatT_bottomHeader2); $ws->write($row_SeqContext6, $colStart_matrixSeqContext+9, $c_tg6_g, $formatT_bottomHeader2);
  		if($total_SBS_genomic != $refH_file->{$k_file}{'TotalSBSGenomic'}) { print STDERR "Error in the calculation of the total number of SBS on the genomic strand!!!!\nFrom hash table $refH_file->{$k_file}{'TotalSBSGenomic'}\tVS\t$total_SBS_genomic\n"; exit; }

		  # Write the total number of SBS per mutation type: PERCENT
		  $ws->write($row_SeqContext6, $colStart_matrixSeqContext+14, $p_ca6_g, $formatT_bottomHeader2); $ws->write($row_SeqContext6, $colStart_matrixSeqContext+15, $p_cg6_g, $formatT_bottomHeader2); $ws->write($row_SeqContext6, $colStart_matrixSeqContext+16, $p_ct6_g, $formatT_bottomHeader2);
		  $ws->write($row_SeqContext6, $colStart_matrixSeqContext+17, $p_ta6_g, $formatT_bottomHeader2); $ws->write($row_SeqContext6, $colStart_matrixSeqContext+18, $p_tc6_g, $formatT_bottomHeader2); $ws->write($row_SeqContext6, $colStart_matrixSeqContext+19, $p_tg6_g, $formatT_bottomHeader2);
		  my $totalPercent_genomic = $p_ca6_g + $p_cg6_g + $p_ct6_g + $p_ta6_g + $p_tc6_g + $p_tg6_g; $totalPercent_genomic = sprintf("%.0f", $totalPercent_genomic);
		  if($totalPercent_genomic != 100) { print STDERR "Error in the calculation of the total percentages on the genomic strand!!!\nThe total is equal to=\t$totalPercent_genomic\n"; exit; }


		  #----------------------------------------------------------------------------------------------------------------------------------------------------------------#
		  # For the input matrix for NMF
		  if($k_file ne "Pool_Data") { push(@{$h_inputNMF{'Sample'}}, $k_file); }


			###########################################################################################################################################################
  		######################################################## SEQUENCE CONTEXT ON CODING STRAND	###############################################################
  		###########################################################################################################################################################
  		my $row_SeqContext12 = $rowStart_SBSdistrBySeg+6; my $row_SeqContext12Percent = $rowStart_SBSdistrBySeg+27;
  		# Reset the total count and percent calculated for the strand bias
  		($ca_NonTr, $ca_Tr, $cg_NonTr, $cg_Tr, $ct_NonTr, $ct_Tr, $ta_NonTr, $ta_Tr, $tc_NonTr, $tc_Tr, $tg_NonTr, $tg_Tr) = (0,0,0, 0,0,0, 0,0,0, 0,0,0);
  		($percent_ca_NonTr, $percent_ca_Tr, $percent_cg_NonTr, $percent_cg_Tr, $percent_ct_NonTr, $percent_ct_Tr, $percent_ta_NonTr, $percent_ta_Tr, $percent_tc_NonTr, $percent_tc_Tr, $percent_tg_NonTr, $percent_tg_Tr) = (0,0,0, 0,0,0, 0,0,0, 0,0,0);

  		# For checking if the total number of SBS is correct
  		my $total_SBS_coding  = 0;

  		open(COUNT, ">", "$folderFigure/Stranded_Analysis/$k_file/$k_file-StrandedSignatureCount.txt") or die "$!: $folderFigure/Stranded_Analysis/$k_file/$k_file-StrandedSignatureCount.txt\n";
		  print COUNT "MutationTypeContext\tStrand\tValue\tSample\n";
		  open(PERCENT, ">", "$folderFigure/Stranded_Analysis/$k_file/$k_file-StrandedSignaturePercent.txt") or die "$!: $folderFigure/Stranded_Analysis/$k_file/$k_file-StrandedSignaturePercent.txt\n";
		  print PERCENT "MutationTypeContext\tStrand\tValue\tSample\n";

  		foreach my $k_context (sort keys $refH_file->{$k_file}{'SeqContextC'})
  		{
  			if( ($k_context =~ /N/) || (length($k_context) != 3) ) { next; }

  			# Write the context: 12 mut type on coding strand
  			$ws->write($row_SeqContext12 , $colStart_matrixSeqContext, $k_context, $formatT_left); $ws->write($row_SeqContext12Percent , $colStart_matrixSeqContext, $k_context, $formatT_left);

  			foreach my $k_mutation (sort keys $refH_file->{$k_file}{'SeqContextC'}{$k_context})
  			{
  				# Percent: 12 mut type on coding strand
  				my ($percent_NonTr, $percent_Tr) = (0, 0);
  				if($refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'NonTr'} == 0) { $percent_NonTr = 0; }
  				else { $percent_NonTr = ( $refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'NonTr'} / $refH_file->{$k_file}{'TotalSBSCoding'} ) * 100 }
  				if($refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'Tr'}    == 0) { $percent_Tr    = 0; }
  				else { $percent_Tr = ( $refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'Tr'} / $refH_file->{$k_file}{'TotalSBSCoding'} ) * 100 }


  				# Calculate the total number for each mutation types
  				if($k_mutation eq "C>A")
  				{
						$ca_NonTr += $refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'NonTr'};
						$ca_Tr    += $refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'Tr'};

	  				# COUNT : 12 mutation type (stranded bar graph)
  					$ws->write($row_SeqContext12, $colStart_matrixSeqContext+1, $refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'NonTr'}, $format_A10);
  					print COUNT "$k_mutation:$k_context\tNonTranscribed\t$refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'NonTr'}\t$k_file\n";
  					$ws->write($row_SeqContext12, $colStart_matrixSeqContext+2, $refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'Tr'}, $format_A10);
  					print COUNT "$k_mutation:$k_context\tTranscribed\t$refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'Tr'}\t$k_file\n";

	  				## PERCENT : 12 mutation type (stranded bar graph)
	  				my $percent_NonTr = ($refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'NonTr'}/$refH_file->{$k_file}{'TotalSBSCoding'})*100;
	  				$percent_NonTr    = sprintf("%.2f", $percent_NonTr); $percent_ca_NonTr += $percent_NonTr;
	  				my $percent_Tr    = ($refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'Tr'}/$refH_file->{$k_file}{'TotalSBSCoding'})*100;
	  				$percent_Tr       = sprintf("%.2f", $percent_Tr); $percent_ca_Tr += $percent_Tr;
	  				print PERCENT "$k_mutation:$k_context\tNonTranscribed\t$percent_NonTr\t$k_file\n";
	  				print PERCENT "$k_mutation:$k_context\tTranscribed\t$percent_Tr\t$k_file\n";

	  				$ws->write($row_SeqContext12Percent, $colStart_matrixSeqContext+1, $percent_NonTr, $format_A10);
	  				$ws->write($row_SeqContext12Percent, $colStart_matrixSeqContext+2, $percent_Tr, $format_A10);
  				}
  				if($k_mutation eq "C>G")
  				{
						$cg_NonTr += $refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'NonTr'};
						$cg_Tr    += $refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'Tr'};

	  				# COUNT : 12 mutation type (stranded bar graph)
  					$ws->write($row_SeqContext12, $colStart_matrixSeqContext+3, $refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'NonTr'}, $format_A10);
  					print COUNT "$k_mutation:$k_context\tNonTranscribed\t$refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'NonTr'}\t$k_file\n";
  					$ws->write($row_SeqContext12, $colStart_matrixSeqContext+4, $refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'Tr'}, $format_A10);
  					print COUNT "$k_mutation:$k_context\tTranscribed\t$refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'Tr'}\t$k_file\n";

	  				## PERCENT : 12 mutation type (stranded bar graph)
	  				my $percent_NonTr = ($refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'NonTr'}/$refH_file->{$k_file}{'TotalSBSCoding'})*100;
	  				$percent_NonTr    = sprintf("%.2f", $percent_NonTr); $percent_cg_NonTr += $percent_NonTr;
	  				my $percent_Tr    = ($refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'Tr'}/$refH_file->{$k_file}{'TotalSBSCoding'})*100;
	  				$percent_Tr       = sprintf("%.2f", $percent_Tr); $percent_cg_Tr += $percent_Tr;
	  				print PERCENT "$k_mutation:$k_context\tNonTranscribed\t$percent_NonTr\t$k_file\n";
	  				print PERCENT "$k_mutation:$k_context\tTranscribed\t$percent_Tr\t$k_file\n";

	  				$ws->write($row_SeqContext12Percent, $colStart_matrixSeqContext+3, $percent_NonTr, $format_A10);
	  				$ws->write($row_SeqContext12Percent, $colStart_matrixSeqContext+4, $percent_Tr, $format_A10);
  				}
  				if($k_mutation eq "C>T")
  				{
						$ct_NonTr += $refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'NonTr'};
						$ct_Tr    += $refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'Tr'};

	  				# COUNT : 12 mutation type (stranded bar graph)
  					$ws->write($row_SeqContext12, $colStart_matrixSeqContext+5, $refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'NonTr'}, $format_A10);
  					print COUNT "$k_mutation:$k_context\tNonTranscribed\t$refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'NonTr'}\t$k_file\n";
  					$ws->write($row_SeqContext12, $colStart_matrixSeqContext+6, $refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'Tr'}, $format_A10);
  					print COUNT "$k_mutation:$k_context\tTranscribed\t$refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'Tr'}\t$k_file\n";

	  				## PERCENT : 12 mutation type (stranded bar graph)
	  				my $percent_NonTr = ($refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'NonTr'}/$refH_file->{$k_file}{'TotalSBSCoding'})*100;
	  				$percent_NonTr    = sprintf("%.2f", $percent_NonTr); $percent_ct_NonTr += $percent_NonTr;
	  				my $percent_Tr    = ($refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'Tr'}/$refH_file->{$k_file}{'TotalSBSCoding'})*100;
	  				$percent_Tr       = sprintf("%.2f", $percent_Tr); $percent_ct_Tr += $percent_Tr;
	  				print PERCENT "$k_mutation:$k_context\tNonTranscribed\t$percent_NonTr\t$k_file\n";
	  				print PERCENT "$k_mutation:$k_context\tTranscribed\t$percent_Tr\t$k_file\n";

	  				$ws->write($row_SeqContext12Percent, $colStart_matrixSeqContext+5, $percent_NonTr, $format_A10);
	  				$ws->write($row_SeqContext12Percent, $colStart_matrixSeqContext+6, $percent_Tr, $format_A10);
  				}
  				if($k_mutation eq "T>A")
  				{
						$ta_NonTr += $refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'NonTr'};
						$ta_Tr    += $refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'Tr'};

	  				# COUNT : 12 mutation type (stranded bar graph)
  					$ws->write($row_SeqContext12, $colStart_matrixSeqContext+7, $refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'NonTr'}, $format_A10);
  					print COUNT "$k_mutation:$k_context\tNonTranscribed\t$refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'NonTr'}\t$k_file\n";
  					$ws->write($row_SeqContext12, $colStart_matrixSeqContext+8, $refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'Tr'}, $format_A10);
  					print COUNT "$k_mutation:$k_context\tTranscribed\t$refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'Tr'}\t$k_file\n";

	  				## PERCENT : 12 mutation type (stranded bar graph)
	  				my $percent_NonTr = ($refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'NonTr'}/$refH_file->{$k_file}{'TotalSBSCoding'})*100;
	  				$percent_NonTr    = sprintf("%.2f", $percent_NonTr); $percent_ta_NonTr += $percent_NonTr;
	  				my $percent_Tr    = ($refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'Tr'}/$refH_file->{$k_file}{'TotalSBSCoding'})*100;
	  				$percent_Tr       = sprintf("%.2f", $percent_Tr); $percent_ta_Tr += $percent_Tr;
	  				print PERCENT "$k_mutation:$k_context\tNonTranscribed\t$percent_NonTr\t$k_file\n";
	  				print PERCENT "$k_mutation:$k_context\tTranscribed\t$percent_Tr\t$k_file\n";

	  				$ws->write($row_SeqContext12Percent, $colStart_matrixSeqContext+7, $percent_NonTr, $format_A10);
	  				$ws->write($row_SeqContext12Percent, $colStart_matrixSeqContext+8, $percent_Tr, $format_A10);
  				}
  				if($k_mutation eq "T>C")
  				{
						$tc_NonTr += $refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'NonTr'};
						$tc_Tr    += $refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'Tr'};

	  				# COUNT : 12 mutation type (stranded bar graph)
  					$ws->write($row_SeqContext12, $colStart_matrixSeqContext+9, $refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'NonTr'}, $format_A10);
  					print COUNT "$k_mutation:$k_context\tNonTranscribed\t$refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'NonTr'}\t$k_file\n";
  					$ws->write($row_SeqContext12, $colStart_matrixSeqContext+10, $refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'Tr'}, $format_A10);
  					print COUNT "$k_mutation:$k_context\tTranscribed\t$refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'Tr'}\t$k_file\n";

	  				## PERCENT : 12 mutation type (stranded bar graph)
	  				my $percent_NonTr = ($refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'NonTr'}/$refH_file->{$k_file}{'TotalSBSCoding'})*100;
	  				$percent_NonTr    = sprintf("%.2f", $percent_NonTr); $percent_tc_NonTr += $percent_NonTr;
	  				my $percent_Tr    = ($refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'Tr'}/$refH_file->{$k_file}{'TotalSBSCoding'})*100;
	  				$percent_Tr       = sprintf("%.2f", $percent_Tr); $percent_tc_Tr += $percent_Tr;
	  				print PERCENT "$k_mutation:$k_context\tNonTranscribed\t$percent_NonTr\t$k_file\n";
	  				print PERCENT "$k_mutation:$k_context\tTranscribed\t$percent_Tr\t$k_file\n";

	  				$ws->write($row_SeqContext12Percent, $colStart_matrixSeqContext+9, $percent_NonTr, $format_A10);
	  				$ws->write($row_SeqContext12Percent, $colStart_matrixSeqContext+10, $percent_Tr, $format_A10);
  				}
  				if($k_mutation eq "T>G")
  				{
						$tg_NonTr += $refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'NonTr'};
						$tg_Tr    += $refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'Tr'};

	  				# COUNT : 12 mutation type (stranded bar graph)
  					$ws->write($row_SeqContext12, $colStart_matrixSeqContext+11, $refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'NonTr'}, $format_A10);
  					print COUNT "$k_mutation:$k_context\tNonTranscribed\t$refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'NonTr'}\t$k_file\n";
  					$ws->write($row_SeqContext12, $colStart_matrixSeqContext+12, $refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'Tr'}, $format_A10);
  					print COUNT "$k_mutation:$k_context\tTranscribed\t$refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'Tr'}\t$k_file\n";

	  				## PERCENT : 12 mutation type (stranded bar graph)
	  				my $percent_NonTr = ($refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'NonTr'}/$refH_file->{$k_file}{'TotalSBSCoding'})*100;
	  				$percent_NonTr    = sprintf("%.2f", $percent_NonTr); $percent_tg_NonTr += $percent_NonTr;
	  				my $percent_Tr    = ($refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'Tr'}/$refH_file->{$k_file}{'TotalSBSCoding'})*100;
	  				$percent_Tr       = sprintf("%.2f", $percent_Tr); $percent_tg_Tr += $percent_Tr;
	  				print PERCENT "$k_mutation:$k_context\tNonTranscribed\t$percent_NonTr\t$k_file\n";
	  				print PERCENT "$k_mutation:$k_context\tTranscribed\t$percent_Tr\t$k_file\n";

	  				$ws->write($row_SeqContext12Percent, $colStart_matrixSeqContext+11, $percent_NonTr, $format_A10);
	  				$ws->write($row_SeqContext12Percent, $colStart_matrixSeqContext+12, $percent_Tr, $format_A10);
  				}

  				# For checking if the total number of SBS is correct
  				$total_SBS_coding  += $refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'NonTr'} + $refH_file->{$k_file}{'SeqContextC'}{$k_context}{$k_mutation}{'Tr'};
  			}
  			$row_SeqContext12++; $row_SeqContext12Percent++;
  		}
  		close COUNT; close PERCENT;

  		## Write the total of each mutation types: 12 mut type on coding strand
			$ws->write($row_SeqContext12, $colStart_matrixSeqContext+1, $ca_NonTr, $formatT_bottomHeader2); $ws->write($row_SeqContext12, $colStart_matrixSeqContext+2, $ca_Tr, $formatT_bottomHeader2);
			$ws->write($row_SeqContext12, $colStart_matrixSeqContext+3, $cg_NonTr, $formatT_bottomHeader2); $ws->write($row_SeqContext12, $colStart_matrixSeqContext+4, $cg_Tr, $formatT_bottomHeader2);
			$ws->write($row_SeqContext12, $colStart_matrixSeqContext+5, $ct_NonTr, $formatT_bottomHeader2); $ws->write($row_SeqContext12, $colStart_matrixSeqContext+6, $ct_Tr, $formatT_bottomHeader2);
			$ws->write($row_SeqContext12, $colStart_matrixSeqContext+7, $ta_NonTr, $formatT_bottomHeader2); $ws->write($row_SeqContext12, $colStart_matrixSeqContext+8, $ta_Tr, $formatT_bottomHeader2);
			$ws->write($row_SeqContext12, $colStart_matrixSeqContext+9, $tc_NonTr, $formatT_bottomHeader2); $ws->write($row_SeqContext12, $colStart_matrixSeqContext+10, $tc_Tr, $formatT_bottomHeader2);
			$ws->write($row_SeqContext12, $colStart_matrixSeqContext+11, $tg_NonTr, $formatT_bottomHeader2); $ws->write($row_SeqContext12, $colStart_matrixSeqContext+12, $tg_Tr, $formatT_bottomHeader2);
  		# Write the total percentages of each mutation types: 12 mut type on coding strand
  		$ws->write($row_SeqContext12Percent, $colStart_matrixSeqContext+1, $percent_ca_NonTr, $formatT_bottomHeader); $ws->write($row_SeqContext12Percent, $colStart_matrixSeqContext+2, $percent_ca_Tr, $formatT_bottomHeader);
  		$ws->write($row_SeqContext12Percent, $colStart_matrixSeqContext+3, $percent_cg_NonTr, $formatT_bottomHeader); $ws->write($row_SeqContext12Percent, $colStart_matrixSeqContext+4, $percent_cg_Tr, $formatT_bottomHeader);
  		$ws->write($row_SeqContext12Percent, $colStart_matrixSeqContext+5, $percent_ct_NonTr, $formatT_bottomHeader); $ws->write($row_SeqContext12Percent, $colStart_matrixSeqContext+6, $percent_ct_Tr, $formatT_bottomHeader);
  		$ws->write($row_SeqContext12Percent, $colStart_matrixSeqContext+7, $percent_ta_NonTr, $formatT_bottomHeader); $ws->write($row_SeqContext12Percent, $colStart_matrixSeqContext+8, $percent_ta_Tr, $formatT_bottomHeader);
  		$ws->write($row_SeqContext12Percent, $colStart_matrixSeqContext+9, $percent_tc_NonTr, $formatT_bottomHeader); $ws->write($row_SeqContext12Percent, $colStart_matrixSeqContext+10, $percent_tc_Tr, $formatT_bottomHeader);
  		$ws->write($row_SeqContext12Percent, $colStart_matrixSeqContext+11, $percent_tg_NonTr, $formatT_bottomHeader); $ws->write($row_SeqContext12Percent, $colStart_matrixSeqContext+12, $percent_tg_Tr, $formatT_bottomHeader);

  		if($total_SBS_coding == $refH_file->{$k_file}{'TotalSBSCoding'}) { $ws->write($row_SeqContext12, $colStart_matrixSeqContext+13, $refH_file->{$k_file}{'TotalSBSCoding'}, $formatT_bottomHeader2) }
  		else { print STDERR "Error in the calculation of the total number of SBS on the coding strand!!!!\nFrom hash table $refH_file->{$k_file}{'TotalSBSCoding'}\tVS\t$total_SBS_coding\n"; exit; }


  		my $totalP_SBS_coding = $percent_ca_NonTr + $percent_ca_Tr + $percent_cg_NonTr + $percent_cg_Tr + $percent_ct_NonTr + $percent_ct_Tr + $percent_ta_NonTr + $percent_ta_Tr + $percent_tc_NonTr + $percent_tc_Tr + $percent_tg_NonTr + $percent_tg_Tr; $totalP_SBS_coding = sprintf("%.0f", $totalP_SBS_coding);
  		if($totalP_SBS_coding != 100) { print STDERR "The percentages for the trinucleotide sequence context on the coding strand for 12 mutation types is not equal to 100!!!\n$totalP_SBS_coding\n"; exit; }


  		###########################################################################################################################################################
  		################################################################### GRAPHS & TABLES	#######################################################################
  		###########################################################################################################################################################
  		Create_Graph($folderFigure, $k_file, $maxValue);

  		## Distribution of SBS into the Excel report (Figure 1 + Table 1)
			$ws->write(0, 0, "Graph 1. SBS distribution", $formatT_graphTitle); $ws->set_row(0, 18);
			$ws->insert_image(1, 0, "$folder_temp/$k_file-SBS_distribution-Report.png", 0, 0, .2, .2);
			$ws->write(29, 0, "Table 1. Frequency and counts of all SBS", $format_A10Boldleft);
			$ws->write(30, 0, "Mutation type", $table_topleft); $ws->write(30, 1, "Percentage", $table_top); $ws->write(30, 2, "Count", $table_topRight);
			$ws->write(31, 0, "C:G>A:T", $table_left); $ws->write(31, 1, $percent_ca, $format_A10); $ws->write(31, 2, $ca_genomique, $table_right);
			$ws->write(32, 0, "C:G>G:C", $table_left); $ws->write(32, 1, $percent_cg, $format_A10); $ws->write(32, 2, $cg_genomique, $table_right);
			$ws->write(33, 0, "C:G>T:A", $table_left); $ws->write(33, 1, $percent_ct, $format_A10); $ws->write(33, 2, $ct_genomique, $table_right);
			$ws->write(34, 0, "T:A>A:T", $table_left); $ws->write(34, 1, $percent_ta, $format_A10); $ws->write(34, 2, $ta_genomique, $table_right);
			$ws->write(35, 0, "T:A>C:G", $table_left); $ws->write(35, 1, $percent_tc, $format_A10); $ws->write(35, 2, $tc_genomique, $table_right);
			$ws->write(36, 0, "T:A>G:C", $table_left); $ws->write(36, 1, $percent_tg, $format_A10); $ws->write(36, 2, $tg_genomique, $table_right);
			$ws->write(37, 0, "", $table_bottomleft); $ws->write(37, 1, "", $table_bottom); $ws->write(37, 2, $refH_file->{$k_file}{'TotalSBSGenomic'}, $table_bottomrightHeader);

			## Impact of the SBS on the protein
			$ws->write(0, 6, "Graph 2. Impact on protein sequence", $formatT_graphTitle);
			$ws->insert_image(1, 6, "$folder_temp/$k_file-DistributionExoFunc-Report.png", 0, 0, .2, .2);

			## Strand Bias
			$ws->write(0, 11, "Graph 3. Stranded distribution of SBS", $formatT_graphTitle);
			$ws->insert_image(1, 11, "$folder_temp/$k_file-StrandBias-Report.png", 0, 0, .2, .2);

			## Stranded signature (Scale the inserted image: width x 0.7, height x 0.8)
			$ws->insert_image($rowStart_SBSdistrBySeg+3, $colStart_matrixSeqContext+15, "$folder_temp/$k_file-StrandedSignatureCount-Report.png", 0, 0, .16, .16);
			$ws->insert_image($rowStart_SBSdistrBySeg+24, $colStart_matrixSeqContext+15, "$folder_temp/$k_file-StrandedSignaturePercent-Report.png", 0, 0, .16, .16);


			# Heatamp for the sequence context on the genomic strand (6 mutation types)
			$ws->insert_image(4, $colStart_matrixSeqContext, "$folder_temp/$k_file-HeatmapCount-Genomic-Report.png");
			$ws->insert_image(4, $colStart_matrixSeqContext+10, "$folder_temp/$k_file-HeatmapPercent-Genomic-Report.png");


			## Bar plot for representing the sequence context (NMF like style)
			`Rscript $pathRScriptMutSpectrum $folderFigure/Trinucleotide_Sequence_Context/$k_file/$k_file-MutationSpectraPercent-Genomic.txt $k_file $folderFigure/Trinucleotide_Sequence_Context/$k_file $folder_temp $c_ca6_g $c_cg6_g $c_ct6_g $c_ta6_g $c_tc6_g $c_tg6_g 2>&1`;

			# Bar plot for the sequence context on the genomic strand (6 mutation types)
			$ws->insert_image(27, $colStart_matrixSeqContext+3, "$folder_temp/$k_file-MutationSpectraPercent-Genomic-Report.png");

			# Next sample
			$row_SumSheet++;
		} # End $k_file

		#----------------------------------------------------------------------------------------------------------------------------------------------------------------#
		# Write the input matrix for NMF
		open(OUTINPUTNMFC, ">", "$folderNMF/Input_NMF_Count.txt") or die "$!: $folderNMF/Input_NMF_Count.txt\n"; # with the count
		open(OUTINPUTNMFP, ">", "$folderNMF/Input_NMF_Frequency.txt") or die "$!: $folderNMF/Input_NMF_Frequency.txt\n"; # With the frequency un-normalized

		foreach my $k_sample (@{$h_inputNMF{'Sample'}}) { print OUTINPUTNMFC "\t$k_sample"; print OUTINPUTNMFP "\t$k_sample"; }
		print OUTINPUTNMFC "\n"; print OUTINPUTNMFP "\n";

		my $row_inputNMF = 1;
		foreach my $k_context (sort keys $h_inputNMF{'Count'})
		{
			$k_context =~ /(\w)_(\w)/; my ($base5, $base3) = ($1, $2);
			foreach my $k_mutation (sort keys $h_inputNMF{'Count'}{$k_context})
			{
				my ($col_inputNMF_Count, $col_inputNMF_Percent) = (1, 1);
			  my $contextNMF = $base5."[$k_mutation]".$base3;
			  $ws_inputNMF_count->write($row_inputNMF, 0, $contextNMF); $ws_inputNMF_percent->write($row_inputNMF, 0, $contextNMF);
			  print OUTINPUTNMFC $contextNMF,"\t"; print OUTINPUTNMFP $contextNMF,"\t";

			  foreach (@{$h_inputNMF{'Count'}{$k_context}{$k_mutation}})   { print OUTINPUTNMFC "$_\t"; } print OUTINPUTNMFC "\n";
			  foreach (@{$h_inputNMF{'Percent'}{$k_context}{$k_mutation}}) { print OUTINPUTNMFP "$_\t"; } print OUTINPUTNMFP "\n";

			  foreach (@{$h_inputNMF{'Count'}{$k_context}{$k_mutation}})
			  {
			  	# print "\t$k_context\t$k_mutation\t";
			  	# print "\t$row_inputNMF\t$col_inputNMF_Count\t$_\n";
			  	$ws_inputNMF_count->write($row_inputNMF, $col_inputNMF_Count, $_); $col_inputNMF_Count++;
			  }
			  foreach (@{$h_inputNMF{'Percent'}{$k_context}{$k_mutation}}) { $ws_inputNMF_percent->write($row_inputNMF, $col_inputNMF_Percent, $_); $col_inputNMF_Percent++; }
			  $row_inputNMF++;
			}
		}
		close OUTINPUTNMFP; close OUTINPUTNMFC;


		# Close the workbook
		$wb->close();
	}
	# Calculate the chi2 for the strand bias
	sub CalculateChi2
	{
		my ($refH_file, $folderChi2) = @_;

		# No value for the chi2
		if(scalar (keys $refH_file) == 0) { print STDERR "No value for calculating the chi2 for the strand bias\n"; exit; }

		# Strand bias for one mutation type for all the samples
		my %h_tempchi2 = ();
		my ($ca_NonTr, $ca_Tr, $cg_NonTr, $cg_Tr, $ct_NonTr, $ct_Tr, $ta_NonTr, $ta_Tr, $tc_NonTr, $tc_Tr, $tg_NonTr, $tg_Tr) = (0,0,0,0,0,0, 0,0,0,0,0,0);

		my $nb_file = 0;

		foreach my $k_file (sort keys $refH_file)
		{
			$nb_file++;
			foreach my $k_func (sort keys $refH_file->{$k_file}{'6mutType'})
			{
				foreach my $k_mutation (sort keys $refH_file->{$k_file}{'6mutType'}{$k_func})
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
			foreach my $k_file (sort keys $h_tempchi2{$k_mutation})
			{
				print CHI2 "$k_file\t$h_tempchi2{$k_mutation}{$k_file}{'NonTr'}\t$h_tempchi2{$k_mutation}{$k_file}{'Tr'}\t$k_mutation\n";
			}
		}
		close CHI2;


		# Open the connection with R
		my $R = Statistics::R->new() or die "Impossible to create a communication bridge with R\n";

		$R->send(qq`## Load the data. There is one column with the mutation type and the sample name but it's just for knowing what is corresponding to each line. The two columns with the number of variant per strand would be sufficient.
					      strBias<-read.delim("$folderChi2/Input_chi2_strandBias.txt", dec=".");`);
		$R->send(q`# Chi2
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
											confInt[k]        <- NA
											proportion[k]     <- NA
											sampleSize[k]     <- NA
											pValChi2_round[k] <- NA
											k <- k+1
								    }
									}
								}
								# Adjust with FDR
								FDR<-p.adjust(pValChi2, method="BH")

								# Rount the p-value
								for(i in 1:nrow(strBias))
								{
									 if( (! is.na(pValChi2[i])) && (pValChi2[i] < 0.0001) )
									 {
									   pValChi2_round[i] <- format(pValChi2[i], scientific=T, digits=3)
									 } else if(! is.na(pValChi2[i]))
									 {
									   pValChi2_round[i] <- as.character(round(pValChi2[i], 3))
									 }
								}

								# The option for the pool is specified
								if(!is.null(pValChi2_PoolData))
								{
									# Round the p-value for Pool_Data
									for(i in 1:6)
									{
									  if( (! is.na(pValChi2_PoolData[i])) && (pValChi2_PoolData[i] < 0.0001) )
									  {
									    pValChi2_PoolData_Round[i] <- format(pValChi2_PoolData[i], scientific=T, digits=3)
									  } else if(! is.na(pValChi2_PoolData[i]))
									  {
									    pValChi2_PoolData_Round[i] <- as.character(round(pValChi2_PoolData[i], 3))
									  }
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
								colnames(matrix)[5] <- "FDR"`);

		$R->send(qq`# Export the file
								# dec=".": Set the separator for the decimal by "."
								write.table(matrix,file="$folderChi2/Output_chi2_strandBias.txt",quote = FALSE,sep="\t",row.names = FALSE,dec=".");`);

		# Stop the connection with R
		$R->stop();
	}
	# Pearson correlation
	sub PearsonCoefficient
	{
		our ($refH_file, $filename) = @_;

		#### Calculate the Pearson coefficient
		my @total_SBS = (); # Pearson for all mutation types

		# Create a 2D array
		foreach my $k_mutation (sort keys $refH_file->{$filename}{'SBSPerChr'})
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
				if( $testZero == keys $refH_file->{$filename}{'SBSPerChr'}{$k_mutation}{'CHR'} ) { $correlation = 0; }
				# Pass the 2D array to the correlation subroutine
				else { $correlation = correlation($x); }

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
		if($testZero == keys $refH_file->{$filename}{'SBSPerChr'}{'TotalPerChr'}) { $correlation = 0; }
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
	# Complement bases (for the sequence context)
	sub complement
	{
		if($_[0] eq "A") { return "T"; }
		if($_[0] eq "C") { return "G"; }
		if($_[0] eq "G") { return "C"; }
		if($_[0] eq "T") { return "A"; }
	}
	# Create and write some graphics
	sub Create_Graph
	{
		our ($folderFigure, $filename, $maxValue) = @_;

		# Open the connection with R
		my $R = Statistics::R->new() or die "Impossible to create a communication bridge with R\n";
		$R->startR() ;
		# Load the Library
		$R->send(q`library(ggplot2)`);
		$R->send(q`library(gplots)`);
		$R->send(q`library(gtable)`);


		$R->send(qq`##########################################
								##		OVERALL MUTATION DISTRIBUTION 		##
								##########################################
								distrMut <- read.table("$folderFigure/Overall_mutation_distribution/$filename/$filename-OverallMutationDistribution.txt", header=T)`);
		$R->send(q`# Add the count of each category in the legend
							 distrMut$Legend[[1]] <- paste0(distrMut$Variant_type[[1]], " (", distrMut$Count[[1]], ")")
							 distrMut$Legend[[2]] <- paste0(distrMut$Variant_type[[2]], " (", distrMut$Count[[2]], ")")
							 distrMut$Legend[[3]] <- paste0(distrMut$Variant_type[[3]], " (", distrMut$Count[[3]], ")")`);

		$R->send(qq`# Base plot
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
								png("$folderFigure/Overall_mutation_distribution/$filename/$filename-OverallMutationDistribution.png", width=700, height=1100, res=300)
								print(pie)
								dev.off()


		         		##########################################
								##			SBS MUTATION DISTRIBUTION 			##
								##########################################
		         		distrSBS <- read.delim("$folderFigure/SBS_distribution/$filename/$filename-SBS_distribution.txt")
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
								options(bitmapType='cairo')
								png("$folderFigure/SBS_distribution/$filename/$filename-SBS_distribution.png", width=1800, height=1500, res=300)
								print(bar);
								dev.off()
								# Save the plot for the report
								bar
								ggsave("$folder_temp/$filename-SBS_distribution-Report.png")


								##########################################
								##					IMPACT ON PROTEIN 					##
								##########################################
								impactProt <- read.table("$folderFigure/Impact_protein_sequence/$filename/$filename-DistributionExoFunc.txt", header=T)
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
								options(bitmapType='cairo')
								png("$folderFigure/Impact_protein_sequence/$filename/$filename-DistributionExoFunc.png", width=1600, height=1800, res=300)
								print(pie)
								dev.off()
								# Save the plot for the report
								pie
								ggsave("$folder_temp/$filename-DistributionExoFunc-Report.png")


								##########################################
								##							STRAND BIAS 					  ##
								##########################################
								cb_palette_SB <- c("#0072B2", "#CC0000")
								file_sb       <- read.table("$folderFigure/Stranded_Analysis/$filename/$filename-StrandBias.txt", header=T);
								p_sb          <- ggplot(file_sb, aes(x=Alteration, y=Count, fill=Strand)) + theme_classic() + geom_bar(stat="identity", position="dodge") + scale_fill_manual(values=cb_palette_SB) + theme(axis.text.x = element_text(angle=60, hjust=1)) + xlab("") + theme(legend.position="bottom")
								# Save the plot for the HTML page (higher resolution)
								options(bitmapType='cairo')
								png("$folderFigure/Stranded_Analysis/$filename/$filename-StrandBias.png", width=1000, height=1200, res=300)
								print(p_sb)
								dev.off()
								# Save the plot for the report
								p_sb
								ggsave("$folder_temp/$filename-StrandBias-Report.png")


								##########################################
								##			HEATMAP SEQUENCE CONTEXT 				##
								##						GENOMIC STRAND 						##
								##########################################
								## COUNT
								heatmap_G <- read.table("$folderFigure/Trinucleotide_Sequence_Context/$filename/$filename-HeatmapCount-Genomic.txt", header=T)
								# Save the plot for the report
								options(bitmapType='cairo')
								png(filename="$folder_temp/$filename-HeatmapCount-Genomic-Report.png", bg="transparent", width=240, height=360)
								# Heatmap with an absolute scale
								heatmap.2(as.matrix(heatmap_G),Rowv=F,Colv=F,col=colorpanel(384,low="yellow",high="red"),dendrogram="none",scale="none",trace="none",key=F,labRow=rownames(as.matrix(heatmap_G)),labCol=colnames(as.matrix(heatmap_G)),lmat=rbind(c(5,1,4),c(3,1,2)), lhei=c(0.75,0.75),lwid=c(0.5,1.5,0.5))
								dev.off()
								# Save the plot for the HTML page (higher resolution)
								options(bitmapType='cairo')
								png(filename="$folderFigure/Trinucleotide_Sequence_Context/$filename/$filename-HeatmapCount-Genomic.png", width=1100, height=1600, res=300)
								heatmap.2(as.matrix(heatmap_G),Rowv=F,Colv=F,col=colorpanel(384,low="yellow",high="red"),dendrogram="none",scale="none",trace="none",key=F,labRow=rownames(as.matrix(heatmap_G)),labCol=colnames(as.matrix(heatmap_G)),lmat=rbind(c(5,1,4),c(3,1,2)), lhei=c(0.75,0.75),lwid=c(0.5,1.5,0.5))
								dev.off()

								## PERCENT
								heatmap_G <- read.table("$folderFigure/Trinucleotide_Sequence_Context/$filename/$filename-HeatmapPercent-Genomic.txt", header=T)
								# Save the plot for the report
								options(bitmapType='cairo')
								png(filename="$folder_temp/$filename-HeatmapPercent-Genomic-Report.png",bg="transparent", width=240, height=360)
								# Heatmap with an absolute scale
								heatmap.2(as.matrix(heatmap_G),Rowv=F,Colv=F,col=colorpanel(384,low="yellow",high="red"),dendrogram="none",scale="none",trace="none",key=F,labRow=rownames(as.matrix(heatmap_G)),labCol=colnames(as.matrix(heatmap_G)),lmat=rbind(c(5,1,4),c(3,1,2)), lhei=c(0.75,0.75),lwid=c(0.5,1.5,0.5))
								dev.off()
								# Save the plot for the HTML page (higher resolution)
								options(bitmapType='cairo')
								png(filename="$folderFigure/Trinucleotide_Sequence_Context/$filename/$filename-HeatmapPercent-Genomic.png", width=1100, height=1600, res=300)
								heatmap.2(as.matrix(heatmap_G),Rowv=F,Colv=F,col=colorpanel(384,low="yellow",high="red"),dendrogram="none",scale="none",trace="none",key=F,labRow=rownames(as.matrix(heatmap_G)),labCol=colnames(as.matrix(heatmap_G)),lmat=rbind(c(5,1,4),c(3,1,2)), lhei=c(0.75,0.75),lwid=c(0.5,1.5,0.5))
								dev.off()`);
		$R->stopR() ;

		## Plot the transcriptional strand bias in mutation signature
		`Rscript $pathRScriptTxnSB $folderFigure/Stranded_Analysis/$filename/$filename-StrandedSignatureCount.txt $folderFigure/Stranded_Analysis/$filename/$filename-StrandedSignatureCount $folder_temp/$filename-StrandedSignatureCount Count 2>&1`;
		`Rscript $pathRScriptTxnSB $folderFigure/Stranded_Analysis/$filename/$filename-StrandedSignaturePercent.txt $folderFigure/Stranded_Analysis/$filename/$filename-StrandedSignaturePercent $folder_temp/$filename-StrandedSignaturePercent Percent 2>&1`;
	}
	# Write the titles of the different sections of the report
	sub WriteBoderSection
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
			# Top-Left
			$ws->write($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+4, $colStart_SBSdistrBySeg, "Table 6. SBS distribution per chromosome", $format_topLeft); $ws->set_row($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+4, 18); # Set the height of the row to 0.25"
			# Top
			for(my $i=1; $i<8; $i++) { $ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+4, $colStart_SBSdistrBySeg+$i, $format_top); }
			# Top-Right
			$ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+4, $colStart_SBSdistrBySeg+8, $format_topRight);
			# Right
			$ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+5, $colStart_SBSdistrBySeg+8, $format_right); $ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+6, $colStart_SBSdistrBySeg+8, $format_right);
			# Bottom-Right
			if($refGenome =~ /hg/) { $ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+33, $colStart_SBSdistrBySeg+8, $format_bottomRight); }
			else { $ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+30, $colStart_SBSdistrBySeg+8, $format_bottomRight); }
			# Bottom
			if($refGenome =~ /hg/)
			{
				$ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+33, $colStart_SBSdistrBySeg+1, $format_bottom);
				for(my $i=3; $i<=7; $i++) { $ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+33, $colStart_SBSdistrBySeg+$i, $format_bottom); }
			}
			else
			{
				$ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+30, $colStart_SBSdistrBySeg+1, $format_bottom);
				for(my $i=3; $i<=7; $i++) { $ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+30, $colStart_SBSdistrBySeg+$i, $format_bottom); }
			}
			# Left
			$ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+5, $colStart_SBSdistrBySeg, $format_left); $ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+6, $colStart_SBSdistrBySeg, $format_left); $ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+7, $colStart_SBSdistrBySeg, $format_left);

			# Bottom-left
			if($refGenome =~ /hg/) { $ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+33, $colStart_SBSdistrBySeg, $format_bottomLeft);  }
			else { $ws->write_blank($rowStart_SBSdistrBySeg+8+$nb_func+(($nb_func+4)*2)+30, $colStart_SBSdistrBySeg, $format_bottomLeft); }
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
	# Create logo for representing the sequence context with n bases
	sub CreateLogo
	{
		my ($refH_file, $folderWebLogo) = @_;

		my $folderSample = "";

		foreach my $k_file (sort keys $refH_file)
		{
			$folderSample      = "$folderWebLogo/$k_file";
			if(!-e $folderSample) { mkdir($folderSample) or die "Can't create the directory $folderSample\n"; }

			my $test_lengthSeqContext = 0;

			foreach my $k_mutation (sort keys $refH_file->{$k_file}{'WebLogo3'})
			{
				open(WEBLOGO, ">", "$folderSample/$k_file-$k_mutation.fa") or die "$!: $folderSample/$k_file-$k_mutation.fa\n";
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

				## Test if there is fasta sequence for the mutation type
				my $nbLigne_temp = `wc -l $fastaFile`;
				my @nbLigne = split(" ", $nbLigne_temp);

				if($nbLigne[0] == 0) { print "WARNING: No sequence for $filename\n"; next; }

				# When length sequence context is lower than 10 the image is to small for adding a title
				if($test_lengthSeqContext == 1) { system("weblogo -c classic -F png -U probability --title $title < $fastaFile > $folderSample/$filename-Probability.png"); }
				else { system("weblogo -c classic -F png -U probability < $fastaFile > $folderSample/$filename-Probability.png"); }
			}
		}
	}


	# Define the format of the worksheet: Arial font size=10
	sub Format_A10
	{
		my ($wb, $format) = @_;
		$$format = $wb->add_format(font=>'Arial', size=>10); $$format->set_align('center');
	}
	# Define the format of the worksheet: Arial font size=11 bold and center
	sub Format_A11Bold
	{
		my ($wb, $format) = @_;
		$$format = $wb->add_format(font=>'Arial', size=>11, bold=>1); $$format->set_align('center');
	}
	# Define the format of the worksheet: Arial font size=10 italic and center
	sub Format_A10Italic
	{
		my ($wb, $format) = @_;
		$$format = $wb->add_format(font=>'Arial', size=>10, italic=>1); $$format->set_align('center');
	}
	# Defile the format of the worksheet: Arialt font size=11 bold and left
	sub Format_A11BoldLeft
	{
		my ($wb, $format) = @_;
		$$format = $wb->add_format(valign =>'left', font=>'Arial', size=>11, bold=>1);
	}
	# Defile the format of the worksheet: Arialt font size=10 bold and left
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
	# Define the mutation type header for the Strand bias by segment
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
	# Define the mutation type header for the trinucleotide sequence context on the coding strand
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
		my ($wb, $table_topleft, $table_topRight, $table_bottomleft, $table_bottomRight, $table_top, $table_right, $table_bottom, $table_bottomItalic, $table_left, $table_bottomrightHeader, $table_left2, $table_middleHeader, $table_middleHeader2) = @_;

		$$table_topleft = $wb->add_format(valign=>'center', bold=>1, font=>'Arial', size=>10); $$table_topleft->set_top(1); $$table_topleft->set_left(1);
		$$table_topRight = $wb->add_format(valign=>'center', bold=>1, font=>'Arial', size=>10); $$table_topRight->set_top(1); $$table_topRight->set_right(1);
		$$table_bottomleft = $wb->add_format(valign=>'center', bold=>1, font=>'Arial', size=>10); $$table_bottomleft->set_bottom(1); $$table_bottomleft->set_left(1);
		$$table_bottomRight = $wb->add_format(valign=>'center', font=>'Arial', size=>10); $$table_bottomRight->set_bottom(1); $$table_bottomRight->set_right(1);

		$$table_top          = $wb->add_format(valign=>'center', bold=>1, font=>'Arial', size=>10);   $$table_top->set_top(1);
		$$table_right        = $wb->add_format(valign=>'center', font=>'Arial', size=>10);            $$table_right->set_right(1);
		$$table_bottom       = $wb->add_format(valign=>'center', font=>'Arial', size=>10);            $$table_bottom->set_bottom(1);
		$$table_bottomItalic = $wb->add_format(valign=>'center', font=>'Arial', size=>10, italic=>1); $$table_bottomItalic->set_bottom(1);
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
  if($name_of_column_NB eq "toto") { print STDERR "Error recoverNumCol(): the column named $name_of_column doesn't exits in the input file $input!!!!!\n"; exit; }
  else                             { return $name_of_column_NB; }
}




=head1 NAME

mutSpec-Stat

=head1 SYNOPSIS

	mutSpecstat.pl [arguments] <query-file>

  <query-file>                                   can be a folder with multiple VCF or a single VCF

  Arguments:
        -h,        --help                        print help message
        -m,        --man                         print complete documentation
        -v,        --verbose                     use verbose output
                   --refGenome                   the reference genome to use (hg19 or mm9)
        -o,        --outfile <string>            output directory for the result. If none is specify the result will be write in the same directory as the input file
        -temp      --pathTemporary <string>      the path for saving the temporary files
                   --pathSeqRefGenome            the path to the fasta reference sequences
                   --poolData                    generate the pool of all the samples (optional)
                   --reportSample                generate a report for each sample (optional)


Function: automatically run a pipeline and calculate various statistics on mutations

 Example: mutSpecstat.pl --refGenome hg19 --outfile output_directory --temp path_to_temporary_directory --pathRscript path_to_R_scripts --pathSeqRefGenome path_fasta_ref_seq --poolData --reportSample input

 Version: 10-2015 (Oct 2015)


=head1 OPTIONS

=over 8

=item B<--help>

print a brief usage message and detailed explanation of options.

=item B<--man>

print the complete manual of the program.

=item B<--verbose>

use verbose output.

=item B<--refGenome>

the reference genome to use, could be hg19 or mm9.

=item B<--outfile>

the directory of output file names. If it is nor specify the same directory as the input file is used.

=item B<--pathTemporary>

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
