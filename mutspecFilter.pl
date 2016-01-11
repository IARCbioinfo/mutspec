# !/usr/bin/perl

#-----------------------------------#
# Author: Maude                     #
# Script: mutspecFilter.pl          #
# Last update: 26/08/15             #
#-----------------------------------#

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Basename; # my ($filename, $directories, $suffix) = fileparse($file, qr/\.[^.]*/);
use File::Path;

################################################################################################################################################################################
#																																		Filter an Annotaed file with Annovar																		  																 #
################################################################################################################################################################################

our ($verbose, $man, $help)             = (0, 0, 0);    # Parse options and print usage if there is a syntax error, or if usage was explicitly requested.
our ($dbSNP_value, $segDup, $esp, $thG) = (0, 0, 0, 0); # For filtering agains the databases dbSNP, genomic duplicate segments, Exome Sequencing Project and 1000 genome.
our ($output, $refGenome)               = ("", "");     # The path for saving the result; The reference genome to use (hg19 or mm9).
our ($listAVDB)                         = "empty";      # Text file with the list Annovar databases.
our ($dir)                              = "";

GetOptions('dir|d=s'=>\$dir,'verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'dbSNP=i'=>\$dbSNP_value, 'segDup'=>\$segDup, 'esp'=>\$esp, 'thG'=>\$thG, 'outfile|o=s' => \$output, 'refGenome=s'=>\$refGenome, 'pathAVDBList=s' => \$listAVDB) or pod2usage(2);

our ($input) = @ARGV;

pod2usage(-verbose=>1, -exitval=>1, -output=>\*STDERR) if ($help);
pod2usage(-verbose=>2, -exitval=>1, -output=>\*STDERR) if ($man);
pod2usage(-verbose=>0, -exitval=>1, -output=>\*STDERR) if(@ARGV == 0); # No argument is pass to the command line print the usage of the script
pod2usage(-verbose=>0, -exitval=>1, -output=>\*STDERR) if(@ARGV == 2); # Only one argument is expected to be pass to @ARGV (the input)


############ Check flags ############
if($listAVDB eq "empty") { $listAVDB = "$dir/${refGenome}_listAVDB.txt" }
# If the dbSNP value is not equal to zero filter using dbSNP DB
our $dbSNP = 0;
if($dbSNP_value > 0) { $dbSNP = 1; }


############ Annovar databases ############
my $protocol = "";
ExtractAVDBName($listAVDB, \$protocol);
my @tab_protocol = split(",", $protocol);
############ Annovar databases ############


############ Process Argument ############
my ($filename, $directories, $suffix) = fileparse($input, qr/\.[^.]*/);

if(scalar(@tab_protocol) == 3)
{
	my ($segDup_name, $espAll_name, $thousandGenome_name) = ("", "", "");

	for(my $i=0; $i<=$#tab_protocol; $i++)
	{
		if($tab_protocol[$i] =~ /genomicSuperDups/) { $segDup_name = $tab_protocol[$i]; }
		elsif($tab_protocol[$i] =~ /1000g/)         { $thousandGenome_name = $tab_protocol[$i]; }
		elsif($tab_protocol[$i] =~ /esp/)           { $espAll_name = $tab_protocol[$i]; }
	}

	FilterAV_human($input, $filename, $segDup_name, $espAll_name, $thousandGenome_name);
}
elsif(scalar(@tab_protocol) == 1)
{
	my $segDup_name = $tab_protocol[0];

	FilterAV_mouse($input, $filename, $segDup_name);
}
else
{
	print STDERR "There is more databases specify than the possibility!!!\nAvailable databases for filtering are dbSNP, segDup (for human and mouse genome) and esp, 1000g (for human genomes only)!!!\nYour have specify @tab_protocol databases\n";
}
############ Process Argument ############

# Filter an annotated file using the databases dbSNP, SegDup, ESP and 1000genome
sub FilterAV_human
{
	my ($inputFile, $filename, $segDup_name, $espAll_name, $thousandGenome_name) = @_;

	my $segDup_value         = recoverNumCol($inputFile, $segDup_name);
	my $espAll_value         = recoverNumCol($inputFile, $espAll_name);
	my $thousandGenome_value = recoverNumCol($inputFile, $thousandGenome_name);

	$filename =~ s/\./TOTO/; my @temp = split("TOTO", $filename); $filename =~ s/TOTO/\./;

	open(FILTER, ">", "$output") or die "$!: $output\n";
	open(F1, $inputFile) or die "$!: $inputFile\n";
	my $header = <F1>; print FILTER $header;

	while(<F1>)
	{
		if ($_ =~ /^#/) { next; }

		$_      =~ s/[\r\n]+$//;
		my @tab = split("\t", $_);

		my $segDupInfo = "";
		if($tab[$segDup_value] ne "NA") # Score=0.907883;Name=chr9:36302931
		{
			my @segDup = split(";", $tab[$segDup_value]);
			$segDup[0] =~ /Score=(.+)/;
			$segDupInfo = $1;
		}
		else { $segDupInfo = $tab[$segDup_value]; }

		# Replace NA by 0 for making test on the same type of variable
		$segDupInfo =~ s/NA/0/; $tab[$espAll_value] =~ s/NA/0/; $tab[$thousandGenome_value] =~ s/NA/0/;

		##############################
		#   			One Filter 				 #
		##############################
		# Remove all the variants present in dbSNP
		if( ($dbSNP == 1) && ($segDup==0) && ($esp==0) && ($thG==0) ) { if($tab[$dbSNP_value-1] eq "NA")           { print FILTER "$_\n"; } }
		# Remove all the variants with a frequency greater than or equal to 0.9  in genomic duplicate segments database
		if( ($dbSNP==0) && ($segDup == 1) && ($esp==0) && ($thG==0) ) { if($segDupInfo < 0.9)                    { print FILTER "$_\n"; } }
		# Remove all the variants with greater than 0.001 in Exome sequencing project
		if( ($dbSNP==0) && ($segDup==0) && ($esp == 1) && ($thG==0) )    { if($tab[$espAll_value] <= 0.001)         { print FILTER "$_\n"; } }
		# Remove all the variants with greater than 0.001 in 1000 genome database
		if( ($dbSNP==0) && ($segDup==0) && ($esp==0) && ($thG == 1) )    { if($tab[$thousandGenome_value] <= 0.001) { print FILTER "$_\n"; } }

		##############################
		#   			Two Filter 				 #
		##############################
		if( ($dbSNP==1) && ($segDup==1) && ($esp==0) && ($thG== 0) ) { if( ($tab[$dbSNP_value-1] eq "NA") && ($segDupInfo < 0.9) )                    { print FILTER "$_\n"; } }
		if( ($dbSNP==1) && ($segDup==0) && ($esp==1) && ($thG==0) )  { if( ($tab[$dbSNP_value-1] eq "NA") && ($tab[$espAll_value] <= 0.001) )         { print FILTER "$_\n"; } }
		if( ($dbSNP==1) && ($segDup==0) && ($esp==0) && ($thG==1) )  { if( ($tab[$dbSNP_value-1] eq "NA") && ($tab[$thousandGenome_value] <= 0.001) ) { print FILTER "$_\n"; } }

		if( ($dbSNP==0) && ($segDup==1) && ($esp==1) && ($thG==0) )   { if( ($segDupInfo < 0.9) && ($tab[$espAll_value] <= 0.001) )                 { print FILTER "$_\n"; } }
		if( ($dbSNP==0) && ($segDup==1) && ($esp==0) && ($thG==1) ) { if( ($segDupInfo < 0.9) && ($tab[$thousandGenome_value] <= 0.001) )         { print FILTER "$_\n"; } }

		if( ($dbSNP==0) && ($segDup==0) && ($esp==1) && ($thG==1) )   { if( ($tab[$espAll_value] <= 0.001) && ($tab[$thousandGenome_value] <= 0.001) ) { print FILTER "$_\n"; } }

		##############################
		#   		Three Filter 				 #
		##############################
		if( ($dbSNP==1) && ($segDup==1) && ($esp==1) && ($thG==0) ) { if( ($tab[$dbSNP_value-1] eq "NA") && ($segDupInfo < 0.9) && ($tab[$espAll_value] <= 0.001) )
		{ print FILTER "$_\n"; } }
		if( ($dbSNP==1) && ($segDup==1) && ($esp==0) && ($thG==1) ) { if( ($tab[$dbSNP_value-1] eq "NA") && ($segDupInfo < 0.9) && ($tab[$thousandGenome_value] <= 0.001) )
		{ print FILTER "$_\n"; } }
		if( ($dbSNP==1) && ($segDup==0) && ($esp==1) && ($thG==1) ) { if( ($tab[$dbSNP_value-1] eq "NA") && ($tab[$espAll_value] <= 0.001) && ($tab[$thousandGenome_value] <= 0.001) )
		{ print FILTER "$_\n"; } }
		if( ($dbSNP==0) && ($segDup==1) && ($esp==1) && ($thG==1) ) { if( ($segDupInfo < 0.9) && ($tab[$espAll_value] <= 0.001) && ($tab[$thousandGenome_value] <= 0.001) )
		{ print FILTER "$_\n"; } }

		##############################
		#   		FOUR Filter 				 #
		##############################
		if( ($dbSNP==1) && ($segDup==1) && ($esp==1) && ($thG==1) ) { if( ($tab[$dbSNP_value-1] eq "NA") && ($segDupInfo < 0.9) && ($tab[$espAll_value] <= 0.001) && ($tab[$thousandGenome_value] <= 0.001) )
		{ print FILTER "$_\n"; } }
	}
	close F1; close FILTER;
}
# Filter an annotated file using the databases dbSNP and segDup
sub FilterAV_mouse
{
	my ($inputFile, $filename, $segDup_name) = @_;

	my $segDup_value         = recoverNumCol($inputFile, $segDup_name);

	$filename =~ s/\./TOTO/; my @temp = split("TOTO", $filename); $filename =~ s/TOTO/\./;

	open(FILTER, ">", "$output") or die "$!: $output\n";
	open(F1, $inputFile) or die "$!: $inputFile\n";
	my $header = <F1>; print FILTER $header;

	while(<F1>)
	{
		$_      =~ s/[\r\n]+$//;
		my @tab = split("\t", $_);

		my $segDupInfo = "";
		if($tab[$segDup_value] ne "NA") # Score=0.907883;Name=chr9:36302931
		{
			my @segDup = split(";", $tab[$segDup_value]);
			$segDup[0] =~ /Score=(.+)/;
			$segDupInfo = $1;
		}
		else { $segDupInfo = $tab[$segDup_value]; }

		# Replace NA by 0 for making test on the same type of variable
		$segDupInfo =~ s/NA/0/;

		##############################
		#   			One Filter 				 #
		##############################
		# Remove all the variants present in dbSNP
		if( ($dbSNP == 1) && ($segDup==0) ) { if($tab[$dbSNP_value-1] eq "NA") { print FILTER "$_\n"; } }
		# Remove all the variants with a frequency greater than or equal to 0.9  in genomic duplicate segments database
		if( ($dbSNP==0) && ($segDup == 1) ) { if($segDupInfo < 0.9)          { print FILTER "$_\n"; } }

		##############################
		#   			Two Filter 				 #
		##############################
		if( ($dbSNP==1) && ($segDup==1) ) { if( ($tab[$dbSNP_value-1] eq "NA") && ($segDupInfo < 0.9) ) { print FILTER "$_\n"; } }
	}
	close F1; close FILTER;
}


sub ExtractAVDBName
{
	my ($listAVDB, $refS_protocol) = @_;

	open(F1, $listAVDB) or die "$!: $listAVDB\n";
	while(<F1>)
	{
		if ($_ =~ /^#/) { next; }

		$_      =~ s/[\r\n]+$//;
		my @tab = split("\t", $_);

		# db name like refGenome_dbName.txt
		if( ($tab[0] =~ /\w+_(\w+)\.txt/) && ($tab[0] !~ /sites/) && ($tab[0] !~ /esp/) && ($tab[0] !~ /sift/) && ($tab[0] !~ /pp2/) )
		{
			my $temp = $1;
			if($temp =~ /genomicSuperDups/) { $$refS_protocol .= $temp.","; }
		}
		# 1000 genome
		if($tab[0] =~ /sites/)
		{
			$tab[0] =~ /\w+_(\w+)\.sites.(\d+)_(\d+)\.txt/;
			my ($dbName, $year, $month) = ($1, $2, $3);
			$dbName =~ tr/A-Z/a-z/;

			# convert the month number into the month name
			ConvertMonth(\$month);

			my $AVdbName_final = "1000g".$year.$month."_".$dbName;

			if($dbName eq "all") { $$refS_protocol .=$AVdbName_final.","; }
		}
		# ESP
		if($tab[0] =~ /esp/)
		{
			$tab[0] =~ /\w+_(\w+)_(\w+)\.txt/;
			my $AVdbName_final = $1."_".$2;

			if($2 eq "all") { $$refS_protocol .=$AVdbName_final.","; }
		}
	}
	close F1;

	sub ConvertMonth
	{
		my ($refS_month) = @_;

		if($$refS_month == 1)  { $$refS_month = "janv"; }
		elsif($$refS_month == 2)  { $$refS_month = "feb"; }
		elsif($$refS_month == 3)  { $$refS_month = "mar"; }
		elsif($$refS_month == 4)  { $$refS_month = "apr"; }
		elsif($$refS_month == 5)  { $$refS_month = "may"; }
		elsif($$refS_month == 6)  { $$refS_month = "jun"; }
		elsif($$refS_month == 7)  { $$refS_month = "jul"; }
		elsif($$refS_month == 8)  { $$refS_month = "aug"; }
		elsif($$refS_month == 9)  { $$refS_month = "sept"; }
		elsif($$refS_month == 10) { $$refS_month = "oct"; }
		elsif($$refS_month == 11) { $$refS_month = "nov"; }
		elsif($$refS_month == 12) { $$refS_month = "dec"; }
		else { print STDERR "Month number don't considered\n"; exit; }
	}
}


sub recoverNumCol
{
	my ($input, $name_of_column) = @_;

	# With Annovar updates the databases name changed and are present in an array
	if( ref($name_of_column) eq "ARRAY" )
	{
		my $test = "";
		my @tab = @$name_of_column;
		foreach (@tab)
		{
			open(F1,$input) or die "$!: $input\n";
		  # For having the name of the columns
		  my $search_header = <F1>; $search_header =~ s/[\r\n]+$//; my @tab_search_header = split("\t",$search_header);
		  close F1;
		  # The number of the column
		  my $name_of_column_NB  = "toto";
		  for(my $i=0; $i<=$#tab_search_header; $i++)
		  {
		    if($tab_search_header[$i] eq $_) { $name_of_column_NB = $i; }
		  }
		  if($name_of_column_NB eq "toto") { next; }
		  else                             { return $name_of_column_NB; }
		}
		if($name_of_column eq "toto") { print "Error recoverNumCol: the column named $name_of_column doesn't exits in the input file $input!!!!!\n"; exit; }
	}
	# Only one name is pass
	else
	{
		open(F1,$input) or die "$!: $input\n";
	  # For having the name of the columns
	  my $search_header = <F1>; $search_header =~ s/[\r\n]+$//; my @tab_search_header = split("\t",$search_header);
	  close F1;
	  # The number of the column
	  my $name_of_column_NB  = "toto";
	  for(my $i=0; $i<=$#tab_search_header; $i++)
	  {
	    if($tab_search_header[$i] eq $name_of_column) { $name_of_column_NB = $i; }
	  }
	  if($name_of_column_NB eq "toto") { print "Error recoverNumCol: the column named $name_of_column doesn't exits in the input file $input!!!!!\n"; exit; }
	  else                        { return $name_of_column_NB; }
	}
}

=head1 NAME

mutspecFilter - Filter a file annotated with MutSpec-Annot tool. Variants present in public databases (dbSNP, SegDup, ESP, 1000 genome obtained from Annovar) will be removed from the input file (with frequency limits described above)

=head1 SYNOPSIS

	mutspecFilter.pl [arguments] <query-file>

  <query-file>                                   an annotated file

  Arguments:
        -h,        --help                        print help message
        -m,        --man                         print complete documentation
        -v,        --verbose                     use verbose output
									 --dbSNP <value>               filter against dbSNP database. Specify the number of the dbSNP column in the file
									 --segDup                      filter against genomic duplicate database
									 --esp                         filter against Exome Sequencing Project database (only for human)
									 --thG                         filter against 1000 genome database (onyl for human)
			  -o,        --outfile <string>            name of output file
			             --refGenome                   reference genome to use
			             --pathAVDBList                path to the list of Annovar databases installed


Function: Filter out variants present in public databases

 Example: # Filter against dbSNP
 					mutspecFilter.pl --dbSNP col_number --refGenome hg19 --pathAVDBList path_to_the_list_of_annovar_DB --outfile output_filename input

 					# Filter against the four databases
 					mutspecFilter.pl --dbSNP col_number --segDup --esp --thG --refGenome hg19 --pathAVDBList path_to_the_list_of_annovar_DB --outfile output_filename input


 Version: 08-2015 (Aug 2015)


=head1 OPTIONS

=over 8

=item B<--help>

print a brief usage message and detailed explanation of options.

=item B<--man>

print the complete manual of the program.

=item B<--verbose>

use verbose output.

=item B<--dbSNP>

Remove all the variants presents in the dbSNP databases
Specify the number of column containing the annotation
For human and mouse genome

=item B<--segDup>

Remove all the variants with a frequency greater or equal to 0.9 in genomic duplicate segments database
For human and mouse genome

=item B<--esp>

Remove all the variants with a frequency greater than 0.001 in Exome sequencing project
For human genome only

=item B<--thG>

Remove all the variants with a frequency greater than 0.001 in 1000 genome database

=item B<--refGenome>

the reference genome to use, could be hg19 or mm9.

=item B<--outfile>

the name of the output file

=item B<--pathAVDBList>

the path to a texte file containing the list of the Annovar databases installed.

=back

=head1 DESCRIPTION

mutspecFilter - Filter a file annotated with MutSpec-Annot tool. Variants present in public databases (dbSNP, SegDup, ESP, 1000 genome obtained from Annovar) will be removed from the input file (with frequency limits described above)

=cut
