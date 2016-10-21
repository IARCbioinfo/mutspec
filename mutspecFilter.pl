# !/usr/bin/perl

#-----------------------------------#
# Author: Maude / Vincent           #
# Script: mutspecFilter.pl          #
# Last update: 21/10/16             #
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
our ($dbSNP_value, $segDup, $esp, $thG, $exac) = (0, 0, 0, 0, 0); # For filtering agains the databases dbSNP, genomic duplicate segments, Exome Sequencing Project and 1000 genome, Exac.
our ($output, $refGenome)               = ("", "");     # The path for saving the result; The reference genome to use.
our ($listAVDB)                         = "empty";      # Text file with the list Annovar databases.
our ($dir)                              = "";           # Path to BED file to filter against
our (@filters);

GetOptions('dir|d=s'=>\$dir,'verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'dbSNP=i'=>\$dbSNP_value, 'segDup'=>\$segDup, 'esp'=>\$esp, 'thG'=>\$thG, 'exac'=>\$exac, 'outfile|o=s' => \$output, 'refGenome=s'=>\$refGenome, 'pathAVDBList=s' => \$listAVDB, 'filter=s'=> \@filters) or pod2usage(2);

our ($input) = @ARGV;

pod2usage(-verbose=>1, -exitval=>1, -output=>\*STDERR) if ($help);
pod2usage(-verbose=>2, -exitval=>1, -output=>\*STDERR) if ($man);
pod2usage(-verbose=>0, -exitval=>1, -output=>\*STDERR) if(@ARGV == 0); # No argument is pass to the command line print the usage of the script
pod2usage(-verbose=>0, -exitval=>1, -output=>\*STDERR) if(@ARGV == 2); # Only one argument is expected to be pass to @ARGV (the input)



# If the dbSNP value is not equal to zero filter using the dbSNP column specify
our $dbSNP = 0;
if($dbSNP_value > 0) { $dbSNP = 1; }


############ Check flags ############
if($listAVDB eq "empty") { $listAVDB = "$dir/${refGenome}_listAVDB.txt" }

# Zero databases is specified
if( ($dbSNP == 0) && ($segDup == 0) && ($esp == 0) && ($thG == 0) && ($exac == 0) && (scalar(@filters) == 0) )
{
	print STDERR "There is no databases selected for filtering against!!!\n";
	print STDERR "Please chose at least one between dbSNP, SegDup, ESP (only for human genome), 1000 genome (only for human genome) or ExAC (only for human genome)\n";
	print STDERR "Or specify a BED file\n";
	exit;
}



############ Recover the name of the databases to filter against ############
my ($segDup_name, $espAll_name, $thousandGenome_name, $exac_name) = ("", "", "", "");
my @tab_protocol = ();

if( ($segDup == 1) || ($esp == 1) || ($thG == 1)  || ($exac == 1))
{
	### Recover the name of the column
	my $protocol = "";
	ExtractAVDBName($listAVDB, \$protocol);
	@tab_protocol = split(",", $protocol);

	for(my $i=0; $i<=$#tab_protocol; $i++)
	{
		if($tab_protocol[$i] =~ /genomicSuperDups/) { $segDup_name = $tab_protocol[$i]; }
		elsif($tab_protocol[$i] =~ /1000g/)         { $thousandGenome_name = $tab_protocol[$i]; }
		elsif($tab_protocol[$i] =~ /esp/)           { $espAll_name = $tab_protocol[$i]; }
		elsif($tab_protocol[$i] =~ /exac/i)          { $exac_name = $tab_protocol[$i]; }
	}
}


############ Filter the file ############
filterAgainstPublicDB();


print STDOUT "\tFilter selected\tdbSNP = ".$dbSNP."\tsegDup = ".$segDup."\tesp = ".$esp."\tthG = ".$thG."\tEXac = ". $exac . "\n";


### Write a message if the input file contains zero variants or if all the variants are filtered out
my ($filename, $directories, $suffix) = fileparse($input, qr/\.[^.]*/);
my $nbVariantsIn = `wc -l $input`;
$nbVariantsIn =~ /(\d+).+/;
my $nbLineIn  = $1;
if($nbLineIn == 1)
{
	print STDOUT "\nThere is no variant to be filtered for $filename\n";
	print STDOUT "Check MutSpecAnnot tool standard output for more informations\n";
}
else
{
	my ($filenameOut, $directoriesOut, $suffixOut) = fileparse($output, qr/\.[^.]*/);
	my $nbVariantsOut = `wc -l $output`;
	$nbVariantsOut =~ /(\d+).+/;
	my $nbLineOut  = $1;
	if($nbLineOut == 1)
	{
		print STDOUT "\nAll the variants were filtered for $filenameOut\n";
	}
}

### filter versus additional VCF files if provided.
if ( scalar(@filters) > 0) { filterAdditionalBED(); }




sub filterAgainstPublicDB
{
	open(FILTER, ">", "$output") or die "$!: $output\n";

	open(F1, $input) or die "$!: $input\n";
	my $header = <F1>; print FILTER $header;
	while(<F1>)
	{
		$_      =~ s/[\r\n]+$//;
		my @tab = split("\t", $_);

		my ($segDupInfo, $espAllInfo, $thgInfo, $exacInfo) = (0, 0 ,0, 0);

		if($segDup == 1)
		{
			my $segDup_value = recoverNumCol($input, $segDup_name);
			$segDupInfo      = formatSegDupInfo($tab[$segDup_value]);
			# Replace NA by 0 for making test on the same type of variable
			$segDupInfo =~ s/NA/0/;
		}
		if($esp == 1)
		{
			my $espAll_value = recoverNumCol($input, $espAll_name);
			$espAllInfo      = $tab[$espAll_value];
			# Replace NA by 0 for making test on the same type of variable
			$espAllInfo      =~ s/NA/0/;
		}
		if($thG == 1)
		{
			my $thousandGenome_value = recoverNumCol($input, $thousandGenome_name);
			# Replace NA by 0 for making test on the same type of variable
			$thgInfo = $tab[$thousandGenome_value];
			$thgInfo =~ s/NA/0/;
		}
		if($exac == 1)
		{
			my $exac_value = recoverNumCol($input, $exac_name);
			# Replace NA by 0 for making test on the same type of variable
			$exacInfo = $tab[$exac_value];
			$exacInfo =~ s/NA/0/;
		}

		my $filter = 0;
		if( $dbSNP  == 1 && $tab[$dbSNP_value-1] ne "NA" ){ $filter = 1; }
		if( $segDup == 1 && $segDupInfo >= 0.9)  		  { $filter = 1; }
		if( $esp    == 1 && $espAllInfo > 0.001)  		  { $filter = 1; }
		if( $thG 	== 1 && $thgInfo > 0.001)  		  	  { $filter = 1; }
		if( $thG 	== 1 && $exacInfo > 0.001)  		  { $filter = 1; }

		if (!$filter) { print FILTER "$_\n"; }

	}
	close F1; close FILTER;
}


sub formatSegDupInfo
{
	my ($segDup_info) = @_;

	if($segDup_info ne "NA") # Score=0.907883;Name=chr9:36302931
	{
		my @segDup = split(";", $segDup_info);
		$segDup[0] =~ /Score=(.+)/;
		return $1;
	}
	else { return $segDup_info; }
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
		if( ($tab[0] =~ /\w+_(\w+)\.txt/) && ($tab[0] !~ /sites/) && ($tab[0] !~ /esp/) && ($tab[0] !~ /sift/) && ($tab[0] !~ /pp2/) && ($tab[0] !~ /exac/i) )
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
		# EXAC
		if($tab[0] =~ /exac/i)
		{
			$tab[0] =~ /\w+_(\w+)_(\w+)\.txt/;
			my $AVdbName_final = "ExAC_ALL";

			$$refS_protocol .= $AVdbName_final.",";
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
		open(FT,$input) or die "$!: $input\n";
	  # For having the name of the columns
	  my $search_header = <FT>; $search_header =~ s/[\r\n]+$//; my @tab_search_header = split("\t",$search_header);
	  close FT;
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


sub filterAdditionalBED{

	#create bed
	open(TABLE, "$output") or die "$!: $output\n";
	open(F1, ">bed") or die "cannot create bed file";
	my $NL=1;
	my $headers=<TABLE>;
	while(<TABLE>)
	{
        $NL++;
        my @line=split("\t", $_);
        print F1 "$line[0]\t$line[1]\t$line[2]\t$NL\n";
	}
	close F1;
	close TABLE;
	#and sort it
	`sort -k1,1 -k2,2n bed > sorted`;

	foreach my $filter (@filters){

		my ($filename, $directories, $suffix) = fileparse($filter, qr/\.[^.]*/);

		print "\tFilter against BED: $filename\n";

		#find intersect
		`sort -k1,1 -k2,2n $filter > ref`;
		`bedtools intersect -a sorted -b ref -v -sorted > bed`;
		`sort -k1,1 -k2,2n bed > sorted`;
	}

	#generate new output
	`sort -k4n sorted > bed`;
	`cp $output table`;

	open(F1, "bed") or die "error no sorted file";
	open(F2, "table") or die "error no table file";
	open(OUT, ">$output") or die "error cannot open output file";
	print OUT $headers;
	$NL=1;
	my $line = <F2>;
	while(<F1>)
	{
        my @NR=split("\t", $_);
        while( $NL < $NR[3]){ $line = <F2>; $NL++; }
        print OUT $line;
	}
	close F1;
	close F2;
	close OUT;

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
									 --dbSNP <value>               filter against dbSNP database. Specify the number of the dbSNP column in the file (start to count from 1)
									 --segDup                      filter against genomic duplicate database
									 --esp                         filter against Exome Sequencing Project database (only for human)
									 --thG                         filter against 1000 genome database (onyl for human)
			  -o,        --outfile <string>            name of output file
			             --refGenome                   reference genome to use
			             --pathAVDBList                path to the list of Annovar databases installed
			             --filter                      path to a bed file


Function: Filter out variants present in public databases

 Example: # Filter against dbSNP
 					mutspecFilter.pl --dbSNP col_number (start to count from 1) --refGenome hg19 --pathAVDBList path_to_the_list_of_annovar_DB --outfile output_filename input

 					# Filter against all Annovar databases
 					mutspecFilter.pl --dbSNP col_number (start to count from 1) --segDup --esp --thG --exac --refGenome hg19 --pathAVDBList path_to_the_list_of_annovar_DB --outfile output_filename input

 					# Filter against additional databases in BED format
 					mutspecFilter.pl --filter path_to_bed --refGenome hg19 --pathAVDBList path_to_the_list_of_annovar_DB --outfile output_filename input


 Version: 10-2016 (October 2016)


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
Specify the number of the dbSNP column in the file (start to count from 1)
For human and mouse genome

=item B<--segDup>

Remove all the variants with a frequency greater or equal to 0.9 in genomic duplicate segments database
For human and mouse genome

=item B<--esp>

Remove all the variants with a frequency greater than 0.001 in Exome sequencing project
For human genome only

=item B<--thG>

Remove all the variants with a frequency greater than 0.001 in 1000 genome database


=item B<--exac>

Remove all the variants with a frequency greater than 0.001 in ExAC database


=item B<--filter>

Remove all variants present in the BED file


=item B<--refGenome>

The reference genome to use.

=item B<--outfile>

the name of the output file

=item B<--pathAVDBList>

the path to a texte file containing the list of the Annovar databases installed.

=back

=head1 DESCRIPTION

mutspecFilter - Filter a file annotated with MutSpec-Annot tool.
Variants present in public databases (dbSNP, SegDup, ESP, 1000 genome, exac obtained from Annovar) will be removed from the input file (with frequency limits described above).
Additionally, using the --filter option, any variants present in a specified bed file will be removed from the input file.

=cut
