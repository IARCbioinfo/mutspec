#!/usr/bin/env perl

#-----------------------------------#
# Author: Maude                     #
# Script: mutspecAnnot.pl           #
# Last update: 28/10/16             #
#-----------------------------------#

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Basename; # my ($filename, $directories, $suffix) = fileparse($file, qr/\.[^.]*/);
use File::Path;
use Parallel::ForkManager;


our ($verbose, $man, $help) = (0, 0, 0); # Parse options and print usage if there is a syntax error, or if usage was explicitly requested.
our ($refGenome, $output, $path_AVDB, $pathAVDBList, $folder_temp) = ("empty", "empty", "empty", "empty", "empty");  # The reference genome to use; The path for saving the result; The path to Annovar database; Text file with the list of the databases for Annovar; the path for saving the temporary files
our ($intervalEnd)          = (10); # Number of bases for the flanking region for the sequence context.
our ($fullAVDB)             = "yes"; # Add an option for using all Annovar databases for the annotation or only refGene + strand + context for having a quicker annotation (for large file with million of lines)

GetOptions('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'refGenome=s'=>\$refGenome, 'interval=i' => \$intervalEnd, 'fullAnnotation=s' => \$fullAVDB, 'outfile|o=s' => \$output, 'pathAnnovarDB|AVDB=s' => \$path_AVDB, 'pathAVDBList=s' => \$pathAVDBList, 'pathTemporary|temp=s' => \$folder_temp) or pod2usage(2);

our ($input) = @ARGV;

pod2usage(-verbose=>1, -exitval=>1, -output=>\*STDERR) if ($help);
pod2usage(-verbose=>2, -exitval=>1, -output=>\*STDERR) if ($man);
pod2usage(-verbose=>0, -exitval=>1, -output=>\*STDERR) if(@ARGV == 0); # No argument is pass to the command line print the usage of the script
pod2usage(-verbose=>0, -exitval=>1, -output=>\*STDERR) if(@ARGV == 2); # Only one argument is expected to be pass to @ARGV (the input)



######################################################################################################################################################
#																																			GLOBAL VARIABLES																															 #
######################################################################################################################################################

#########################################
###     SPECIFY THE NUMBER OF CPU     ###
#########################################
our $max_cpu = 12; # Max number of CPU to use for the annotation


# Recover the current path
our $pwd = `pwd`;
chomp($pwd);

# Input file path
our @pathInput = split("/", $input);
# Output directories
our ($folderMutAnalysis, $folderAnnovar) = ("", "");
# File with the list of Annovar databases to use
our $listAVDB = "";
# Initialisation of chromosome, position, ref and alt values
our ($chrValue, $positionValue, $refValue, $altValue) = ("c", "s", "r", "a");


######################################################################################################################################################
#																																								MAIN 																																 #
######################################################################################################################################################
## Check the presence of the flags and create the output and temp directories
CheckFlags();

## Format the file in the correct format if they are vcf or MuTect output and recover the column positions
FormatingInputFile();

# Annotate the file with Annovar, add the strand orientation and the sequence context
FullAnnotation();

######################################################################################################################################################
#																																							FUNCTIONS																															 #
######################################################################################################################################################

## Check the presence of the flags and create the output and temp directories
sub CheckFlags
{
	# Check the reference genome
	if($refGenome eq "empty")   { print STDERR "You forget to specify the name for the reference genome!!!\nPlease specify it with the flag --refGenome\n"; exit 2; }
	if($intervalEnd eq "empty") { print STDERR "You forget to specify the length for the sequence context!!!\nPlease specify it with the flag --intervalEnd\n"; exit 2; }
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
	# Create the output folder for Annovar
	$folderAnnovar         = "$folderMutAnalysis/Annovar";
	if(!-e $folderAnnovar) { mkdir($folderAnnovar) or die "$!: $folderAnnovar\n"; }

	# Verify the access to Annovar databases
	if($path_AVDB eq "empty") { print STDERR "You forget to specify the path to Annovar databases!!!\nPlease specify it with the flag --pathAnnovarDB\n"; exit 2; }
	elsif(!-e $path_AVDB) { print STDERR"\nCan't access Annovar databases!\nPlease check the access to the disk\n"; exit 3; }

	# Check the file list AV DB
	if($pathAVDBList eq "empty") { print STDERR "You forget to specify the path to the list of Annovar databases!!!\nPlease specify it with the flag --pathAVDBList\n"; exit 2; }
	else { $listAVDB = "$pathAVDBList/${refGenome}_listAVDB.txt" }

	# If no temp folder is specified write the result in the current path
	if($folder_temp eq "empty") { $folder_temp   = "$pwd/TEMP_MutationalAnalysis_$pathInput[$#pathInput]"; }
	if(!-e $folder_temp)        { mkdir($folder_temp) or die "$!: $folder_temp\n"; }
}

## Format the file in the correct format if they are vcf or MuTect output and recover the column positions
sub FormatingInputFile
{
	# The input is a folder
	if(-d $input)
	{
		foreach my $file (`ls $input`)
		{
			my $headerOriginalFile = "";
			chomp($file);
			my ($filename, $directories, $suffix) = fileparse("$input/$file", qr/\.[^.]*/);

			CheckLengthFilename("$input/$file");

			#################################################
			###						Recover the input format 				###
			#################################################
			RecoverInputFormat("$input/$file", \$headerOriginalFile);
		}
	}
	# The input is one file
	else
	{
		my $headerOriginalFile = "";

		CheckLengthFilename($input);
		my ($filename, $directories, $suffix) = fileparse($input, qr/\.[^.]*/);

		#################################################
		###						Recover the input format 				###
		#################################################
		RecoverInputFormat($input, \$headerOriginalFile);
	}
}

# The name for the Excel sheet can't be longer than 31 characters
sub CheckLengthFilename
{
	my ($inputFile) = @_;

	## Verify the name of file, must be <= 31 chars for the sheet name
	my ($filename, $directories, $suffix) = fileparse($inputFile, qr/\.[^.]*/);

	if(length($filename) > 32) { print STDERR "The file: $inputFile must be <= 31 chars\nPlease modify it before running the script\n"; exit 4; }
}

# Recover the input format (vcf or txt) and depending on the format convert the input file in a suitable format for Annovar (ex: for MuTect files keep only the confident variants)
sub RecoverInputFormat
{
	my ($file, $refS_headerOriginalFile) = @_;

	my ($filename, $directories, $suffix) = fileparse($file, qr/\.[^.]*/);

	my $inputFormat = "";

	open(F1, $file) or die "$!: $file\n";
	my $header = <F1>;
	close F1;

	### VCF and MuTect files have in their first line the type of the file
	if($header =~ /fileformat=VCF/i) { $inputFormat = "vcf"; }
	elsif($header =~ /mutect/i)      { $inputFormat = "mutect"; }
	else                             { $inputFormat = "unknown"; }


	### VCF files
	if($inputFormat eq "vcf")
	{
		### MuTect2 output VCFs
		my $testVC = `grep MuTect2 $file`;
		if($testVC =~ /MuTect2/)
		{
			# Keep only the variants passing MuTect2 filters
			`grep PASS $file > $folder_temp/$filename-PASS.txt`;

			# Recover the header
			$$refS_headerOriginalFile = `grep '#CHROM' $file`;

			# Add the header
			# Sed command doesn't work... sed 's/^/some text\n/' file > res
			open(OUT, ">", "$folder_temp/$filename-KEEP.txt") or die "$!: $folder_temp/$filename-KEEP.txt\n";
			print OUT $$refS_headerOriginalFile;
			open(F1, "$folder_temp/$filename-PASS.txt") or die "$!: $folder_temp/$filename-PASS.txt\n";
			while(<F1>) { print OUT $_; }
			close F1; close OUT;

			`rm $folder_temp/$filename-PASS.txt`;

			# Check if there if no empty column
			CheckEmptyColumn("$folder_temp/$filename-KEEP.txt");
			`rm $folder_temp/$filename-KEEP.txt`;

			# Set the col number for the chr,start,ref and alt
			($chrValue, $positionValue, $refValue, $altValue) = (0, 1, 3, 4);
		}
		else
		{
			open(F1, $file) or die "$!: $file\n";
			open(OUT, ">", "$folder_temp/$filename.txt") or die "$!: $folder_temp/$filename.txt\n";
			while (<F1>)
			{
				$_      =~ s/[\r\n]+$//;
				my @tab = split("\t", $_);
				# Print the VCF header
				if($tab[0] eq "#CHROM")
				{
					$tab[0] =~ /#(.+)/;
					print OUT "$1";
					for(my $i=1; $i<=$#tab; $i++) { print OUT "\t$tab[$i]"; }
					print OUT "\n";
				}
				elsif($tab[0] !~ /##/)
				{
					# Don't consider chromosome random, GL and MT
					if( ($tab[0] =~ /random/) || ($tab[0] =~ /GL/i) ) { next; }
					print OUT "$_\n";
				}
			}
			close F1; close OUT;

			## Recover the header
			open(F1, "$folder_temp/$filename.txt") or die "$!: $folder_temp/$filename.txt\n";
			$$refS_headerOriginalFile = <F1>;
			close F1;

			# Check if there if no empty column
			CheckEmptyColumn("$folder_temp/$filename.txt");
			`rm $folder_temp/$filename.txt`;


			# Set the col number for the chr,start,ref and alt
			($chrValue, $positionValue, $refValue, $altValue) = (0, 1, 3, 4);
		}
	}
	### MuTect files
	elsif($inputFormat eq "mutect")
	{
		`sed '1d' $file > $folder_temp/$filename-HeaderOK`;
		# Keep only the SNVs of good quality
		`grep -v REJECT $folder_temp/$filename-HeaderOK > $folder_temp/$filename-KEEP.txt`;
		`rm $folder_temp/$filename-HeaderOK`;

		# Recover the header
		open(F1, "$folder_temp/$filename-KEEP.txt") or die "$!: $folder_temp/$filename-KEEP.txt\n";
		$$refS_headerOriginalFile = <F1>;
		close F1;

		# Check if there if no empty column
		CheckEmptyColumn("$folder_temp/$filename-KEEP.txt");
		`rm $folder_temp/$filename-KEEP.txt`;

		# Recover the name and the number of the columns that contain the chromosome number, the start position, the ref and alt alleles.
		# Use the dictionary for recover the names and the position of the columns
		RecoverColNameAuto("$folder_temp/$filename-KEEPColumnCorrect.txt", $$refS_headerOriginalFile, \$chrValue, \$positionValue, \$refValue, \$altValue);
	}
	### Unknown type
	else
	{
		## Recover the header
		open(F1, $file) or die "$!: $file\n";
		$$refS_headerOriginalFile = <F1>;
		close F1;

		# Check if there if no empty column
		CheckEmptyColumn($file);

		## Recover the name and the number of the columns that contain the chromosome number, the start position, the ref and alt alleles.
		# Use the dictionary for recover the names and the position of the columns
		RecoverColNameAuto($file, $$refS_headerOriginalFile, \$chrValue, \$positionValue, \$refValue, \$altValue);
	}
}

# Some files can have empty column with no information
sub CheckEmptyColumn
{
	my ($inputFile) = @_;

	my ($filename, $directories, $suffix) = fileparse($inputFile, qr/\.[^.]*/);

	if($filename =~ /(.+)-KEEP/) { $filename = $1; }

	open(OUT, ">", "$folder_temp/$filename-ColumnCorrect.txt") or die "$!: $folder_temp/$filename-ColumnCorrect.txt\n";

	open(F1, $inputFile) or die "$!: $inputFile\n";
	my $header = <F1>; $header =~ s/[\r\n]+$//; my @tabHeader = split("\t", $header);
	print OUT $header, "\n";
	while(<F1>)
	{
		$_      =~ s/[\r\n]+$//;
		my @tab = split("\t", $_);

		if(scalar(@tab) != scalar(@tabHeader))
		{
			print OUT $tab[0];
			for(my $i=1; $i<=$#tabHeader; $i++)
			{
				if(defined($tab[$i])) { print OUT "\t$tab[$i]"; }
				else { print OUT "\tNA"; }
			}
			print OUT "\n";
		}
		else { print OUT "$_\n"; }
	}
	close F1; close OUT;
}

# Dictionnary for extracting the name and number of columns for the chromosome, start position, ref and alt alleles.
sub RecoverColNameAuto
{
	our ($inputFile, $header, $ref_chrValue, $ref_positionValue, $ref_refValue, $ref_altValue) = @_;

	$header      =~ s/[\r\n]+$//;

	## Name of the columns
	my @mutect     = qw(contig position ref_allele alt_allele);
	my @vcf        = qw(CHROM POS REF ALT);
	my @cosmic     = qw(Mutation_GRCh37_chromosome_number Mutation_GRCh37_genome_position Description_Ref_Genomic Description_Alt_Genomic);
	my @icgc       = qw(chromosome chromosome_start reference_genome_allele mutated_to_allele);
	my @tcga       = qw(Chromosome Start_position Reference_Allele Tumor_Seq_Allele2);
	my @ionTorrent = qw(chr Position Ref Alt);
	my @proton     = qw(Chrom Position Ref Variant);
	my @varScan2   = qw(Chrom Position Ref VarAllele);
	my @annovar    = qw(Chr Start Ref Obs);
	my @custom     = qw(Chromosome Start Wild_Type Mutant);

	my @allTab = (\@mutect, \@vcf, \@cosmic, \@icgc, \@tcga, \@ionTorrent, \@proton, \@varScan2, \@annovar, \@custom);
	my $timer  = 0; # For controlling if the names are present on the dictionnary or not

	foreach my $refTab (@allTab)
	{
		my @tab = @$refTab;

		SearchCol(\@tab);

		# The columns names were find
		if( ($$ref_chrValue ne "c") && ($$ref_positionValue ne "s") && ($$ref_refValue ne "r") && ($$ref_altValue ne "a") ) { last; }
		# The names of the columns are not present in the dictionnary
		else { $timer++; }
	}

	if($timer == scalar(@allTab))
	{
		print STDERR "The columns name are not in the dictionnary please change them before running the tool again\nFile concerning: $inputFile\n";
		print STDERR "TIP: Use one of the columns names proposed in the section Input formats of the tool\n";
		exit 4;
	}

	# Extract the number of the column that contain the information
	sub SearchCol
	{
		my ($refTab) = @_;

		my @tabNames          = @$refTab;
		my @tabHeader         = split("\t", $header);

		# For VCF
		if($tabHeader[0] eq "#CHROM") { ($$ref_chrValue, $$ref_positionValue, $$ref_refValue, $$ref_altValue) = (0, 1, 3, 4); }
		# For tabular files
		else
		{
			for(my $i=0; $i<=$#tabNames; $i++)
			{
				for(my $j=0; $j<=$#tabHeader; $j++)
				{
					if($tabHeader[$j] eq $tabNames[$i])
					{
						if($i == 0) { $$ref_chrValue = $j; }
						if($i == 1) { $$ref_positionValue = $j; }
						if($i == 2) { $$ref_refValue = $j; }
						if($i == 3) { $$ref_altValue = $j; }
						last; # Once find pass to the next name
					}
				}
			}
		}
	}
}

# Annotate the file with Annovar, add the strand orientation and the sequence context
sub FullAnnotation
{
	print "-----------------------------------------------------------------\n";
	print "---------------------------Annotation----------------------------\n";
	print "-----------------------------------------------------------------\n";


	# If the input is a folder
	if(-d $input)
	{
		foreach my $file (`ls $folder_temp/*.txt`)
		{
			chomp($file);

			# For recover the name of the file without extension, the directory where the file is and the extension of the file
			my ($filename, $directories, $suffix)    = fileparse("$folder_temp/$file", qr/\.[^.]*/);
			my $filenameOK = "";
			# For removing the ColumnCorrect for txt files
			if($filename =~ /(.+)-ColumnCorrect/)
			{
				if($filename =~ /(.+)-VariantListVCF-ColumnCorrect/) { $filenameOK = $1; }
				else { $filenameOK = $1; }
			}
			else { print STDERR "Case not considered for $filename!!!\n"; exit 4; }


			#################################################
			###						 Cut the files in n part  		  ###
			#################################################
			# Recover the number of variants in the file for deciding the number of CPU to use
			my $cpu = 0;
			my $nbVariants = `wc -l $file`;
			$nbVariants =~ /(\d+).+/;
			my $nbLine  = $1;

			if($nbLine-1 <= 5000)      { $cpu = 1; }
			elsif( ($nbLine-1 > 5000) && ($nbLine-1 < 25000) )    { $cpu = 2; }
			elsif( ($nbLine-1 >= 25000) && ($nbLine-1 < 100000) ) { $cpu = 8; }
			else { $cpu = $max_cpu; }

			# If the number predefined can't be used on the machine use the maximum number specify by the administrator
			if($cpu > $max_cpu) { $cpu = $max_cpu }

			## Recover the header
			open(F1, $file) or die "$!: $file\n";
			my $headerOriginalFile = <F1>;
			close F1;

			## Remove the first line of the file
			my $fileNoHeader = "$folder_temp/${filenameOK}-NoHeader";
			`sed 1d $file > $fileNoHeader`;

			if(!-e "$folder_temp/$filenameOK") { mkdir("$folder_temp/$filenameOK") or die "Can't create the directory $folder_temp/$filenameOK\n"; }
			my $lines_per_temp = int(1+($1 / $cpu)); # +1 in case of the div == 0
			`split -l $lines_per_temp $fileNoHeader $folder_temp/$filenameOK/$filenameOK-`;

			if($headerOriginalFile eq "") { print STDERR "No header for the file $file!!!\nPlease check the format of your file\n"; exit 4; }
			my @files = <$folder_temp/$filenameOK/$filenameOK-*>;

			#################################################
			###							Annotate the n part  		 		  ###
			#################################################
			my $pm = Parallel::ForkManager->new($cpu);

			foreach my $tempFile (@files)
			{
				# Forks and returns the pid for the child:
			 	my $pid = $pm->start and next;

			 	# Convert the file in a correct format for Annovar: Chr Start End Ref Alt Otherinfo
			 	my ($filename, $directories, $suffix) = fileparse($tempFile, qr/\-[^.]*/);
				my $outFilenameTemp = $filename.$suffix;
				Convert2AV($tempFile, $chrValue, $positionValue, $refValue, $altValue, "$folder_temp/$outFilenameTemp-AVInput");

				# Annotate the file with Annovar
				my $tempFileName_AVOutput = $filename.$suffix.".".${refGenome}."_multianno.txt";
				if($fullAVDB eq "yes") { AnnotateAV("$folder_temp/$outFilenameTemp-AVInput", "$folder_temp/$outFilenameTemp"); }
				else { annotateAV_min("$folder_temp/$outFilenameTemp-AVInput", "$folder_temp/$outFilenameTemp"); }

				# Check if the annotations worked
				open(F1, "$folderMutAnalysis/log_annovar.txt") or die "$!: $folderMutAnalysis/log_annovar.txt\n";
				while(<F1>)
				{
					if($_ =~ /ERROR/i)
					{
						print STDERR "\n\n\t\tANNOVAR LOG FILE\n\n";
						print STDERR $_;
						print STDERR "\n\n\t\tANNOVAR LOG FILE\n\n\n";
						exit 5;
					}
				}
				close F1;

				# Recover the strand orientation
				my $length_AVheader = 0;
				RecoverStrand("$folder_temp/$tempFileName_AVOutput", $headerOriginalFile, $path_AVDB, $refGenome, "$folder_temp/$outFilenameTemp-Strand", \$length_AVheader);

				# Recover the sequence context
				RecoverGenomicSequence("$folder_temp/$outFilenameTemp-Strand", $length_AVheader, $intervalEnd, $refGenome, $path_AVDB, "$folder_temp/$filenameOK/$outFilenameTemp".".".${refGenome}."_multianno.txt");

				$pm->finish; # Terminates the child process
			}
			# Wait all the child process
			$pm->wait_all_children;



			#################################################
			###					Paste the file together 		 		  ###
			#################################################
			## For MuTect and MuTect2 calling only variants passing MuTect filters are kept and sometines there is no variant passing these filters making an error in Galaxy when using "collection".
			if($nbLine == 1)
			{
				print STDOUT "\nThe sample $filenameOK didn't pass MuTect filters\n";

				### Print Annovar minimal header + the original header of the input file
				my $outputFile = "$folderAnnovar/$filenameOK".".".${refGenome}."_multianno.txt";
				open(OUT, ">", $outputFile) or die "$!: $outputFile\n";

				if($fullAVDB eq "no")
				{
					print OUT "Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tGene.refGene\tGeneDetail.refGene\tExonicFunc.refGene\tAAChange.refGene\tStrand\tcontext";
					print OUT "\t".$headerOriginalFile;
				}
				### Print complete Annovar header (using the database name present in the file listAVDB) + the original header of the input file
				else
				{
					print OUT "Chr\tStart\tEnd\tRef\tAlt";
					open(F1, $listAVDB) or die "$!: $listAVDB\n";
					while(<F1>)
					{
						if($_ =~ /^#/) { next; }

						my @tab = split("\t", $_);
						$tab[0] =~ /$refGenome\_(.+)\.txt/;
						my $dbName = $1;

						if($dbName =~ /refGene|knownGene|ensGene/)
						{
							print OUT "\t"."Func.$dbName\tGene.$dbName\tGeneDetail.$dbName\tExonicFunc.$dbName\tAAChange.$dbName";
						}
						else
						{
							print OUT "\t".$dbName;
						}
					}
					print OUT "\tStrand\tcontext\t".$headerOriginalFile;

					close F1;
				}
				close OUT;
			}
			else
			{
				CombinedTempFile("$folder_temp/$filenameOK", "$folderAnnovar/$filenameOK".".".${refGenome}."_multianno.txt");
			}
		}
	}
	# The input file is one file
	else
	{
		my ($filenameO, $directoriesO, $suffixO) = fileparse($input, qr/\.[^.]*/);

		#################################################
		###						 Cut the files in n part  		  ###
		#################################################
		# Recover the number of variants in the file for deciding the number of CPU to use
		my $cpu = 0;
		my $nbVariants = `wc -l $folder_temp/$filenameO-ColumnCorrect.txt`;
		$nbVariants =~ /(\d+).+/;
		my $nbLine  = $1;

		if($nbLine-1 <= 5000)      { $cpu = 1; }
		elsif( ($nbLine-1 > 5000) && ($nbLine-1 < 25000) )    { $cpu = 2; }
		elsif( ($nbLine-1 >= 25000) && ($nbLine-1 < 100000) ) { $cpu = 8; }
		else { $cpu = $max_cpu; }

		# If the number predefined can't be used on the machine use the maximum number specify by the administrator
		if($cpu > $max_cpu) { $cpu = $max_cpu }

		## Recover the header
		open(F1, "$folder_temp/$filenameO-ColumnCorrect.txt") or die "$!: $folder_temp/$filenameO-ColumnCorrect.txt\n";
		my $headerOriginalFile = <F1>;
		close F1;

		## Remove the first line of the file
		my $fileNoHeader = "$folder_temp/$filenameO-NoHeader";
		`sed 1d $folder_temp/$filenameO-ColumnCorrect.txt > $fileNoHeader`;

		if(!-e "$folder_temp/$filenameO") { mkdir("$folder_temp/$filenameO") or die "Can't create the directory $folder_temp/$filenameO\n"; }
		my $lines_per_temp = int(1+($1 / $cpu)); # +1 in case of the div == 0
		`split -l $lines_per_temp $fileNoHeader $folder_temp/$filenameO/$filenameO-`;

		if($headerOriginalFile eq "") { print STDERR "No header for the file $input!!!\nPlease check the format of your file\n"; exit 4; }
		my @files = <$folder_temp/$filenameO/$filenameO-*>;

		#################################################
		###							Annotate the n part  		 		  ###
		#################################################
		my $pm = Parallel::ForkManager->new($cpu);
		foreach my $tempFile (@files)
		{
			# Forks and returns the pid for the child:
			my $pid = $pm->start and next;

			# Convert the file in a correct format for Annovar: Chr Start End Ref Alt Otherinfo
			# For recover the name of the file without extension, the directory were the file is and the extension of the file
			my ($filename, $directories, $suffix) = fileparse($tempFile, qr/\.[^.]*/);
			my $outFilenameTemp = $filename.$suffix;
			Convert2AV($tempFile, $chrValue, $positionValue, $refValue, $altValue, "$folder_temp/$outFilenameTemp-AVInput");

			# Annotate the file with Annovar
			my $tempFileName_AVOutput = $outFilenameTemp.".".${refGenome}."_multianno.txt";
			if($fullAVDB eq "yes") { AnnotateAV("$folder_temp/$outFilenameTemp-AVInput", "$folder_temp/$outFilenameTemp"); }
			else { annotateAV_min("$folder_temp/$outFilenameTemp-AVInput", "$folder_temp/$outFilenameTemp"); }

			# Check if the annotations worked
				open(F1, "$folderMutAnalysis/log_annovar.txt") or die "$!: $folderMutAnalysis/log_annovar.txt\n";
				while(<F1>)
				{
					if($_ =~ /ERROR/i)
					{
						print STDERR "\n\n\t\tANNOVAR LOG FILE\n\n";
						print STDERR $_;
						print STDERR "\n\n\t\tANNOVAR LOG FILE\n\n\n";
						exit 5;
					}
				}
				close F1;

			# Recover the strand orientation
			my $length_AVheader = 0;
			RecoverStrand("$folder_temp/$tempFileName_AVOutput",  $headerOriginalFile, $path_AVDB, $refGenome, "$folder_temp/$outFilenameTemp-Strand", \$length_AVheader);

			# Recover the sequence context
			RecoverGenomicSequence("$folder_temp/$outFilenameTemp-Strand", $length_AVheader, $intervalEnd, $refGenome, $path_AVDB, "$folder_temp/$filenameO/$tempFileName_AVOutput");

			$pm->finish; # Terminates the child process
		}
		# Wait all the child process
		$pm->wait_all_children;

		#################################################
		###					Paste the file together 		 		  ###
		#################################################
		## For MuTect and MuTect2 calling only variants passing MuTect filters are kept and sometines there is no variant passing these filters making an error in Galaxy for the next tool when using "Collection".
		if($nbLine == 1)
		{
			print STDOUT "\nThe sample $filenameO didn't pass MuTect filters\n";

			### Print Annovar minimal header + the original header of the input file
			my $outputFile = "$folderAnnovar/$filenameO".".".${refGenome}."_multianno.txt";
			open(OUT, ">", $outputFile) or die "$!: $outputFile\n";

			if($fullAVDB eq "no")
			{
				print OUT "Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tGene.refGene\tGeneDetail.refGene\tExonicFunc.refGene\tAAChange.refGene\tStrand\tcontext";
				print OUT "\t".$headerOriginalFile;
			}
			### Print complete Annovar header (using the database name present in the file listAVDB) + the original header of the input file
			else
			{
				print OUT "Chr\tStart\tEnd\tRef\tAlt";
				open(F1, $listAVDB) or die "$!: $listAVDB\n";
				while(<F1>)
				{
					if($_ =~ /^#/) { next; }

					my @tab = split("\t", $_);
					$tab[0] =~ /$refGenome\_(.+)\.txt/;
					my $dbName = $1;

					if($dbName =~ /refGene|knownGene|ensGene/)
					{
						print OUT "\t"."Func.$dbName\tGene.$dbName\tGeneDetail.$dbName\tExonicFunc.$dbName\tAAChange.$dbName";
					}
					else
					{
						print OUT "\t".$dbName;
					}
				}
				print OUT "\tStrand\tcontext\t".$headerOriginalFile;

				close F1;
			}
			close OUT;
		}
		else
		{
			CombinedTempFile("$folder_temp/$filenameO", "$folderAnnovar/$filenameO".".".${refGenome}."_multianno.txt");
		}
	}
	# Remove the temporary directory
	rmtree($folder_temp);
}

sub Convert2AV
{
	my ($inputFile, $chr_value, $start_value, $ref_value, $alt_value, $output) = @_;

	my ($filename, $directories, $suffix) = fileparse($inputFile, qr/\.[^.]*/);

	open(F1, $inputFile) or die "$!: $inputFile\n";

	open(OUT, ">", $output) or die "$!: $output\n";
	while(<F1>)
	{
		$_ =~ s/[\r\n]+$//;
		my @tab = split("\t", $_);
		my $chr = "";

		# Don't consider chrM and GL
		if($tab[$chr_value] =~ /M|GL/i) { next; }

		# Replace chr23 or chr24 by X or Y
		if($tab[$chr_value] =~ /23/)     { $chr = "chrX"; }
		elsif($tab[$chr_value] =~ /24/)  { $chr = "chrY"; }
		elsif($tab[$chr_value] =~ /chr/) { $chr = $tab[$chr_value]; }
		else                             { $chr = "chr".$tab[$chr_value]; }

		### Reformat the Indels for Annovar
			# chr1	85642631	C	    CT  => chr1	85642631	85642631	-	  T   (mm10)
			# chr5	26085724	ACTT	A   => chr5	26085725	26085727	CTT	-   (mm10)
			if( ((length($tab[$ref_value]) != 1) || (length($tab[$alt_value]) != 1)) || (($tab[$ref_value] eq "-") || ($tab[$alt_value] eq "-") ) )
			{
				### First check if the indels in the file are not already correctly formated
				if( ($tab[$ref_value] eq "-") || ($tab[$alt_value] eq "-") )
				{
					# For indels count the number of bases deleted or inserted for modifying the end position (if start + end is the same the annotations are not retrieved for indels)
					# Insertion: start = start & end = start
					if($tab[$ref_value] =~ /\-/)
					{
						print OUT "$chr\t$tab[$start_value]\t$tab[$start_value]\t$tab[$ref_value]\t$tab[$alt_value]";
					}
					## Deletion: start = start & end = start + length(del) -1
					else
					{
						my $end = $tab[$start_value] + (length($tab[$ref_value]) - 1);
						print OUT "$chr\t$tab[$start_value]\t$end\t$tab[$ref_value]\t$tab[$alt_value]";
					}
				}
				### Indels not correctly formated for Annovar
				else
				{
					my @tabRef = split("", $tab[$ref_value]);
					my @tabAlt = split("", $tab[$alt_value]);

					# Remove the first base
					my $ref2 = join("", @tabRef[1 .. $#tabRef]);
					my $alt2 = join("", @tabAlt[1 .. $#tabAlt]);

					if(length($alt2) == 0)
					{
						my $altOK   = "-";
						my $startOK = $tab[$start_value] + 1;
						my $stopOK  = $startOK + length($ref2) - 1;
						print OUT $chr."\t".$startOK."\t".$stopOK."\t".$ref2."\t".$altOK;
					}

					if(length($ref2) == 0)
					{
						my $refOK = "-";
						print OUT $chr."\t".$tab[$start_value]."\t".$tab[$start_value]."\t".$refOK."\t".$alt2;
					}
				}
			}
			### SBS
			else
			{
				print OUT $chr."\t".$tab[$start_value]."\t".$tab[$start_value]."\t".$tab[$ref_value]."\t".$tab[$alt_value];
			}

			## Print the original file at the end
			foreach  (@tab) {  print OUT "\t$_"; }
			print OUT "\n";
	}
	close F1; close OUT;
}

sub AnnotateAV
{
	my ($inputFile, $output) = @_;

	if(!-e $path_AVDB) { print STDERR "The Annovar database doesn't exists for the reference genome $refGenome!!!\n"; print STDERR "Please install the database for this genome before running Annovar\n"; exit 4; }

	# Extract the name of the databases
	my $protocol = ""; my $operation = "";
	ExtractAVDBName($listAVDB, \$protocol, \$operation);

	`table_annovar.pl $inputFile $path_AVDB -buildver $refGenome -protocol $protocol -operation $operation -remove -nastring NA -otherinfo -outfile $output > $folderMutAnalysis/log_annovar.txt 2>&1`;

	sub ExtractAVDBName
	{
		my ($listAVDB, $refS_protocol, $refS_operation) = @_;

		open(F1, $listAVDB) or die "$!: $listAVDB\n";
		while(<F1>)
		{
			if ($_ =~ /^#/) { next; }

			$_      =~ s/[\r\n]+$//;
			my @tab = split("\t", $_);

			# db name like refGenome_dbName.txt
			if( ($tab[0] =~ /\w+_(\w+)\.txt/) && ($tab[0] !~ /sites/) && ($tab[0] !~ /esp/) && ($tab[0] !~ /ljb26/) )
			{
				$$refS_protocol .= $1.","; $$refS_operation .= $tab[1].",";
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
				$$refS_protocol .=$AVdbName_final.","; $$refS_operation .= $tab[1].",";
			}
			# ESP
			if( ($tab[0] =~ /esp/) || ($tab[0] =~ /ljb26/) )
			{
				$tab[0] =~ /\w+_(\w+)_(\w+)\.txt/;
				my $AVdbName_final = $1."_".$2;
				$$refS_protocol .=$AVdbName_final.","; $$refS_operation .= $tab[1].",";
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
		}
	}
}

### Add the minimum of annotations (refGene + strand + context)
sub annotateAV_min
{
	my ($inputFile, $output) = @_;

	if(!-e $path_AVDB) { print STDERR "The Annovar database doesn't exists for the reference genome $refGenome!!!\n"; print STDERR "Please install the database for this genome before running Annovar\n"; exit 4; }

	# Extract the name of the databases
	my ($protocol, $operation) = ("refGene", "g");

	`table_annovar.pl $inputFile $path_AVDB -buildver $refGenome -protocol $protocol -operation $operation -remove -nastring NA -otherinfo -outfile $output > $folderMutAnalysis/log_annovar.txt 2>&1`;
}

sub RecoverStrand
{
	my ($input, $headerOriginalFile, $pathDB, $refGenome, $output, $refS_lengthAVheader) = @_;

	my ($chr_value, $start_value, $ref_value, $alt_value, $func_value, $geneSymbol_value) = ("", "", "", "", "", "", "", "");

	$chr_value        = recoverNumCol($input, "Chr");
	$start_value      = recoverNumCol($input, "Start");
	$ref_value        = recoverNumCol($input, "Ref");
	$alt_value        = recoverNumCol($input, "Alt");
	$func_value       = recoverNumCol($input, "Func.refGene");
	$geneSymbol_value = recoverNumCol($input, "Gene.refGene");

	#################### Convert the input file into a hash table
	my %h_inputFile = ();
	open(F1, $input) or die "$!: $input\n";
	my $annovar_header  = <F1>;

	while(<F1>)
	{
		$_ =~ s/[\r\n]+$//;
		my @tab = split("\t", $_);

		# In COSMIC the chromosome X and Y are annotated 23 and 24
		my $chr = "";
		if($tab[$chr_value] eq "chr23")    { $chr = "chrX"; }
		elsif($tab[$chr_value] eq "chr24") { $chr = "chrY"; }
		elsif($tab[$chr_value] eq "chr25") { $chr = "chrM"; }
		else { $chr = $tab[$chr_value]; }

		# Verify if the element exists
		if($chr eq "")                       { print STDERR "Error RecoverStrand: The chromosome value is nor defined for $_\n"; exit 4; }
		if(! exists $tab[$start_value])      { print STDERR "Error RecoverStrand: The start value is nor defined for $_\n"; exit 4; }
		if(! exists $tab[$ref_value])        { print STDERR "Error RecoverStrand: The reference value is nor defined for $_\n"; exit 4; }
		if(! exists $tab[$alt_value])        { print STDERR "Error RecoverStrand: The alternate value is nor defined for $_\n"; exit 4; }
		if(! exists $tab[$func_value])       { print STDERR "Error RecoverStrand: The functional value is nor defined for $_\n"; exit 4; }
		if(! exists $tab[$geneSymbol_value]) { print STDERR "Error RecoverStrand: The gene symbol value is nor defined for $_\n"; exit 4; }

		my $geneSymbol = "";
		######## For the splicing annotation we separate the gene symbol from the aa change
		if($tab[$func_value] eq "splicing")
		{
			if($tab[$geneSymbol_value] =~ /(.+)\((.+)\)/) { $geneSymbol = $1; }
			else { $geneSymbol = $tab[$geneSymbol_value]; }
		}
		else { $geneSymbol = $tab[$geneSymbol_value]; }

		push(@{$h_inputFile{"$chr:$tab[$start_value]:$tab[$start_value]:$tab[$ref_value]:$tab[$alt_value]:$geneSymbol"}}, $_);
	}
	close F1;

	# print "\t\tRecoverStrand: $input\n";

	#################### Convert the database file into a hash table
	my %h_database  = ();
	my ($db_geneSymbol_value, $db_strandInfo_value, $db_chr_value) = (12, 3, 2);

	my $folderNameDB = $refGenome."db";
	my $fileNameDB   = $refGenome."_refGene.txt";

	open(F1, "$pathDB/$fileNameDB") or die "$!: $pathDB/$fileNameDB\n";
	while(<F1>)
	{
		$_ =~ s/[\r\n]+$//;
		my @tab = split("\t", $_);
		my $strand = "";
		$strand = $tab[$db_strandInfo_value];
		if($strand eq "") { print STDERR "Error: the strand orientation is not specify in the database refGene\n$_\n"; exit 4; }
		else
		{
			# Some genes have several strand orientation, keep the first in the database
			if(! exists $h_database{"$tab[$db_geneSymbol_value]:$tab[$db_chr_value]"}) { $h_database{"$tab[$db_geneSymbol_value]:$tab[$db_chr_value]"} = $strand; }
		}
	}
	close F1;

	#################### Parse the two hash tables for recover the strand information
	open(OUT, ">", $output) or die "$!: $output\n";


	## Add the header only for the firts part of the files
	if($input =~ /\-aa/)
	{
		my @tabHeaderInput  = "";
		$annovar_header =~ s/[\r\n]+$//; @tabHeaderInput = split("\t", $annovar_header);
		# Save the length of the Annovar header for the next function (RecoverGenomicSequence)
		$$refS_lengthAVheader = $#tabHeaderInput;

		# Print the Annovar header until the column before OtherInfo
		print OUT "$tabHeaderInput[0]";
		for(my $i=1; $i<$#tabHeaderInput; $i++) { print OUT "\t$tabHeaderInput[$i]"; }
		print OUT "\tStrand";
		print OUT "\t",$headerOriginalFile;
	}


	# Timer for comparing the number of SNVs present in the hash table
	my $timerUniqueSNVs = 0;
	# Timer for comparing the number of SNVs with the strand orientation
	my $timerSNVsStrand = 0;

	foreach my $kFile (sort keys %h_inputFile)
	{
		my $test = 0;
		my @tab = split(":", $kFile);

		# Sometimes the line is not printed correctely !!!!! :@
		my @tHeaderInput        = split("\t", $annovar_header); my @lengthLine = split("\t", $h_inputFile{$kFile}[0]);
		my @tHeaderOriginalFile = split("\t", $headerOriginalFile);
		my $lengthHeader = @tHeaderInput + (scalar(@tHeaderOriginalFile)-1) ; my $lengthLine = @lengthLine;

		# Save the length of the Annovar header for the next function (RecoverGenomicSequence)
		$$refS_lengthAVheader = $#tHeaderInput;

		foreach my $kDB (sort keys %h_database)
		{
			if("$tab[5]:$tab[0]" eq $kDB)
			{
				if($lengthHeader != $lengthLine) { print STDERR "Error Recover Strand the length of the current line is not valid!!!!!\nExpected length: $lengthHeader\tlength of the line: $lengthLine\n$h_inputFile{$kFile}[0]\n"; exit 4; }

				foreach my $line (@{$h_inputFile{$kFile}})
				{
					my @tab = split("\t", $line);
					my $j = 0;

					for(my $i=0; $i<$#tHeaderInput; $i++) { print OUT $tab[$i],"\t"; $j=$i }
					print OUT $h_database{$kDB};
					for(my $i=$j+1; $i<=$#tab; $i++) { print OUT "\t$tab[$i]"; }
					print OUT "\n";
				}
				$timerSNVsStrand++;
				$test = 1; last;
			}
		}
		# The strand orientation isn't defined
		if($test == 0)
		{
			my @tHeaderInput = split("\t", $annovar_header);
			foreach my $line (@{$h_inputFile{$kFile}})
			{
				my @tab = split("\t", $line);
				my $j = 0;
				for(my $i=0; $i<$#tHeaderInput; $i++) { print OUT $tab[$i],"\t"; $j=$i }
				print OUT "NA";
				for(my $i=$j+1; $i<=$#tab; $i++) { print OUT "\t$tab[$i]"; }
				print OUT "\n";
			}
			$timerSNVsStrand++;
		}
		$timerUniqueSNVs++;
	}
	close OUT;

	# print "Strand orientation recover for $timerSNVsStrand SNVs out of $timerUniqueSNVs uniques\n";
}

sub RecoverGenomicSequence
{
	my ($inputFile, $length_AVheader, $intervalEnd, $referenceGenome, $pathToRefSeq, $output) = @_;

	############ 1) Transform the input file in a hash table: one for recover the sequence context and one for keeping the original file
	my %h_inputFileForSeqContext = (); my %h_inputFile = ();
	my $header                   = "";
	CreateHashTable_from_InputFile($inputFile, $length_AVheader, \$header, $intervalEnd, \%h_inputFileForSeqContext, \%h_inputFile);

	sub CreateHashTable_from_InputFile
	{
		my ($input, $length_AVheader, $refS_header, $intervalEnd, $refH_inputFileForSeqContext, $refH_inputFile) = @_;

		my ($chr_value, $start_value, $strand_value) = (0, 1, $length_AVheader);

		my $countregion   = 0;
		my %allchr        = ();

		open(F1, $input) or die "$!: $input\n";
		if($input =~ /\-aa/) { $$refS_header = <F1>; }

		while(<F1>)
		{
			$_ =~ s/[\r\n]+$//;
			my @tab = split("\t", $_);

			my $name  = "$tab[$chr_value]:$tab[$start_value]";
			my $start = $tab[$start_value] - $intervalEnd;
			my $end   = $tab[$start_value] + $intervalEnd;

			$start--;		#make zero-start coordinate, to be consistent with UCSC
			my $exonpos = "$tab[$chr_value]:$start";

			push @{$refH_inputFileForSeqContext->{$tab[$chr_value]}}, [$name, $start, $end, $tab[$strand_value], $exonpos];
			push(@{$refH_inputFile->{"$tab[$chr_value]\t$start\t$end"}}, $_);
			$countregion++;
			$allchr{$tab[$chr_value]}++;
		}
		close F1;
	}

	############ 2) Extract the sequence context from the hash table
	my %h_allRegionSeqContext = ();
	my $refSeq = $pathToRefSeq;
	Extract_SequenceContext(\%h_inputFileForSeqContext, $referenceGenome, $refSeq, \%h_allRegionSeqContext);

	sub Extract_SequenceContext
	{
		my ($refH_allRegion, $referenceGenome, $refSeq, $refH_allRegionSeqContext) = @_;

		my $folderDB  = $referenceGenome."db";
		my $folderSeq = $referenceGenome."_seq";
		my $seqdir = "$refSeq/$folderSeq";

		my %seqhash  = (); #database sequence for each chromosome
		my %name_seq = (); #sequence for each region
		my (%seqlen, %discordlen, %badorf);	#store the length of each sequence, and the ID of sequences with discordant length, ORF contains stop codon
		my ($count_success, @failure) = (0);

		for my $curchr (sort keys $refH_allRegion)
		{
			my ($seqid, $curseq) = ('', '');
			my $fastafile        = "";
			if ($curchr =~ m/^chr/)
			{
				%seqhash   = (); #clear the seqhash storage
				$fastafile = "$seqdir/$curchr.fa"; #by default, all FASTA files should be saved at fastadir, with the same name
			}
			else
			{
				%seqhash   = ();		#clear the seqhash storage
				$fastafile = "$seqdir/chr$curchr.fa"; #by default, all FASTA files should be saved at fastadir, with the same name
			}
			if (not -e $fastafile) {			#to handle cases where no "chr" prefix is given
				print "WARNING: the FASTA file $curchr.fa cannot be retrieved from the specified directory $seqdir. Sequences in this chromosome will not be processed\n";
				next;
			}

			if (not %seqhash)
			{
				open (FASTA, $fastafile) or print "WARNING: cannot read from FASTA file $fastafile so sequences in $curchr will not be processed: $!\n" and next;
				while (<FASTA>)
				{
					if (m/^>(\S+)/)
					{
						$seqid and $seqhash{$seqid} = $curseq; #finish reading the sequence for seqid and save it
						$seqid = $1;
						$curseq = '';
					}
					else
					{
						s/[\r\n]+$//;
						$curseq .= uc $_; #only use upper case characters
					}
				}
				close FASTA;
				$seqhash{$seqid} = $curseq;
			}
			if (not $seqhash{$curchr})
			{
				#this chromosome just do not have FASTA sequences (maybe users used a wrong seqdir
				print "WARNING: Unable to retrieve regions at $curchr due to lack of sequence information\n";
				next;
			}

			for my $i (0 .. @{$refH_allRegion->{$curchr}}-1)
			{
				my ($name, $start, $end, $strand, $exonpos) = @{$refH_allRegion->{$curchr}[$i]};
				my @start = split (/,/, $start);
				my @end   = split (/,/, $end);
				my $seq;
				for my $i (0..@start-1)
				{
					if ($start[$i] >= length ($seqhash{$curchr}))
					{
						#here there must be an annotation error in user-specified gene/region definition file
						print "WARNING: Ignoring the start position start=$start[$i] since it is longer than the $curchr sequence (length=" , length($seqhash{$curchr}), ")\n";
						undef $seq;
						last;
					}
					$seq .= substr ($seqhash{$curchr}, $start[$i], $end[$i]-$start[$i]);
				}

				if (defined $seq)
				{
					if (defined $seqlen{$name})
					{
						$seqlen{$name} != length ($seq) and warn "WARNING: the sequence $name was found more than once with different sequence lengths\n";
						$seqlen{$name} != length ($seq) and $discordlen{$name}++;
					}
					else { $seqlen{$name} = length ($seq); }

					$name_seq{$name, $exonpos} = $seq;
					$count_success++;

					# Put the sequence context in a hash table for Write the result after
					## Some sequence context are NNNNNN or empty
					if( ($seq ne "NA") && ($seq =~ /N/i) ) { $refH_allRegionSeqContext->{"$curchr\t$start\t$end"} = "NA"; }
					else { $refH_allRegionSeqContext->{"$curchr\t$start\t$end"} = $seq; }
				}
				else
				{
					print "WARNING: DNA sequence for $name cannot be inferred\n";
					push @failure, $name;
				}
			}
		} # End for $curchr
	}

	############ 3) Create a file with the sequence context
	WriteFile_SeqContext($inputFile, $length_AVheader, \%h_inputFile, $header, \%h_allRegionSeqContext, $output);

	sub WriteFile_SeqContext
	{
		my ($inputFile, $length_AVheader, $refH_InputFile, $header, $refH_allRegionSeqContext, $output) = @_;

		open(OUT, ">", $output) or die "$!: $output\n";

		## Add the header only for the firts part of the files
		if($inputFile =~ /\-aa/)
		{
			my @tabHeaderInput  = "";

			$header =~ s/[\r\n]+$//; @tabHeaderInput = split("\t", $header);
			# Print the Annovar header until the column before OtherInfo
			print OUT "$tabHeaderInput[0]";
			my $j = 0;
			for(my $i=1; $i<$length_AVheader+1; $i++) { print OUT "\t$tabHeaderInput[$i]"; $j=$i; }
			print OUT "\tcontext";
			for(my $i=$j+1; $i<=$#tabHeaderInput; $i++) { print OUT "\t$tabHeaderInput[$i]"; }
			print OUT "\n";
		}

		foreach my $k_hFile (sort keys $refH_InputFile)
		{
			foreach my $k_allRegonSeqContext (sort keys $refH_allRegionSeqContext)
			{
				if($k_hFile eq $k_allRegonSeqContext)
				{
					my $j=0;

					for(my $k=0; $k<=$#{$refH_InputFile->{$k_hFile}};$k++)
					{
						my @tab = split("\t", ${$refH_InputFile->{$k_hFile}}[$k]);

						for(my $i=0; $i<$length_AVheader+1; $i++) { print OUT $tab[$i],"\t"; $j=$i; }
						print OUT $refH_allRegionSeqContext->{$k_allRegonSeqContext};
						for(my $i=$j+1; $i<=$#tab; $i++) { print OUT "\t$tab[$i]"; }
						print OUT "\n";
					}
					last;
				}
			}
		}
		close OUT;
	}
}

sub CombinedTempFile
{
	my ($folderTempFile, $output) = @_;

	my $cmd_cat_mt_results = "cat ";

	foreach my $file (`ls $folderTempFile/*.txt`)
	{
		chomp($file);
		$cmd_cat_mt_results = $cmd_cat_mt_results." $file";
	}
	$cmd_cat_mt_results = $cmd_cat_mt_results." > $output";
	`$cmd_cat_mt_results`;
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
  if($name_of_column_NB eq "toto") { print STDERR "Error recoverNumCol(): the column named $name_of_column doesn't exits in the input file $input!!!!!\n"; exit 4; }
  else                             { return $name_of_column_NB; }
}




=head1 NAME

mutspec-Annot

=head1 SYNOPSIS

	mutspecannot.pl [arguments] <query-file>

  <query-file>                                   can be a folder with multiple VCF or a single VCF

  Arguments:
        -h,        --help                        print help message
        -m,        --man                         print complete documentation
        -v,        --verbose                     use verbose output
                   --refGenome                   the reference genome to use
                   --interval <interger>         the number of bases for the sequence context
        -o,        --outfile <string>            output directory for the result. If none is specify the result will be write in the same directory as the input file
        -AVDB      --pathAnnovarDB <string>      the path to Annovar database and the files with the chromosome size
                   --pathAVDBList                the path to the list of AV databases installed
        -temp      --pathTemporary <string>      the path for saving the temporary files
                   --fullAnnotation <string>     recover all Annovar annotations (yes) or only the minimum for MutSpec-Stat (no)


Function: automatically run a pipeline on a list of variants and annote them using Annovar

 Example: # Annotation only
          mutspecannot.pl --refGenome hg19 --interval 10 --outfile output_directory --pathAnnovarDB path_to_annovar_database --pathAVDBList path_to_the_list_of_annovar_DB --temp path_to_temporary_directory --fullAnnotation yes|no input


 Version: 10-2016 (October 2016)


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

=item B<--interval>

the number of bases surrounding the mutated bases, for the sequence context analysis.

=item B<--outfile>

the directory of output file names. If it is nor specify the same directory as the input file is used.

=item B<--pathAnnovarDB>

the path to the directory containing the Annovar databases and the files with the chromosome size.

=item B<--pathAVDBList>

the path to a texte file containing the list of the Annovar databases installed.

=item B<--pathTemporary>

the path for saving temporary files generated by the script.
If any is specify a temporary folder is created in the same directory where the script is running.
Deleted when the script is finish

=item B<--fullAnnotation>

Use all Annovar databases for the annotation (set to yes) or only refGene + strand + context (set to no) for having a quicker annotation (for large file with million of lines)

=head1 DESCRIPTION

MutSpec-Annot is a perl script for added annotations on a list of genetic variants generated with NGS.
Functional annotations are added using ANNOVAR software. Strand transcript orientation is added using RefSeq database and the sequence context for x bases flanking the variant positions is also added.
A text tab delimited file is produced.

=cut
