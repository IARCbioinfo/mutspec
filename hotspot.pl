#/usr/bin/perl

#	~~~~	HOTSPOT tool	~~~~
#	13/06/2016
#	Alexis ROBITAILLE
#	robitaillea@students.iarc.fr
#	Version : 1.0

#################
#				#
#	HOTSPOT		#
#				#
#################

#	- The goal of this tool is to compute variant frequency in a define datasets
#	- If information provided on duplicate sample, infers somatics status of the variants in this dataset

########################################################################################################################################
#Input :	1) Directory of VCF/tabular variants files (same format)
#			2) Tabular file containing the name of the input files
#			The variants file for a same individual mlust be in the same line
#				2a) One column : not differenciated files (ex : all tumor sample)
#				2b) Two column : differenciated files (ex : normal-tumor sample, same individual)
#				2c) Three or more column : differenciated files + duplicates (ex : normal-tumor-duplicates, same individual)
#
#Output :1) Variants_Summary.vcf = One uniqu variant per line + annotation on the frequence of this variant in this dataset
#			If input 2a : frequency in all the dataset
#			If input 2b : frequency in the two separate series	
#		 2) Collection of the same input 2) files + annotation as describe below (frequency of variant in the dataset of differenciated or not series)
#			If input 2c : Add annotation on somatic status of the variants
########################################################################################################################################

#	Library
use strict;
use warnings;
use Getopt::Std;
use Math::Round ':all';

#	Options
my %opts = ();
getopts( 'd:o:i:s:p:h', \%opts ) or print_usage();
if ( defined( $opts{d} ) ) {
	#print "-d $opts{d}\n";
}
else {
	do_help();
	exit;
}
if ( defined( $opts{i} ) ) {
	#print "-i $opts{i}\n";
}
else {
	do_help();
	exit;
}
if ( defined( $opts{o} ) ) {
	#print "-o $opts{o}\n";
}
else {
	do_help();
	exit;
}
if ( defined( $opts{s} ) ) {
	#print "-s $opts{s}\n";
}
else {
	do_help();
	exit;
}
if ( defined( $opts{p} ) ) {
	#print "-p $opts{p}\n";
}
else {
	do_help();
	exit;
}
if ( $opts{h} ) {
	do_help();
}

sub do_help {
	printf "
	HOTSPOT tool :
The goal of this tool is to compute variant frequency in a define datasets.
If information provided on duplicate sample, infers somatics status of the variants in this dataset.

	Option -i : INPUT : tabular file with sample name
	Option -d : INPUT : directory name where VCF/tabular file are located
	Option -o : OUTPUT: directory name
	Option -s : OUTPUT: directory name for variant_summary
	Option -p : INPUT: Paires analysis (Y/N)
	Option -h : describe the help\n\n"
}


my $fileformat;		#Format of the input VCF/tabular variants files
my @Files;			#Table of path of the input VCF/tabular variants files
my $InfoFile;		#Path of tabular file containing name of input VCF/tabular variants files
my $check; 			#Provide the information of why the infofile is wrong, or in case of rigth, it's a hash of the series contain in the file --> become $serie variable
my $cat;			#Ref to hashtable K=categorie, V=ref of table containing sample name that are in this categorie
my $uniqs;			#Ref to table containing the name of the sample in InfoFile (uniqu)
my @filesformat;	#Table of format of vcf input files
my $series;			#Ref to hashtable K=integer, V=ref of table containing sample name (string sample name in case of normal only)

######################
##	Accepted format	## of the input VCF/tabular file in the directory opts{d}
######################
my @mutect2=("#CHROM", "POS","REF", "ALT");
my @mutect=("contig", "position", "ref_allele", "alt_allele");
my @vcf=("CHROM", "POS", "REF", "ALT");
my @cosmic=("Mutation_GRCh37_chromosome_number", "Mutation_GRCh37_genome_position", "Description_Ref_Genomic", "Description_Alt_Genomic");
my @icgc=("chromosome", "chromosome_start", "reference_genome_allele", "mutated_to_allele");
my @tcga=("Chromosome", "Start_position", "Reference_Allele", "Tumor_Seq_Allele2");
my @ionTorrent=("chr", "Position", "Ref", "Alt");
my @proton=("Chrom", "Position", "Ref", "Variant");
my @varScan2=("Chrom", "Position", "Ref", "VarAllele");
my @annovar=("Chr", "Start", "Ref", "Obs");
my @custom=("Chromosome", "Start", "Wild_Type", "Mutant");
my @GATK16=("Chromosome", "Position", "Reference Allele", "Variant Allele");
my @annotated=("Chr", "Start", "Ref", "Alt");

my %acceptedformat=(
"mutect2" => \@mutect2,
"mutect" => \@mutect,
"vcf" => \@vcf,
"cosmic" => \@vcf,
"icgc" => \@icgc,
"tcga" => \@tcga,
"ionTorrent" => \@ionTorrent,
"proton" => \@proton,
"varScan2" => \@varScan2,
"annovar" => \@annovar,
"custom" => \@custom,
#"GATK16" => \@GATK16,
"annotated" => \@annotated
);
######################


$InfoFile=$opts{i};						#Get the path to InfoFile

if($InfoFile eq "None"){
	if($opts{p} eq "pair"){
		print "Program STOP - You have to provide an InfoFile if you want to do a paired analysis\n";
		exit;
	}
}



###################################################################################################
#	This part concern the verification of the integrity of the InfoFile, and the fact that it fit with the sample file name present in input directory
#
my $type="";														#Type Normal-Tumor-Dup or Traitement
if($InfoFile ne "None"){
	($check,$cat,$uniqs)=CheckInfoFile($InfoFile, $opts{p});
	if($check=~/^Err/){													#InfoFile does not respect the format
		print "Program STOP - InfoFile does not respect the format\n";
		print $check;													#check correspond to error message
		exit;
	}
	else{
		@Files = GetFilesList( $opts{d});				#Get path of file present in input directory
		$series=$check;
		my @files; 										#Get file name
		foreach my $f (@Files){
			my @tmp=split('/',$f);
			push(@files,renam($tmp[$#tmp]));
		}
		foreach my $k (keys %$series){					#Check sample name in InfoFile are present in input directory
			foreach my $sample (@{$$series{$k}}){
				if($sample ne "NA"){
					if (!(grep {$_ eq $sample} @files)){
						print "Program STOP - InfoFile contains name of sample not present in the input collection of samples\n";
						print "File not found in the input collection : $sample\n";
						exit;
					}
				}
			}
		}
	}
}
else{
	$type="SIMPLE";
	@Files = GetAllFilesList( $opts{d} );	
	my @files;	
	foreach my $f (@Files){
		my @tmp=split('/',$f);
		push(@files,renam($tmp[$#tmp]));
	}
	for(my $i=0; $i<=$#files;$i++){
		$$series{$i}=$files[$i];
		push(@{$$cat{"All"}},$files[$i]);
	}
}
###################################################################################################



#####################################################################################
#	This part display some informations about the type of analysis
#	It also get the variant of the differents file present in the input directory
#
if($type eq "Other"){
	#print "You're working with several category\n";
}
elsif($type eq "Normal"){
	#print "You're working with only one category of sample\n";
}
elsif($type eq "Normal-Tumor"){
	#print "You're working with two category of sample : Normal-Tumor\n";
}
else{
	#print "You're working with three category of sample : Normal-Tumor-Duplicates\n";
}
my ($head, 	#Header line of the input VCF/tabular file
%h, 		#hashtable(K=path,V=chr|pos|ref|alt)
%hline);	#hastable(K=name.chr|pos|ref|alt, V=line)

#print "GETTING THE VARIANT OF FILE : \n";
foreach my $f (@Files){
	#print "\t";
	$head=GetVar($f,\%h, \%hline);
}

#############################################################
#	Check if input file format is the same for all files	#
#############################################################
my $officialformat=$filesformat[0];
foreach my $format (@filesformat){
	if($format ne $officialformat){
		print "Program STOP - You need to provide all input file in the SAME format\n";
		print "Please see the documentation\n";
		exit;
	}
}


#####################################################
#	Get all the variant one time only in a table	#
#####################################################
my @uniquevariant;				#Table for store all the variants in an uniqu way
my @uniquevariant2;				#Table for store all the variants in an uniqu way, the one with chr_random...
my %uniqueline;					#line choose from a random file for a variant
foreach my $name (keys %h){
	foreach my $var (@{$h{$name}}){
		$uniqueline{$var}=$hline{$name.$var};
		if (!(grep {$_ eq $var} @uniquevariant)){
			my @tmp=split('\|',$var);
			if($tmp[0]=~/chr\d+$/){
				push(@uniquevariant,$var);
			}
			else{
				push(@uniquevariant2,$var);
			}
		}
	}
}
my $nbsample=keys %h;	#Get the number of sample corresponding to number of input file

my @ordcat=sort (keys %$cat);				#Ordered the categories name

#####################################
#	Get the header of input files	#
#####################################
my $headInfoFile;
if($type ne "SIMPLE"){
	open(IN, $InfoFile) or die ("Unable to open the file containing name of input VCF/tabular variants files, option -i.\n");
	$headInfoFile=<IN>;		#Get header infoFile
	chomp($headInfoFile);
	close(IN);
}


if($opts{p} eq "pair"){
	my @field=split("\t",$headInfoFile);
	if($field[0]!~/^Normal$/){
		print "Err : You must provide an InfoFile with first column name \"Normal\" if you want to do a paired analysis\n";
		exit;
	}
	elsif($#field>2){
		print "Err : You can't have more than 3 column (Normal, Tumor and Duplicates) if ou want to do a paired analysis\n";
		exit;
	}
}

#########################
#	Apply HotSpot Tool	#
#########################

my ($count,$nom)=HotSpotv(\%h,$type);

#count is a HASH : K=chr|pos|ref|alt.categorie, V=number of time we see this variant in this categorie
#nom is a HASH : K=chr|pos|ref|alt.categorie, V=table containing the name of the sample file in this categorie where a variant was find



##################################################################""
#	Count variant frequency in the population
#		input :	1)hash variant
#				2)type of analysis
#				3)hash categorie
sub HotSpotv{
	my($h,$t)=@_;
	my(%count,%nom);
	if($t!~/^Normal-Tumor-Duplicates$/){				#Case several categories, or only normal, or only normal-tumor
		foreach my $cate (keys %$cat){						#Foreach categories
			foreach my $name (@{$$cat{$cate}}){				#Foreach file in this categorie
				if($name ne "NA"){
					foreach my $var (@{$h{$name}}){				#Foreach variant in this file of this categorie
						$count{$var.$cate}++;					#increment number of time you see this variant in this categorie
						push(@{$nom{$var.$cate}},$name);		#insert in a table the name of file where this variants was find in this categorie
					}
				}
			}
		}
	}
	else{
		if($opts{p} eq "pair"){												#Case Normal-Tumor-Duplicates
			my @header=split("\t",$headInfoFile);				#split the header of the InfoFile to have the categories name
			foreach my $cate (keys %$cat){						#Foreach categorie describe in InfoFile
				foreach my $name (@{$$cat{$cate}}){				#For each sample of this category
					if($cate eq $header[$#header]){				#If name of the categorie is duplicates
						$cate=$header[$#header-1];				#So change the categorie as tumor (to considers duplicates as tumor sample)
						foreach my $var (@{$h{$name}}){			#Foreach variant
							$count{$var.$cate}++;				#Count number of time we see this variants in this category
							push(@{$nom{$var.$cate}},$name);	#save the sample name of this categorie where the variants was find
						}
					}
					else{										#If not duplicates for this line
						foreach my $var (@{$h{$name}}){			#for each variant of the sample
							$count{$var.$cate}++;				#Count number of time we see this variants in this category
							push(@{$nom{$var.$cate}},$name);	#save the sample name of this categorie where the variants was find
						}
					}
				}
			}
		}
		else{											#Case catagories but they are name Normal Tumor...
			foreach my $cate (keys %$cat){						#Foreach categories
				foreach my $name (@{$$cat{$cate}}){				#Foreach file in this categorie
					if($name ne "NA"){
						foreach my $var (@{$h{$name}}){				#Foreach variant in this file of this categorie
							$count{$var.$cate}++;					#increment number of time you see this variant in this categorie
							push(@{$nom{$var.$cate}},$name);		#insert in a table the name of file where this variants was find in this categorie
						}
					}
				}
			}
		}
		
	}
	return(\%count,\%nom);			#return the 2 HASHTABLE
}

#########################
#	Header construction	#
#########################		
my $headfinal;							#if mutect2, remove # of header, to obtain homogeneity with MutSpec Annot
if($officialformat eq "mutect2"){
	$headfinal=$head;
	$headfinal=~s/^.//;
}
else{
	$headfinal=$head;
}
if($type!~/^Normal-Tumor-Duplicates$/){							#All the case without case with duplicates
	foreach my $c (@ordcat){										#Foreach ordered categries name
		$headfinal.="\t".$c."_count\t".$c."_freq\t".$c."_Sample";	#Add the future column in the header
	}
}
else{															#Case normal-tumor-duplicate
	foreach my $c (@ordcat){
		my @header=split("\t",$headInfoFile);
		if($c eq $header[$#header]){				#$header[$#header] eq "Duplicates" --> If duplicates next because considers as the tumor
			next;
		}
		$headfinal.="\t".$c."_count\t".$c."_freq\t".$c."_Sample";
	}
}

#########################
#	Output Directory	#
#########################
#print "Create Output Directory\n";
if(! -e $opts{s}){
	mkdir($opts{s}) or die("Erreur creation repertoire $opts{s}\n");
}
chdir($opts{s}) or die("Erreur chdir pour aller dans le répertoire $opts{s} \n");

#########################
#	Sample_Summary.vcf	#
#########################
#Foreach variants
#	Foreach categories
#		If variants in this category
#			1 Number of time the variant was see in this categorie
#			2 Number of time the variant was see in this categorie/Number of sample in this categorie
#			3 Name of the sample of the categorie where this variant was find
#		else
#			1 0
#			2 0%
#			3 NA
#	End
#End
#print "Writing variants_summary.vcf...\n";
open(OUT, ">"."variants_summary.vcf") or die ("Unable to open output writing file variants_summary.vcf");
print OUT $headfinal."\n";
my %varsampleannotation;						#HASH K=var, V=information of the frequency of this variants in the differents categories
my @uniquevariantsort= sort { return (split(/r/,(split('\|',$a))[0]))[1] <=> (split(/r/,(split('\|',$b))[0]))[1] || (split('\|',$b))[1] <=> (split('\|',$b))[1] } @uniquevariant;
push(@uniquevariantsort,@uniquevariant2);
foreach my $var (@uniquevariantsort){					#Foreach variants
	my $aecrire="";									#annotation on frequency
	if($type!~/^Normal-Tumor-Duplicates$/){			#Case without duplicates
		foreach my $c (@ordcat){					#Foreach categories
			if($#{$$nom{$var.$c}}>-1){					#If there is at least one sample of this categorie where a variant was find
				# $$count{$var.$c};						#Number of time the variant was see in this categorie
				# $$count{$var.$c}/$#{$$cat{$c}};		#Number of time the variant was see in this categorie/Number of sample in this categorie
				my $tot=0;
				foreach my $tp (@{$$cat{$c}}){
					if($tp ne "NA"){
						$tot++;
					}
				}
				my $perc=nearest(.001,($$count{$var.$c}/($tot))*100);
				my $concernsample;						#Name of the sample of the categorie where this variant was find
				foreach my $n (sort @{$$nom{$var.$c}}){
					$concernsample.=$n.',';				#Concatenate
				}
				chop($concernsample);
				$aecrire.="\t".$$count{$var.$c}."\t".$perc." %\t".$concernsample;
			}
			else{
				$aecrire.="\t0\t0 %\tNA";				#Else no sample of this categorie where the variant is find
			}
		}
	}
	else{											#Case Normal-tumor-duplicate
		my @header=split("\t",$headInfoFile);
		foreach my $c (@ordcat){						#Foreach categorie
			if($c eq $header[$#header]){				#If categorie duplicates, so pass car considers as categorie tumor
				next;
			}
			elsif($c eq $header[$#header-1]){			#Case categorie Tumor
				if($#{$$nom{$var.$c}}>-1){
					# $$count{$var.$c};						#Number of time the variant was see in this categorie
					# $$count{$var.$c}/$#{$$cat{$c}};		#Number of time the variant was see in this categorie/Number of sample in this categorie
					my $perc=nearest(.001,($$count{$var.$c}/((($#{$$cat{$c}})+1)+(($#{$$cat{$header[$#header]}})+1)))*100);
					my $concernsample;						#Name of the sample of the categorie where this variant was find
					foreach my $n (sort @{$$nom{$var.$c}}){
						$concernsample.=$n.',';				#Concatenate
					}
					chop($concernsample);
					$aecrire.="\t".$$count{$var.$c}."\t".$perc." %\t".$concernsample;
				}
				else{
					$aecrire.="\t0\t0 %\tNA";
				}
			}
			else{										#Case categorie Normal
				if($#{$$nom{$var.$c}}>-1){
					# $$count{$var.$c};						#Number of time the variant was see in this categorie
					# $$count{$var.$c}/$#{$$cat{$c}};		#Number of time the variant was see in this categorie/Number of sample in this categorie
					my $perc=nearest(.001,($$count{$var.$c}/(($#{$$cat{$c}})+1))*100);
					my $concernsample;						#name of the sample of the categorie where this variant was find
					foreach my $n (sort @{$$nom{$var.$c}}){
						$concernsample.=$n.',';				#Concatenate
					}
					chop($concernsample);
					$aecrire.="\t".$$count{$var.$c}."\t".$perc." %\t".$concernsample;
				}
				else{
					$aecrire.="\t0\t0 %\tNA";
				}
			}
		}
	}
	print OUT $uniqueline{$var}.$aecrire."\n";	#write in the file
	$varsampleannotation{$var}=$aecrire;		#store the new annotation on the variant frequency
}
close(OUT);

chdir("..") or die("Erreur chdir pour sortir du répertoire $opts{s}\n");
#################
#	Annotation	#
#################
#Foreach individual	$k key %series
#	Foreach sample of this individual	my for my $i=0;$i<=@{$$series{$k}};$i++
#		if $i==0 Case sample in normal categorie
#			Foreach variant of this normal sample of this individual	my $var @{$h{${$$series{$k}}[$i]}}
#				push(@{$var.$individu},normal);
#		if $i==1 #case tumeur
#			Foreach variant of this tumor sample of this individual		my $var @{$h{${$$series{$k}}[$i]}}
#				push(@{$var.$individu},tumeur);
#		if $i>1
#			Foreach variant of this duplicates sample of this individual	my $var @{$h{${$$series{$k}}[$i]}}
#				push(@{$var.$individu},dup);
#	end
#end
my %hannot;						#HASHTABLE K=integer identifier of individual.chr|pos|ref|alt , V=table storing the type of categorie sample
if($type ne "SIMPLE"){
	foreach my $k (keys %{$series}){						#Foreach individual
		my $i=0;
		for ($i=0;$i<=$#{$$series{$k}};$i++){			#Foreach sample of this individual
			if($i==0){							#if $i==0, case sample in normal categorie
				if(${$$series{$k}}[$i] ne "NA"){
					#print ${$$series{$k}}[$i]."\n";
					foreach my $var (@{$h{${$$series{$k}}[$i]}}){	#Foreach variant of this normal sample of this individual
						push(@{$hannot{$k.$var}},"normal");
					}
				}
				else{
					#print ${$$series{$k}}[$i]."\n";
					push(@{$hannot{$k}},"Nabsent");
				}
			}
			elsif($i==1){			
				if(${$$series{$k}}[$i] ne "NA"){				#if $i==1, case tumeur
					foreach my $var (@{$h{${$$series{$k}}[$i]}}){	#Foreach variant of this tumor sample of this individual
						push(@{$hannot{$k.$var}},"tumeur");
					}
				}
				else{
					push(@{$hannot{$k}},"Tabsent");
				}
			}
			elsif($i>1){							#if $i>1, case duplicates
				if(${$$series{$k}}[$i] ne "NA"){	
					foreach my $var (@{$h{${$$series{$k}}[$i]}}){	#Foreach variant of this duplicates sample of this individual
						push(@{$hannot{$k.$var}},"dup");
					}
				}
				else{										#Si pas d'echantillons pour cette category pour cet individu	
					push(@{$hannot{$k}},"Dabsent");
				}
			}	
		}
	}
}



## To test if it's fonctionnal


#########################
#	Output Directory	#
#########################
#print "Create Output Directory\n";
if(! -e $opts{o}){
	mkdir($opts{o}) or die("Erreur creation repertoire $opts{o}\n");
}
chdir($opts{o}) or die("Erreur chdir pour aller dans le répertoire $opts{o} \n");


#################
#	Collection	#
#################
#
#	Create the collection of the input file with the annotation on frequency
foreach my $n (keys %h){				#Foreach sample
	#print "\t".$n."_hotspot.vcf\n";
	open(OUT, ">".$n."_hotspot.vcf") or die ("Unable to create and open ".$n."_hotspot.vcf");
	if($type!~/^Normal-Tumor/){		#Case no duplicates
		print OUT $headfinal."\n";
		my @temp;					#variant that have a number for chromosome, not a letter
		my @temp2;					#variant that have a letter for chromosome, not a number
		foreach my $tmp (@{$h{$n}}){
			my @tp=split('\|',$tmp);
			if($tp[0]=~/chr\d+$/){
				push(@temp,$tmp);
			}
			else{
				push(@temp2,$tmp);
			}
		}							#It permit to class them by chromosome number without having warnings
		my @sortvar= sort { return (split(/r/,(split('\|',$a))[0]))[1] <=> (split(/r/,(split('\|',$b))[0]))[1] || (split('\|',$b))[1] <=> (split('\|',$b))[1] } @temp;
		push(@sortvar,@temp2);
		foreach my $var (@sortvar){
			print OUT $hline{$n.$var}.$varsampleannotation{$var}."\n";
		}
	}
	else{
		print OUT $headfinal."\tAnnotation\tStatus\n";
		my @temp;					#variant that have a number for chromosome, not a letter
		my @temp2;					#variant that have a letter for chromosome, not a number
		foreach my $tmp (@{$h{$n}}){
			my @tp=split('\|',$tmp);
			if($tp[0]=~/chr\d+$/){
				push(@temp,$tmp);
			}
			else{
				push(@temp2,$tmp);
			}
		}							#It permit to class them by chromosome number without having warnings
		my @sortvar= sort { return (split(/r/,(split('\|',$a))[0]))[1] <=> (split(/r/,(split('\|',$b))[0]))[1] || (split('\|',$b))[1] <=> (split('\|',$b))[1] } @temp;
		push(@sortvar,@temp2);
		foreach my $var (@sortvar){
			my $individu="";
			my $annotation="";
			my $confirmation="";
			
			my $normalsample="";
			my $tumorsample="";
			my $dupsample="";
			
			foreach my $tmp (keys %$series){		#Foreach individual
				if ((grep {$_ eq $n} @{$$series{$tmp}})){	#Search individual where come from the file
					$individu=$tmp;
				}
			}
			
			if((grep {$_=~/Nabsent/} @{$hannot{$individu}})){	
				$normalsample="false";
			}
			else{
				$normalsample="true";
			}
			if((grep {$_=~/Tabsent/} @{$hannot{$individu}})){	
				$tumorsample="false";
			}
			else{
				$tumorsample="true";
			}
			if((grep {$_=~/Tabsent/} @{$hannot{$individu}})){	
				$dupsample="false";
			}
			else{
				$dupsample="true";
			}
			
			
			
			
			if (!(grep {$_=~/normal/} @{$hannot{$individu.$var}})){		#Not find in the normal file
				if ((grep {$_=~/tumeur/} @{$hannot{$individu.$var}})){		#Find in the tumor file
					if ((grep {$_=~/dup/} @{$hannot{$individu.$var}})){		#Find in the duplicate
						$annotation="Somatic";
						if($normalsample eq "true"){
							$confirmation="confirmed";
						}
						else{
							$confirmation="NA";
						}
						
					}
					else{													#Not find in the duplicates
						$annotation="Somatic";
						if($normalsample eq "true"){
							$confirmation="not confirmed";
						}
						else{
							$confirmation="NA";
						}
					}
				}
				else{
					if ((grep {$_=~/dup/} @{$hannot{$individu.$var}})){		#Find in the duplicates but not in tumor sample
						$annotation="Somatic";
						if($normalsample eq "true"){
							$confirmation="not confirmed";
						}
						else{
							$confirmation="NA";
						}
					}
				}
			}
			else{					#Find in the normal file
				$annotation="Germline";			
				$confirmation="";
			}
			
			
			
			
			
			
			
			
			
			print OUT $hline{$n.$var}.$varsampleannotation{$var}."\t".$annotation."\t".$confirmation."\n"; #write for each file the annotation 
		}
	}
	
	close(OUT);
}

chdir("..") or die("Erreur chdir pour sortir du répertoire $opts{o}\n");




########################################################################################################################
#	Function that check the InfoFiles
#	Input : Name of the Infofile
#	Output : 	If error, the Error message that begin with "Err"
#				If correct, a hash table containing the information of Infofile
#					1) Case category Other, keys=integer, value=ref on table containing sample name
#					2) Case category "Normal", keys=integer, value=sample name
#					3) Case category "Normal-Tumor", keys=integer, value=ref on table containing sample name
#					4) Case category "Normal-Tumor-Duplicate", keys=integer, value=ref on table containing sample name
#
sub CheckInfoFile{
	my $file=shift;
	my $p=shift;
	my $error_message="";	#Error message return by the function
	my %series;				#The series of sample for individual n
	my %category;			#The category of each input file
	open(IN, $file) or die ("Unable to open the file containing name of input VCF/tabular variants files, option -i.\n");
	my $head= <IN>;
	chomp($head);
	my @col=split("\t",$head);		#column if the header
	my $nb=0;						#Number of line, corresponding to number of individuals
	my @categories=split("\t",$head);	#Several categories presents in this file
	if($col[0]=~/^Normal$/i && $p eq "pair"){			#Case Normal-Tumor-Duplicates
		if($#col eq 0){					#Case normal
			$type="Normal";						#Not differenciated serie
			while ( defined( my $line = <IN> ) ) {
				chomp($line);
				my @infos = split("\t",$line);
				if($#col<$#infos){		#More column in the differents lines than in the header, impossible because only one category, so not any line with 2 sample or more
					$type="error";
					$error_message="Err : You can not have lines with more column than the header, if not the name of the sample will not be assign to a category\n";
				}
				else{
					$nb++;
					if(defined($infos[0]) && $infos[0] ne ""){
						$series{$nb}=renamtab(@infos);		#Get the name of input sample for the different individuals (each line = one individual)
					}
					else{
						$type="error";
						$error_message="Err : You can not leave a blank in your InfoFile. If sample missing, please name it \"NA\"\n";
						last;
					}
					for(my $i=0; $i<=$#categories; $i++){
						if($infos[$i] ne "NA"){
							push(@{$category{$categories[$i]}},renam($infos[$i]));	#This file is associated with this categorie
						}
					}
				}
			}
		}
		elsif($#col eq 1){				#Case normal-tumor
			$type="Normal-Tumor";		#differenciated serie
			while ( defined( my $line = <IN> ) ) {
				chomp($line);
				my @infos = split("\t",$line);
				if($#col<$#infos){		#More column in the differents lines than in the header, impossible because only two category, so not any line with 3 sample or more
					$type="error";
					$error_message="Err : You can not have lines with more column than the header, if not the name of the sample will not be assign to a category\n";
				}
				else{
					$nb++;
					if(defined($infos[0]) && $infos[0] ne "" && defined($infos[1]) && $infos[1] ne ""){
						$series{$nb}=renamtab(@infos);		#Get the name of input sample for the different individuals (each line = one individual)
					}
					else{
						$type="error";
						$error_message="Err : You can not leave a blank in your InfoFile. If sample missing, please name it \"NA\"\n";
						last;
					}
					for(my $i=0; $i<=$#categories; $i++){
						if($infos[$i] ne "NA"){
							push(@{$category{$categories[$i]}},renam($infos[$i]));	#This file is associated with this categorie
						}
					}
				}
			}
		}
		elsif($#col eq 2){						#Case normal-tumor-duplicate
			$type="Normal-Tumor-Duplicates";	#differenciated serie + duplicates
			while ( defined( my $line = <IN> ) ) {
				chomp($line);
				my @infos = split("\t",$line);
				$nb++;
				foreach my $tmp (@infos){
					if(!(defined($tmp)) or $tmp eq ""){
						$type="error";
						$error_message="Err : You can not leave a blank in your InfoFile. If sample missing, please name it \"NA\"\n";
						last;
					}
				}
				if($type eq "error"){
					last;
				}
				if(defined($infos[0]) && $infos[0] ne "" && defined($infos[1]) && $infos[1] ne ""){
						$series{$nb}=renamtab(@infos);		#Get the name of input sample for the different individuals (each line = one individual)
				}
				else{
					$type="error";
					$error_message="Err : You can not leave a blank in your InfoFile. If sample missing, please name it \"NA\"\n";
					last;
				}
				for(my $i=0; $i<=$#categories; $i++){
					if(defined($infos[$i]) && $infos[$i] ne "NA"){			#If duplicate, put it in the same category as tumor
						push(@{$category{$categories[$i]}},renam($infos[$i]));		#This file is associated with this categorie
						#print $categories[$i]."\t".renam($infos[$i])."\n";
					}
				}
				if($#infos>$#categories){			#Case several duplicates
					for(my $i=$#categories+1; $i<=$#infos; $i++){
						#print $categories[$#categories]."\t".renam($infos[$i])."\n";
						push(@{$category{$categories[$#categories]}},renam($infos[$i]));
					}
				}
			}
		}
		else{
			$type="error";
			$error_message="Err : If you have an header with first column \"Normal\", you can provide 3 column maximum (Normal-Tumor-Duplicates).\n";
			$error_message.="Note : It concern the header, so you can have more than 3 column after the header, like the case of several duplicates for one individual, but must be tabular separate.\n";
		}
	}
	else{						#Case with several column for several treatment
		$type="Other";
		while ( defined( my $line = <IN> ) ) {
			my @tabtmp=();
			chomp($line);
			my @infos = split("\t",$line);
			
			
#			if($#col<$#infos){		#More column in the differents lines than in the header, impossible because defined number of categories, so not any line with pmore sample than categories in the header
#				$type="error";
#				$error_message="Err : You can not have lines with more column than the header, if not the name of the sample will not be assign to a category\n";
#				last;
#			}
#			else{
#				$nb++;
#				if($#infos==$#categories){
#					foreach my $tmp (@infos){
#						if((!(defined($tmp))) or $tmp eq ""){
#							$type="error";
#							$error_message="Err : You can not leave a blank in your InfoFile. If sample missing, please name it \"NA\"\n";
#							last;
#						}
#					}
#				}
#				else{
#					$type="error";
#					$error_message="Err : You can not leave a blank in your InfoFile. If sample missing, please name it \"NA\"\n";
#					last;
#				}
#				if($type eq "error"){
#					last;
#				}
#				$series{$nb}=renamtab(@infos);		#Get the name of input samples for the different individuals (each line = one individual)
#				for(my $i=0; $i<=$#categories; $i++){
#					if($infos[$i] ne "NA"){
#						push(@{$category{$categories[$i]}},renam($infos[$i]));	#This file is associated with this category
#					}
#				}
#			}
			
			
			
			
			
			
			if($#categories<$#infos){		#More column in the differents lines than in the header, impossible because defined number of categories, so not any line with pmore sample than categories in the header
				$type="error";
				$error_message="Err : You can not have lines with more column than the header, if not the name of the sample will not be assign to a category\n";
				last;
			}
			elsif($#infos==$#categories){		#Same amount of column between header and the treated line
				$nb++;
				for(my $i=0; $i<=$#categories; $i++){
					if(defined($infos[$i]) && $infos[$i] ne ""){
						push(@{$category{$categories[$i]}},renam($infos[$i]));	#This file is associated with this category
						push(@tabtmp,$infos[$i]);
					}
					else{
						push(@{$category{$categories[$i]}},"NA");	#This file is associated with this category
						push(@tabtmp,"NA");
					}
				}
			}
			else{
				$nb++;
				$series{$nb}=renamtab(@infos);
				for(my $i=0; $i<=$#categories; $i++){
					if(defined($infos[$i]) && $infos[$i] ne ""){
						push(@{$category{$categories[$i]}},renam($infos[$i]));	#This file is associated with this category
						push(@tabtmp,$infos[$i]);
					}
					else{
						push(@{$category{$categories[$i]}},"NA");	#This file is associated with this category
						push(@tabtmp,"NA");
					}
				}
			}
			$series{$nb}=renamtab(@tabtmp);
		}
	}
	close(IN);
	my @uniquesample=();
	if($type ne "Other"){
	#Check not 2 time the same sample name in InfoFile
		foreach my $k (keys %series){
			foreach my $name (@{$series{$k}}){
				if($name ne "NA"){
					if (!(grep {$_ eq $name} @uniquesample)){
						push(@uniquesample,$name);
					}
					else{
						$type="error";
						$error_message="Err : You have a sample name present twice or more in your InfoFile but you are in a paired analysis\n";
						$error_message.="Name : $name \n";
						last;
					}
				}
			}
			if($type eq "error"){
				last;
			}
		}	
	}
	else{
		foreach my $k (keys %series){
			foreach my $name (@{$series{$k}}){
				if($name ne "NA"){
					if (!(grep {$_ eq $name} @uniquesample)){
						push(@uniquesample,$name);
					}
				}
			}
		}	
	}
	if($type eq "Other"){
		foreach my $cat (keys %category){
			my @uniquesampletmp=();
			
			#print $cat."\n";
			foreach my $name (@{$category{$cat}}){
				#print "\t".$name."\n";
				
				if($name ne "NA"){
					if (!(grep {$_ eq $name} @uniquesampletmp)){
						push(@uniquesampletmp,$name);
					}
					else{
						$type="error";
						$error_message="Err : You have a sample name present twice or more in the same category in your InfoFile\n";
						$error_message.="Category : $cat ; Name : $name \n";
						last;
					}
				}
			}
			if($type eq "error"){
				last;
			}
			print "\n";
		}
	}
	
	if($type eq "error"){
		return($error_message, '0');
	}
	else{
		return(\%series,\%category, \@uniquesample);
	}
}
########################################################################################################################



########################################
# Input : directory path
# Output : Table of file path @FilesList
sub GetFilesList {
	##
	##	This script get all file name of a specify directory and return an array of this file list
	##
	my $Path = $_[0];
	my $FileFound;
	my @FilesList = ();

	# Read file list
	opendir( my $FhRep, $Path )
	  or die "Impossible d'ouvrir le repertoire $Path\n";
	my @Contenu = grep { !/^\.\.?$/ } readdir($FhRep);
	closedir($FhRep);

	foreach my $FileFound (@Contenu) {			#Foreach file found
		my $interest="false";					#Are we interested in this file?
		my $renamfile=renam($FileFound);
		foreach my $k (keys %$cat){					#Foreach categorie
			foreach my $interestFile (@{$$cat{$k}}){	#Foreach file of this categorie
				if($interestFile eq $renamfile){		#If the file found match with one of the file store in a category
					$interest="true";					#So yes, we are interested in
				}
			}
		}
		if($interest eq "true"){
			# Treatment of file
			if ( -f "$Path/$FileFound" ) {
				push( @FilesList, "$Path/$FileFound" );
			}
	
			# Treatment of directory
			elsif ( -d "$Path/$FileFound" ) {
	
				# Recursive reshearch
				push( @FilesList, GetFilesList("$Path/$FileFound") );
			}
		}

	}
	return @FilesList;
}
########################################

##############################################################################################################
#	Function the Get variant from a VCF/tabular file
#		Input : Path of the file; hashtable(K=path,V=chr|pos|ref|alt); hastable2(K=name.chr|pos|ref|alt, V=line)
#		Output : the header line of the file
#	Check if the format is respected
#
sub GetVar{
	my ($File,$h,$hline) =@_;
	my ($header);
	my $name=renam($File);
	#print "$name - Format : ";
	open( IN, $File ) or die("Fail to open in txt file $File");
	my @data;
	my $vcfheader="";			#Correspond to the header if vcf input file
	my $som =-1;				#Correspond to column judgment of somatic or not by the variant caller in case of mutect and mutect 2
	my $first = -1;
	my ( $chr)="";
	my ($pos)="";
	my ($ref)="";
	my ($alt)="";
	while ( defined( my $line = <IN> ) ) {
		chomp($line);
		if ( $line !~ /^##/ ) {
			@data = split( '\t', $line );
			if ( $first < 0 ) {
				$first = 1;
				$header=$line;
				chomp($header);
				foreach my $k (keys %acceptedformat){
					for ( my $i = 0 ; $i <= $#data ; $i++ ) {
						if($data[$i]=~/^${$acceptedformat{$k}}[0]$/){
							$chr=$i;
						}
						elsif($data[$i]=~/^${$acceptedformat{$k}}[1]$/){
							$pos=$i;
						}
						elsif($data[$i]=~/^${$acceptedformat{$k}}[2]$/){
							$ref=$i;
						}
						elsif($data[$i]=~/^${$acceptedformat{$k}}[3]$/){
							$alt=$i;
						}
						if($vcfheader=~/VCFv4.[1|2]/ && $vcfheader=~/MuTect2/){
							if($data[$i]=~/^FILTER$/){
								$som=$i;
							}
						}
						elsif($vcfheader=~/muTector/){
							if($data[$i]=~/judgement/){
								$som=$i;
							}
						}
					}
					if($chr ne "" && $pos ne "" && $ref ne "" && $alt ne ""){
						#print $k."\n";
						$fileformat=$k;
						push(@filesformat,$fileformat);
						last;
					}
				}
				if(!(defined($fileformat))){
					print "Err : Incorrect input file format. Please refer to the documentation\n";
					exit;
				}
			}			
			else {
				 if($data[$chr]!~/^chr/i){
				 	$data[$chr]="chr".$data[$chr];
				 }
				 my $var =
					    $data[$chr] . "|"
					  . $data[$pos] . "|"
					  . $data[$ref] . "|"
					  . $data[$alt];
				if($som>-1){
					if($data[$som] eq "PASS" || $data[$som] eq "KEEP" || $data[$som] eq "Somatic"){
						push( @{ $$h{ $name } },         $var );
						$$hline{ $name . $var }= $line ;
					}
				}
				else{
					push( @{ $$h{ $name } },         $var );
					$$hline{ $name . $var }= $line ;
				}
			}
		}
		else{				#Take in account header if vcf file
			$vcfheader.=$line;
		}
	}
	close(IN);
	return($header);
}

##############################################################################################################
sub renam{
	#Input : sample name | Output : sample name reformat
	my $name=shift;
	my $newname="";
	if($name=~/\/+/){					#if path
		my @names=split('\/',$name);
		my @n=split('\.',$names[$#names]);
		if($#n>0){
			for(my $i=0; $i<$#n; $i++){
				$newname.=$n[$i]."_";
			}
			chop($newname);
		}
		else{
			$newname=$n[0];
		}
	}
	else{							#if name
		my @n=split('\.',$name);	#remove extension
		if($#n>0){
			for(my $i=0; $i<$#n; $i++){
				$newname.=$n[$i]."_";
			}
			chop($newname);
		}
		else{
			$newname=$n[0];
		}
	}
	return ($newname);
}

##############################################################################################################
sub renamtab{
	#Input : table of sample name | Output : table of sample name reformat
	my @tab=@_;
	my @tabnewname=();
	foreach my $name (@tab){
		my $newname="";
		if($name=~/\/+/){					#if path
			my @names=split('\/',$name);
			my @n=split('\.',$names[$#names]);
			if($#n>0){
				for(my $i=0; $i<$#n; $i++){
					$newname.=$n[$i]."_";
				}
				chop($newname);
			}
			else{
				$newname=$n[0];
			}
		}
		else{							#if name
			my @n=split('\.',$name);	#remove extension
			if($#n>0){
				for(my $i=0; $i<$#n; $i++){
					$newname.=$n[$i]."_";
				}
				chop($newname);
			}
			else{
				$newname=$n[0];
			}
		}
		push(@tabnewname,$newname);
	}
	return(\@tabnewname);
}

sub GetAllFilesList {
	##
	##	This script get all file name of a specify directory and return an array of this file list
	##
	my $Path = $_[0];
	my $FileFound;
	my @FilesList = ();

	# Lecture de la liste des fichiers
	opendir( my $FhRep, $Path )
	  or die "Impossible d'ouvrir le repertoire $Path\n";
	my @Contenu = grep { !/^\.\.?$/ } readdir($FhRep);
	closedir($FhRep);

	foreach my $FileFound (@Contenu) {

		# Traitement des fichiers
		# Traitement des fichiers
		if ( -f "$Path/$FileFound" ) {
			push( @FilesList, "$Path/$FileFound" );
		}

		# Traitement des repertoires
		elsif ( -d "$Path/$FileFound" ) {

			# Boucle pour lancer la recherche en mode recursif
			push( @FilesList, GetFilesList("$Path/$FileFound") );
		}

	}
	return @FilesList;
}