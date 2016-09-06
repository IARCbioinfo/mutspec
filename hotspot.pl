#/usr/bin/perl

# ~~~~ HOTSPOT tool ~~~~
# 13/06/2016
# Alexis ROBITAILLE
# robitaillea@students.iarc.fr
# Version : 1.0

#################
# #
# HOTSPOT #
# #
#################

# - The goal of this tool is to compute variant frequency in a define datasets
# - If information provided on duplicate sample, infers somatics status of the variants in this dataset

########################################################################################################################################
#Input : 1) Directory of VCF/tabular variants files (same format)
# 2) Tabular file containing the name of the input files
# The variants file for a same individual mlust be in the same line
# 2a) One column : not differenciated files (ex : all tumor sample)
# 2b) Two column : differenciated files (ex : normal-tumor sample, same individual)
# 2c) Three or more column : differenciated files + duplicates (ex : normal-tumor-duplicates, same individual)
#
#Output :1) Variants_Summary.vcf = One uniqu variant per line + annotation on the frequence of this variant in this dataset
# If input 2a : frequency in all the dataset
# If input 2b : frequency in the two separate series
#  2) Collection of the same input 2) files + annotation as describe below (frequency of variant in the dataset of differenciated or not series)
# If input 2c : Add annotation on somatic status of the variants
########################################################################################################################################

# Library
use strict;
use warnings;
use Getopt::Std;
use Math::Round ':all';

# Options
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


my $fileformat; #Format of the input VCF/tabular variants files
my @Files; #Table of path of the input VCF/tabular variants files
my $InfoFile; #Path of tabular file containing name of input VCF/tabular variants files
my $check;  #Provide the information of why the infofile is wrong, or in case of rigth, it's a hash of the series contain in the file --> become $serie variable
my $cat; #Ref to hashtable K=categorie, V=ref of table containing sample name that are in this categorie
my $uniqs; #Ref to table containing the name of the sample in InfoFile (uniqu)
my @filesformat; #Table of format of vcf input files
my $series; #Ref to hashtable K=integer, V=ref of table containing sample name (string sample name in case of normal only)

######################
## Accepted format ## of the input VCF/tabular file in the directory opts{d}
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


$InfoFile=$opts{i}; #Get the path to InfoFile

if($InfoFile eq "None"){
if($opts{p} eq "pair"){
print "Program STOP - You have to provide an InfoFile if you want to do a paired analysis\n";
exit;
}
}



###################################################################################################
# This part concern the verification of the integrity of the InfoFile, and the fact that it fit with the sample file name present in input directory
#
my $valid=1
my $type=""; #Type Normal-Tumor-Dup or Traitement
if($InfoFile ne "None"){
($check,$cat,$uniqs)=CheckInfoFile($InfoFile, $opts{p});
if($check=~/^Err/){ #InfoFile does not respect the format
print "Program STOP - InfoFile does not respect the format\n";
print $check; #check correspond to error message
exit;
}
else{
@Files = GetFilesList( $opts{d}); #Get path of file present in input directory
$series=$check;
my @files;  #Get file name
foreach my $f (@Files){
my @tmp=split('/',$f);
push(@files,renam($tmp[$#tmp]));
}
foreach my $k (keys %$series){ #Check sample name in InfoFile are present in input directory
foreach my $sample (@{$$series{$k}}){
if($sample ne "NA"){
if (!(grep {$_ eq $sample} @files)){
print "Program STOP - InfoFile contains name of sample not present in the input collection of samples\n";
print "File not found in the input collection : $sample\n";
$valid=0;
}
}
}
if (! $valid) {exit;}
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
# This part display some informations about the type of analysis
# It also get the variant of the differents file present in the input directory
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
my ($head,  #Header line of the input VCF/tabular file
%h,  #hashtable(K=path,V=chr|pos|ref|alt)
%hline); #hastable(K=name.chr|pos|ref|alt, V=line)

#print "GETTING THE VARIANT OF FILE : \n";
foreach my $f (@Files){
#print "\t";
$head=GetVar($f,\%h, \%hline);
}

#############################################################
# Check if input file format is the same for all files #
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
# Get all the variant one time only in a table #
#####################################################
my @uniquevariant; #Table for store all the variants in an uniqu way
my @uniquevariant2; #Table for store all the variants in an uniqu way, the one with chr_random...
my %uniqueline; #line choose from a random file for a variant
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
my $nbsample=keys %h; #Get the number of sample corresponding to number of input file

my @ordcat=sort (keys %$cat); #Ordered the categories name

#####################################
# Get the header of input files #
#####################################
my $headInfoFile;
if($type ne "SIMPLE"){
open(IN, $InfoFile) or die ("Unable to open the file containing name of input VCF/tabular variants files, option -i.\n");
$headInfoFile=<IN>; #Get header infoFile
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
# Apply HotSpot Tool #
#########################

my ($count,$nom)=HotSpotv(\%h,$type);

#count is a HASH : K=chr|pos|ref|alt.categorie, V=number 

...some text eliminated here because file size is larger than maximum viewing size of 32.0 KB...

name of the sample will not be assign to a category\n";
# last;
# }
# else{
# $nb++;
# if($#infos==$#categories){
# foreach my $tmp (@infos){
# if((!(defined($tmp))) or $tmp eq ""){
# $type="error";
# $error_message="Err : You can not leave a blank in your InfoFile. If sample missing, please name it \"NA\"\n";
# last;
# }
# }
# }
# else{
# $type="error";
# $error_message="Err : You can not leave a blank in your InfoFile. If sample missing, please name it \"NA\"\n";
# last;
# }
# if($type eq "error"){
# last;
# }
# $series{$nb}=renamtab(@infos); #Get the name of input samples for the different individuals (each line = one individual)
# for(my $i=0; $i<=$#categories; $i++){
# if($infos[$i] ne "NA"){
# push(@{$category{$categories[$i]}},renam($infos[$i])); #This file is associated with this category
# }
# }
# }






if($#categories<$#infos){ #More column in the differents lines than in the header, impossible because defined number of categories, so not any line with pmore sample than categories in the header
$type="error";
$error_message="Err : You can not have lines with more column than the header, if not the name of the sample will not be assign to a category\n";
last;
}
elsif($#infos==$#categories){ #Same amount of column between header and the treated line
$nb++;
for(my $i=0; $i<=$#categories; $i++){
if(defined($infos[$i]) && $infos[$i] ne ""){
push(@{$category{$categories[$i]}},renam($infos[$i])); #This file is associated with this category
push(@tabtmp,$infos[$i]);
}
else{
push(@{$category{$categories[$i]}},"NA"); #This file is associated with this category
push(@tabtmp,"NA");
}
}
}
else{
$nb++;
$series{$nb}=renamtab(@infos);
for(my $i=0; $i<=$#categories; $i++){
if(defined($infos[$i]) && $infos[$i] ne ""){
push(@{$category{$categories[$i]}},renam($infos[$i])); #This file is associated with this category
push(@tabtmp,$infos[$i]);
}
else{
push(@{$category{$categories[$i]}},"NA"); #This file is associated with this category
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
## This script get all file name of a specify directory and return an array of this file list
##
my $Path = $_[0];
my $FileFound;
my @FilesList = ();

# Read file list
opendir( my $FhRep, $Path )
  or die "Impossible d'ouvrir le repertoire $Path\n";
my @Contenu = grep { !/^\.\.?$/ } readdir($FhRep);
closedir($FhRep);

foreach my $FileFound (@Contenu) { #Foreach file found
my $interest="false"; #Are we interested in this file?
my $renamfile=renam($FileFound);
foreach my $k (keys %$cat){ #Foreach categorie
foreach my $interestFile (@{$$cat{$k}}){ #Foreach file of this categorie
if($interestFile eq $renamfile){ #If the file found match with one of the file store in a category
$interest="true"; #So yes, we are interested in
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
# Function the Get variant from a VCF/tabular file
# Input : Path of the file; hashtable(K=path,V=chr|pos|ref|alt); hastable2(K=name.chr|pos|ref|alt, V=line)
# Output : the header line of the file
# Check if the format is respected
#
sub GetVar{
my ($File,$h,$hline) =@_;
my ($header);
my $name=renam($File);
#print "$name - Format : ";
open( IN, $File ) or die("Fail to open in txt file $File");
my @data;
my $vcfheader=""; #Correspond to the header if vcf input file
my $som =-1; #Correspond to column judgment of somatic or not by the variant caller in case of mutect and mutect 2
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
else{ #Take in account header if vcf file
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
if($name=~/\/+/){ #if path
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
else{ #if name
my @n=split('\.',$name); #remove extension
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
if($name=~/\/+/){ #if path
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
else{ #if name
my @n=split('\.',$name); #remove extension
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
## This script get all file name of a specify directory and return an array of this file list
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