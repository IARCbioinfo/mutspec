# !/usr/bin/perl

#-----------------------------------#
# Author: Vincent                   #
# Script: mutspecSplit.pl           #
# Last update: 01/07/14             #
#-----------------------------------#


use strict;
use warnings;
use Getopt::Long;

our $file="";
our $column="";
our $path="";
our $key="";


GetOptions('file|f=s'		=>\$file,
		   'key|k=s'      	=>\$key,
		   'column|i=s'		=>\$column,
		   'path|p=s'     	=>\$path);


mkdir ("outputFiles") or die ("Erreur creation repertoire\n");
# print $file,"\n", $key,"\n", $column,"\n", $path,"\n"; exit;

my %tab;
if ($column==0) {$column++;}
$column--;

open(FILE, "$file") or die "cannot open $file\n";

$_=<FILE>; #skip headers
chomp;
my @line = split(/\t/,$_);
my $headers = join("\t", @line[0..($column-1),($column+1)..$#line]);

while(<FILE>){
	chomp;
	my @line = split(/\t/,$_);
	#if (!exists($tab{$line[$column]})) { $tab{$line[$column]}=[]; }
	#push( @{ $tab{$line[$column]} }, join("\t", @line[0..($column-1),($column+1)..$#line]) );
	my $tmp = join("\t", @line[0..($column-1),($column+1)..$#line]) ;
	my $id = $line[$column];
	push( @{ $tab{$id} }, $tmp);
}


while( my ($name,$lines) = each(%tab) ) {
	my $output="outputFiles/$name";
	#my $output="primary_$key" . "_$name" . "_visible_tabular";
	# my $output=$name;
	open(FILE, ">$output") or die "cannot create file $output \n";
	print FILE $headers."\n";
	foreach my $line (@{$lines}){
		print FILE "$line\n";
	}
	close FILE;
}

my $list=`ls outputFiles/*`;
print ($list);
