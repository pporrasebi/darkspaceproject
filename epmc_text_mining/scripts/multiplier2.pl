#!/usr/bin/perl -w
# Program to create a single line for each duplicated id given in fields separated by ;
# v1.0, 27.03.2014

use strict; 
use warnings;

# First, let's declare initial variables

my $infile = $ARGV[0];
my $outfile = $ARGV[1];
my @lines = ();

my @geneid_1=();
my @geneid_2=();
my @genename_1=();
my @genename_2=();
my @nof_docs=();
my @nof_co_occr_in_title_abs=();
my @nof_co_occr_in_body=();
my @pmcids=();

my @pmids = ();

# Start the program

print("Thinking...\n");

# Read in the data

open( INFILE, "<$infile" ) or die( "Couldn't open $infile: $!\n" );
@lines = <INFILE>;
close( INFILE );

open( OUTFILE, ">$outfile" ) or die( "Couldn't open $outfile: $!\n" );

chomp @lines;

# Prepare the lines to be readable 
# (split the columns). 

for my $i ( 0 .. $#lines ) {
	( $geneid_1[$i],$geneid_2[$i],$genename_1[$i],$genename_2[$i],$nof_docs[$i],$nof_co_occr_in_title_abs[$i],$nof_co_occr_in_body[$i],$pmcids[$i] ) = split( /\t/, $lines[$i] );
	
	if ($pmcids[$i] =~ "###") {
	
		@pmids = split(/###/, $pmcids[$i]);
		for my $j ( 0 .. $#pmids) {
			print(OUTFILE "$geneid_1[$i]\t$geneid_2[$i]\t$genename_1[$i]\t$genename_2[$i]\t$nof_docs[$i]\t$nof_co_occr_in_title_abs[$i]\t$nof_co_occr_in_body[$i]\t$pmids[$j]\n");
		}
	}
	
	else {
		print(OUTFILE "$geneid_1[$i]\t$geneid_2[$i]\t$genename_1[$i]\t$genename_2[$i]\t$nof_docs[$i]\t$nof_co_occr_in_title_abs[$i]\t$nof_co_occr_in_body[$i]\t$pmcids[$i]\n");
	}
}

print("Done!\n");

close OUTFILE;
exit;