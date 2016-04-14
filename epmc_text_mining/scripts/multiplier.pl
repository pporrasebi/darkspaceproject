#!/usr/bin/perl -w
# Program to create a single line for each duplicated id given in fields separated by ,
# v1.0, 06.06.2016

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

my @up_ac_1 = ();
my @up_ac_2 = ();

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
	
	if ( ($geneid_1[$i] =~ ",") && ($geneid_2[$i] =~ ",")) {
		@up_ac_1 = split(/,/, $geneid_1[$i]);
		@up_ac_2 = split(/,/, $geneid_2[$i]);
		for my $j ( 0 .. $#up_ac_1) {
			for my $k ( 0 .. $#up_ac_2) {
				print(OUTFILE "$up_ac_1[$j]\t$up_ac_2[$k]\t$genename_1[$i]\t$genename_2[$i]\t$nof_docs[$i]\t$nof_co_occr_in_title_abs[$i]\t$nof_co_occr_in_body[$i]\t$pmcids[$i]\n");
			}
		}
	}
	elsif ($geneid_1[$i] =~ ",") {
		@up_ac_1 = split(/,/, $geneid_1[$i]);
		for my $j ( 0 .. $#up_ac_1) {
			print(OUTFILE "$up_ac_1[$j]\t$geneid_2[$i]\t$genename_1[$i]\t$genename_2[$i]\t$nof_docs[$i]\t$nof_co_occr_in_title_abs[$i]\t$nof_co_occr_in_body[$i]\t$pmcids[$i]\n");
		}
	}
	elsif ($geneid_2[$i] =~ ",") {
		@up_ac_2 = split(/,/, $geneid_2[$i]);
		for my $j ( 0 .. $#up_ac_2) {
			print(OUTFILE "$geneid_1[$i]\t$up_ac_2[$j]\t$genename_1[$i]\t$genename_2[$i]\t$nof_docs[$i]\t$nof_co_occr_in_title_abs[$i]\t$nof_co_occr_in_body[$i]\t$pmcids[$i]\n");
		}
	}
	else {
		print(OUTFILE "$geneid_1[$i]\t$geneid_2[$i]\t$genename_1[$i]\t$genename_2[$i]\t$nof_docs[$i]\t$nof_co_occr_in_title_abs[$i]\t$nof_co_occr_in_body[$i]\t$pmcids[$i]\n");
	}
}

print("Done!\n");

close OUTFILE;
exit;