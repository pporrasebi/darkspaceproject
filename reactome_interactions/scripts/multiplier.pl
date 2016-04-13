#!/usr/bin/perl -w
# Program to create a single line for each duplicated id given in fields separated by ;
# v1.0, 27.03.2014

use strict; 
use warnings;

# First, let's declare initial variables

my $infile = $ARGV[0];
my $outfile = $ARGV[1];
my @lines = ();

my @id = ();
my @redup_ac = ();
my @taxid = ();
my @up_ac = ();


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
	( $id[$i], $redup_ac[$i], $taxid[$i] ) = split( /\t/, $lines[ $i ] );
	@up_ac = split(/,/, $redup_ac[$i]);
	for my $j ( 0 .. $#up_ac) {
	print(OUTFILE "$up_ac[$j]\t$taxid[$i]\n");
	}
}

print("Done!\n");

close OUTFILE;
exit;