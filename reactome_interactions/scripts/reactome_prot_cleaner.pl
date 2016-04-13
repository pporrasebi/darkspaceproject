#!/usr/bin/perl -w
# Program to clean up protein identifiers given in the reactome proteins file.
# v1.0, 12.04.2016


use strict; 

print "Thinking...\n\n";

# First, let's declare initial variables

my @lines = ();


my @prot = ();	
my @react = ();	
my @type = ();	


my @prot_clean = ();

# Ask for the input file and read in the data

my $infile = $ARGV[0];
my $outfile = $ARGV[1];

open( INFILE, "<$infile" ) or die( "Couldn't open $infile: $!\n" );
@lines = <INFILE>;
close( INFILE );

# Prepare the lines to be readable 
# (split the columns). 

chomp @lines;

for my $i ( 0 .. $#lines ) {
	( $prot[$i], $react[$i], $type[$i]) = split( /\t/, $lines[ $i ] );
}

open ( OUTFILE, ">$outfile" ) or die( "Could't create $outfile: $!\n" );

for my $i ( 0 .. $#lines ) {

# First let's obtain the clean IDs
	
	if ( $prot[ $i ] =~ m/
			(([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}))
			-([0-9]{1,2})	# new uniprot regular expression, plus isoforms
			/x ) {
			$prot_clean[ $i ] = $1;
		}
	elsif ( $prot[ $i ] =~ m/
			(([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}))
			-(PRO_\d{10})	# new uniprot regular expression, plus PRO identifiers
			/x ) {
			$prot_clean[ $i ] = $1;
		}
	elsif ( $prot[ $i ] =~ m/
			([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})	# new uniprot regular expression
			/x ) {
			$prot_clean[ $i ] = $1;
		}
	else {
		$prot_clean[ $i ] = $prot[$i];
	}
	
    print( OUTFILE  "$prot[$i]\t$prot_clean[$i]\t$react[$i]\t$type[$i]\n");
}


close( OUTFILE );


print( "Finished!\nYour results can be found in $outfile.\n\n");
