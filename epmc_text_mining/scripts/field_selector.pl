#!usr/bin/perl
use strict;
use warnings;

print "Thinking...\n\n"; 

my $infile = $ARGV[0];
open (INFILE, $infile) or die "Couldn't open input obo file: $!";;

my $outfile = $ARGV[1];

open (OUTFILE, ">$outfile") or die( "Couldn't open $outfile: $!\n" );


# Read through the file one line at a time and split the lines, so each col is stored as a different variable

print ("This will print in an output file cols 1 to 7 and 11 to 14 from an input file with 14 cols.\n\n");

while (<INFILE>) {

	my $col1=();
	my $col2=();
	my $col3=();
	my $col4=();
	my $col5=();
	my $col6=();
	my $col7=();
	my $col8=();
	my $col9=();

	( $col1, $col2, $col3, $col4, $col5, $col6, $col7, $col8, $col9 ) = split( /\t/ );


	print OUTFILE ("$col1\t$col2\t$col3\t$col4\t$col5\t$col6\t$col7\t$col8\n");

}
close (INFILE);
close (OUTFILE);
print "Done!\n";
exit;