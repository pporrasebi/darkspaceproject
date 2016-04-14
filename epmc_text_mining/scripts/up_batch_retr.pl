#!/usr/bin/perl

use strict;
use warnings;
use LWP::UserAgent;


my $file = $ARGV[0];

open(FILE, $file) || die "cannot open $file";

while(<FILE>) {
        chomp;
        my $query_term = $_;
        # query UniProt for the given identifier and retrieve tab-delimited format
        system ("wget -O output -q  \"http://www.uniprot.org/uniprot/?query=id:$query_term&format=tab&columns=id,organism-id\"");
        open(WGET, "output") || die "cannot open output";
        while(<WGET>) {
                print $_ if (/[A-Z]/ && !/Entry/); # do not print header
        }
        close(WGET);
}
close(FILE);
