#!/usr/bin/perl -w
# Program to extract pairIDs, protein pairs, interaction ACs, database and PMIDs from a MITAB25 file. 
# v1.0, 03.07.2015


use strict; 

print "Thinking...\n\n";

# First, let's declare initial variables

my @mitab_lines = ();


my @id_a = (); #ID(s) interactor A
my @id_b = (); #ID(s) interactor B
my @altid_a = (); #Alt. ID(s) interactor A
my @altid_b = (); #Alt. ID(s) interactor B
my @alias_a = (); #Alias(es) interactor A
my @alias_b = (); #Alias(es) interactor B
my @intdet = (); #Interaction detection method(s)
my @firstauth = (); #Publication 1st author(s)
my @pubid = (); #Publication Identifier(s)
my @taxid_a = (); #Taxid interactor A
my @taxid_b = (); #Taxid interactor B
my @inttype = (); #Interaction type(s)
my @sourcedb = (); #Source database(s)
my @int_id = (); #Interaction identifier(s)
my @score = (); #Confidence value(s)

my @pair_id = ();
my @pair_id_clean = ();
my @id_a_clean = ();
my @id_b_clean = ();


# Ask for the input file and read in the data

my $mitabfile = $ARGV[0];
my $outfile = $ARGV[1];

open( MITAB, "<$mitabfile" ) or die( "Couldn't open $mitabfile: $!\n" );
@mitab_lines = <MITAB>;
close( MITAB );

# Prepare the lines to be readable 
# (split the columns). 

chomp @mitab_lines;

for my $i ( 0 .. $#mitab_lines ) {
	( $id_a[ $i ], $id_b[ $i ], $altid_a[ $i ], $altid_b[ $i ], $alias_a[ $i ], $alias_b[ $i ], $intdet[ $i ], $firstauth[ $i ], $pubid[ $i ], $taxid_a[ $i ], $taxid_b[ $i ], $inttype[ $i ], $sourcedb[ $i ], $int_id[ $i ], $score[ $i ]) = split( /\t/, $mitab_lines[ $i ] );
}

open ( OUTFILE, ">$outfile" ) or die( "Could't create $outfile: $!\n" );

print( OUTFILE "pair_id\tid_a\tid_b\tpair_id_clean\tid_a_clean\tid_b_clean\ttaxid_a\ttaxid_b\tpubid\n");

for my $i ( 1 .. $#mitab_lines ) {

# First let's obtain the clean IDs
	
	if ( $id_a[ $i ] =~ m/uniprotkb:	# uniprotkb:
			([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})	# new uniprot regular expression
			/x ) {
			$id_a_clean[ $i ] = $1;
		}
	elsif ( $id_a[ $i ] !~ m/uniprotkb:	# no uniprotkb found:
            ([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})	# new uniprot regular expression
            /x ) {
                if ( $altid_a[ $i ] =~ m/uniprotkb:	# uniprotkb:
                    ([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})	# new uniprot regular expression
                /x ) {
                    $id_a_clean[ $i ] = $1;
                }
                
                elsif ( $id_a[ $i ] =~ m/intact:	#
                (EBI-\S{1,7})
                /x ) {
                    $id_a_clean[ $i ] = $1; 
                }
                
        
                elsif ( $id_a[ $i ] =~ m/intact:	#
                (\S{1,})
                /x ) {
                    $id_a_clean[ $i ] = $1;
                }
        
                elsif ( $id_a[ $i ] =~ m/(ddbj\S{1,}:)	#
                (\S{1,})
                /x ) {
                    $id_a_clean[ $i ] = $2;
                }
                        
                elsif ( $id_a[ $i ] =~ m/chebi:	# chebi:
                \"
                (CHEBI:\S{1,})
                \"
                /x ) {
                    $id_a_clean[ $i ] = $1; 
                }

                elsif ( $id_a[ $i ] =~ m/irefindex:	#
                (\S{1,})
                /x ) {
                    $id_a_clean[ $i ] = $1;
                }
                
                elsif ( $id_a[ $i ] =~ m/hgnc:	#
                (\S{1,})
                /x ) {
                    $id_a_clean[ $i ] = $1;
                }
                
                elsif ( $id_a[ $i ] =~ m/ensembl:	#
                (ENS\S{1,})	#
                /x ) {
                    $id_a_clean[ $i ] = $1; 
                }

		elsif ( $id_a[ $i ] =~ m/refseq:	#
                (\S{1,})
                /x ) {
                    $id_a_clean[ $i ] = $1; 
                }
		
                elsif ( $id_a[ $i ] =~ m/dip:	#
                (DIP-\S{1,})
                /x ) {
                    $id_a_clean[ $i ] = $1; 
                }
                
            }

	else {
		$id_a_clean[ $i ] = $id_a[ $i ];
	}
	
	if ($id_b[ $i ] eq "-") {
	$id_b[ $i ] = $id_a[$i];
	}
    
	elsif ( $id_b[ $i ] =~ m/uniprotkb:	# uniprotkb:
		([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})	# new uniprot regular expression	
			/x ) {
			$id_b_clean[ $i ] = $1;
		}
	elsif ( $id_b[ $i ] !~ m/uniprotkb:	# no uniprotkb found:
    	([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})	# new uniprot regular expression
    /x ) {
		if ( $altid_b[ $i ] =~ m/uniprotkb:	# uniprotkb:
		    ([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})	# new uniprot regular expression
		/x ) {
		    $id_b_clean[ $i ] = $1;
		}
		
		elsif ( $id_b[ $i ] =~ m/intact:	#
		(EBI-\S{1,7})
		/x ) {
		    $id_b_clean[ $i ] = $1;
		}
        
		elsif ( $id_b[ $i ] =~ m/intact:	#
		(\S{1,})
		/x ) {
		    $id_b_clean[ $i ] = $1;
		}
        
		elsif ( $id_b[ $i ] =~ m/(ddbj\S{1,}:)	#
		(\S{1,})
		/x ) {
		    $id_b_clean[ $i ] = $2;
		}
        
		elsif ( $id_b[ $i ] =~ m/chebi:	# chebi:
		\"
		(CHEBI:\S{1,})
		\"
		/x ) {
		    $id_b_clean[ $i ] = $1;
		}
        
		elsif ( $id_b[ $i ] =~ m/irefindex:	#
		(\S{1,})
		/x ) {
		    $id_b_clean[ $i ] = $1;
		}
        
		elsif ( $id_b[ $i ] =~ m/hgnc:	#
		(\S{1,})
		/x ) {
		    $id_b_clean[ $i ] = $1;
		}
        
		elsif ( $id_b[ $i ] =~ m/ensembl:	#
		(ENS\S{1,})	#
		/x ) {
		    $id_b_clean[ $i ] = $1;
		}
        
		elsif ( $id_b[ $i ] =~ m/refseq:	#
		(\S{1,})
		/x ) {
		    $id_b_clean[ $i ] = $1;
		}
	
		elsif ( $id_b[ $i ] =~ m/dip:	#
		(DIP-\S{1,})
		/x ) {
		    $id_b_clean[ $i ] = $1;
		}
        
	}
 
	else {
		$id_b_clean[ $i ] = $id_b[ $i ];
	}

# Now let's obtain the uniprot accessions including isoforms and PRO identifiers

	if ( $id_a[ $i ] =~ m/uniprotkb:	# uniprotkb:
			(([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})-([0-9]{1,2}))	# new uniprot regular expression, plus isoforms
			/x ) {
			$id_a[ $i ] = $1;
		}
	elsif ( $id_a[ $i ] =~ m/uniprotkb:	# uniprotkb:
			(([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})-(PRO_\d{10}))	# new uniprot regular expression, plus PRO identifiers
			/x ) {
			$id_a[ $i ] = $1;
		}
	elsif ( $id_a[ $i ] =~ m/uniprotkb:	# uniprotkb:
			([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})	# new uniprot regular expression
			/x ) {
			$id_a[ $i ] = $1;
		}
	elsif ( $id_a[ $i ] !~ m/uniprotkb:	# no uniprotkb found:
            ([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})	# new uniprot regular expression
            /x ) {
                if ( $altid_a[ $i ] =~ m/uniprotkb:	# uniprotkb:
                    ([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})	# new uniprot regular expression
                /x ) {
                    $id_a[ $i ] = $1;
                }
                
                elsif ( $id_a[ $i ] =~ m/intact:	#
                (EBI-\S{1,7})
                /x ) {
                    $id_a[ $i ] = $1; 
                }
                
                elsif ( $id_a[ $i ] =~ m/intact:	#
                (\S{1,})
                /x ) {
                    $id_a[ $i ] = $1;
                }
		
                elsif ( $id_a[ $i ] =~ m/(ddbj\S{1,}:)	#
                (\S{1,})
                /x ) {
                    $id_a[ $i ] = $2;
                }
                
                elsif ( $id_a[ $i ] =~ m/chebi:	# chebi:
                \"
                (CHEBI:\S{1,})
                \"
                /x ) {
                    $id_a[ $i ] = $1; 
                }

                elsif ( $id_a[ $i ] =~ m/irefindex:	#
                (\S{1,})
                /x ) {
                    $id_a[ $i ] = $1;
                }
                
                elsif ( $id_a[ $i ] =~ m/hgnc:	#
                (\S{1,})
                /x ) {
                    $id_a[ $i ] = $1;
                }
                
                elsif ( $id_a[ $i ] =~ m/ensembl:	#
                (ENS\S{1,})	#
                /x ) {
                    $id_a[ $i ] = $1; 
                }

                elsif ( $id_a[ $i ] =~ m/refseq:	#
                (\S{1,})
                /x ) {
                    $id_a[ $i ] = $1; 
                }
		
                elsif ( $id_a[ $i ] =~ m/dip:	#
                (DIP-\S{1,})
                /x ) {
                    $id_a[ $i ] = $1; 
                }
                
            }

	else {
		$id_a[ $i ] = $id_a[ $i ];
	}
	
	
	if ( $id_b[ $i ] =~ m/uniprotkb:	# uniprotkb:
		(([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})-([0-9]{1,2}))	# new uniprot regular expression, plus isoforms
		/x ) {
		$id_b[ $i ] = $1;
		}
	elsif ( $id_b[ $i ] =~ m/uniprotkb:	# uniprotkb:
			(([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})-(PRO_\d{10}))	# new uniprot regular expression, plus PRO identifiers
			/x ) {
			$id_b[ $i ] = $1;
		}
	elsif ( $id_b[ $i ] =~ m/uniprotkb:	# uniprotkb:
			([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})	# new uniprot regular expression
			/x ) {
			$id_b[ $i ] = $1;
		}
	elsif ( $id_b[ $i ] !~ m/uniprotkb:	# no uniprotkb found:
            ([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})	# new uniprot regular expression
            /x ) {
                if ( $altid_b[ $i ] =~ m/uniprotkb:	# uniprotkb:
                    ([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})	# new uniprot regular expression
                /x ) {
                    $id_b[ $i ] = $1;
                }
                
                elsif ( $id_b[ $i ] =~ m/intact:	#
                (EBI-\S{1,7})
                /x ) {
                    $id_b[ $i ] = $1; 
                }
                
                elsif ( $id_b[ $i ] =~ m/intact:	#
                (\S{1,})
                /x ) {
                    $id_b[ $i ] = $1;
                }
		
                elsif ( $id_b[ $i ] =~ m/(ddbj\S{1,}:)	#
                (\S{1,})
                /x ) {
                    $id_b[ $i ] = $2;
                }
                
                elsif ( $id_b[ $i ] =~ m/chebi:	# chebi:
                \"
                (CHEBI:\S{1,})
                \"
                /x ) {
                    $id_b[ $i ] = $1; 
                }

                elsif ( $id_b[ $i ] =~ m/irefindex:	#
                (\S{1,})
                /x ) {
                    $id_b[ $i ] = $1;
                }
                
                elsif ( $id_b[ $i ] =~ m/hgnc:	#
                (\S{1,})
                /x ) {
                    $id_b[ $i ] = $1;
                }
                
                elsif ( $id_b[ $i ] =~ m/ensembl:	#
                (ENS\S{1,})	#
                /x ) {
                    $id_b[ $i ] = $1; 
                }

                elsif ( $id_b[ $i ] =~ m/refseq:	#
                (\S{1,})
                /x ) {
                    $id_b[ $i ] = $1; 
                }
		
                elsif ( $id_b[ $i ] =~ m/dip:	#
                (DIP-\S{1,})
                /x ) {
                    $id_b[ $i ] = $1; 
                }
                
            }

	else {
		$id_b[ $i ] = $id_b[ $i ];
	}
	

# Clean up the taxid field 

	if ( $taxid_a[ $i ] =~ m/taxid:	#
		(\-?\d{1,})
		\(
		/x ) {
		$taxid_a[$i] = $1;
	}
	
	if ( $taxid_b[ $i ] =~ m/taxid:	#
		(\-?\d{1,})
		\(
		/x ) {
		$taxid_b[$i] = $1;
	}
	
	
# Obtain a clean PMID
	
	if ( $pubid[ $i ] =~ m/pubmed:	#
		(\d{1,})
		/x ) {
		$pubid[$i] = $1;
	}
	elsif ( $pubid[ $i ] =~ m/pubmed:	#
		(unassigned\d{1,})
		/x ) {
		$pubid[$i] = $1;
	}
	else {
		$pubid[$i] = $pubid[$i];
	}



# Now let's create a unique id for each pair taking the two identifiers and sorting them out. 
	
	my @id_array = ();
	
	push (@id_array, $id_a[$i]);
	push (@id_array, $id_b[$i]);
	@id_array = sort @id_array;
	
	$pair_id[$i] = "$id_array[0]_$id_array[1]";
        
        my @id_array_clean = ();
	
	push (@id_array_clean, $id_a_clean[$i]);
	push (@id_array_clean, $id_b_clean[$i]);
	@id_array_clean = sort @id_array_clean;
	
	$pair_id_clean[$i] = "$id_array_clean[0]_$id_array_clean[1]";
		
	
	
    print( OUTFILE  "$pair_id[$i]\t$id_a[$i]\t$id_b[$i]\t$pair_id_clean[$i]\t$id_a_clean[$i]\t$id_b_clean[$i]\t$taxid_a[$i]\t$taxid_b[$i]\t$pubid[$i]\n");
}


close( OUTFILE );


print( "Finished!\nYour results can be found in $outfile.\n\n");
