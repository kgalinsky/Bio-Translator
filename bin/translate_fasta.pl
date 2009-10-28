#!/usr/bin/env perl

# $Author$
# $Date$
# $Revision$
# $HeadURL$

use strict;
#use warnings;

use JCVI::Translator;
use Getopt::Euclid ':vars';

# Instantiate the translator
my $t = JCVI::Translator->new($ARGV_translation_table);

# Build the list of file handles or standard input

my @handles;

foreach my $file (@ARGV_fasta_files) {
    # For "-" use standard input
    if ( $file eq '-' ) {
        push @handles, \*STDIN;
        next;
    }

    # Open the file or die
    open FH, $file
      or die qq{Can't open file "$file" for reading: $!};
    push @handles, \*FH;
}

# Open standard input if no file names provided
@handles = ( \*STDIN ) unless (@handles);

# Set the input record separator
local $/ = '>';

foreach my $fh (@handles) {
    # Throw away the first "record" (i.e. it is just ">")
    $_ = <$fh>;

    while (<$fh>) {
        chomp;

        next unless ( $_ =~ /\S/ );

        # Extract the sequence and translate it
        my ( $header, $sequence ) = split /\n/, $_, 2;
        my $pep_ref = $t->translate( \$sequence );
        
        # Format the peptide and print out the record
        $$pep_ref =~ s/(.{1,60})/$1\n/g;
        print ">$header\n$$pep_ref";
    }

    close $fh;
}

=head1 NAME

translate_fasta.pl - translate fasta files

=head1 OPTIONS

=over

=item <fasta_files>...

Translate list of fasta files. Input "-" for standard input. If no files
given, translate from standard input.

=item [-]-t[ranslation[_table]] <table_id>

Translation table ID. Default is 1.

=for Euclid:
    table_id.type:      +int
    table_id.default:   1

=back

=cut

