# Translator
#
# $Author$
# $Date$
# $Revision$
# $HeadURL$

=head1 NAME

JCVI::Translator - JCVI Translator object

=head1 SYNOPSES

 use JCVI::Translator;

 my $translator = new Translator(
			   id => $id,
                           name => $name,
                           tableRef => $tableRef
			   );

=head1 DESCRIPTION

JCVI::Translator tries to be a robust translator object
featuring translation tables based off the the ones provided
by NCBI
(http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).
Key features include the ability to handle degenerate
nucleotides and to translate to ambiguous amino acids.

The way to work with JCVI::Translator is you create a new
translator using an internal translation table or a provided
one, and from there you can perform translations and other
functions requiring knowledge of the translation table.

Translator uses interbase numbering. See below for the
difference between interbase numbering and traditional
numbering methods:

 Traditional   1 2 3 4
               A C G T ...
 Interbase    0 1 2 3 4

Conversion methods between the two methods can depend upon
what you are trying to do, but in general, just add 1 to the
start base for interbase numbering to get the boundaries
for traditional numbering (i.e. 0-4 in interbase numbering
corresponds to bases 1-4).

For logging, it uses Log::Log4Perl. This needs to be
initialized to work. See http://log4perl.sourceforge.net/.

For parameter validation, uses Params::Validate. This
introduces a bit of overhead, however, for scripts that are
fully tested, validation can be disabled. See the
Params::Validate documentation.

=head1 AUTHOR

Kevin Galinsky, <kgalinsk@jcvi.org>

=head1 FUNCTIONS

=over

=cut

package JCVI::Translator;

use strict;
use warnings;

our $VERSION = '0.3.1';

use Log::Log4perl qw(:easy);
use Params::Validate qw(:all);

use JCVI::DNATools qw(%degenerateMap
    $degenMatch
    $nucs
    $nucMatch
    cleanDNA
    reverseComplement);

use JCVI::AATools qw(%ambiguousForward);

my $DEFAULT_ID = 1;

=item new()

=item $translator = new Translator(%params);

Creates a new translator using the given parameters. If no
parameters are given, the standard translation table will be
used. The parameters are:

 id - id of an internal translation table
 name - name or part of the name of an internal table
 tableRef - a reference to a table string
 complete - does the table string contain all degenerate nucs?

table takes precedence, and the format of the table string should
reflect that of the internal tables.

 name "(\w+(; )?)+"
 name "(\w+(; )?)+"
 id \d+
 ncbieaa "\w+"
 sncbieaa "[M-]+"
 base1 [ACGT]+
 base2 [ACGT]+
 base3 [ACGT]+

Examples:

 $translator = new Translator(); # Default translator
 $translator = new Translator('id' => 4);
 $translator = new Translator('name' => 'mitochondrial');
 $translator = new Translator('table' =>
               'name "All Alanines; All the Time"
                id 9000
                ncbieaa  "AAAAAAAA"
                sncbieaa "----M---"
                base1 AAAAAAAA
                base2 AACCGGTT
                base3 ACACACAC'
               );

=cut

sub new {
    my $type = shift;

    my %params = validate( @_,
                           {  id => { default => $DEFAULT_ID,
                                      regex   => qr/^\d+$/
                              },
                              name => { optional => 1,
                                        type     => SCALAR
                              },
                              tableRef => { optional => 1,
                                            type     => SCALARREF
                              },
                              complete => { default => 0,
                                            regex   => qr/^[01]$/
                              }
                           }
    );

    my $self = { seqRef => undef };
    bless $self, $type;

    DEBUG("New $type");

    TRACE("Params:");
    map { TRACE("$_ => $params{$_}") } keys %params;

    my $error;

    ########################################
    # If a tableref is passed, go directly
    # to loadTable. Else, go to
    # loadInternalTable which will then call
    # loadTable.

    ########################################
    # Maybe this should be split into two
    # constructors?

    if ( $params{tableRef} ) {
        DEBUG('Loading custom table');
        $error = $self->_loadTable( @params{qw(tableRef complete)} );
    }
    else {
        DEBUG('Loading internal table');
        $error = $self->_loadInternalTable( @params{qw(id name)} );
    }

    if ($error) {
        ERROR("Error: $error");
        return $error;
    }
    else {
        DEBUG("New $type successful");
        return $self;
    }
}

=item _loadInternalTable()

=item $err = $translator->_loadInternalTable($id, $name);

This method loads a table from the internal list. Gets
called from "new" if no table string is provided. Passed
either $id or $name. Not recommend to pass both, but if both
are passed, then whichever matches first is used.

=cut

sub _loadInternalTable {
    my $self = shift;

    my ( $id, $name ) = validate_pos( @_, { regex => qr/^\d+$/ }, 0 );

    ########################################
    # Set up regular expression match for
    # searching.

    my $match = $name ? qr/name ".*$name.*"/i : qr/id $id/i;

    DEBUG("_loadInternalTable called");
    TRACE("ID $id")     if ( defined $id );
    TRACE("Name $name") if ( defined $name );

    my $found = 0;

    # Get the beginning DATA input
    my $startPos = tell DATA;

    ########################################
    # Go through every internal table until
    # it matches on id or name.

    local $/ = "}";
    local $_;
    while (<DATA>) {
        if ( $_ =~ $match ) {
            $found = 1;
            last;
        }
    }

    # Reset DATA input
    seek DATA, $startPos, 0;

    ########################################
    # Call loadTable with internal table.
    # Complete is set to 1.
    return $found
        ? $self->_loadTable( \$_, 1 )
        : 'Translation table not found';
}

=item _loadTable()

=item $err = $translator->loadTable($tableRef, $complete);

Loads a table based off a passed table reference for custom
translation tables. Gets called from "new" if a table string
is provided. Loads degenerate nucs if $complete isn't set.

=cut

sub _loadTable {
    my $self = shift;

    my ( $tableRef, $complete )
        = validate_pos( @_,
                        { type => SCALARREF },
                        { default => 0,
                          regex   => qr/^[01]$/
                        }
        );

    DEBUG("_loadTable called");

    ( $$self{info}{id} ) = $$tableRef =~ /id\s+(\d+)/i;

    ########################################
    # Extract each name, massage, and push
    # it onto names array
    while ( $$tableRef =~ /name\s+"(.+?)"/gis ) {
        my @names = split( /;/, $1 );
        foreach (@names) {
            s/^\s+//;
            s/\s+$//;
            s/\n/ /g;
            s/\s{2,}/ /g;
            push @{ $$self{info}{names} }, $_ if $_;
        }
    }

    ########################################
    # Pull each string to be used for
    # translation table generation.

    my ($residues) = $$tableRef =~ /ncbieaa.+?([a-z*]+)/i;
    my ($starts)   = $$tableRef =~ /sncbieaa.+?([a-z-]+)/i;
    my ($base1)    = $$tableRef =~ /base1.+?([a-z]+)/i;
    my ($base2)    = $$tableRef =~ /base2.+?([a-z]+)/i;
    my ($base3)    = $$tableRef =~ /base3.+?([a-z]+)/i;

    ########################################
    # Chop is used to efficiently get the
    # last character from each string; like
    # pop, but for strings.

    while ( my $residue = uc( chop $residues ) ) {
        my $start = uc( chop $starts );
        my $codon = uc( chop($base1) . chop($base2) . chop($base3) );

        my $rc_codon_ref = reverseComplement( \$codon );

        if ( $residue ne 'X' ) {
            $$self{forward}{$codon}            = $residue;
            $$self{rc_forward}{$$rc_codon_ref} = $residue;
        }
        if ( ( $start ne '-' ) ) {
            $$self{starts}{$codon}            = $start;
            $$self{rc_starts}{$$rc_codon_ref} = $start;
        }

        push @{ $$self{reverse}{$residue} },    $codon;
        push @{ $$self{rc_reverse}{$residue} }, $$rc_codon_ref;
    }

    ########################################
    # Use the printTable command to fill in
    # the gaps in the translation table
    # unless the table has been marked as
    # complete.

    $self->printTable() unless $complete;

    return $$self{forward} ? 0 : 'Translation table could not be loaded';
}

=item loadSequence

=item $translator->loadSequence(\$seqRef);

Translator can cache sequences internally to be translated
from later. This can save time when translating from the
same sequence multiple times, but not all the ranges are
available immediately.

=cut

sub loadSequence {
    my $self = shift;

    my ( $seqRef, $sanitized ) = validate_pos(
        @_,
        {  type     => SCALARREF,
           callback => {
               'Sequence contains invalid nucleotides' =>
                   sub { ${ $_[0] } !~ /[^$nucMatch]/ }
           }
        },
        0
    );

    cleanDNA($seqRef) unless ($sanitized);

    DEBUG('loadSequence called');
    TRACE( 'Sequence starts with ' . substr $$seqRef, 0, 5 );
    TRACE( 'Sequence length is ' . length $$seqRef );

    $$self{seqRef} = $seqRef;
}

=item clearSequence

=item $translator->clearSequence(\$seqRef);

Clear cached sequences.

=cut

sub clearSequence {
    my $self = shift;

    DEBUG('clearSequence called');

    undef $$self{seqRef};
}

=item printTable()

=item $tableString = $translator->printTable();

Returns the complete, absolute version of the table string.
Unrolls all degenerates and everything. Due to the caching
nature of translateCodon, this routine will also store all
possibilities for translation.

=cut

sub printTable {
    my $self = shift;

    DEBUG('printTable called');

    my $names = join( '; ', @{ $$self{info}{names} } );
    my ( $residues, $starts, @base, @b );

    my @NUCS = split '', $nucs;

    ########################################
    # Rigorous loop unrolls the nucleotide
    # table. This causes the subroutine to
    # take a while to run.

    foreach (@NUCS) {
        $b[0] = $_;
        foreach (@NUCS) {
            $b[1] = $_;
            foreach (@NUCS) {
                $b[2] = $_;
                my $aa = $self->translateCodon( join( '', @b ), 1 );
                my $st = $self->translateCodon( join( '', @b ), 1, 1 );

                unless (    ( $aa eq 'X' )
                         && ( $st eq '-' ) )
                {
                    $residues .= $aa;
                    $starts   .= $st;
                    $base[$_] .= $b[$_] foreach ( 0 .. 2 );
                }
            }
        }
    }

    return
        join( "\n",
              '{',
              qq(name "$names" ,),
              qq(id $$self{info}{id} ,),
              qq(ncbieaa  "$residues",),
              qq(sncbieaa "$starts"),
              map( {"-- Base$_  $base[$_ - 1]"} ( 1 .. 3 ) ),
              '}' );
}

=item translate()

=item $pepRef = $translator->translate(%params);

The basic function of this module. Translate the specified
region of the sequence and return a reference to the
translated string. The parameters are:

 strand - 1 or -1; mandatory
 lower  - integer between 0 and length; optional
          defaults to 0
 upper  - integer between 0 and length; optional
          defaults to length
 seqRef - reference to a string; optional
 sequence - string; optional
 partial - 0 or 1; optional

Translator uses interbase coordinates. lower and upper are
optional parameters such that:

 0 <= lower <= upper <= length

Translator will log and die if those conditions are not
satisfied. strand is the only mandatory parameter. sequence
and seqRef are both optional. If both are provided, sequence
takes priority. If neither is provided, then Translator will
use a previously loaded sequence. If no sequence has been
loaded, translator will log and die.

Partial sets whether or not the sequence is a 5' partial.By
default, partial is taken to be false, and the translator
will try to translate the first codon as if it is a start
codon. If you specify partial, the translator will skip
that step.


To translate the following:

 0 1 2 3 4 5 6 7 8 9
  C G C G C A G G A
    ---------->

 $pepRef = $translator->translate(seqRef => \$sequence,
                                  strand => 1,
                                  lower  => 1,
                                  upper  => 7);

 0 1 2 3 4 5 6 7 8 9
  C G C G C A G G A
      <----------

 $pepRef = $translator->translate(seqRef => \$sequence,
                                  strand => -1,
                                  lower  => 2,
                                  upper  => 8);

Examples:

 $pepRef = $translator->translate(strand => -1);

 $pepRef = $translator->translate(strand => 1,
                                  seqRef => \'acttgacgt');

 $pepRef = $translator->translate(strand => -1,
                                  sequence => 'acttgacgt',
                                  lower => 2,
                                  upper => 5);

 $pepRef = $translator->translate(strand => +1,
                                  seqRef => \'acttgacgt',
                                  lower => 0,
                                  upper => 8,
                                  partial => 1);

=cut

sub translate {
    my $self = shift;

    DEBUG('translate called');

    my %params = validate(
        @_,
        {  strand => { default => 1,
                       regex   => qr/^[+-]?1$/,
                       type    => SCALAR
           },
           upper  => 0,
           lower  => { default => 0 },
           seqRef => {
               optional => 1,
               type     => SCALARREF,
               callback => {
                   'Sequence contains invalid nucleotides' => sub {
                       ${ $_[0] } !~ /[^$nucMatch]/;
                       }
               }
           },
           partial   => 0,
           sanitized => 0
        }
    );

    unless ( defined $params{upper} ) {
        $params{upper} =
            length( defined $params{seqRef}
                    ? ${ $params{seqRef} }
                    : $self->{seqRef}
            );
    }

    $params{exons} = [ [ $params{lower}, $params{upper} ] ];
    delete $params{lower};
    delete $params{upper};

    return $self->translateExons(%params);
}

=item translate6()

=item $pepRefs = $translator->translate6($seqRef);

Translate the sequence in every possible way. Returns an
array reference of all the translations. The
structure of the array is as follows:

 0: ---------->
 1:  --------->
 2:   -------->
    NNNN...NNNN
 3: <----------
 4: <---------
 5: <--------

I should stress that $$pepRefs[x] is not a sequence, but a
reference to a sequence. So to access a sequence, you need
to do the following:

 $sequence = ${$$pepRefs[x]}

Example:

 $pepRefs = $translator->translate6(\'acttgacgt');

Output:

 $pepRefs = [$pepRef0,
             $pepRef1,
             $pepRef2,
             $reversePepRef0,
             $reversePepRef1,
             $reversePepRef2]

=cut

sub translate6 {
    my ( $self, $seqRef, $sanitized ) = @_;

    DEBUG("translate6 called");

    cleanDNA($seqRef) unless ($sanitized);

    my @pepRefs;

    push @pepRefs,
        $self->translate( seqRef    => $seqRef,
                          lower     => $_,
                          strand    => 1,
                          sanitized => 1
        ) foreach ( 0 .. 2 );
    push @pepRefs,
        $self->translate( seqRef    => $seqRef,
                          upper     => length($$seqRef) - $_,
                          strand    => -1,
                          sanitized => 1
        ) foreach ( 0 .. 2 );

    return \@pepRefs;
}

=item translateExons()

=item $pepRef = translateExons(%params);

Translate a gene spanning multiple exons. Paramters are:

 strand: 1 or -1; mandatory
 seqRef: reference to sequence; optional if one is loaded
 partial: '0' or '1'; optional, defaults to '0'

Input:

 $exonRangesRef = [
                   [$start0, $stop0],
                   [$start1, $stop1],
                    ...
                  ]

Example:

 $pepRef = translateExons(\'actgcat', [ [0,2], [3,7] ], 1);

=cut

sub translateExons {
    my $self = shift;

    my %params = validate(
        @_,
        {  strand => { regex => qr/^[+-]?1$/,
                       type  => SCALAR
           },
           seqRef => {
               default  => $$self{seqRef},
               type     => SCALARREF,
               callback => {
                   'Sequence contains invalid nucleotides' => sub {
                       ${ $_[0] } !~ /[^$nucMatch]/;
                       }
               }
           },
           exons => {
               type      => ARRAYREF,
               callbacks => {
                   'Bound not an integer' => sub {
                       foreach my $bounds ( @{ $_[0] } ) {
                           foreach my $bound (@$bounds) {
                               return 0 unless ( $bound =~ /^\d+$/ );
                           }
                       }
                       return 1;
                   },
                   'Bound out of range' => sub {
                       foreach my $bounds ( @{ $_[0] } ) {
                           foreach my $bound (@$bounds) {
                               return 0
                                   unless ( ( $bound >= 0 )
                                    && ( $bound <= length ${ $_[1]{seqRef} } )
                                   );
                           }
                       }
                       return 1;
                       }
               }
           },
           partial   => 0,
           sanitized => 0
        }
    );

    my $seqRef;

    ########################################
    # Do some further validation.

VALIDATION: {
        if ( $params{seqRef} ) {
            $seqRef = $params{seqRef};
            cleanDNA($seqRef) unless ( $params{sanitized} );
        }
        else {
            $seqRef = $$self{seqRef};
        }

        unless ( defined $seqRef ) {
            my $logger = get_logger();
            $logger->logcroak('Sequence undefined');
        }
    }

    DEBUG('translateExons called');
    TRACE("Strand is $params{strand}");
    TRACE("5' Partial") if ( $params{partial} );
    TRACE( 'Sequence starts with ' . substr ${ $params{seqRef} }, 0, 5 );
    TRACE( 'Sequence length is ' . length ${ $params{seqRef} } );

    my $increment;
    my $prefix;
    my $offset;

    if ( $params{strand} == 1 ) {
        $increment = 3;
        $prefix    = '';
        $offset    = 0;
    }
    else {
        $increment = -3;
        $prefix    = 'rc_';
        $offset    = -3;
    }

    my $peptide  = '';
    my $leftover = '';

EXON: foreach my $i ( 0 .. $#{ $params{exons} } ) {
        TRACE("Exon $i");

        my ( $lower, $upper );

    VALIDATE: {
            my $exon
                = $params{strand} == 1
                ? $params{exons}[$i]
                : $params{exons}[ -( $i + 1 ) ];

            ( $lower, $upper ) = validate_pos(
                @$exon,
                {  'Lower out of range' => sub {
                       ( ( 0 <= $_[0] ) && ( $_[0] <= $_[1][1] ) );
                       }
                },
                {  'Upper out of range' => sub {
                       (     ( $_[1][0] <= $_[0] )
                          && ( $_[0] <= length $$seqRef ) );
                       }
                }
            );

            TRACE("Lower $lower");
            TRACE("Upper $upper");
        }

    LEFTOVER: {

            ########################################
            # Deal with leftovers (exons that cut
            # codons) and exons that are short.
            # If the exon has fewer nucleotides than
            # what is required to complete the
            # codon, dump the nucleotides into the
            # growing codon and then continue.
            # Otherwise, complete the exon and
            # increment the start index.

            my $togo = 3 - length $leftover;

            if ( ( my $length = $upper - $lower ) <= $togo ) {
                WARN("Exon very short: $length <= $togo");

                ########################################
                # For forward direction, append to
                # $leftover, otherwise, prepend.

                unless ($prefix) {
                    $leftover .= substr( $$seqRef, $lower, $length );
                }
                else {
                    $leftover
                        = substr( $$seqRef, $lower, $length ) . $leftover;
                }

                next EXON;
            }
            else {
                unless ($prefix) {
                    $leftover .= substr( $$seqRef, $lower, $togo );
                    $lower += $togo;
                }
                else {
                    $upper -= $togo;
                    $leftover = substr( $$seqRef, $upper, $togo ) . $leftover;
                }
            }
        }

    PARTIAL: {

            ########################################
            # Handle 5' partials. The first exon may
            # be the actual start of the gene, so
            # the option to have the start codon
            # translated is left in there. Otherwise
            # just translate the leftover codon like
            # a regular codon.

            if ( $params{partial} ) {
                $peptide .= $$self{ $prefix . 'forward' }{$leftover};
            }
            else {
                $peptide 
                    = $$self{ $prefix . 'starts' }{$leftover}
                    || $$self{ $prefix . 'forward' }{$leftover}
                    || 'X';
                $params{partial} = 1;
            }
        }

        my $start;
        my $stop;

    BOUNDS: {
            ########################################
            # Similar to translate. Set up looping
            # variables.

            my $phaseDiff = ( $upper - $lower ) % 3;

            unless ($prefix) {
                $start    = $lower;
                $stop     = $upper - $phaseDiff;
                $leftover = substr( $$seqRef, $stop, $phaseDiff );
            }
            else {
                $start    = $upper - 3;
                $stop     = $lower + $phaseDiff - 3;
                $leftover = substr( $$seqRef, $lower, $phaseDiff );
            }
        }

        ########################################
        # The guts of the translation routine.
        # This routine is very fast and is what
        # should be running over the course of
        # most of the gene. Dealing with where
        # the exons start and stop slows down
        # the execution.

        for ( $start = $start; $start != $stop; $start += $increment ) {
            $peptide .= $$self{ $prefix . 'forward' }
                { substr( $$seqRef, $start, 3 ) } || 'X';
        }
    }

    return \$peptide;
}

=item translateCodon()

=item $residue = $translator->translateCodon($codon, $start);

Translate a codon. Return 'X' or '-' if it isn't in the
codon table. Handles degenerate nucleotides, so if all
possible codons for an ambiguity map to the same residue,
return that residue. Will also handle ambiguous amino acids.
$start dictates whether or not to translate this as a start
codon. Will also cache any new translations it finds.

For those looking for the translateStart routine, it has
been merged into translateCodon.

Example:

 $residue = $translator->translateCodon('atg');
 $residue = $translator->translateCodon('tty', 1);
 $residue = $translator->translateCodon('cat', -1, 1);

=cut

sub translateCodon {
    my $self = shift;

    my ( $codon, $strand, $start )
        = validate_pos( @_,
                        { regex => qr/^${nucMatch}{3}$/ },
                        { default => 1,
                          regex   => qr/^[+-]?1$/,
                          type    => SCALAR
                        },
                        { default => 0,
                          regex   => qr/^[01]$/,
                          type    => SCALAR
                        }
        );
    $codon = uc $codon;

    my $prefix = $strand == 1 ? '' : 'rc_';

    DEBUG("translateCodon called");
    TRACE("Codon $codon");

    my $table;
    my $notFound;
    unless ($start) {
        $table    = $prefix . 'forward';
        $notFound = 'X';
    }
    else {
        $table    = $prefix . 'starts';
        $notFound = '-';
    }
    TRACE("Using $table table");

    return $$self{$table}{$codon} if ( defined $$self{$table}{$codon} );

    ########################################
    # Handles codons with degenerate
    # nucleotides: [RYMKWS] [BDHV] or N
    # Several codons may map to the same
    # amino acid. If all possible codons for
    # an amibguity map to the same residue,
    # return that residue rather than X

    if ( my ($nuc) = $codon =~ /($degenMatch)/ ) {
        my $consensus;

        ########################################
        # Replace the nucleotide with every
        # possiblity from degenerate map hash.

        foreach ( @{ $degenerateMap{$nuc} } ) {
            my $newCodon = $codon;
            $newCodon =~ s/$nuc/$_/;
            my $residue = $self->translateCodon( $newCodon, $strand, $start );

            ########################################
            # If consensus isn't set, set it to the
            # current residue.

            $consensus = $residue unless $consensus;

            ########################################
            # If the returned residue was an 'X' or
            # '-', we return that.

            return $notFound if ( $residue eq $notFound );

            ########################################
            # This is an interesting step. If the
            # residue isn't the same as the
            # consensus, check to see if they map to
            # the same ambiguous amino acid. If
            # true, then change the consensus to
            # that ambiguous acid and proceed.
            # Otherwise, return $notFound.

            if ( $residue ne $consensus ) {
                if (    ( defined $ambiguousForward{$residue} )
                     && ( defined $ambiguousForward{$consensus} )
                     && ( $ambiguousForward{$residue} eq
                          $ambiguousForward{$consensus} )
                    )
                {
                    $consensus = $ambiguousForward{$consensus};
                }
                else {
                    return $notFound;
                }
            }
        }

        ########################################
        # Cache the residue in the translation
        # tables.

        DEBUG("New codon translation found: $codon => $consensus");

        my $rc_codon_ref = reverseComplement( \$codon );

        $$self{$table}{$codon} = $consensus;
        $$self{ 'rc_' . $table }{$$rc_codon_ref} = $consensus;

        ########################################
        # In the case of regular codons, push
        # that codon onto the reverse
        # translation array.

        unless ($start) {
            push @{ $$self{reverse}{$consensus} },    $codon;
            push @{ $$self{rc_reverse}{$consensus} }, $$rc_codon_ref;
        }

        return $consensus;
    }

    return $notFound;
}

1;    #end of module

=back

=cut

=head1 MISC

These are the original translation tables. The translation
tables below have been unrolled - degenerate nucleotides are
in place as well as ambiguous amino acids.

{
name "Standard" ,
name "SGC0" ,
id 1 ,
ncbieaa  "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
sncbieaa "---M---------------M---------------M----------------------------"
-- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
-- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
-- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
},
{
name "Vertebrate Mitochondrial" ,
name "SGC1" ,
id 2 ,
ncbieaa  "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
sncbieaa "--------------------------------MMMM---------------M------------"
-- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
-- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
-- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
},
{
name "Yeast Mitochondrial" ,
name "SGC2" ,
id 3 ,
ncbieaa  "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
sncbieaa "----------------------------------MM----------------------------"
-- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
-- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
-- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
},
{
name "Mold Mitochondrial; Protozoan Mitochondrial;"
name "Coelenterate Mitochondrial; Mycoplasma; Spiroplasma" ,
name "SGC3" ,
id 4 ,
ncbieaa  "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
sncbieaa "--MM---------------M------------MMMM---------------M------------"
-- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
-- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
-- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
},
{
name "Invertebrate Mitochondrial" ,
name "SGC4" ,
id 5 ,
ncbieaa  "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
sncbieaa "---M----------------------------MMMM---------------M------------"
-- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
-- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
-- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
},
{
name "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear" ,
name "SGC5" ,
id 6 ,
ncbieaa  "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
sncbieaa "-----------------------------------M----------------------------"
-- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
-- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
-- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
},
{
name "Echinoderm Mitochondrial; Flatworm Mitochondrial" ,
name "SGC8" ,
id 9 ,
ncbieaa  "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
sncbieaa "-----------------------------------M---------------M------------"
-- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
-- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
-- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
},
{
name "Euplotid Nuclear" ,
name "SGC9" ,
id 10 ,
ncbieaa  "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
sncbieaa "-----------------------------------M----------------------------"
-- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
-- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
-- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
},
{
name "Bacterial and Plant Plastid" ,
id 11 ,
ncbieaa  "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
sncbieaa "---M---------------M------------MMMM---------------M------------"
-- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
-- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
-- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
},
{
name "Alternative Yeast Nuclear" ,
id 12 ,
ncbieaa  "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
sncbieaa "-------------------M---------------M----------------------------"
-- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
-- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
-- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
},
{
name "Ascidian Mitochondrial" ,
id 13 ,
ncbieaa  "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
sncbieaa "---M------------------------------MM---------------M------------"
-- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
-- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
-- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
},
{
name "Alternative Flatworm Mitochondrial" ,
id 14 ,
ncbieaa  "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
sncbieaa "-----------------------------------M----------------------------"
-- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
-- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
-- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
} ,
{
name "Blepharisma Macronuclear" ,
id 15 ,
ncbieaa  "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
sncbieaa "-----------------------------------M----------------------------"
-- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
-- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
-- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
} ,
{
name "Chlorophycean Mitochondrial" ,
id 16 ,
ncbieaa  "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
sncbieaa "-----------------------------------M----------------------------"
-- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
-- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
-- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
} ,
{
name "Trematode Mitochondrial" ,
id 21 ,
ncbieaa  "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
sncbieaa "-----------------------------------M---------------M------------"
-- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
-- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
-- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
} ,
{
name "Scenedesmus obliquus Mitochondrial" ,
id 22 ,
ncbieaa  "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
sncbieaa "-----------------------------------M----------------------------"
-- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
-- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
-- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
} ,
{
name "Thraustochytrium Mitochondrial" ,
id 23 ,
ncbieaa  "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
sncbieaa "--------------------------------M--M---------------M------------"
-- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
-- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
-- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
}

=cut

__DATA__

{
name "Standard; SGC0" ,
id 1 ,
ncbieaa  "KNKNKNTTTTTTTTTTTTTTTRSRSRSIIMIIIIIQHQHQHPPPPPPPPPPPPPPPRRRRRRRRRRRRRRRLLLLLLLLLLLLLLLEDEDEDAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGVVVVVVVVVVVVVVVBBB*Y*Y*YSSSSSSSSSSSSSSS*CWCCLFLFLF*JXRRRJJXJJJJJZZZJXLLL",
sncbieaa "-----------------------------M-------------------------------------------M----------------------------------------------------------------------------------------------M-----M-----M---------M-M-"
-- Base1  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTHHMMMMMMMMMMMSSSWWYYY
-- Base2  AAAAAACCCCCCCCCCCCCCCGGGGGGTTTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTMMMAAAAAACCCCCCCCCCCCCCCGGGGGTTTTTTRTTGGGTTTTTTTTAAATTTTT
-- Base3  ACGTRYACGTBDHKMNRSVWYACGTRYACGTHMWYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYCTYACGTRYACGTBDHKMNRSVWYACGTYACGTRYAAGAGRACGTHMWYAGRAGAGR
}
{
name "Vertebrate Mitochondrial; SGC1" ,
id 2 ,
ncbieaa  "KNKNKNTTTTTTTTTTTTTTT*S*S*SMIMIXXXXXXMXXXIQHQHQHPPPPPPPPPPPPPPPRRRRRRRRRRRRRRRLLLLLLLLLLLLLLLEDEDEDAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGVVVVVVVVVVVVVVVBBB*Y*Y*YSSSSSSSSSSSSSSSWCWCWCLFLFLFJJJXZZZLLL",
sncbieaa "---------------------------MMMMMMMMMMMMMMM-----------------------------------------------------------------------------------------M---------------------------------------------------M------"
-- Base1  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTMMMRSSSYYY
-- Base2  AAAAAACCCCCCCCCCCCCCCGGGGGGTTTTTTTTTTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTMMMAAAAAACCCCCCCCCCCCCCCGGGGGGTTTTTTTTTTAAATTT
-- Base3  ACGTRYACGTBDHKMNRSVWYACGTRYACGTBDHKMNRSVWYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYCTYACGTRYACGTBDHKMNRSVWYACGTRYACGTRYCTYGAGRAGR
}
{
name "Yeast Mitochondrial; SGC2" ,
id 3 ,
ncbieaa  "KNKNKNTTTTTTTTTTTTTTTRSRSRSMIMIMIQHQHQHPPPPPPPPPPPPPPPRRRRRRRRRRRRRRRTTTTTTTTTTTTTTTEDEDEDAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGVVVVVVVVVVVVVVVBBB*Y*Y*YSSSSSSSSSSSSSSSWCWCWCLFLFLFRRRZZZ",
sncbieaa "---------------------------M-M-M-------------------------------------------------------------------------------------------------------------------------------------------------"
-- Base1  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTMMMSSS
-- Base2  AAAAAACCCCCCCCCCCCCCCGGGGGGTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTMMMAAAAAACCCCCCCCCCCCCCCGGGGGGTTTTTTGGGAAA
-- Base3  ACGTRYACGTBDHKMNRSVWYACGTRYACGTRYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYCTYACGTRYACGTBDHKMNRSVWYACGTRYACGTRYAGRAGR
}
{
name "Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma; SGC3" ,
id 4 ,
ncbieaa  "KNKNKNTTTTTTTTTTTTTTTRSRSRSIIMIXXIXIXXXXIIQHQHQHPPPPPPPPPPPPPPPRRRRRRRRRRRRRRRLLLLLLLLLLLLLLLEDEDEDAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGVVVVVVVVVVVVVVVBBB*Y*Y*YSSSSSSSSSSSSSSSWCWCWCLFLFLFXXJXXRRRJJXJJJJJXXZZZXXJXXLLL",
sncbieaa "---------------------------MMMMMMMMMMMMMMM--------------------------------------M--------------------------------------------------M------------------------------------------M-M-M-MM-MM-----M-----MM---MMMMM-M-"
-- Base1  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTBDHHKMMMMMMMMMMMNRSSSSVWWWYYY
-- Base2  AAAAAACCCCCCCCCCCCCCCGGGGGGTTTTTTTTTTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTMMMAAAAAACCCCCCCCCCCCCCCGGGGGGTTTTTTTTTTTGGGTTTTTTTTTTAAATTTTTTTT
-- Base3  ACGTRYACGTBDHKMNRSVWYACGTRYACGTBDHKMNRSVWYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYCTYACGTRYACGTBDHKMNRSVWYACGTRYACGTRYGGAGGAGRACGTHMWYGGAGRGGAGRAGR
}
{
name "Invertebrate Mitochondrial; SGC4" ,
id 5 ,
ncbieaa  "KNKNKNTTTTTTTTTTTTTTTSSSSSSSSSSSSSSSMIMIXXXXXXMXXXIQHQHQHPPPPPPPPPPPPPPPRRRRRRRRRRRRRRRLLLLLLLLLLLLLLLEDEDEDAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGVVVVVVVVVVVVVVVBBB*Y*Y*YSSSSSSSSSSSSSSSWCWCWCLFLFLFXXJJJXZZZXLLL",
sncbieaa "------------------------------------MMMMMMMMMMMMMMM-----------------------------------------------------------------------------------------M--------------------------------------------M---MM---M---M---"
-- Base1  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTDKMMMRSSSWYYY
-- Base2  AAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTMMMAAAAAACCCCCCCCCCCCCCCGGGGGGTTTTTTTTTTTTAAATTTT
-- Base3  ACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYCTYACGTRYACGTBDHKMNRSVWYACGTRYACGTRYGGCTYGAGRGAGR
}
{
name "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear; SGC5" ,
id 6 ,
ncbieaa  "KNKNKNTTTTTTTTTTTTTTTRSRSRSIIMIIIIIQHQHQHPPPPPPPPPPPPPPPRRRRRRRRRRRRRRRLLLLLLLLLLLLLLLEDEDEDAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGVVVVVVVVVVVVVVVBBBQYQYQYSSSSSSSSSSSSSSS*CWCCLFLFLFZZZJZZZRRRJJJJJJJZZZJQQQLLL",
sncbieaa "-----------------------------M-------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
-- Base1  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTBBBHKKKMMMMMMMMMMSSSWYYYYYY
-- Base2  AAAAAACCCCCCCCCCCCCCCGGGGGGTTTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTMMMAAAAAACCCCCCCCCCCCCCCGGGGGTTTTTTAAATAAAGGGTTTTTTTAAATAAATTT
-- Base3  ACGTRYACGTBDHKMNRSVWYACGTRYACGTHMWYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYCTYACGTRYACGTBDHKMNRSVWYACGTYACGTRYAGRAAGRAGRACTHMWYAGRAAGRAGR
}
{
name "Echinoderm Mitochondrial; Flatworm Mitochondrial; SGC8" ,
id 9 ,
ncbieaa  "NNKNNNNNTTTTTTTTTTTTTTTSSSSSSSSSSSSSSSIIMIIIIIQHQHQHPPPPPPPPPPPPPPPRRRRRRRRRRRRRRRLLLLLLLLLLLLLLLEDEDEDAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGVVVVVVVVVVVVVVVBBB*Y*Y*YSSSSSSSSSSSSSSSWCWCWCLFLFLFJJJJJJJJXZZZJLLL",
sncbieaa "----------------------------------------M----------------------------------------------------------------------------------------------M--------------------------------------------------------M-------"
-- Base1  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTHMMMMMMMRSSSWYYY
-- Base2  AAAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTMMMAAAAAACCCCCCCCCCCCCCCGGGGGGTTTTTTTTTTTTTTTAAATTTT
-- Base3  ACGTHMWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTHMWYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYCTYACGTRYACGTBDHKMNRSVWYACGTRYACGTRYAACTHMWYGAGRAAGR
}
{
name "Euplotid Nuclear; SGC9" ,
id 10 ,
ncbieaa  "KNKNKNTTTTTTTTTTTTTTTRSRSRSIIMIIIIIQHQHQHPPPPPPPPPPPPPPPRRRRRRRRRRRRRRRLLLLLLLLLLLLLLLEDEDEDAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGVVVVVVVVVVVVVVVBBB*Y*Y*YSSSSSSSSSSSSSSSCCWCCCCCLFLFLFJRRRJJJJJJJZZZJLLL",
sncbieaa "-----------------------------M-------------------------------------------------------------------------------------------------------------------------------------------------------------------"
-- Base1  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTHMMMMMMMMMMSSSWYYY
-- Base2  AAAAAACCCCCCCCCCCCCCCGGGGGGTTTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTMMMAAAAAACCCCCCCCCCCCCCCGGGGGGGGTTTTTTTGGGTTTTTTTAAATTTT
-- Base3  ACGTRYACGTBDHKMNRSVWYACGTRYACGTHMWYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYCTYACGTRYACGTBDHKMNRSVWYACGTHMWYACGTRYAAGRACTHMWYAGRAAGR
}
{
name "Bacterial and Plant Plastid" ,
id 11 ,
ncbieaa  "KNKNKNTTTTTTTTTTTTTTTRSRSRSIIMIXXIXIXXXXIIQHQHQHPPPPPPPPPPPPPPPRRRRRRRRRRRRRRRLLLLLLLLLLLLLLLEDEDEDAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGVVVVVVVVVVVVVVVBBB*Y*Y*YSSSSSSSSSSSSSSS*CWCCLFLFLF*XXJXXRRRJJXJJJJJXXZZZXXJXLLL",
sncbieaa "---------------------------MMMMMMMMMMMMMMM--------------------------------------M--------------------------------------------------M-------------------------------------------M----MM-MM-----M-----MM---MM-M-M-"
-- Base1  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTBDHHKMMMMMMMMMMMNRSSSSVWWYYY
-- Base2  AAAAAACCCCCCCCCCCCCCCGGGGGGTTTTTTTTTTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTMMMAAAAAACCCCCCCCCCCCCCCGGGGGTTTTTTRTTTTTGGGTTTTTTTTTTAAATTTTTTT
-- Base3  ACGTRYACGTBDHKMNRSVWYACGTRYACGTBDHKMNRSVWYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYCTYACGTRYACGTBDHKMNRSVWYACGTYACGTRYAGGAGGAGRACGTHMWYGGAGRGGAGAGR
}
{
name "Alternative Yeast Nuclear" ,
id 12 ,
ncbieaa  "KNKNKNTTTTTTTTTTTTTTTRSRSRSIIMIIIIIQHQHQHPPPPPPPPPPPPPPPRRRRRRRRRRRRRRRLLSLLLLLEDEDEDAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGVVVVVVVVVVVVVVVBBB*Y*Y*YSSSSSSSSSSSSSSS*CWCCLFLFLF*JRRRJJXJJJJJZZZJL",
sncbieaa "-----------------------------M-------------------------------------------M--------------------------------------------------------------------------------------------------M----------"
-- Base1  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTHMMMMMMMMMMMSSSWY
-- Base2  AAAAAACCCCCCCCCCCCCCCGGGGGGTTTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTMMMAAAAAACCCCCCCCCCCCCCCGGGGGTTTTTTRTGGGTTTTTTTTAAATT
-- Base3  ACGTRYACGTBDHKMNRSVWYACGTRYACGTHMWYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTHMWYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYCTYACGTRYACGTBDHKMNRSVWYACGTYACGTRYAAAGRACGTHMWYAGRAA
}
{
name "Ascidian Mitochondrial" ,
id 13 ,
ncbieaa  "KNKNKNTTTTTTTTTTTTTTTGSGSGSMIMIMIQHQHQHPPPPPPPPPPPPPPPRRRRRRRRRRRRRRRLLLLLLLLLLLLLLLEDEDEDAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGVVVVVVVVVVVVVVVBBB*Y*Y*YSSSSSSSSSSSSSSSWCWCWCLFLFLFXXJJJGGGXZZZXLLL",
sncbieaa "---------------------------M-M-M------------------------------------------------------------------------------------------M--------------------------------------------M---MM------M---M---"
-- Base1  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTDKMMMRRRRSSSWYYY
-- Base2  AAAAAACCCCCCCCCCCCCCCGGGGGGTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTMMMAAAAAACCCCCCCCCCCCCCCGGGGGGTTTTTTTTTTTGGGTAAATTTT
-- Base3  ACGTRYACGTBDHKMNRSVWYACGTRYACGTRYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYCTYACGTRYACGTBDHKMNRSVWYACGTRYACGTRYGGCTYAGRGAGRGAGR
}
{
name "Alternative Flatworm Mitochondrial" ,
id 14 ,
ncbieaa  "NNKNNNNNTTTTTTTTTTTTTTTSSSSSSSSSSSSSSSIIMIIIIIQHQHQHPPPPPPPPPPPPPPPRRRRRRRRRRRRRRRLLLLLLLLLLLLLLLEDEDEDAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGVVVVVVVVVVVVVVVBBBYY*YYYYYSSSSSSSSSSSSSSSWCWCWCLFLFLFJJJJJJJJZZZJLLL",
sncbieaa "----------------------------------------M----------------------------------------------------------------------------------------------------------------------------------------------------------------"
-- Base1  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTHMMMMMMMSSSWYYY
-- Base2  AAAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTMMMAAAAAAAACCCCCCCCCCCCCCCGGGGGGTTTTTTTTTTTTTTAAATTTT
-- Base3  ACGTHMWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTHMWYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYCTYACGTHMWYACGTBDHKMNRSVWYACGTRYACGTRYAACTHMWYAGRAAGR
}
{
name "Blepharisma Macronuclear" ,
id 15 ,
ncbieaa  "KNKNKNTTTTTTTTTTTTTTTRSRSRSIIMIIIIIQHQHQHPPPPPPPPPPPPPPPRRRRRRRRRRRRRRRLLLLLLLLLLLLLLLEDEDEDAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGVVVVVVVVVVVVVVVBBB*YQYYSSSSSSSSSSSSSSS*CWCCLFLFLF*ZJZRRRJJJJJJJZZZJQLLL",
sncbieaa "-----------------------------M-------------------------------------------------------------------------------------------------------------------------------------------------------------------"
-- Base1  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTBHKMMMMMMMMMMSSSWYYYY
-- Base2  AAAAAACCCCCCCCCCCCCCCGGGGGGTTTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTMMMAAAAACCCCCCCCCCCCCCCGGGGGTTTTTTRATAGGGTTTTTTTAAATATTT
-- Base3  ACGTRYACGTBDHKMNRSVWYACGTRYACGTHMWYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYCTYACGTYACGTBDHKMNRSVWYACGTYACGTRYAGAGAGRACTHMWYAGRAGAGR
}
{
name "Chlorophycean Mitochondrial" ,
id 16 ,
ncbieaa  "KNKNKNTTTTTTTTTTTTTTTRSRSRSIIMIIIIIQHQHQHPPPPPPPPPPPPPPPRRRRRRRRRRRRRRRLLLLLLLLLLLLLLLEDEDEDAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGVVVVVVVVVVVVVVVBBB*YLYYSSSSSSSSSSSSSSS*CWCCLFLFLF*LJRRRJJJJJJJZZZJLLL",
sncbieaa "-----------------------------M-----------------------------------------------------------------------------------------------------------------------------------------------------------------"
-- Base1  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTHMMMMMMMMMMSSSWYYY
-- Base2  AAAAAACCCCCCCCCCCCCCCGGGGGGTTTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTMMMAAAAACCCCCCCCCCCCCCCGGGGGTTTTTTRWTGGGTTTTTTTAAATTTT
-- Base3  ACGTRYACGTBDHKMNRSVWYACGTRYACGTHMWYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYCTYACGTYACGTBDHKMNRSVWYACGTYACGTRYAGAAGRACTHMWYAGRAAGR
}
{
name "Trematode Mitochondrial" ,
id 21 ,
ncbieaa  "NNKNNNNNTTTTTTTTTTTTTTTSSSSSSSSSSSSSSSMIMIMIQHQHQHPPPPPPPPPPPPPPPRRRRRRRRRRRRRRRLLLLLLLLLLLLLLLEDEDEDAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGVVVVVVVVVVVVVVVBBB*Y*Y*YSSSSSSSSSSSSSSSWCWCWCLFLFLFJJJXZZZLLL",
sncbieaa "----------------------------------------M--------------------------------------------------------------------------------------------M---------------------------------------------------M------"
-- Base1  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTMMMRSSSYYY
-- Base2  AAAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTMMMAAAAAACCCCCCCCCCCCCCCGGGGGGTTTTTTTTTTAAATTT
-- Base3  ACGTHMWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTRYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYCTYACGTRYACGTBDHKMNRSVWYACGTRYACGTRYCTYGAGRAGR
}
{
name "Scenedesmus obliquus Mitochondrial" ,
id 22 ,
ncbieaa  "KNKNKNTTTTTTTTTTTTTTTRSRSRSIIMIIIIIQHQHQHPPPPPPPPPPPPPPPRRRRRRRRRRRRRRRLLLLLLLLLLLLLLLEDEDEDAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGVVVVVVVVVVVVVVVBBB*YLYY*SSSSSSS*CWCCLFLFLF****LJRRRJJJJJJJZZZJLLL",
sncbieaa "-----------------------------M-------------------------------------------------------------------------------------------------------------------------------------------------------------"
-- Base1  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTHMMMMMMMMMMSSSWYYY
-- Base2  AAAAAACCCCCCCCCCCCCCCGGGGGGTTTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTMMMAAAAACCCCCCCCGGGGGTTTTTTMRSVWTGGGTTTTTTTAAATTTT
-- Base3  ACGTRYACGTBDHKMNRSVWYACGTRYACGTHMWYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYCTYACGTYACGTBKSYACGTYACGTRYAAAAGAAGRACTHMWYAGRAAGR
}
{
name "Thraustochytrium Mitochondrial" ,
id 23 ,
ncbieaa  "KNKNKNTTTTTTTTTTTTTTTRSRSRSIIMIIXIIIQHQHQHPPPPPPPPPPPPPPPRRRRRRRRRRRRRRRLLLLLLLLLLLLLLLEDEDEDAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGVVVVVVVVVVVVVVVBBB*Y*Y*YSSSSSSSSSSSSSSS*CWCC*FLFF****RRRJJJJJJJXZZZL",
sncbieaa "-----------------------------MM-M--------------------------------------------------------------------------------------------M------------------------------------------------------------M----"
-- Base1  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTMMMMMMMMMMRSSSY
-- Base2  AAAAAACCCCCCCCCCCCCCCGGGGGGTTTTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTMMMAAAAAACCCCCCCCCCCCCCCGGGGGTTTTTDKRWGGGTTTTTTTTAAAT
-- Base3  ACGTRYACGTBDHKMNRSVWYACGTRYACGTHKMWYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTRYACGTBDHKMNRSVWYACGTBDHKMNRSVWYACGTBDHKMNRSVWYCTYACGTRYACGTBDHKMNRSVWYACGTYACGTYAAAAAGRACTHMWYGAGRG
}
{
name "Strict Standard" ,
ncbieaa  "KNKKNNTTTTTTTTTTTTTTTRSRRSSIIMIIIIIQHQQHHPPPPPPPPPPPPPPPRRRRRRRRRRRRRRRLLLLLLLLLLLLLLLEDEEDDAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGBBBVVVVVVVVVVVVVVVJRRRJJJJJJJZZZ*Y**YYSSSSSSSSSSSSSSS*CWCC*LFLLFFJLLL",
sncbieaa "-----------------------------M-----------------------------------------------------------------------------------------------------------------------------------------------------------------"
-- Base1  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGHMMMMMMMMMMSSSTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTWYYY
-- Base2  AAAAAACCCCCCCCCCCCCCCGGGGGGTTTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTAAAAAACCCCCCCCCCCCCCCGGGGGGGGGGGGGGGMMMTTTTTTTTTTTTTTTTGGGTTTTTTTAAAAAAAAACCCCCCCCCCCCCCCGGGGGRTTTTTTTTTT
-- Base3  ACGRTYABCDGHKMNRSTVWYACGRTYACGHMTWYACGRTYABCDGHKMNRSTVWYABCDGHKMNRSTVWYABCDGHKMNRSTVWYACGRTYABCDGHKMNRSTVWYABCDGHKMNRSTVWYCTYABCDGHKMNRSTVWYAAGRACHMTWYAGRACGRTYABCDGHKMNRSTVWYACGTYAACGRTYAAGR
}