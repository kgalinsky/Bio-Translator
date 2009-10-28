# JCVI::Translator
#
# $Author$
# $Date$
# $Revision$
# $HeadURL$

=head1 NAME

JCVI::Translator - Translate DNA sequences

=head1 SYNOPSIS

use JCVI::Translator;

    my $translator = new JCVI::Translator();
    my $translator = new JCVI::Translator(11);
    my $translator = new JCVI::Translator( 12, 'id' );
    my $translator = new JCVI::Translator( 'Yeast Mitochondrial', 'name' );
    my $translator = new JCVI::Translator( 'mito', 'name' );

    my $translator = custom JCVI::Translator( \$custom_table );
    my $translator = custom JCVI::Translator( \$custom_table, 1 );

    $translator->translate( \$seq );
    $translator->translate( \$seq, { strand => 1 } );
    $translator->translate( \$seq, { strand => -1 } );

=head1 DESCRIPTION

JCVI::Translator tries to be a robust translator object featuring translation
tables based off the the ones provided by NCBI
(http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).
Key features include the ability to handle degenerate nucleotides and to
translate to ambiguous amino acids.

The way to work with JCVI::Translator is you create a new translator using an
internal translation table or a provided one, this module will translate DNA
sequences for you.

Translator uses interbase coordinates. See below for the difference between
interbase coordinates and traditional numbering methods:

    Traditional   1 2 3 4
                  A C G T ...
    Interbase    0 1 2 3 4

Conversion methods between the two methods can depend upon what you are trying
to do, but the simple way to do this is:

    strand = 3' end <=> 5' end
    lower  = min( 5' end, 3' end ) - 1
    upper  = max( 5' end, 3' end )

For logging, it uses Log::Log4Perl. This needs to be initialized to work.

For parameter validation, uses Params::Validate. This
introduces a bit of overhead, however, for scripts that are
fully tested, validation can be disabled. See the
Params::Validate documentation.

=cut

package JCVI::Translator;

use strict;
use warnings;

use version;
our $VERSION = qv('0.4.2');

use base qw(Class::Accessor::Fast);
__PACKAGE__->mk_accessors(qw(id names _table _starts _reverse));

use Log::Log4perl qw(:easy);
use Params::Validate;

use JCVI::DNATools qw(
  %degenerate_map
  $degen_match
  @nucs
  $nuc_match
  cleanDNA
  reverse_complement
);

use JCVI::AATools qw(
  %ambiguous_forward
  $aa_match
);

# Defaults. Used by the validation functions.
our $DEFAULT_ID        = 1;
our $DEFAULT_TYPE      = 'id';
our $DEFAULT_COMPLETE  = 0;
our $DEFAULT_BOOTSTRAP = 1;
our $DEFAULT_STRAND    = 1;
our $DEFAULT_PARTIAL   = 0;
our $DEFAULT_SANITIZED = 0;

=head1 CONSTRUCTORS

=cut

=head2 new

    my $translator = new JCVI::Translator();
    my $translator = new JCVI::Translator( $id );
    my $translator = new JCVI::Translator( $id, $type );

This method creates a translator by loading a translation table from the
internal list. Pass an ID and the type of ID. By default, it will load the
tranlation table with id 1. The type of ID may be "id" or "name," which
correspond to the numeric id of the translation table or the long name of the
translation table. For instance, below are the headers for the first 3
translation tables.

    {
    name "Standard" ,
    name "SGC0" ,
    id 1 ,
    ...
    },
    {
    name "Vertebrate Mitochondrial" ,
    name "SGC1" ,
    id 2 ,
    ...
    },
    {
    name "Yeast Mitochondrial" ,
    name "SGC2" ,
    id 3 ,
    ...
    },
    ...

By default, the "Standard" translation table will be loaded. You may create a
translator with this table by calling any of the following:

    my $t = new JCVI::Translator();                     # default table
    my $t = new JCVI::Translator(1);                    # explicitly set id
    my $t = new JCVI::Translator( 1, 'id' );            # set id and type
    my $t = new JCVI::Translator( 'Standard', 'name' ); # set name
    my $t = new JCVI::Translator( 'SGC0', 'name' );     # alternate name
    my $t = new JCVI::Translator( 'standard', 'name' ); # not case-sensitive
    my $t = new JCVI::Translator( 'stan', 'name' );     # partial match ok

For partial matches, JCVI::Translator will use the first matching translation
table.

    my $t = new JCVI::Translator( 'mitochondrial', 'name' );

This will use translation table with ID 2, "Vertebrate Mitochondrial," because
that is the first match (even though "Yeast Mitochondrial" would also match).

=cut

sub new {
    TRACE('new called');

    my $class = shift;

    my ( $id, $type ) = validate_pos(
        @_,
        { default => $DEFAULT_ID },
        { default => $DEFAULT_TYPE, regex => qr/id|name/ }
    );

    TRACE( uc($type) . ': ' . $id );

    # Get the beginning DATA so that we can seek back to it
    my $start_pos = tell DATA;

    # Set up regular expression for searching.
    my $match = ( $type eq 'id' ) ? qr/id $id\b/ : qr/name ".*$id.*"/i;

    # Go through every internal table until it matches on id or name.
    my $found = 0;
    local $/ = "}";
    local $_;
    while (<DATA>) {
        if ( $_ =~ $match ) {
            $found = 1;
            last;
        }
    }

    # Reset DATA
    seek DATA, $start_pos, 0;

    # Call custom with internal table. Complete is set to 1.
    return $class->custom( \$_, 1 ) if ($found);

    # Internal table not matched.
    ERROR("Table with $type of $id not found");
    return undef;
}

=head2 custom()

    my $translator = $translator->custom( $table_ref );
    my $translator = $translator->custom( $table_ref, $complete );

Create a translator table based off a passed table reference for custom
translation tables. Loads degenerate nucleotides if $complete isn't set (this
can take a little time). The format of the translation table should reflect
those of the internal tables:

    name "Names separated; by semicolons"
    name "May have multiple lines"
    id 99
    ncbieaa  "AMINOACIDS...",
    sncbieaa "-M--------..."
    -- Base1  AAAAAAAAAA...
    -- Base2  AAAACCCCGG...
    -- Base3  ACGTACTGAC...

JCVI::Translator is a bit more permissive, see the $TABLE_REGEX regular
expression to see that actual format.

Examples:

    $translator = new Translator(
        table_ref => \'name "All Alanines; All the Time"
                       id 9000
                       ncbieaa  "AAAAAAAA"
                       sncbieaa "----M---"
                       base1     AAAAAAAA
                       base2     AACCGGTT
                       base3     ACACACAC'
    );

    $translator = new Translator(
        table_ref => \$table,
        complete  => 1
    );

=cut

# Regular expression which should match translation tables and also extracts
# relevant information.
our $TABLE_REGEX = qr/
                        ( (?:name\s+".+?".*?) + )
                        id\s+(\d+).*
                        ncbieaa\s+"([a-z*]+)".*
                        sncbieaa\s+"([a-z-]+)".*
                        base1\s+([a-z]+).*
                        base2\s+([a-z]+).*
                        base3\s+([a-z]+).*
                     /isx;

sub custom {
    TRACE('custom called');

    my $class = shift;

    my ( $table_ref, $complete ) = validate_pos(
        @_,
        { type => Params::Validate::SCALARREF },
        {
            default => $DEFAULT_COMPLETE,
            regex   => qr/^[01]$/
        }
    );

    # Match the table or return undef.
    unless ( $$table_ref =~ $TABLE_REGEX ) {
        ERROR( 'Translation table is in invalid format', $$table_ref );
        return undef;
    }

    # Store the data that has been stripped using descriptive names;
    my $names    = $1;
    my $id       = $2;
    my $residues = $3;
    my $starts   = $4;
    my $base1    = $5;
    my $base2    = $6;
    my $base3    = $7;

    my $self = $class->_new();

    $self->id($id);

    # Extract each name, massage, and push it onto names array
    while ( $names =~ /"(.+?)"/gis ) {
        my @names = split( /;/, $1 );
        local $_;
        foreach (@names) {
            s/^\s+//;
            s/\s+$//;
            s/\n/ /g;
            s/\s{2,}/ /g;
            push @{ $self->names }, $_ if $_;
        }
    }

    # Store all the hashes in $self so we don't have to keep using accessors
    my $forward_hash    = $self->_table->[0];
    my $rc_forward_hash = $self->_table->[1];

    my $starts_hash    = $self->_starts->[0];
    my $rc_starts_hash = $self->_starts->[1];

    my $reverse_hash    = $self->_reverse->[0];
    my $rc_reverse_hash = $self->_reverse->[1];

    # Chop is used to efficiently get the last character from each string
    while ( my $residue = uc( chop $residues ) ) {
        my $start = uc( chop $starts );
        my $codon = uc( chop($base1) . chop($base2) . chop($base3) );

        my $rc_codon_ref = reverse_complement( \$codon );

        # If the residue is valid, store it
        if ( $residue ne 'X' ) {
            $forward_hash->{$codon}            = $residue;
            $rc_forward_hash->{$$rc_codon_ref} = $residue;

            push @{ $reverse_hash->{$residue} },    $codon;
            push @{ $rc_reverse_hash->{$residue} }, $$rc_codon_ref;
        }

        # If the start is valid, store it
        if ( ( $start ne '-' ) ) {
            $starts_hash->{$codon}            = $start;
            $rc_starts_hash->{$$rc_codon_ref} = $start;

            push @{ $reverse_hash->{start} },    $codon;
            push @{ $rc_reverse_hash->{start} }, $$rc_codon_ref;
        }
    }

    # Unroll the translation table unless it has been marked complete
    $self->bootstrap() unless ($complete);

    return $self;
}

# Helper constructor. Instantiates the object with arrayrefs and hashrefs in
# the right places
sub _new {
    my $self  = shift->SUPER::new(
        {
            names    => [],
            _table   => [],
            _starts  => [],
            _reverse => []
        }
    );

    foreach my $func (qw( _table _starts _reverse )) {
        foreach my $rc ( 0 .. 1 ) {
            $self->$func->[$rc] = {};
        }
    }

    return $self;
}

=head1 METHODS

=cut

=head2 add_translation

    $translator->add_translation( $codon, $residue );
    $translator->add_translation( $codon, $residue, \%params );

Add a codon-to-residue translation to the translation table. $start inidicates
if this is a start codon.

Examples:

    # THESE AREN'T REAL!!!
    $translator->add_translation( 'ABA', 'G' );
    $translator->add_translation( 'ABA', 'M', 1 );

=cut

sub add_translation {
    TRACE('add_translation called');

    my $self = shift;

    my ( $codon, $residue, @p );

    ( $codon, $residue, $p[0] ) = validate_pos(
        @_,
        { regex => qr/^${nuc_match}{3}$/ },
        { regex => qr/^$aa_match$/ },
        { type  => Params::Validate::HASHREF, default => {} }
    );

    my %p = validate(
        @p,
        {
            strand => {
                default => 1,
                regex   => qr/^[+-]?1$/,
                type    => Params::Validate::SCALAR
            },
            start => {
                default => 0,
                regex   => qr/^[01]$/,
                type    => Params::Validate::SCALAR
            }
        }
    );

    my $codon_ref;
    my $rc_codon_ref;

    if ( $p{strand} == 1 ) {
        $codon_ref    = \$codon;
        $rc_codon_ref = reverse_complement( \$codon );
    }
    else {
        $rc_codon_ref = \$codon;
        $codon_ref    = reverse_complement( \$codon );
    }

    # Store residue in the starts or regular translation table.
    my $table = $p{start} ? '_starts' : '_table';
    $self->$table->[0]->{$$codon_ref}    = $residue;
    $self->$table->[1]->{$$rc_codon_ref} = $residue;

    # Store the reverse lookup
    $residue = 'start' if ( $p{start} );
    push @{ $self->_reverse->[0]->{$residue} }, $$codon_ref;
    push @{ $self->_reverse->[1]->{$residue} }, $$rc_codon_ref;
}

=head2 bootstrap

    $translator->bootstrap();

Bootstrap the translation table. Find every possible translation, even those
that involve degenerate nucleotides or ambiguous amino acids.

=cut

sub bootstrap {
    TRACE('bootstrap called');

    my $self = shift;

    # Loop through every nucleotide combination and run _translate_codon on
    # each.
    foreach my $n1 (@nucs) {
        foreach my $n2 (@nucs) {
            foreach my $n3 (@nucs) {
                $self->_translate_codon( $n1 . $n2 . $n3, $self->_table->[0] );
                $self->_translate_codon(
                    $n1 . $n2 . $n3,
                    $self->_starts->[0],
                    { start => 1 }
                );
            }
        }
    }
}

=head2 table_string

    my $table_string_ref = $translator->_table_string();
    my $table_string_ref = $translator->_table_string( $bootstrap );

Returns the table string. $bootstrap specifies whether or not this table should
try to bootstrap itself using the bootstrap function above. By default, it is
1.

Examples:

    my $table_string_ref = $translator->_table_string();
    my $table_string_ref = $translator->_table_string(0); # To not bootstrap

=cut

sub table_string {
    TRACE('table_string called');

    my $self = shift;

    my $bootstrap =
      validate_pos( @_,
        { default => $DEFAULT_BOOTSTRAP, regex => qr/^[01]$/ } );

    # Bootstrap if necessary
    $self->bootstrap() if ($bootstrap);

    # Generate the names string
    my $names = join( '; ', @{ $self->names } );

    my ( $residues, $starts );    # starts/residues string
    my @base = (undef) x 3;       # this will store the base strings

    # Loop over all stored codons. Sort the codons in the translation table and
    # starts table, then use grep to get the unique ones with the help of $prev
    # which stores the previous value
    my $prev = '';
    foreach my $codon (
        grep ( ( $_ ne $prev ) && ( $prev = $_ ),
            sort { $a cmp $b } (
                keys( %{ $self->_table->[0] } ),
                keys( %{ $self->_starts->[0] } )
              ) )
      )
    {
        $residues .= $self->_table->[0]->{$codon}  || 'X';
        $starts   .= $self->_starts->[0]->{$codon} || '-';

        # Chop up the codon because the bases are stored on separate lines
        $base[ -$_ ] .= chop $codon foreach ( 1 .. 3 );
    }

    # Generate the string
    my $string = join( "\n",
        '{',
        qq(name "$names" ,),
        qq(id $self->{id} ,),
        qq(ncbieaa  "$residues",),
        qq(sncbieaa "$starts"),
        map( {"-- Base$_  $base[$_ - 1]"} ( 1 .. 3 ) ),
        '}' );

    return \$string;
}

=head2 translate

    $pep_ref = $translator->translate( $seq_ref, \%params );

The basic function of this module. Translate the specified region of the
sequence (passed as $seq_ref) and return a reference to the translated string.
The parameters are:

    strand:     1 or -1; default = 1
    lower:      integer between 0 and length; default = 0
    upper:      integer between 0 and length; default = length
    partial:    0 or 1; default = 0
    sanitized:  0 or 1; default = 0

Translator uses interbase coordinates. lower and upper are optional parameters
such that:

    0 <= lower <= upper <= length

Translator will log and die if those conditions are not satisfied.

partial sets whether or not the sequence is a 5' partial. By default, partial
is taken to be false  and the translator will try to translate the first codon
as if it were a start codon. You can specify that the sequence is 5' partial
and the translator will skip that step.

sanitized is a flag translator know that this sequence has been stripped of
whitespace and that all the codons are capitalized. Otherwise, translator will
do that in order to speed up the translation process (see
JCVI::DNATools::cleanDNA).

To translate the following:

 0 1 2 3 4 5 6 7 8 9
  C G C G C A G G A
    ---------->

    $pep_ref = $translator->translate(
        \$sequence,
        {
            strand => 1,
            lower  => 1,
            upper  => 7
        }
    );

 0 1 2 3 4 5 6 7 8 9
  C G C G C A G G A
      <----------

    $pep_ref = $translator->translate(
        \$sequence,
        {
            strand => -1,
            lower  => 2,
            upper  => 8
        }
    );

Examples:

    my $pep_ref = $translator->translate( \'acttgacgt' );

    my $pep_ref = $translator->translate( \'acttgacgt', { strand => -1 } );

    my $pep_ref = $translator->translate(
        \'acttgacgt',
        {
            strand => -1,
            lower  => 2,
            upper  => 5
        }
    );

    my $pep_ref = $translator->translate(
        \'acttgacgt',
        {
            strand  => 1,
            lower   => 0,
            upper   => 8,
            partial => 0
        }
    );

=cut

sub translate {
    TRACE('translate called');

    my $self = shift;

    my ( $seq_ref, @p );
    ( $seq_ref, $p[0] ) = validate_pos(
        @_,
        { type => Params::Validate::SCALARREF, },
        { type => Params::Validate::HASHREF, default => {} }
    );

    my %p = validate(
        @p,
        {
            strand => {
                default => $DEFAULT_STRAND,
                regex   => qr/^[+-]?1$/,
                type    => Params::Validate::SCALAR
            },
            lower => {
                default   => 0,
                regex     => qr/^[0-9]+$/,
                type      => Params::Validate::SCALAR,
                callbacks => {
                    'lower >= 0'          => sub { $_[0] >= 0 },
                    'lower <= seq_length' => sub { $_[0] <= length($$seq_ref) }
                }
            },
            upper => {
                default   => length($$seq_ref),
                regex     => qr/^[0-9]+$/,
                type      => Params::Validate::SCALAR,
                callbacks => {
                    'upper >= 0'          => sub { $_[0] >= 0 },
                    'upper <= seq_length' => sub { $_[0] <= length($$seq_ref) }
                }
            },
            partial   => { default => $DEFAULT_PARTIAL },
            sanitized => { default => $DEFAULT_SANITIZED }
        }
    );

    # Die if upper < lower
    if ( $p{upper} < $p{lower} ) {
        FATAL "Upper $p{upper} < Lower $p{lower}";
        die "Upper $p{upper} < Lower $p{lower}";
    }

    $seq_ref = cleanDNA($seq_ref) unless ( $p{sanitized} );

    # These are necessary for the _translate function
    my $prep = $self->_prepare( $p{strand} );
    my $ends = $self->_endpoints( @p{qw(strand lower upper)} );

    my $peptide = '';

    $self->_start( $seq_ref, \$peptide, $ends, $prep ) unless ( $p{partial} );
    $self->_translate( $seq_ref, \$peptide, $ends, $prep );

    return \$peptide;
}

=head2 translate6

    my $pep_refs = $translator->translate6( $seq_ref );
    my $pep_refs = $translator->translate6( $seq_ref, \%params );

Translate the sequence in every possible way. Returns an array reference of all
the translations. The structure of the array is as follows:

    0: ---------->
    1:  --------->
    2:   -------->
       NNNN...NNNN
    3: <----------
    4: <---------
    5: <--------

The parameters are similar to those use in translate:

    lower:      integer between 0 and length; default = 0
    upper:      integer between 0 and length; default = length
    partial:    0 or 1; default = 0
    sanitized:  0 or 1; default = 0

Example:

    $pep_refs = $translator->translate6(\'acttgacgt');

Output:

    $pep_refs = [
                    $pep1,
                    $pep2,
                    $pep3,
                    $reverse_pep1,
                    $reverse_pep2,
                    $reverse_pep3
                ]

=cut

sub translate6 {
    TRACE('translate6 called');

    my $self = shift;

    my ( $seq_ref, @p );
    ( $seq_ref, $p[0] ) = validate_pos(
        @_,
        { type => Params::Validate::SCALARREF },
        { type => Params::Validate::HASHREF, default => {} }
    );

    my %p = validate(
        @p,
        {
            lower => {
                default   => 0,
                regex     => qr/^[0-9]+$/,
                type      => Params::Validate::SCALAR,
                callbacks => {
                    'lower >= 0'          => sub { $_[0] >= 0 },
                    'lower <= seq_length' => sub { $_[0] <= length($$seq_ref) }
                }
            },
            upper => {
                default   => length($$seq_ref),
                regex     => qr/^[0-9]+$/,
                type      => Params::Validate::SCALAR,
                callbacks => {
                    'upper >= 0'          => sub { $_[0] >= 0 },
                    'upper <= seq_length' => sub { $_[0] <= length($$seq_ref) }
                }
            },
            partial   => { default => $DEFAULT_PARTIAL },
            sanitized => { default => $DEFAULT_SANITIZED }
        }
    );

    $seq_ref = cleanDNA($seq_ref) unless ( $p{sanitized} );

    my @peptides;

    foreach my $strand ( -1, 1 ) {

        # We only need to calculate prep once for a given strand
        my $prep = $self->_prepare($strand);
        my $rc   = $prep->[0];                 # True if reverse complement
        my $fw   = ( $rc + 1 ) % 2;            # True if forward strand
        foreach ( 0 .. 2 ) {

            # Calculate endpoints and translate
            my $ends = $self->_endpoints(
                $strand,
                $p{lower} + $fw * $_,
                $p{upper} - $rc * $_
            );
            $self->_start( $seq_ref, \$peptides[ $rc * 3 + $_ ], $ends, $prep )
              unless ( $p{partial} );
            $self->_translate( $seq_ref, \$peptides[ $rc * 3 + $_ ],
                $ends, $prep );
        }
    }

    return \@peptides;
}

=head2 translate_exons

    my $pep_ref = translate_exons( $str_ref, $exons_array_ref );
    my $pep_ref = translate_exons( $str_ref, $exons_array_ref, \%params );

Translate a gene spanning multiple exons. Paramters are:

    strand:     1 or -1; default = 1
    partial:    0 or 1;  default = 0
    sanitized:  0 or 1;  default = 0

Input:

    $exons_array_ref = [
                            [$start0, $stop0],
                            [$start1, $stop1],
                            ...
                       ];

The order of the exons in the array doesn't matter. translate_exons will sort
the exons.

Example:

    $pep_ref = translate_exons(\'actgcat', [ [0,2], [3,7] ]);
    $pep_ref = translate_exons(\'actgcat', [ [0,2], [3,7] ], { strand => -1});

=cut

sub translate_exons {
    TRACE('translate_exons called');

    my $self = shift;

    my ( $seq_ref, $exons, @p );
    ( $seq_ref, $exons, $p[0] ) = validate_pos(
        @_,
        { type => Params::Validate::SCALARREF },
        { type => Params::Validate::ARRAYREF },
        { type => Params::Validate::HASHREF, default => {} }
    );

    validate_pos(
        @$exons,
        (
            {
                type      => Params::Validate::ARRAYREF,
                callbacks => {
                    'Bound not an integer' => sub {
                        foreach my $bound ( @{ $_[0] } ) {
                            return 0 unless ( $bound =~ /^\d+$/ );
                        }
                        return 1;
                    },
                    'Bound out of range' => sub {
                        foreach my $bound ( @{ $_[0] } ) {
                            return 0
                              unless ( ( $bound >= 0 )
                                && ( $bound <= length $$seq_ref ) );
                        }
                        return 1;
                    },
                    'lower <= upper' => sub {
                        return $_[0][0] <= $_[0][1];
                      }
                }
            }
          ) x @$exons
    );

    my %p = validate(
        @p,
        {
            strand => {
                default => $DEFAULT_STRAND,
                regex   => qr/^[+-]?1$/,
                type    => Params::Validate::SCALAR
            },
            partial   => { default => $DEFAULT_PARTIAL },
            sanitized => { default => $DEFAULT_SANITIZED }
        }
    );

    my @exons =
      sort { ( $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] ) * $p{strand} }
      @$exons;

    my $prep     = $self->_prepare( $p{strand} );
    my $leftover = '';
    my $peptide;

  EXON: foreach my $exon (@exons) {
        my ( $lower, $upper ) = @$exon;

      LEFTOVER: {

            # Deal with leftovers. These are codons that have been cut by
            # splicing. In the event that no codon has been cut, the leftover
            # will be the first codon of the exon.

            my $to_go = 3 - length($leftover);

            # If the exon has fewer nucleotides than what is required to
            # complete the codon, set $to_go to be the length of that exon.
            if ( ( my $length = $upper - $lower ) < $to_go ) {
                $to_go = $length;
            }

            #  Complete the leftover and increment the start index.
            unless ( $prep->[0] ) {
                $leftover .= substr( $$seq_ref, $lower, $to_go );
                $lower += $to_go;
            }
            else {
                $upper -= $to_go;
                $leftover = substr( $$seq_ref, $upper, $to_go ) . $leftover;
            }

            # If leftover isn't long enough, then move to the next exon.
            next EXON if ( length($leftover) < 3 );
        }

      START: {

            # Handle the start codon. After the start codon has been
            # translated, set the partial flag so we don't try it again.

            my $ends = [ 0, $prep->[1] ];
            unless ( $p{partial} ) {
                $self->_start( \$leftover, \$peptide, $ends, $prep );
                $p{partial} = 1;
            }

            $self->_translate( \$leftover, \$peptide, $ends, $prep );
        }

        my $ends = $self->_endpoints( $p{strand}, $lower, $upper );
      BOUNDS: {
            my $phase_diff = ( $upper - $lower ) % 3;
            $leftover =
              $prep->[0]
              ? substr( $$seq_ref, $lower,     $phase_diff )
              : substr( $$seq_ref, $ends->[1], $phase_diff );
        }

        $self->_translate( $seq_ref, \$peptide, $ends, $prep );
    }

    return \$peptide;
}

# Returns [ $is_this_a_reverse_complement, $increment_for_loop ]
sub _prepare {
    my $self = shift;
    my ($strand) = @_;
    return [ ( $strand == 1 ? 0 : 1 ), $strand * 3 ];
}

# Convert (lower, upper) into endpoints for a loop. For the + strand, we just
# adjust upper so that it is in phase with lower. However, for the - strand,
# we not only adjust lower for phase, but also subtract 3 from everything so
# that we can take the substring properly. Here is a picture that might make
# sense of this:

# Positions:             0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4
# Region of interest (-): . . . .4- - - - - - - - - - - - - - - -20 . . .
# For + strand:                  4- - -|- - -|- - -|- - -|- - >19
# For - strand:              2. . .|< - -|- - -|- - -|- - -17 - -|

# For a + strand, it returns [4, 19]. For a - strand, it returns [17, 2]. To
# get the first codon on the - strand, we take the substring starting at 17,
# and the loop will end once we decrement our counter to 2.

sub _endpoints {
    my $self = shift;
    my ( $strand, $lower, $upper ) = @_;

    return $strand == 1
      ? [ $lower, $upper - ( ( $upper - $lower ) % 3 ) ]
      : [ $upper - 3, $lower - 3 + ( ( $upper - $lower ) % 3 ) ];
}

# The actual translation function. Goes from start to stop, appends to the
# peptide sequence using the translation table provided.

sub _translate {
    my $self = shift;
    my ( $seq_ref, $pep_ref, $ends, $prep ) = @_;

    my $table = $self->_table->[ $prep->[0] ];
    while ( $ends->[0] != $ends->[1] ) {
        $$pep_ref .= $table->{ substr( $$seq_ref, $ends->[0], 3 ) } || 'X';
        $ends->[0] += $prep->[1];
    }
}

# Perform translation for only one frame and adjusts the start only if it finds
# a codon in the translation table.

sub _start {
    my $self = shift;
    my ( $seq_ref, $pep_ref, $ends, $prep ) = @_;

    return if ( $ends->[0] == $ends->[1] );

    my $start =
      $self->_starts->[ $prep->[0] ]->{ substr( $$seq_ref, $ends->[0], 3 ) };
    if ($start) {
        $$pep_ref = $start;
        $ends->[0] += $prep->[1];
    }
}

=head2 translate_codon

    my $residue = $translator->translate_codon( $codon );
    my $residue = $translator->translate_codon( $codon, \%params );

Translate a codon. Return 'X' or '-' if it isn't in the
codon table. Handles degenerate nucleotides, so if all
possible codons for an ambiguity map to the same residue,
return that residue. Will also handle ambiguous amino acids.
start dictates whether or not to translate this as a start
codon. Will also cache any new translations it finds.

Example:

    $residue = $translator->translate_codon('atg');
    $residue = $translator->translate_codon( 'tty', { strand => -1 } );
    $residue = $translator->translate_codon( 'cat', { start => 1 } );

=cut

sub translate_codon {
    TRACE("translate_codon called");

    my $self = shift;

    my ( $codon, @p );

    ( $codon, $p[0] ) = validate_pos(
        @_,
        { regex => qr/^${nuc_match}{3}$/ },
        { type  => Params::Validate::HASHREF, default => {} }

    );

    my %p = validate(
        @p,
        {
            strand => {
                default => 1,
                regex   => qr/^[+-]?1$/,
                type    => Params::Validate::SCALAR
            },
            start => {
                default => 0,
                regex   => qr/^[01]$/,
                type    => Params::Validate::SCALAR
            }
        }
    );

    $codon = uc $codon;
    my $rc = $p{strand} == 1 ? 0 : 1;

    my ( $table, $not_found );
    unless ( $p{start} ) {
        $table     = $self->_table->[$rc];
        $not_found = 'X';
    }
    else {
        $table     = $self->_starts->[$rc];
        $not_found = '-';
    }

    #    return $table->{$codon} if ( defined $table->{$codon} );
    return $self->_translate_codon( $codon, $table, \%p ) || $not_found;
}

# This is the helper function for translate_codon. It is designed to speed
# things up because it doesn't perform validation or try to figure out which
# tables to use, which can slow things down since this is a recursive function.
# Handles codons with degenerate nucleotides: [RYMKWS] [BDHV] or N. Several
# codons may map to the same amino acid. If all possible codons for an
# amibguity map to the same residue, return that residue rather than X.

sub _translate_codon {
    my $self  = shift;
    my $codon = shift;
    my $table = shift;

    # Check for base case: no degenerate nucleotides; we can't unroll further.
    unless ( $codon =~ /($degen_match)/ ) {
        return $table->{$codon};
    }

    # Check to see if this degenerate-containing codon has been computed
    return $table->{$codon} if ( $table->{$codon} );

    my $consensus;
    my $nuc = $1;

    # Replace the nucleotide with every possiblity from degenerate map hash.
    foreach ( @{ $degenerate_map{$nuc} } ) {
        my $new_codon = $codon;
        $new_codon =~ s/$nuc/$_/;

        # Recursively call this function
        my $residue = $self->_translate_codon( $new_codon, $table, @_ );

        # If the new_codon didn't come to a consensus, or if the translation
        # isn't defined for new_codon in a custom translation table, return
        # undef.
        return undef unless ( defined $residue );

        # If consensus isn't set, set it to the current residue.
        $consensus = $residue unless ($consensus);

        # This is an interesting step. If the residue isn't the same as the
        # consensus, check to see if they map to the same ambiguous amino acid.
        # If true, then change the consensus to that ambiguous acid and proceed.
        # Otherwise, return undef (consensus could not be reached).
        if ( $residue ne $consensus ) {
            if (
                   ( defined $ambiguous_forward{$residue} )
                && ( defined $ambiguous_forward{$consensus} )
                && ( $ambiguous_forward{$residue} eq
                    $ambiguous_forward{$consensus} )
              )
            {
                $consensus = $ambiguous_forward{$consensus};
            }
            else {
                return undef;
            }
        }
    }

    # If we got this far, it means that we have a valid consensus sequence for
    # a degenerate-nucleotide-containing codon. Cache and return results.
    DEBUG("New codon translation found: $codon => $consensus");
    $self->add_translation( $codon, $consensus, @_ );
    return $consensus;
}

1;

=head1 MISC

These are the original translation tables. The translation tables used by this
module have been boostrapped - they include translations for degenerate
nucleotides and allow ambiguous amino acids to be the targets of translation
(e.g. every effort has been made to give a translation that isn't "X").

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

=head1 AUTHOR

Kevin Galinsky, C<< <kgalinsk at jcvi.org> >>

=head1 BUGS

Please report any bugs or feature requests to
C<bug-jcvi-translator at rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=JCVI-Translator>.
I will be notified, and then you'll automatically be notified of progress on
your bug as I make changes.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc JCVI::Translator

You can also look for information at:

=over 4

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/JCVI-Translator>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/JCVI-Translator>

=item * RT: CPAN's request tracker

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=JCVI-Translator>

=item * Search CPAN

L<http://search.cpan.org/dist/JCVI-Translator>

=back

=head1 ACKNOWLEDGEMENTS

Log::Log4perl
Params::Validate

=head1 COPYRIGHT & LICENSE

Copyright 2008-2009 J. Craig Venter Institute, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

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
