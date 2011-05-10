# Bio::Tiny::Translator::Table
#
# $Author$
# $Date$
# $Revision$
# $HeadURL$

=head1 NAME

Bio::Tiny::Translator::Table - Translation table for Bio::Tiny::Translator

=head1 SYNOPSIS

    use Bio::Tiny::Translator::Table;
    
    my $table = new Bio::Tiny::Translator();
    my $table = new Bio::Tiny::Translator(11);
    my $table = new Bio::Tiny::Translator( 12, { type => 'id' } );
    my $table = new Bio::Tiny::Translator( 'Yeast Mitochondrial', { type => 'name' } );
    my $table = new Bio::Tiny::Translator( 'mito', { type => 'name' } );

    my $table = custom Bio::Tiny::Translator( \$custom_table );
    my $tale = custom Bio::Tiny::Translator( \$custom_table, { bootstrap => 0 } );


=cut

package Bio::Tiny::Translator::Table;

use strict;
use warnings;

use base qw(Class::Accessor::Fast);
__PACKAGE__->mk_accessors(qw(id names codon2aa codon2start aa2codons));

use Params::Validate;

use Bio::Tiny::Util::DNA qw(
  %degenerate2nucleotides
  %nucleotides2degenerate
  %degenerate_hierarchy

  $degenerate_match

  @DNAs

  @all_nucleotides
  $all_nucleotide_match

  unrollDNA
  reverse_complement
);

use Bio::Tiny::Util::AA qw(
  %ambiguous_forward
  %ambiguous_map
  $aa_match
);

our $DEFAULT_ID        = 1;
our $DEFAULT_TYPE      = 'id';
our $DEFAULT_BOOTSTRAP = 1;

# Helper constructor. Instantiates the object with arrayrefs and hashrefs in
# the right places
sub _new {
    shift->SUPER::new(
        {
            names        => [],
            codon2aa    => Bio::Tiny::Translator::Table::Pair->new(),
            codon2start => Bio::Tiny::Translator::Table::Pair->new(),
            aa2codons   => Bio::Tiny::Translator::Table::Pair->new()
        }
    );
}

=head1 CONSTRUCTORS

=cut

=head2 new

    my $table = Bio::Tiny::Translator::Table->new();
    my $table = Bio::Tiny::Translator::Table->new( $id );
    my $table = Bio::Tiny::Translator::Table->new( $id, \%params );

This method creates a translation table by loading a table string from the
internal list. Pass an ID and the type of ID. By default, it will load the
translation table with id 1. The type of ID may be "id" or "name," which
correspond to the numeric id of the translation table or the long name of the
translation table. For instance, below are the headers for the first 3 table
strings.

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

By default, the "Standard" translation table will be loaded. You may instantiate
this translation table by calling any of the following:

    my $t = Bio::Tiny::Translator::Table->new();
    my $t = Bio::Tiny::Translator::Table->new(1);
    my $t = Bio::Tiny::Translator::Table->new( 1,          { type => 'id' } );
    my $t = Bio::Tiny::Translator::Table->new( 'Standard', { type => 'name' } );
    my $t = Bio::Tiny::Translator::Table->new( 'SGC0',     { type => 'name' } );
    my $t = Bio::Tiny::Translator::Table->new( 'standard', { type => 'name' } );
    my $t = Bio::Tiny::Translator::Table->new( 'stan',     { type => 'name' } );

For partial matches, this module will use the first matching translation
table.

    my $t = Bio::Tiny::Translator::Table->new( 'mitochondrial', { type => 'name' } );

This will use translation table with ID 2, "Vertebrate Mitochondrial," because
that is the first match (even though "Yeast Mitochondrial" would also match).

=cut

sub new {
    my $class = shift;

    my ( $id, @p );

    # id has a default, but if supplied, must be a scalar
    ( $id, $p[0] ) = validate_pos(
        @_,
        { type => Params::Validate::SCALAR,  default => $DEFAULT_ID },
        { type => Params::Validate::HASHREF, default => {} }
    );

    # type must be either id or name
    my %p = validate(
        @p,
        {
            type => {
                default => $DEFAULT_TYPE,
                regex   => qr/id|name/
            }
        }
    );
    # Get the beginning DATA so that we can seek back to it
    my $start_pos = tell DATA;

    # Set up regular expression for searching.
    my $match = ( $p{type} eq 'id' ) ? qr/id $id\b/ : qr/name ".*$id.*"/i;

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

    # Call custom with internal table. We don't want to bootstrap.
    return $class->custom( \$_, { bootstrap => 0 } ) if ($found);

    # Internal table not matched.
    ERROR("Table with $p{type} of $id not found");
    return undef;
}

=head2 custom()

    my $table = Bio::Tiny::Translator::Table->custom( $table_ref );
    my $table = Bio::Tiny::Translator::Table->custom( $table_ref, \%params );

Create a translation table based off a passed table reference for custom
translation tables. Loads degenerate nucleotides if bootstrap isn't set (this
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

This module is a bit more permissive than that; see the $TABLE_REGEX regular
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
        bootstrap  => 0
    );

=cut

# Regular expression which should match translation tables and also extracts
# relevant information.
my $TABLE_REGEX = qr/
                        ( (?:name\s+".+?".*?) + )
                        (?:id\s+(\d+).*)?
                        ncbieaa\s+"([a-z*]+)".*
                        sncbieaa\s+"([a-z-]+)".*
                        base1\s+([a-z]+).*
                        base2\s+([a-z]+).*
                        base3\s+([a-z]+).*
                     /isx;

sub custom {
    my $class = shift;

    my ( $table_ref, @p );

    # table_ref is required and must be a refrerence to a scalar
    ( $table_ref, $p[0] ) = validate_pos(
        @_,
        { type => Params::Validate::SCALARREF },
        { type => Params::Validate::HASHREF, default => {} }
    );

    # get the bootstrap parameter
    my %p = validate(
        @p,
        {
            bootstrap => {
                default => $DEFAULT_BOOTSTRAP,
                regex   => qr/^[01]$/
            }
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
    my @bases    = ( $5, $6, $7 );

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

    # Get all the table pairs so we don't have to keep using accessors
    my $codon2aa    = $self->codon2aa;
    my $codon2start = $self->codon2start;
    my $aa2codons   = $self->aa2codons;

    # Chop is used to efficiently get the last character from each string
    while ( my $residue = uc( chop $residues ) ) {
        my $start = uc( chop $starts );

        my @nucleotides;
        foreach my $i ( 0 .. 2 ) {

            # Get the base
            my $base = chop( $bases[$i] );
            $nucleotides[$i][0] = $base;

            # If the base is a degenerate
            if ( $base =~ m/$degenerate_match/ ) {

                # Store the nucleotides it can be
                push @{ $nucleotides[$i] }, @{ $degenerate2nucleotides{$base} };

                # Store the other degenerates it can be
                if ( my $hierarchy = $degenerate_hierarchy{$base} ) {
                    push @{ $nucleotides[$i] }, @$hierarchy;
                }
            }
        }

        # Add each potential codon to the translation table
        foreach my $base1 ( @{ $nucleotides[0] } ) {
            foreach my $base2 ( @{ $nucleotides[1] } ) {
                foreach my $base3 ( @{ $nucleotides[2] } ) {
                    my $codon = join( '', $base1, $base2, $base3 );

                    my $reverse = ${ reverse_complement( \$codon ) };

                    # If the residue is valid, store it
                    if ( $residue ne 'X' ) {
                        $codon2aa->store( $residue, $codon, $reverse );
                        $aa2codons->push( $residue, $codon, $reverse );
                    }

                    # If the start is valid, store it
                    if ( ( $start ne '-' ) ) {
                        $codon2start->store( $start, $codon, $reverse );
                        $aa2codons->push( '+', $codon, $reverse );
                    }
                }
            }
        }
    }

    # Bootstrap the translation table
    $self->bootstrap() if ( $p{bootstrap} );

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
    my $self = shift;

    my ( $codon, $residue, @p );

    ( $codon, $residue, $p[0] ) = validate_pos(
        @_,
        { regex => qr/^${all_nucleotide_match}{3}$/ },
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
    my $table = $p{start} ? 'codon2start' : 'codon2aa';
    $table = $self->$table;

    $table->store( $residue, $$codon_ref, $$rc_codon_ref );

    # Store the reverse lookup
    $residue = '+' if ( $p{start} );
    $self->aa2codons->push( $residue, $$codon_ref, $$rc_codon_ref );
}

=head2 bootstrap

    $translator->bootstrap();

Bootstrap the translation table. Find every possible translation, even those
that involve degenerate nucleotides or ambiguous amino acids.

=cut

sub bootstrap {
    my $self = shift;

    # Loop through every nucleotide combination and run _translate_codon on
    # each.
    foreach my $n1 (@all_nucleotides) {
        foreach my $n2 (@all_nucleotides) {
            foreach my $n3 (@all_nucleotides) {
                $self->_unroll( $n1 . $n2 . $n3, $self->codon2aa->[0] );
                $self->_unroll(
                    $n1 . $n2 . $n3,
                    $self->codon2start->[0],
                    { start => 1 }
                );
            }
        }
    }
}

# This is the helper function for bootstrap. Handles codons with degenerate
# nucleotides: [RYMKWS] [BDHV] or N. Several codons may map to the same amino
# acid. If all possible codons for an amibguity map to the same residue, store
# that residue.

sub _unroll {
    my $self  = shift;
    my $codon = shift;
    my $table = shift;

    # Return the codon if we have it
    return $table->{$codon} if ( $table->{$codon} );

    # Check for base case: no degenerate nucleotides; we can't unroll further.
    unless ( $codon =~ /($degenerate_match)/ ) {
        return undef;
    }

    my $consensus;
    my $nuc = $1;

    # Replace the nucleotide with every possiblity from degenerate map hash.
    foreach ( @{ $degenerate2nucleotides{$nuc} } ) {
        my $new_codon = $codon;
        $new_codon =~ s/$nuc/$_/;

        # Recursively call this function
        my $residue = $self->_unroll( $new_codon, $table, @_ );

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
    $self->add_translation( $codon, $consensus, @_ );
    return $consensus;
}

=head2 string

    my $table_string_ref = $translator->string();
    my $table_string_ref = $translator->string( \%params );

Returns the table string. %params can specify whether or not this table should
try to bootstrap itself using the bootstrap function above. By default, it will
try to.

Examples:

    my $table_string_ref = $translator->string();
    my $table_string_ref = $translator->string( { bootstrap => 0 } );

=cut

sub string {
    my $self = shift;

    my $bootstrap =
      validate_pos( @_,
        { default => $DEFAULT_BOOTSTRAP, regex => qr/^[01]$/ } );

    # Bootstrap if necessary
    $self->bootstrap() if ($bootstrap);

    # Generate the names string
    my $names = join( '; ', @{ $self->names } );


    # Make the hashes of amino acid to codons and start amino acids to codons
    my %aa2codons = %{ $self->aa2codons->forward };

    my $codon2start = $self->codon2start->forward;
    my %start2codons;

    foreach my $start_codon ( @{ $aa2codons{'+'} } ) {
        my $start_aa = $codon2start->{$start_codon};
        push @{ $start2codons{$start_aa} }, $start_codon;
    }

    delete $aa2codons{'+'};

    # Minimize the codons in each group by removing those that are implied by
    # a degenerate codon
    foreach my $group ( values(%aa2codons), values(%start2codons) ) {
        my %group_hash = map { $_ => undef } @$group;
        foreach my $codon (@$group) {
            my $possibilities = unrollDNA( \$codon );
            shift(@$possibilities);
            delete @group_hash{@$possibilities};
        }
        @$group = sort keys %group_hash;
    }

    # Create the arrays that will be used to generate the string
    my ( @residues, @starts );
    my @bases = map { [] } (undef) x 3;

    foreach my $start ( sort keys %start2codons ) {
        foreach my $codon ( @{ $start2codons{$start} } ) {
            push @residues, 'X';
            push @starts,   $start;
            for my $i ( 1 .. 3 ) { push @{ $bases[ -$i ] }, chop $codon }
        }
    }
    foreach my $aa ( sort keys %aa2codons ) {
        next if ($ambiguous_map{$aa});
        foreach my $codon ( @{ $aa2codons{$aa} } ) {
            push @residues, $aa;
            push @starts,   '-';
            for my $i ( 1 .. 3 ) { push @{ $bases[ -$i ] }, chop $codon }
        }
    }
    foreach my $aa ( sort keys %ambiguous_map ) {
        my $group = $aa2codons{$aa} or next;
        foreach my $codon ( @$group ) {
            push @residues, $aa;
            push @starts,   '-';
            for my $i ( 1 .. 3 ) { push @{ $bases[ -$i ] }, chop $codon }
        }
    }

    # Generate the string
    my $string = join(
        "\n", '{',
        qq{name "$names" ,},
        'id ' . $self->id . ' ,',
        'ncbieaa  "' . join( '', @residues ) . '",',
        'sncbieaa "' . join( '', @starts ) . '",',
        map( { "-- Base$_  " . join( '', @{ $bases[ $_ - 1 ] } ) . '"' }
            ( 1 .. 3 )
        ),
        '}'
    );

    return \$string;
}

{

    package Bio::Tiny::Translator::Table::Pair;

    use strict;
    use warnings;

    use Bio::Tiny::Util::DNA qw(reverse_complement);

    sub new {
        my $class = shift;
        my $self = [ {}, {} ];
        bless $self, $class;
    }

    sub forward { return $_[0][0] }
    sub reverse { return $_[0][1] }

    sub store {
        my ( $self, $residue, $codon, $reverse ) = @_;

        $reverse ||= ${ reverse_complement($codon) };

        $self->[0]->{$codon}   = $residue;
        $self->[1]->{$reverse} = $residue;
    }

    sub push {
        my ( $self, $residue, $codon, $reverse ) = @_;

        $reverse ||= ${ reverse_complement($codon) };

        $self->[0]->{$residue} ||= [];
        $self->[1]->{$residue} ||= [];

        push @{ $self->[0]->{$residue} }, $codon;
        push @{ $self->[1]->{$residue} }, $reverse;

        foreach my $i ( 0 .. 1 ) {
            my %seen = map { $_ => undef } @{ $self->[$i]->{$residue} };
            $self->[$i]->{$residue} = [ sort { $a cmp $b } keys %seen ];
        }
    }

    1;
}

1;

=head1 MISC

These are the original translation tables. The translation tables used by this
module have been boostrapped and compacted. They were first expanded to include
translations for degenerate nucleotides and allow ambiguous amino acids to be
the targets of translation (e.g. every effort has been made to give a
translation that isn't "X"). Then, the tables were had reduntant columns
removed; any codon that was implied by the presence of another codon containing
degenerate nucleotides was removed.

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

Kevin Galinsky, <kgalinsky plus cpan at gmail dot com>

=head1 COPYRIGHT & LICENSE

Copyright 2008-2009 J. Craig Venter Institute, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

__DATA__

{
name "Standard; SGC0" ,
id 1 ,
ncbieaa  "X**ACDEFGHIKLLMNPQRRSSTVWYBJJZ",
sncbieaa "M-----------------------------",
-- Base1  HTTGTGGTGCAACYAACCCMATAGTTGHMS"
-- Base2  TARCGAATGATATTTACAGGGCCTGAMTTA"
-- Base3  GRANYYRYNYHRNRGYNRNRYNNNGYYAHR"
}
{
name "Vertebrate Mitochondrial; SGC1" ,
id 2 ,
ncbieaa  "XX**ACDEFGHIKLLMNPQRSSTVWYBJZ",
sncbieaa "MM---------------------------",
-- Base1  ARATGTGGTGCAACYAACCCATAGTTGMS"
-- Base2  TTGACGAATGATATTTACAGGCCTGAMTA"
-- Base3  NGRRNYYRYNYYRNRRYNRNYNNNRYYYR"
}
{
name "Yeast Mitochondrial; SGC2" ,
id 3 ,
ncbieaa  "X*ACDEFGHIKLMNPQRRSSTTVWYBZ",
sncbieaa "M--------------------------",
-- Base1  ATGTGGTGCAATAACCCMATACGTTGS"
-- Base2  TACGAATGATATTACAGGGCCTTGAMA"
-- Base3  RRNYYRYNYYRRRYNRNRYNNNNRYYR"
}
{
name "Mold Mitochondrial; Protozoan Mitochondrial;" ,
name "Coelenterate Mitochondrial; Mycoplasma; Spiroplasma;" ,
name "SGC3" ,
id 4 ,
ncbieaa  "XXX*ACDEFGHIKLLMNPQRRSSTVWYBJJZ",
sncbieaa "MMM----------------------------",
-- Base1  ANWTGTGGTGCAACYAACCCMATAGTTGHMS"
-- Base2  TTTACGAATGATATTTACAGGGCCTGAMTTA"
-- Base3  NGRRNYYRYNYHRNRGYNRNRYNNNRYYAHR"
}
{
name "Invertebrate Mitochondrial; SGC4" ,
id 5 ,
ncbieaa  "XX*ACDEFGHIKLLMNPQRSSTVWYBJZ",
sncbieaa "MM--------------------------",
-- Base1  ADTGTGGTGCAACYAACCCATAGTTGMS"
-- Base2  TTACGAATGATATTTACAGGCCTGAMTA"
-- Base3  NGRNYYRYNYYRNRRYNRNNNNNRYYYR"
}
{
name "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear; SGC5" ,
id 6 ,
ncbieaa  "X*ACDEFGHIKLLMNPQRRSSTVWYBJJZ",
sncbieaa "M----------------------------",
-- Base1  ATGTGGTGCAACYAACYCMATAGTTGHMB"
-- Base2  TGCGAATGATATTTACAGGGCCTGAMTTA"
-- Base3  GANYYRYNYHRNRGYNRNRYNNNGYYAHR"
}
{
name "Echinoderm Mitochondrial; Flatworm Mitochondrial; SGC8" ,
id 9 ,
ncbieaa  "X*ACDEFGHIKLLMNPQRSSTVWYBJJZ",
sncbieaa "M---------------------------",
-- Base1  RTGTGGTGCAACYAACCCATAGTTGHMS"
-- Base2  TACGAATGATATTTACAGGCCTGAMTTA"
-- Base3  GRNYYRYNYHGNRGHNRNNNNNRYYAHR"
}
{
name "Euplotid Nuclear; SGC9" ,
id 10 ,
ncbieaa  "X*ACDEFGHIKLLMNPQRRSSTVWYBJJZ",
sncbieaa "M----------------------------",
-- Base1  ATGTGGTGCAACYAACCCMATAGTTGHMS"
-- Base2  TACGAATGATATTTACAGGGCCTGAMTTA"
-- Base3  GRNHYRYNYHRNRGYNRNRYNNNGYYAHR"
}
{
name "Bacterial and Plant Plastid" ,
id 11 ,
ncbieaa  "XX**ACDEFGHIKLLMNPQRRSSTVWYBJJZ",
sncbieaa "MM-----------------------------",
-- Base1  ANTTGTGGTGCAACYAACCCMATAGTTGHMS"
-- Base2  TTARCGAATGATATTTACAGGGCCTGAMTTA"
-- Base3  NGRANYYRYNYHRNRGYNRNRYNNNGYYAHR"
}
{
name "Alternative Yeast Nuclear" ,
id 12 ,
ncbieaa  "X**ACDEFGHIKLLLMNPQRRSSSTVWYBJJZ",
sncbieaa "M-------------------------------",
-- Base1  MTTGTGGTGCAACTYAACCCMACTAGTTGHMS"
-- Base2  TARCGAATGATATTTTACAGGGTCCTGAMTTA"
-- Base3  GRANYYRYNYHRHRAGYNRNRYGNNNGYYAHR"
}
{
name "Ascidian Mitochondrial" ,
id 13 ,
ncbieaa  "XX*ACDEFGGHIKLLMNPQRSSTVWYBJZ",
sncbieaa "MM---------------------------",
-- Base1  ADTGTGGTGRCAACYAACCCATAGTTGMS"
-- Base2  TTACGAATGGATATTTACAGGCCTGAMTA"
-- Base3  RGRNYYRYNRYYRNRRYNRNYNNNRYYYR"
}
{
name "Alternative Flatworm Mitochondrial" ,
id 14 ,
ncbieaa  "X*ACDEFGHIKLLMNPQRSSTVWYBJJZ",
sncbieaa "M---------------------------",
-- Base1  ATGTGGTGCAACYAACCCATAGTTGHMS"
-- Base2  TACGAATGATATTTACAGGCCTGAMTTA"
-- Base3  GGNYYRYNYHGNRGHNRNNNNNRHYAHR"
}
{
name "Blepharisma Macronuclear" ,
id 15 ,
ncbieaa  "X*ACDEFGHIKLLMNPQQRRSSTVWYBJJZZ",
sncbieaa "M------------------------------",
-- Base1  ATGTGGTGCAACYAACCYCMATAGTTGHMBS"
-- Base2  TRCGAATGATATTTACAAGGGCCTGAMTTAA"
-- Base3  GANYYRYNYHRNRGYNRGNRYNNNGYYAHGR"
}
{
name "Chlorophycean Mitochondrial" ,
id 16 ,
ncbieaa  "X*ACDEFGHIKLLLMNPQRRSSTVWYBJJZ",
sncbieaa "M-----------------------------",
-- Base1  ATGTGGTGCAACTYAACCCMATAGTTGHMS"
-- Base2  TRCGAATGATATWTTACAGGGCCTGAMTTA"
-- Base3  GANYYRYNYHRNGRGYNRNRYNNNGYYAHR"
}
{
name "Trematode Mitochondrial" ,
id 21 ,
ncbieaa  "X*ACDEFGHIKLLMNPQRSSTVWYBJZ",
sncbieaa "M--------------------------",
-- Base1  RTGTGGTGCAACYAACCCATAGTTGMS"
-- Base2  TACGAATGATATTTACAGGCCTGAMTA"
-- Base3  GRNYYRYNYYGNRRHNRNNNNNRYYYR"
}
{
name "Scenedesmus obliquus Mitochondrial" ,
id 22 ,
ncbieaa  "X*ACDEFGHIKLLLMNPQRRSSTVWYBJJZ",
sncbieaa "M-----------------------------",
-- Base1  ATGTGGTGCAACTYAACCCMATAGTTGHMS"
-- Base2  TVCGAATGATATWTTACAGGGCCTGAMTTA"
-- Base3  GANYYRYNYHRNGRGYNRNRYBNNGYYAHR"
}
{
name "Thraustochytrium Mitochondrial" ,
id 23 ,
ncbieaa  "XX**ACDEFGHIKLLMNPQRRSSTVWYBJZ",
sncbieaa "MM----------------------------",
-- Base1  ARTTGTGGTGCAACYAACCCMATAGTTGMS"
-- Base2  TTADCGAATGATATTTACAGGGCCTGAMTA"
-- Base3  KGRANYYRYNYHRNGGYNRNRYNNNGYYHR"
}