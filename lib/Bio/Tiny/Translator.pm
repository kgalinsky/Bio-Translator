package Bio::Tiny::Translator;

use strict;
use warnings;

use version; our $VERSION = qv('0.6.0');

=head1 NAME

Bio::Tiny::Translator - Translate DNA sequences

=head1 SYNOPSIS

    use Bio::Tiny::Translator;

    my $translator = new Bio::Tiny::Translator();
    my $translator = new Bio::Tiny::Translator(11);
    my $translator = new Bio::Tiny::Translator( 12, 'id' );
    my $translator = new Bio::Tiny::Translator( 'Yeast Mitochondrial', 'name' );
    my $translator = new Bio::Tiny::Translator( 'mito', 'name' );

    my $translator = custom Bio::Tiny::Translator( \$custom_table );
    my $translator = custom Bio::Tiny::Translator( \$custom_table, 1 );

    $translator->translate( \$seq );
    $translator->translate( \$seq, { strand => 1 } );
    $translator->translate( \$seq, { strand => -1 } );

=head1 DESCRIPTION

C<Bio::Tiny::Translator> tries to be a robust translator object featuring
translation tables based off the the ones provided by
L<NCBI|http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>.
Key features include the ability to handle degenerate nucleotides and to
translate to ambiguous amino acids.

First, create a new translator object using one of the included tables or a
custom one (see C<Bio::Tiny::Translator::Table> for table formats), and then
passing your DNA sequences to your translator object.

The translator uses interbase coordinates. See below for the difference between
interbase coordinates and traditional numbering methods:

    Traditional   1 2 3 4
                  A C G T ...
    Interbase    0 1 2 3 4

Conversion methods between the two methods can depend upon what you are trying
to do, but the simple way to do this is:

    strand = 3' end <=> 5' end # that's the spaceship operator!
    lower  = min( 5' end, 3' end ) - 1
    upper  = max( 5' end, 3' end )

Parameter validation uses L<Params::Validate> which introduces overhead but can
be disabled. See the C<Params::Validate> documentation for more information.

=cut

use base 'Class::Accessor::Faster';
__PACKAGE__->mk_accessors('table');

use Carp;
use Params::Validate;

use Bio::Tiny::Translator::Table;

#use Bio::Tiny::Translator::Base;
#use Bio::Tiny::Translator::Validations qw(:validations);

use Bio::Tiny::Util::DNA qw/ $all_nucleotide_match cleanDNA /;

=head1 CONSTRUCTORS

=cut

sub _new {
    shift->SUPER::new( { table => shift } );
}

=head2 new

    my $translator = new Bio::Tiny::Translator();
    my $translator = new Bio::Tiny::Translator( $id );
    my $translator = new Bio::Tiny::Translator( $id, \%params );

Create a translator with a translation table provided by $id. Please see
Bio::Tiny::Translator::Table for the full list of options.

=cut

sub new {
    my $class = shift;
    my $table = Bio::Tiny::Translator::Table->new(@_) or return;
    $class->_new($table);
}

=head2 custom()

    my $translator = $translator->custom( $table_ref );
    my $translator = $translator->custom( $table_ref, \%params );

Create a translator with a custom translation table. Please see
Bio::Tiny::Translator::Table for the full list of options.

=cut

sub custom {
    my $class = shift;
    my $table = Bio::Tiny::Translator::Table->custom(@_) or return;
    $class->_new($table);
}

=head1 METHODS

=cut

=head2 translate

    $pep_ref = $translator->translate( $seq_ref, \%params );

The basic function of this module. Translate the specified region of the
sequence (passed as $seq_ref) and return a reference to the translated string.
The parameters are:

    strand: [+-]?1; default = 1
    lower:  integer between 0 and seq_length; default = 0
    upper:  integer between 0 and seq_length; default = seq_length
    start:  boolean
    offset: [012]

Translator uses interbase coordinates. "lower" and "upper" are optional
parameters such that:

    0 <= lower <= upper <= seq_length

Translator will croak if those conditions are not satisfied.

"start" sets whether or not to try translating the first codon as a start
codon. By default, translator will try to do this. "offset" allows you to
specify an offset in addition to the lower and upper abounds and have
Translator figure out the correct bound to offset from.

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
            start   => 0
        }
    );

=cut

sub translate {
    my $self = shift;

    my ( $seq_ref, @p ) = validate_pos(
        @_,
        { type    => Params::Validate::SCALARREF | Params::Validate::SCALAR },
        { default => {}, type => Params::Validate::HASHREF }
    );

    $seq_ref = \$seq_ref unless ( ref $seq_ref );

    # perform basic parameter validation
    my %p = validate(
        @p,
        {
            lower => {
                default   => 0,
                type      => Params::Validate::SCALAR,
                regex     => qr/^\d+$/,
                callbacks => {
                    '0 <= lower' => sub { 0 <= $_[0] }
                },
            },
            upper => {
                default   => length($$seq_ref),
                type      => Params::Validate::SCALAR,
                regex     => qr/^\d+$/,
                callbacks => {
                    'upper <= seq_length' => sub { $_[0] <= length($$seq_ref) }
                },
            },
            strand => {
                default => 1,
                type    => Params::Validate::SCALAR,
                regex   => qr/^[+-]?1$/,
            },
            start => {
                default => 1,
                type    => Params::Validate::SCALAR,
            },
            offset => {
                default => 0,
                type    => Params::Validate::SCALAR,
                regex   => qr/^[012]$/,
            }
        }
    );

    my ( $lower, $upper, $strand, $start, $offset ) =
      @p{qw/ lower upper strand start offset /};

    # do we need to do reverse_complement?
    my $rc = $strand == -1 ? 1 : 0;

    # do some final checks on lower upper bound and adjust lower bound
    croak 'lower > upper' if ( $lower > $upper );
    $lower += $rc ? ( $upper - $offset - $lower ) % 3 : $offset;

    # get a list of codon start locations
    my @codon_starts =
      map { $lower + 3 * $_ } ( 0 .. ( int( ( $upper - $lower ) / 3 ) - 1 ) );
    @codon_starts = reverse @codon_starts if ($rc);

    return \'' unless (@codon_starts);

    # try to translate the start codon
    my @start_peptide;
    if ($start) {
        my $start_aa =
          $self->table->codon2start->[$rc]
          { substr( $$seq_ref, $codon_starts[0], 3 ) };
        if ($start_aa) { push @start_peptide, $start_aa; shift @codon_starts }
    }

    # translate the rest of the peptide
    my $codon2aa = $self->table->codon2aa->[$rc];
    my $peptide = join '', @start_peptide, map { $_ || 'X' }
      @$codon2aa{ map { substr $$seq_ref, $_, 3 } @codon_starts };

    return \$peptide;
}

=head2 translate_codon

    my $residue = $translator->translate_codon( $codon );
    my $residue = $translator->translate_codon( $codon, \%params );

Translate a codon. Return 'X' or '-' if it isn't in the
codon table. Handles degenerate nucleotides, so if all
possible codons for an ambiguity map to the same residue,
return that residue.

Example:

    $residue = $translator->translate_codon('atg');
    $residue = $translator->translate_codon( 'tty', { strand => -1 } );
    $residue = $translator->translate_codon( 'cat', { start => 1 } );

=cut

sub translate_codon {
    my $self = shift;

    my ( $codon, @p ) = validate_pos(
        @_,
        { regex => qr/^${all_nucleotide_match}{3}$/ },
        { type  => Params::Validate::HASHREF, default => {} }

    );

    my %p = validate(
        @p,
        {

            # Make sure strand is 1 or -1
            strand => {
                default => 1,
                regex   => qr/^[+-]?1$/,
                type    => Params::Validate::SCALAR
            },

            # Make sure it is a boolean value
            start => {
                default => 0,
                regex   => qr/^[01]$/,
                type    => Params::Validate::SCALAR
            }
        }
    );

    $codon = uc $codon;

    # Set up the translation table given the strand and whether this is
    # searching for stop codons. Set up the not_found string by whether this
    # is a start or not.
    my $rc = $p{strand} == 1 ? 0 : 1;
    my ( $table, $not_found );
    unless ( $p{start} ) {
        $table     = $self->table->codon2aa->[$rc];
        $not_found = 'X';
    }
    else {
        $table     = $self->table->codon2start->[$rc];
        $not_found = '-';
    }

    return $self->table->_unroll( $codon, $table, { start => $p{start} } )
      || $not_found;
}

1;

=head1 AUTHOR

Kevin Galinsky, C<kgalinsky plus cpan at gmail dot com>

=head1 BUGS

Please report any bugs or feature requests to
C<bug-bio-tiny-translator at rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Bio-Tiny-Translator>.
I will be notified, and then you'll automatically be notified of progress on
your bug as I make changes.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Bio::Tiny::Translator

You can also look for information at:

=over 4

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Bio-Tiny-Translator>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/Bio-Tiny-Translator>

=item * RT: CPAN's request tracker

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Bio-Tiny-Translator>

=item * Search CPAN

L<http://search.cpan.org/dist/Bio-Tiny-Translator>

=back

=head1 ACKNOWLEDGEMENTS

JCVI/Paolo Amedeo

=head1 COPYRIGHT & LICENSE

Copyright 2008-2009 J. Craig Venter Institute, 2011 Kevin Galinsky.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut
