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

use version; our $VERSION = qv('0.5.7');

use base qw(Class::Accessor::Fast);
__PACKAGE__->mk_accessors(qw(table base));

use Carp;
use Log::Log4perl qw(:easy);
use Params::Validate;

use JCVI::Translator::Table;
use JCVI::Translator::Base;
use JCVI::Translator::Validations qw(:validations);

use JCVI::DNATools qw( $nuc_match cleanDNA );

=head1 CONSTRUCTORS

=cut

sub _new {
    my $class = shift;
    $class->SUPER::new(
        { table => shift, base => JCVI::Translator::Base->new() } );
}

=head2 new

    my $translator = new JCVI::Translator();
    my $translator = new JCVI::Translator( $id );
    my $translator = new JCVI::Translator( $id, \%params );

Create a translator with a translation table provided by $id. Please see
JCVI::Translator::Table for the full list of options.

=cut

sub new {
    TRACE('new called');

    my $class = shift;
    my $table = JCVI::Translator::Table->new(@_);

    return undef unless ($table);

    $class->_new($table);
}

=head2 custom()

    my $translator = $translator->custom( $table_ref );
    my $translator = $translator->custom( $table_ref, \%params );

Create a translator with a custom translation table. Please see
JCVI::Translator::Table for the full list of options.

=cut

sub custom {
    TRACE('custom called');

    my $class = shift;
    my $table = JCVI::Translator::Table->custom(@_);

    return undef unless ($table);

    $class->_new($table);
}

=head1 METHODS

=cut

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

    my ( $seq_ref, @p ) = validate_seq_params(@_);

    # Check the sanitized value separately
    my $sanitized = validate_sanitized(@p);

    # Clean the sequence and cache it
    $seq_ref = cleanDNA($seq_ref) unless ($sanitized);
    $self->base->set_seq($seq_ref);

    # Get lower and upper bounds
    my ( $lower, $upper ) = validate_lower_upper( $seq_ref, @p );

    my %p = validate( @p, { %VAL_STRAND, %VAL_PARTIAL, %VAL_OFFSET } );

    # Return undef if the offset is bigger than the space between bounds
    return undef if ( $upper <= $lower + $p{offset} );

    # Set the partial status
    $self->base->set_partial( $p{partial} );

    # Prepare for translation
    $self->base->prepare( $p{strand}, $self->table );
    $self->base->endpoints( $lower, $upper, $p{offset} );

    # Translate and convert the resulting arrayref to a string
    my $peptide = join( '', @{ $self->base->translate() } );
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

    my ( $seq_ref, @p ) = validate_seq_params(@_);

    # Check the sanitized value separately
    my $sanitized = validate_sanitized(@p);

    # Clean the sequence and cache it
    $seq_ref = cleanDNA($seq_ref) unless ($sanitized);
    $self->base->set_seq($seq_ref);

    # Get lower and upper bounds
    my @lu = validate_lower_upper( $seq_ref, @p );

    my %p = validate( @p, {%VAL_PARTIAL} );

    $self->base->set_partial( $p{partial} );

    my @peptides;

    foreach my $strand ( 1, -1 ) {

        # We only need to prepare once for a given strand
        $self->base->prepare( $strand, $self->table );

        foreach my $offset ( 0 .. 2 ) {

            # Calculate endpoints and translate
            $self->base->endpoints( @lu, $offset );

            # Translate and push onto array
            push @peptides, join( '', @{ $self->base->translate() } );
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

    my ( $seq_ref, $exons, @p ) = validate_seq_exons_params(@_);

    # Check the sanitized value separately
    my $sanitized = validate_sanitized(@p);

    # Clean the sequence and cache it
    $seq_ref = cleanDNA($seq_ref) unless ($sanitized);
    $self->base->set_seq($seq_ref);

    # Validate optional arguments
    my %p = validate( @p, { %VAL_STRAND, %VAL_PARTIAL, %VAL_SANITIZED } );

    # Validate the exons and sort in the proper order for translation
    validate_exons( $seq_ref, $exons );
    my @exons = sort {
        ( ( $a->[0] <=> $b->[0] ) || ( $a->[1] <=> $b->[1] ) ) * $p{strand}
    } @$exons;

    $self->base->set_partial( $p{partial} );

    my $prep = $self->base->prepare( $p{strand}, $self->table );
    my @peptides;

  EXON: foreach my $exon (@exons) {
        $self->base->endpoints(@$exon);

        next EXON unless ( $self->base->finish_leftover() );

        push @peptides, [ $self->base->translate_leftover() ];
        push @peptides, $self->base->translate();
        $self->base->store_leftover();
    }

    my $peptide = join( '', map { @$_ } @peptides );

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
    TRACE("translate_codon called");

    my $self = shift;

    my ( $codon, @p ) = validate_pos(
        @_,
        { regex => qr/^${nuc_match}{3}$/ },
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
        $table     = $self->table->_forward->[$rc];
        $not_found = 'X';
    }
    else {
        $table     = $self->table->_starts->[$rc];
        $not_found = '-';
    }

    return $self->table->_unroll( $codon, $table, { start => $p{start} } )
      || $not_found;
}

1;

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
