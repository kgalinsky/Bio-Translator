# Bio::Tiny::Translator::Validations
#
# $Author$
# $Date$
# $Revision$
# $HeadURL$

package Bio::Tiny::Translator::Validations;

use strict;
use warnings;

use Log::Log4perl qw(:easy);
use Carp;
use Params::Validate;
use Exporter qw(import);

our %EXPORT_TAGS = (
    defaults => [
        qw(
          $DEFAULT_STRAND
          $DEFAULT_PARTIAL
          $DEFAULT_SANITIZED
          $DEFAULT_OFFSET
          )
    ],
    regexes => [
        qw(
          $REGEX_BOOLEAN
          $REGEX_NON_NEG_INT
          $REGEX_STRAND
          $REGEX_SEARCH_STRAND
          $REGEX_012
          )
    ],
    validations => [
        qw(
          %VAL_STRAND
          %VAL_SEARCH_STRAND
          %VAL_PARTIAL
          %VAL_SANITIZED
          %VAL_OFFSET
          validate_seq_params
          validate_seq_exons_params
          validate_sanitized
          validate_lower_upper
          validate_exons
          )
    ]
);

our @EXPORT_OK = map { @$_ } values %EXPORT_TAGS;

=head1 DEFAULTS

=cut

our $DEFAULT_STRAND        = 1;
our $DEFAULT_SEARCH_STRAND = 0;
our $DEFAULT_PARTIAL       = 0;
our $DEFAULT_SANITIZED     = 0;
our $DEFAULT_OFFSET        = 0;

=head1 REGULAR EXPRESSIONS

=cut

our $REGEX_BOOLEAN       = qr/^[01]$/;
our $REGEX_NON_NEG_INT   = qr/^\d+$/;
our $REGEX_STRAND        = qr/^[+-]?1$/;
our $REGEX_SEARCH_STRAND = qr/^[+-]?[01]$/;
our $REGEX_012           = qr/^[012]$/;

=head1 VALIDATIONS

=cut

# Make sure strand is 1 or -1 and set default
our %VAL_STRAND = (
    strand => {
        default => $DEFAULT_STRAND,
        regex   => $REGEX_STRAND,
        type    => Params::Validate::SCALAR
    }
);

# Make sure strand is 0, 1 or -1 and set default
our %VAL_SEARCH_STRAND = (
    strand => {
        default => $DEFAULT_SEARCH_STRAND,
        regex   => $REGEX_SEARCH_STRAND,
        type    => Params::Validate::SCALAR
    }
);

# Make sure partial is boolean and set default
our %VAL_PARTIAL = (
    partial => {
        default => $DEFAULT_PARTIAL,
        regex   => $REGEX_BOOLEAN,
        type    => Params::Validate::SCALAR
    }
);

# Make sure sanitized is boolean and set default
our %VAL_SANITIZED = (
    sanitized => {
        default => $DEFAULT_SANITIZED,
        regex   => $REGEX_BOOLEAN,
        type    => Params::Validate::SCALAR
    }
);

# Make sure offset is 0, 1 or 2 and set default
our %VAL_OFFSET = (
    offset => {
        default => $DEFAULT_OFFSET,
        regex   => $REGEX_012,
        type    => Params::Validate::SCALAR
    }
);

=head1 VALIDATION METHODS

=cut

=head2 validate_seq_params

    my ( $seq_ref, @p ) = validate_seq_params(@_);

Do validations for methods expecting to be called as:

    method( \$sequence, \%params );

=cut

sub validate_seq_params (\@) {
    validate_pos(
        @{ $_[0] },
        { type => Params::Validate::SCALARREF },
        { type => Params::Validate::HASHREF, default => {} }
    );
}

=head2 validate_seq_exons_params

    my ( $seq_ref, $exons, @p ) = validate_seq_exons_params(@_);

Do validations for methods expecting to be called as:

    method( \$sequence, \@exons \%params );

=cut

sub validate_seq_exons_params (\@) {
    validate_pos(
        @{ $_[0] },
        { type => Params::Validate::SCALARREF },
        { type => Params::Validate::ARRAYREF },
        { type => Params::Validate::HASHREF, default => {} }
    );
}

=head2 validate_sanitized

    my $sanitized = validate_sanitized(@p);

Validate sanitized from the validation hash and delete it so it doesn't
interfere later.

=cut

sub validate_sanitized {
    my $p = $_[0];

    return $DEFAULT_SANITIZED unless ( exists $p->{sanitized} );

    my $s = $p->{sanitized};
    croak qq{Invalid value for sanitized "$s" (must be 0 or 1) }
      unless ( $s =~ m/^[01]$/ );
    delete $p->{sanitized};

    return $s;
}

=head2 validate_lower_upper

    my ( $lower, $upper ) = validate_lower_upper( $seq_ref, @p )

Validate lower and upper bounds from the validation hash and delete them.

=cut

sub validate_lower_upper {
    my ( $seq_ref, $p ) = @_;

    my $upper;
    if ( exists $p->{upper} ) {
        $upper = $p->{upper};
        delete $p->{upper};

        croak 'upper bound is not a non-negative integer'
          unless ( $upper =~ m/$REGEX_NON_NEG_INT/ );

        croak 'upper bound is out range'
          if ( $upper > length($$seq_ref) );
    }
    else { $upper = length($$seq_ref) }

    my $lower;
    if ( exists $p->{lower} ) {
        $lower = $p->{lower};
        delete $p->{lower};

        croak 'lower bound is not a non-negative integer'
          unless ( $lower =~ m/$REGEX_NON_NEG_INT/ );

        croak 'lower bound is greater than upper bound'
          if ( $lower > $upper );
    }
    else { $lower = 0 }

    return ( $lower, $upper );
}

=head2 validate_exons

    validate_exons( $seq_ref, $exons );

=cut

sub validate_exons {
    my ( $seq_ref, $exons ) = @_;

    foreach my $exon (@$exons) {
        my ( $lower, $upper ) = @$exon;

        # Make sure upper and lower bounds are integers
        get_logger()->logcroak("Lower $lower not an integer")
          if ( $lower !~ m/^\d+$/ );
        get_logger()->logcroak("Upper $upper not an integer")
          if ( $upper !~ m/^\d+$/ );

        # Make sure upper >= lower
        get_logger()->logcroak("Upper $upper < Lower $lower")
          if ( $upper < $lower );

        # Make sure upper is within the sequence
        get_logger()->logcroak("Upper $upper not in the sequence")
          if ( $upper > length($$seq_ref) );
    }
}

1;
