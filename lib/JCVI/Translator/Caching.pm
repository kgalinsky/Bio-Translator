package JCVI::Translator::Caching;

use strict;
use warnings;

use base qw(JCVI::Translator);

use Carp;

our $DEFAULT_PARTIAL = 0;

=head2 load_sequence

=cut

sub load_sequence {
    TRACE('load_sequence called');

    my $self = shift;

    my ( $seq_ref, @p );
    ( $seq_ref, $p[0] ) = validate_pos(
        @_,
        { type => Params::Validate::SCALARREF, },
        { type => Params::Validate::HASHREF, default => {} }
    );

    # Check the sanitized value separately
    my $sanitized = $self->_is_sanitized(@p);

    # Clean the sequence and cache it
    $seq_ref = cleanDNA($seq_ref) unless ($sanitized);
    $self->base->set_seq($seq_ref);

    TRACE('load_sequence exiting');
}

=head2 translate

=cut

sub translate {
    TRACE('translate called');

    my $self = shift;

    my $seq_ref = $self->base->{seq_ref};
    croak 'No sequence loaded' unless ( defined $seq_ref );

    my %p = validate(
        @_,

        # Make sure lower is an integer within the sequence
        lower => {
            regex     => qr/^\d+$/,
            type      => Params::Validate::SCALAR,
            callbacks => {
                'lower >= 0'          => sub { $_[0] >= 0 },
                'lower <= seq_length' => sub { $_[0] <= length($$seq_ref) }
            }
        },

        # Make sure upper is an integer within the sequence
        upper => {
            regex     => qr/^[0-9]+$/,
            type      => Params::Validate::SCALAR,
            callbacks => {
                'upper >= 0'          => sub { $_[0] >= 0 },
                'upper <= seq_length' => sub { $_[0] <= length($$seq_ref) }
            }
        },

        # Make sure the offset is 0, 1 or 2.
        offset => {
            default => 0,
            regex   => qr/^[012]$/,
            type    => Params::Validate::SCALAR
        },

        # Make sure they are boolean values
        partial => {
            default => $DEFAULT_PARTIAL,
            regex   => qr/^[01]$/,
            type    => Params::Validate::SCALAR
        }
    );

    TRACE('translate exiting');
}

1;
