# JCVI::Translator:Utils
#
# $Author$
# $Date$
# $Revision$
# $HeadURL$

=head1 NAME

package Utils

=head1 SYNOPSES

 use JCVI::Translator::Utils;

 my $translator = new JCVI::Translator::Utils(
                           id => $id,
                           name => $name,
                           table => $table,
                           tableRef => $tableRef
                           );

=head1 DESCRIPTION

See Translator for more info. Utils extends Translator and
adds a few more functions that are normally not used.

=head1 AUTHOR

Kevin Galinsky, <kgalinsk@jcvi.org>

=head1 FUNCTIONS

=over

=cut

package JCVI::Translator::Utils;
use base JCVI::Translator;

use strict;
use warnings;

our $VERSION = '0.2.2';

use Log::Log4perl qw(:easy);
use Params::Validate qw(:all);

=item getORF()

=item [$start, $stop] = $translator->getORF($seqRef, $strand);

This will get the longest region between stops and return
the strand, lower and upper bounds, inclusive:

 0 1 2 3 4 5 6 7 8 9 10
  T A A A T C T A A G
  *****       *****
        <--------->

Will return [1, 3, 9]. You can also specify which strand
you are looking for the ORF to be on.

For ORFs starting at the very beginning of the strand or
trailing off the end, but not in phase with the start or
ends, this method will cut at the last complete codon.

 Eg:

 0 1 2 3 4 5 6 7 8 9 10
  A C G T A G T T T A
                *****
    <--------->

Will return [-1, 1, 7]. The distance between lower and
upper will always be a multiple of 3. This is to make it
clear which frame the ORF is in.

Example:

 $ref = $translator->getORF(\'TAGAAATAG');

Output:

 $ref = [$strand, $lower, $upper]

=cut

sub getORF {
    my $self = shift;

    my ( $seqRef, $strand )
        = validate_pos( @_,
                        { type => SCALARREF },
                        { default => 0,
                          regex   => qr/^[+-]?[01]$/
                        }
        );

    DEBUG('getORF called');

    my ( $best_strand, $lower, $upper ) = ( 1, 0, 0 );

    foreach my $cur_strand ( $strand == 0 ? ( -1, 1 ) : ($strand) ) {
        my @lowers = ( 0 .. 2 );
        my $stopRegex = $self->regex( 'stop', $cur_strand );

        ########################################
        # Rather than using a regular expression
        # to find regions between stops, it
        # should be  more computationally
        # efficient to find all the stops and
        # compute from there. However, Perl's
        # regular expression engine may be
        # faster than code execution, so this
        # may not be the case

        ########################################
        # A lookahead is used for two reasons:
        # the main one is to get every position
        # within two bases of the end of the
        # sequence, and also to cope with the
        # possibility of overlapping stop
        # codons.

        while (
            $$seqRef =~ /(?=
				($stopRegex)|.{0,2}$
			    )/gx
            )
        {
            my $curUpper = pos $$seqRef;
            my $frame    = $curUpper % 3;

            $curUpper += length $1 if ( $1 && ( $cur_strand == 1 ) );

            ########################################
            # If the current distance between start
            # and stop is greater than the distance
            # between the stored start and stop,
            # change the stored start and stop to be
            # the current one.

            if ( $upper - $lower < $curUpper - $lowers[$frame] ) {
                $best_strand = $cur_strand;
                $lower       = $lowers[$frame];
                $upper       = $curUpper;
            }

            $lowers[$frame] = $curUpper;
        }

    }

    return [ $best_strand, $lower, $upper ];
}

=item getCDS()

=item [$start, $stop] = $translator->getCDS($seqRef, $strand, $strict);

This will return the strand and boundaries of the longest
CDS.

 0 1 2 3 4 5 6 7 8 9 10
  A T G A A A T A A G
  >>>>>       *****
  <--------------->

Will return [1, 0, 9].

Strict controls how strictly getCDS functions.
There are 3 levels of strictness, enumerated 0, 1 and 2. 2 is the most strict,
and in that mode, a region will only be considered a CDS if both the start and
stop is found. In strict level 1, if a start is found, but no stop is present
before the end of the sequence, the CDS will run until the end of the sequence.
Strict level 0 assumes that start codon is present in each frame just before
the start of the molecule.

Example:

 $ref = $translator->getCDSs(\'ATGAAATAG');
 $ref = $translator->getCDSs(\'ATGAAATAG', -1);

Output:

 $ref = [$strand, $lower, $upper]

=cut

sub getCDS {
    my $self = shift;

    my ( $seqRef, $strand, $strict )
        = validate_pos( @_,
                        { type => SCALARREF },
                        { default => 0,
                          regex   => qr/^[+-]?[01]$/
                        },
                        { default => 1,
                          regex   => qr/^[012]$/
                        }
        );

    DEBUG('getCDS called');

    my ( $best_strand, $lower, $upper ) = ( 1, 0, 0 );

    foreach my $cur_strand ( $strand == 0 ? ( -1, 1 ) : ($strand) ) {
        my $lowerRegex = $self->regex( 'lower', $cur_strand );
        my $upperRegex = $self->regex( 'upper', $cur_strand );

        ########################################
        # Initialize
        my @lowers;
        if ( $cur_strand == 1 ) {
            @lowers = ( $strict != 0 ? map {undef} ( 0 .. 2 ) : ( 0 .. 2 ) );
        }
        else {
            @lowers = ( $strict == 2
                        ? map {undef} ( 0 .. 2 )
                        : ( 0 .. 2 )
            );
        }

        ########################################
        # Similar to getORF, rather than
        # using a regular expression to find
        # entire regions, instead find
        # individual starts and stops and react
        # accordingly. It captures the starts
        # and stops separately ($1 vs $2) so
        # that it is easy to tell if a start or
        # a stop was matched.
        #
        # If strict mode is at level 2, we don't
        # tes is a newly added feature which t CDSs trailing off the end
        # of the molecule

        my $regex = qr/(?=($lowerRegex)|($upperRegex))/;
        $regex = qr/$regex|(?=.{0,2}$)/ unless ( $strict == 2 );

        while ( $$seqRef =~ /$regex/g ) {

            my $position = pos $$seqRef;
            my $frame    = $position % 3;

            ########################################
            # If we match the lower regex we:
            #
            # In the case that we are on the '-'
            # strand, that means we found a stop,
            # so we update the lower bound.
            #
            # Otherwise, we are on the positive
            # strand, meaning we have found the
            # start, so only set the lower bound if
            # it is not already set (don't want to
            # overwrite the location of a previous
            # start codon).

            if ( $1
                 && (    ( $cur_strand eq '-' )
                      || ( !defined $lowers[$frame] ) )
                )
            {
                $lowers[$frame] = $position;
            }

            ########################################
            # If we don't match the lower regex:
            #
            # If this is the positive strand, that
            # means that this is a valid stop -
            # either a stop codon or the end of the
            # string. Reset the lower bound in this
            # case.
            #
            # On the negative strand, we only care
            # if we matched a start. In that case,
            # do the compute and update.
            #
            # Another option would be to mark where
            # the start is, and only do the compute
            # when we find a stop.

            elsif ( ( $cur_strand == 1 ) || $2 ) {

                # Move on if the lower is unset
                next unless ( defined $lowers[$frame] );

                $position += length $2 if ($2);

                if ( $upper - $lower < $position - $lowers[$frame] ) {
                    $best_strand = $cur_strand;
                    $lower       = $lowers[$frame];
                    $upper       = $position;
                }

                # Reset lower if we found a stop
                undef $lowers[$frame] if ( $cur_strand == 1 );
            }
        }
    }

    return [ $best_strand, $lower, $upper ];
}

=item find()

=item $positions = $translator->find( $seqRef, $type, $strand )

Find codons of given type (i.e. start or stop) in the sequence. Note, the
second two parameters are passed directly to regex, so please look there for
more information. Returns an arrayref containing all the locations of all those
codons.

=cut

sub find {
    my $self = shift;

    my @seqRef = splice @_, 0, 1;
    my $seqRef = validate_pos( @seqRef, { type => SCALARREF } );

    my $regex = $self->regex(@_);

    my @positions;

    while ( $$seqRef =~ /(?=($regex))/g ) {
        push @positions, pos $$seqRef;
    }

    return \@positions;
}

=item regex()

=item $regex = $translator->regex( $type, $strand )

Returns a regular expression of a certain type for a given strand. These are
'start', 'stop', 'lower' and 'upper.' Lower and upper match the lower and upper
end of a CDS for a given strand (i.e. on the positive strand, lower matches
the start, and upper matches the stop).

=cut

sub regex {
    my $self = shift;
    my ( $type, $strand )
        = validate_pos( @_,
                        { regex => qr/^(?:start|stop|lower|upper)$/ },
                        { default => 1,
                          regex   => qr/^[+-]?1$/
                        }
        );

    my $prefix = $strand == 1 ? '' : 'rc_';

    if    ( $type eq 'lower' ) { $type = $strand == 1  ? 'start' : 'stop' }
    elsif ( $type eq 'upper' ) { $type = $strand == -1 ? 'start' : 'stop' }

    unless ( defined $self->{"${prefix}${type}Regex"} ) {
        my $regex = join '|',
            ( $type eq 'start'
              ? keys %{ $self->{"${prefix}starts"} }
              : @{ $$self{"${prefix}reverse"}{'*'} }
            );
        $self->{"${prefix}${type}Regex"} = qr/$regex/;
    }

    return $self->{"${prefix}${type}Regex"};
}

=item nonstop

=item $frames = $translator->nonstop( $seqRef, $strand )

Returns the frames that contain no stop codons for the sequence. $strand is
optional and defaults to 0. Frames are 1, 2, 3, -1, -2, -3.

 3    ---->
 2   ----->
 1  ------>
    -------
 -1 <------
 -2 <-----
 -3 <----

Example:

 $frames = $translator->nonstop(\'TACGTTGGTTAAGTT');     # [-1, -3, 2, 3]
 $frames = $translator->nonstop(\'TACGTTGGTTAAGTT', 1);  # [2, 3]
 $frames = $translator->nonstop(\'TACGTTGGTTAAGTT', -1); # [-1, -3]

=cut

sub nonstop {
    my $self = shift;
    my ( $seqRef, $strand )
        = validate_pos( @_,
                        { type => SCALARREF },
                        { default => 0,
                          regex   => qr/^[+-]?1$/
                        }
        );

    my @frames;
    foreach my $cur_strand ( $strand == 0 ? ( -1, 1 ) : ($strand) ) {
        my $stop = $self->regex( 'stop', $cur_strand );

        foreach my $frame ( 0 .. 2 ) {
            my $regex = $cur_strand == 1
                ? qr/^.{$frame}(?:.{3})*$stop/
                : qr/$stop(?:.{3})*.{$frame}$/;

            push @frames, ( $frame + 1 ) * $cur_strand
                unless ( $$seqRef =~ m/$regex/ );
        }
    }
    
    return \@frames;
}

1;
