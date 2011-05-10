# Bio::Tiny::Translator:Utils
#
# $Author$
# $Date$
# $Revision$
# $HeadURL$

=head1 NAME

Bio::Tiny::Translator::Utils - Utilities that requrie a translation table

=head1 SYNOPSIS

    use Bio::Tiny::Translator::Utils;

    # Same constructor as Bio::Tiny::Translator
    my $utils = new Bio::Tiny::Translator::Utils();
    my $utils = custom Bio::Tiny::Translator( \$custom_table );

    my $codons = $utils->codons( $residue );
    my $regex  = $utils->regex( $residue );
    my $indices = $utils->find( $residue );

    my $orf = $utils->getORF( $seq_ref );
    my $cds = $utils->getCDS( $seq_ref );

    my $frames = $utils->nonstop( $seq_ref );

=head1 DESCRIPTION

See Bio::Tiny::Translator for more info. Utils contains utilites that require
knowledge of the translation table.

=cut

package Bio::Tiny::Translator::Utils;

use strict;
use warnings;

use base qw(Bio::Tiny::Translator);
__PACKAGE__->mk_accessors(qw( _regexes ));

use Log::Log4perl qw(:easy);
use Params::Validate;

use Bio::Tiny::Translator::Validations qw(:validations :regexes);

use Bio::Tiny::Util::DNA qw( cleanDNA );
use Bio::Tiny::Util::AA qw( $aa_match );

# Default values
our $DEFAULT_STRAND        = 1;
our $DEFAULT_SEARCH_STRAND = 0;
our $DEFAULT_SANITIZED     = 0;

# Precompiled regular expressions for SPOT rule and to save time
our $BOOLEAN_REGEX       = qr/^[01]$/;
our $INTEGER_REGEX       = qr/^\d+$/;
our $STRAND_REGEX        = qr/^[+-]?1$/;
our $SEARCH_STRAND_REGEX = qr/^[+-]?[01]$/;
our $RESIDUE_REGEX       = qr/^(?:$aa_match|\+|start|stop|lower|upper)$/;
our $STRICT_REGEX        = qr/^[012]$/;

sub _new {
    my $self = shift->SUPER::_new(@_);
    $self->_regexes( [ {}, {} ] );
    return $self;
}

=head1 METHODS

=cut

=head2 codons

    my $codon_array = $translator->codons( $residue);
    my $codon_array = $translator->codons( $residue, \%params );

Returns a list of codons for a particular residue or start codon. In addition
to the one-letter codes for amino acids, the following are valid inputs for the
residue:

    start:  Start codons (you may also use "+" which is what the translator
            uses as the 1-letter code for start codons)
    stop:   Stop codons (you may also use "*" which is the 1-letter code)
    lower:  Start or stop codons, depending up on strand
    upper:  Start or stop codons, depending up on strand

"lower" and "upper" match the respective ends of a CDS for a given strand (i.e.
on the positive strand, lower matches the start, and upper matches them stop).
Valid options for the params hash are:

    strand:     1 or -1; default = 1

=cut

sub codons {
    TRACE 'entering codons';
    my $self = shift;

    # Get the residue and the optional validation hash
    my ( $residue, @p ) = validate_pos(
        @_,
        {
            type  => Params::Validate::SCALAR,
            regex => $RESIDUE_REGEX
        },
        { type => Params::Validate::HASHREF, default => {} }
    );

    my %p = validate( @p, {%VAL_STRAND} );

    # Set the reverse comlement variable
    my $rc = $p{strand} == 1 ? 0 : 1;

    # Format start/stop to be '+' and '*' which is how translator stores them
    if    ( $residue eq 'stop' )  { $residue = '*' }
    elsif ( $residue eq 'start' ) { $residue = '+' }

    # Lower bound is stop on the - strand, start on the + strand. Upper bound
    # is the reverse
    elsif ( $residue eq 'lower' ) { $residue = $rc ? '*' : '+' }
    elsif ( $residue eq 'upper' ) { $residue = $rc ? '+' : '*' }

    # Capitalize all other residues
    else { $residue = uc $residue }

    # Get the codons array or set it to the empty array
    my $codons = $self->table->aa2codons->[$rc]->{$residue} || [];

    # Return a copy of the arrayref so that the internal array can't get
    # modified
    TRACE 'leaving codons';
    return [@$codons];
}

=head2 regex

    my $regex = $translator->regex( $residue );
    my $regex = $translator->regex( $residue, \%params );

Returns a regular expression matching codons for a particular amino acid
residue. In addition to the one-letter codes for amino acids, the following are
valid inputs for the residue:

    start:  Start codons (you may also use "+" which is what the translator
            uses as the 1-letter code for start codons)
    stop:   Stop codons (you may also use "*" which is the 1 letter code)
    lower:  Start or stop codons, depending up on strand
    upper:  Start or stop codons, depending up on strand

"lower" and "upper" match the respective ends of a CDS for a given strand (i.e.
on the positive strand, lower matches the start, and upper matches the stop).
Valid options for the params hash are:

    strand: 1 or -1; default = 1

=cut

sub regex {
    TRACE 'entering regex';
    my $self = shift;

    my ( $residue, @p ) = validate_pos(
        @_,
        { type => Params::Validate::SCALAR },
        { type => Params::Validate::HASHREF, default => {} }
    );

    my %p = validate(
        @p,
        {

            # Make sure strand is 1 or -1 and set default
            strand => {
                default => $DEFAULT_STRAND,
                regex   => qr/^[+-]?1$/,
                type    => Params::Validate::SCALAR
            }
        }
    );

    # Get the index for the regex array
    my $rc = $p{strand} == 1 ? 0 : 1;

    # Get the regex, and if it is defined, return it
    my $regex = $self->_regexes->[$rc]->{$residue};
    if ( defined $regex ) {
        TRACE 'regular expression defined, leaving regex';
        return $regex;
    }

    # If the regex wasn't defined, build it by calling codons
    $regex = join '|', @{ $self->codons( $residue, \%p ) };
    $regex = qr/$regex/;

    # Cache the regex and return it
    $self->_regexes->[$rc]->{$residue} = $regex;

    TRACE 'cached regular expression, leaving regex';
    return $regex;
}

=head2 find

    my $locations = $translator->find( $seq_ref, $residue );
    my $locations = $translator->find( $seq_ref, $residue, \%params );

Find the indexes of a given residue in a sequence. In addition to the
one-letter codes for amino acids, the following are valid inputs for the
residue:

    start:  Start codons (you may also use "+" which is what the translator
            uses as the 1-letter code for start codons)
    stop:   Stop codons (you may also use "*" which is the 1 letter code)
    lower:  Start or stop codons, depending up on strand
    upper:  Start or stop codons, depending up on strand

"lower" and "upper" match the respective ends of a CDS for a given strand (i.e.
on the positive strand, lower matches the start, and upper matches the stop).
Valid options for the params hash are:

    strand:     1 or -1; default = 1

=cut

sub find {
    TRACE 'entering find';
    my $self = shift;

    my ( $seq_ref, $residue, @p ) = validate_pos(
        @_,
        { type => Params::Validate::SCALARREF },
        { type => Params::Validate::SCALAR },
        { type => Params::Validate::HASHREF, default => {} }
    );

    # Strand is unnecessary for now. Uncomment this section when other options
    # are added to find.
    #    my %p = validate(
    #        @p,
    #        {
    #            strand => {
    #                default => $DEFAULT_STRAND,
    #                regex   => $STRAND_REGEX,
    #                type    => Params::Validate::SCALAR
    #            }
    #        }
    #    );
    #
    #    my $regex = $self->regex( $residue, \%p );

    my $regex = $self->regex( $residue, @p );

    # Use a look-ahead in the regular expression. For instance, if the amino
    # acid has the codon AAA, and you have a poly-A region, the match will be
    # at every single base, not every 3 bases.
    my @positions;
    while ( $$seq_ref =~ m/(?=$regex)/ig ) {
        push @positions, pos($$seq_ref);
    }

    TRACE 'leaving find';
    return \@positions;
}

=head2 getORF

    my $orf_hash = $translator->getORF( $seq_ref );
    my $orf_hash = $translator->getORF( $seq_ref, \%params );

This will get the longest region between stops and return the strand, lower and
upper bounds, inclusive. Valid options for the params hash are:

    strand:     0, 1 or -1; default = 0 (meaning search both strands)
    lower:      integer between 0 and length; default = 0
    upper:      integer between 0 and length; default = length
    sanitized:  0 or 1; default = 0

Lower and upper are used to specify bounds between which you are searching.
Suppose the following was the longest ORF:

 0 1 2 3 4 5 6 7 8 9 10
  T A A A T C T A A G
  *****       *****
        <--------->

This will return:

    {
        strand => 1,
        lower  => 3,
        upper  => 9
    }

You can also specify which strand you are looking for the ORF to be on.

For ORFs starting at the very beginning of the strand or trailing off the end,
but not in phase with the start or ends, this method will cut at the last
complete codon. For example, if the following was the longest ORF:

    0 1 2 3 4 5 6 7 8 9 10
     A C G T A G T T T A
                   *****
       <--------------->

getORF will return:

    {
        strand => 1,
        lower  => 1,
        upper  => 10
    }

The distance between lower and upper will always be a multiple of 3. This is to
make it clear which frame the ORF is in. The resulting hash may be passed to
the translate method.

Example:

    my $orf_ref = $translator->getORF( \'TAGAAATAG' );
    my $orf_ref = $translator->getORF( \$seq, { strand => -1 } );
    my $orf_ref = $translator->getORF(
        \$seq,
        {
            lower => $lower,
            upper => $upper
        }
    );

=cut

sub getORF {
    TRACE 'entering getORF';

    my $self = shift;

    my ( $seq_ref, @p ) = validate_seq_params(@_);

    my ($sanitized) = validate_sanitized(@p);
    $seq_ref = cleanDNA($seq_ref) unless ($sanitized);

    my ( $lower, $upper ) = validate_lower_upper( $seq_ref, @p );

    my %p = validate( @p, {%VAL_SEARCH_STRAND} );

    # Initialize the longest ORF.
    my %ORF = (
        strand => 0,
        lower  => $lower,
        upper  => $lower
    );

    # Go through each strand which we are looking in
    foreach my $strand ( $p{strand} == 0 ? ( -1, 1 ) : $p{strand} ) {

        # Initialize lower bounds and regular expression for stop
        my @lowers = map { $_ + $lower } ( 0 .. 2 );
        my $stop_regex = $self->regex( '*', { strand => $strand } );

        # Look for all the stops in our sequence using a regular expression. A
        # lookahead is used to cope with the possibility of overlapping stop
        # codons

        pos($$seq_ref) = $lower;

        while ( $$seq_ref =~ /(?=stop_regex)/gx ) {

            # Get the location of the upper bound. Add 3 for the length of the
            # stop codon if we are on the + strand.
            my $cur_upper = pos($$seq_ref) + ( $strand == 1 ? 3 : 0 );

            # End the iteration if we are out of range
            last if ( $cur_upper > $upper );

            # Call our helper function
            $self->_getORF( $strand, \@lowers, $cur_upper, $lower, \%ORF );
        }

        # Now evaluate for the last three ORFS
        foreach my $i ( 0 .. 2 ) {
            my $cur_upper = $upper - $i;
            $self->_getORF( $strand, \@lowers, $cur_upper, $lower, \%ORF );
        }

        # NOTE: Perl's regular expression engine could be faster than code
        # execution, so it may be faster to find ORFS using regular expression
        # matching an entire ORF.
        # m/(?=(^|$stop)((.{3})*)($stop|$))/g
    }

    return \%ORF;
}

# Helper function for getORF above.
sub _getORF {
    my $self = shift;
    my ( $strand, $lowers, $upper, $offset, $longest ) = @_;

    # Calculate the frame relative to the starting offset
    my $frame = ( $upper - $offset ) % 3;

    # Compare if this is better than the longest ORF
    $self->_compare_regions(
        $longest,
        {
            strand => $strand,
            lower  => $lowers->[$frame],
            upper  => $upper
        }
    );

    # Mark the lower bound for this frame
    $lowers->[$frame] = $upper;
}

=head2 getCDS

    my $cds_ref = $translator->getCDS( $seq_ref );
    my $cds_ref = $translator->getCDS( $seq_ref, \%params );

This will return the strand and boundaries of the longest CDS.

 0 1 2 3 4 5 6 7 8 9 10
  A T G A A A T A A G
  >>>>>       *****
  <--------------->

Will return:

    {
        strand => 1,
        lower  => 0,
        upper  => 9
    }

Valid options for the params hash are:

    strand:     0, 1 or -1; default = 0 (meaning search both strands)
    lower:      integer between 0 and length; default = 0
    upper:      integer between 0 and length; default = length
    strict:     0, 1 or 2;  default = 1
    sanitized:  0 or 1; default = 0

Strict controls how strictly getCDS functions. There are 3 levels of
strictness, enumerated 0, 1 and 2. 2 is the most strict, and in that mode, a
region will only be considered a CDS if both the start and stop is found. In
strict level 1, if a start is found, but no stop is present before the end of
the sequence, the CDS will run until the end of the sequence. Strict level 0
assumes that start codon is present in each frame just before the start of the
molecule. Level 1 is a pretty safe bet, so that is the default.

Example:

    my $cds_ref = $translator->getCDS(\'ATGAAATAG');
    my $cds_ref = $translator->getCDS(\$seq, { strand => -1 } );
    my $cds_ref = $translator->getCDS(\$seq, { strict => 2 } );

=cut

sub getCDS {
    TRACE('getCDS called');

    my $self = shift;

    my ( $seq_ref, @p ) = validate_seq_params(@_);

    my $sanitized = validate_sanitized(@p);
    $seq_ref = cleanDNA($seq_ref) unless ($sanitized);

    my ( $lower, $upper ) = validate_lower_upper( $seq_ref, @p );

    my %p = validate(
        @p,
        {
            %VAL_SEARCH_STRAND,

            # Verify that strict is 0, 1 or 2, default to 1
            strict => {
                default => 1,
                regex   => $REGEX_012,
                type    => Params::Validate::SCALAR
            }
        }
    );

    # Initialize the longest CDS. Length is -1.
    my %CDS = (
        strand => 0,
        lower  => 0,
        upper  => -1
    );

    foreach my $strand ( $p{strand} == 0 ? ( -1, 1 ) : $p{strand} ) {
        my $lower_regex = $self->regex( 'lower', { strand => $strand } );
        my $upper_regex = $self->regex( 'upper', { strand => $strand } );

        # Initialize lowers. On the + strand, we don't set the lower bounds
        # unless strict is 0. On the - strand, we don't set the lower bounds if
        # strict is 2. Otherwise, set the lower boudns to be the first bases.
        my @lowers =
          (      ( ( $strand == 1 ) && ( $p{strict} != 0 ) )
              || ( ( $strand == -1 ) && ( $p{strict} == 2 ) ) )
          ? (undef) x 3
          : map { $lower + $_ } ( 0 .. 2 );

        # Similar to getORF, rather than using a regular expression to find
        # entire coding regions, instead find individual starts and stops and
        # react accordingly.
        # The regular expression captures the starts and stops separately
        # ($1 vs $2) so that it is easy to tell if a start or a stop was
        # matched.

        pos($$seq_ref) = $lower;

        while ( $$seq_ref =~ /(?=($lower_regex)|($upper_regex))/g ) {
            my $position = pos $$seq_ref;
            last if ( $position > $upper );

            my $frame = $position % 3;

            # If the lower regex matches:
            #
            # In the case that it is on the '-' strand, that means a stop was
            # found. CDSs always end on stops, so update the lower bound.
            #
            # Otherwise, it is on the positive strand, meaning a start was
            # found. Internal start codons are allowed, so only set the lower
            # bound if it is not already set.
            if ($1) {
                if (   ( $strand == -1 )
                    || ( !defined $lowers[$frame] ) )

                {
                    $lowers[$frame] = $position;
                }
            }

            # If the lower regex wasn't matched, the the upper one was.
            #
            # If this is the positive strand, that means that this is a stop
            # codon. Compute the CDS, update if necessary, and reset the lower
            # bound in this case.
            #
            # On the negative strand, that means that a start was matched.
            # Compute the CDS, update if necessary, but don't reset the lower
            # bound.

            else {
                $position += 3;
                last if ( $position > $upper );

                $self->_getCDS( $strand, \@lowers, $position, $lower, \%CDS );
            }
        }

        # If strict mode is at level 2, we don't allow CDSs to trail off the
        # end of the molecule. We also don't allow the end to trail off if we
        # are on the - strand and strict isn't 0.

        next
          if ( ( $p{strict} == 2 )
            || ( ( $strand == -1 ) && ( $p{strict} != 0 ) ) );

        foreach my $i ( 0 .. 2 ) {
            my $end_upper = $upper - $i;
            $self->_getCDS( $strand, \@lowers, $end_upper, $lower, \%CDS );
        }
    }

    return \%CDS;
}

# Helper function for getCDS above.
sub _getCDS {
    my $self = shift;
    my ( $strand, $lowers, $upper, $offset, $longest ) = @_;

    # Calculate the frame relative to the starting offset
    my $frame = ( $upper - $offset ) % 3;

    # Do nothing if lower bound wasn't defined
    return unless ( defined $lowers->[$frame] );

    # Compare if this is better than the longest ORF
    $self->_compare_regions(
        $longest,
        {
            strand => $strand,
            lower  => $lowers->[$frame],
            upper  => $upper
        }
    );

    # Mark the lower bound for this frame
    undef $lowers->[$frame] if ( $strand == 1 );
}

# If the current ORF is longer than the previously stored longest bounds, store
# the current ORF
sub _compare_regions {
    my $self = shift;
    my ( $longest, $current ) = @_;
    %$longest = %$current
      if ( $longest->{upper} - $longest->{lower} <
        $current->{upper} - $current->{lower} );
}

=head2 nonstop

    my $frames = $translator->nonstop( $seq_ref );
    my $frames = $translator->nonstop( $seq_ref, \%params );

Returns the frames that contain no stop codons for the sequence. Frames are
numbered -3, -2, -1, 1, 2 and 3.

     3   ---->
     2  ----->
     1 ------>
       -------
    -1 <------
    -2 <-----
    -3 <----

Valid options for the params hash are:

    strand:     0, 1 or -1; default = 0 (meaning search both strands)
    sanitized:  0 or 1; default = 0

Example:

    my $frames = $translator->nonstop(\'TACGTTGGTTAAGTT'); # [ 2, 3, -1, -3 ]
    my $frames = $translator->nonstop(\$seq, { strand => 1 }  ); # [ 2, 3 ]
    my $frames = $translator->nonstop(\$seq, { strand => -1 } ); # [ -1, -3 ]

=cut

sub nonstop {
    TRACE('nonstop called');

    my $self = shift;

    my ( $seq_ref, @p ) = validate_seq_params(@_);

    my $sanitized = validate_sanitized(@p);
    $seq_ref = cleanDNA($seq_ref) unless ($sanitized);

    my %p = validate( @p, {%VAL_SEARCH_STRAND} );

    # Go through both strands
    my @frames;
    foreach my $strand ( $p{strand} == 0 ? ( 1, -1 ) : $p{strand} ) {
        my $stop = $self->regex( '*', { strand => $strand } );

        # Go through each frame
        foreach my $frame ( 0 .. 2 ) {
            my $regex =
              $strand == 1
              ? qr/^.{$frame}(?:.{3})*$stop/
              : qr/$stop(?:.{3})*.{$frame}$/;

            # Convert strand = +/-1, frame = [0,1,2] into 1 value +/-[1,2,3]
            push @frames, ( $frame + 1 ) * $strand
              unless ( $$seq_ref =~ m/$regex/ );
        }
    }

    return \@frames;
}

1;

=head1 AUTHOR

Kevin Galinsky, <kgalinsky plus cpan at gmail dot com>

=head1 COPYRIGHT & LICENSE

Copyright 2008-2009 J. Craig Venter Institute, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut
