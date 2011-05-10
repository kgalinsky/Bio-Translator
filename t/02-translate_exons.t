#!perl

use Test::More 'no_plan';

use Bio::Tiny::Translator;

my $seq = 'CTGATATCATGCATGCCATTCTCGACCGCTATGCGCCTCCTGTTCCTCGTGGGCCCAAAA';

my $translator = new Bio::Tiny::Translator();

# Simple cases
is( ${ $translator->translate_exons( \$seq, [ [ 0, 60 ] ] ) },
    'MISCMPFSTAMRLLFLVGPK', 'Translate frame 1' );
is(
    ${
        $translator->translate_exons( \$seq, [ [ 0, 60 ] ], { partial => 1 } )
      },
    'LISCMPFSTAMRLLFLVGPK',
    'Translate frame 1 with partial flag'
);

is(
    ${
        $translator->translate_exons( \$seq, [ [ 0, 60 ] ], { strand => -1 } )
      },
    'FWAHEEQEAHSGREWHA*YQ',
    'Translate frame -1'
);

# Break in the middle
is( ${ $translator->translate_exons( \$seq, [ [ 0, 30 ], [ 30, 60 ] ] ) },
    'MISCMPFSTAMRLLFLVGPK', 'Translate frame 1 with break' );

is(
    ${
        $translator->translate_exons(
            \$seq,
            [ [ 0, 30 ], [ 30, 60 ] ],
            { strand => -1 }
        )
      },
    'FWAHEEQEAHSGREWHA*YQ',
    'Translate frame -1 with break'
);

# Small codons
is(
    ${
        $translator->translate_exons( \$seq,
            [ map { [ $_, $_ + 1 ] } ( 0 .. 59 ) ] )
      },
    'MISCMPFSTAMRLLFLVGPK',
    'Translate frame 1 with small breaks'
);

is(
    ${
        $translator->translate_exons(
            \$seq,
            [ map { [ $_, $_ + 1 ] } ( 0 .. 59 ) ],
            { strand => -1 }
        )
      },
    'FWAHEEQEAHSGREWHA*YQ',
    'Translate frame -1 with small breaks'
);

# Gaps
is( ${ $translator->translate_exons( \$seq, [ [ 0, 18 ], [ 21, 60 ] ] ) },
    'MISCMPSTAMRLLFLVGPK', 'Translate frame 1 with gap' );

is(
    ${
        $translator->translate_exons(
            \$seq,
            [ [ 0, 6 ], [ 9, 60 ] ],
            { strand => -1 }
        )
      },
    'FWAHEEQEAHSGREWHAYQ',
    'Translate frame -1 with gap'
);

# Traverse reading frames
is( ${ $translator->translate_exons( \$seq, [ [ 0, 18 ], [ 28, 60 ] ] ) },
    'MISCMPLCASCSSWAQ', 'Translate with gap going from frame 1 to 2' );

is(
    ${
        $translator->translate_exons(
            \$seq,
            [ [ 6, 21 ], [ 31, 58 ] ],
            { strand => -1 }
        )
      },
    'MGPRGTGGAEWHA*',
    'Translate with gap going from frame -3 to -1'
);

# Translate with random gaps
is(
    ${
        $translator->translate_exons( \$seq,
            [ [ 6, 23 ], [ 26, 32 ], [ 36, 54 ] ] )
      },
    'SCMPFSAISCSSW',
    'Translate with random gaps on + strand'
);

is(
    ${
        $translator->translate_exons(
            \$seq,
            [ [ 6, 23 ], [ 26, 32 ], [ 36, 54 ] ],
            { strand => -1 }
        )
      },
    'AHEEQEIAENGMH',
    'Translate with random gaps on - strand'
);
