#!perl

use Test::More 'no_plan';
use Bio::Tiny::Util::DNA qw(randomDNA);

use Bio::Tiny::Translator;

my $translator = new Bio::Tiny::Translator;

my $dna     = randomDNA();
my $peptide = $translator->translate6($dna);

ok( $peptide, 'translate6 returned something' );

foreach my $strand ( 1, -1 ) {
    foreach my $offset ( 0 .. 2 ) {
        my $reference =
          $translator->translate( $dna,
            { strand => $strand, offset => $offset } );
        is( $peptide->[ $offset + ( ( $strand == 1 ? 0 : 1 ) * 3 ) ],
            $$reference, 'result of translate6 matches translate' );
    }
}
