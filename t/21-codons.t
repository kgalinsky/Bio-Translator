#!perl

use Test::More 'no_plan';
use JCVI::Translator::Utils;
use List::Compare;

my $utils = new JCVI::Translator::Utils();

eval { $utils->codons() };
ok( $@, 'codons died with no parameters' );

eval { $utils->codons('') };
ok( $@, 'codons died with empty string' );

eval { $utils->codons('foo') };
ok( $@, 'codons died on invalid codon' );

eval { $utils->codons('F') };
ok( !$@, 'codons ran with just codon' );

eval { $utils->codons( 'F', { strand => 1 } ) };
ok( !$@, 'codons ran with strand = 1' );

eval { $utils->codons( 'F', { strand => -1 } ) };
ok( !$@, 'codons ran with strand = -1' );

eval { $utils->codons( 'F', { strand => 2 } ) };
ok( $@, 'codons died with strand = 2' );

my @expected = ( [qw(TTT TTC TTY)], [qw(AAA GAA RAA)] );

foreach my $strand ( 1, -1 ) {
    my $codons = $utils->codons( 'F', { strand => $strand } );
    my $expected = shift @expected;
    is( scalar(@$codons), scalar(@$expected), 'Expected number of codons' );

    my $lc = List::Compare->new( $expected, $codons );

    is( scalar( $lc->get_symdiff ), 0, '0 differences between lists' );
}
