#!perl

use Test::More 'no_plan';
use JCVI::DNATools qw(randomDNA);

use Bio::Tiny::Translator::Utils;

my $utils = new Bio::Tiny::Translator::Utils;

ok( $utils->getORF( randomDNA() ), 'getORF ran' );
