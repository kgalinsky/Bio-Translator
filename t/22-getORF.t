#!perl

use Test::More 'no_plan';
use JCVI::DNATools qw(randomDNA);

use JCVI::Translator::Utils;

my $utils = new JCVI::Translator::Utils;

ok( $utils->getORF( randomDNA() ), 'getORF ran' );
