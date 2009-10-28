#!perl

use Test::More 'no_plan';
use JCVI::DNATools qw(randomDNA);

use JCVI::Translator;

my $translator = new JCVI::Translator;

ok( $translator->translate6( randomDNA() ), 'translate6 ran' );
