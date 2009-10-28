#!perl

use Test::More tests => 4;

BEGIN {
    use_ok('JCVI::Translator');
    use_ok('JCVI::Translator::Utils');
    use_ok('JCVI::Translator::Base');
    use_ok('JCVI::Translator::Table');
}

diag("Testing JCVI::Translator $JCVI::Translator::VERSION, Perl $], $^X");
