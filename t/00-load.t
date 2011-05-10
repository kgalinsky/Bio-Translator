#!perl

use Test::More tests => 4;

BEGIN {
    use_ok('Bio::Tiny::Translator');
    use_ok('Bio::Tiny::Translator::Utils');
    use_ok('Bio::Tiny::Translator::Base');
    use_ok('Bio::Tiny::Translator::Table');
}

diag("Testing Bio::Tiny::Translator $Bio::Tiny::Translator::VERSION, Perl $], $^X");
