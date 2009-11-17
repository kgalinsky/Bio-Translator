use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'JCVI::Translator',
    license             => 'perl',
    dist_author         => 'Kevin Galinsky <kgalinsk@jcvi.org>',
    dist_version_from   => 'lib/JCVI/Translator.pm',
    build_requires => {
        'Test::More' => 0,
    },
    requires => {
        'JCVI::AATools' => '>= 0.1.0',
        'JCVI::DNATools' => '>= 0.2.0'
    },
    add_to_cleanup      => [ 'JCVI-Translator-*' ],
    create_makefile_pl => 'traditional',
);

$builder->create_build_script();
