use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'Bio::Tiny::Translator',
    license             => 'perl',
    dist_author         => 'Kevin Galinsky <kgalinsky plus cpan at gmail dot com>',
    dist_version_from   => 'lib/Bio/Tiny/Translator.pm',
    build_requires => {
        'Test::More' => 0,
    },
    requires => {
        'JCVI::AATools' => '>= 0.1.0',
        'JCVI::DNATools' => '>= 0.2.0'
    },
    add_to_cleanup      => [ 'Bio-Tiny-Translator-*' ],
    create_makefile_pl => 'traditional',
);

$builder->create_build_script();
