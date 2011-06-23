use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'Bio::Translator',
    license             => 'perl',
    dist_author         => 'Kevin Galinsky <kgalinsky plus cpan at gmail dot com>',
    dist_version_from   => 'lib/Bio/Translator.pm',
    build_requires => {
        'Test::More' => 0,
    },
    requires => {
        'Bio::Util::AA' => '>= 0.1.0',
        'Bio::Util::DNA' => '>= 0.2.0'
    },
    add_to_cleanup      => [ 'Bio-Translator-*' ],
    create_makefile_pl => 'traditional',
);

$builder->create_build_script();
