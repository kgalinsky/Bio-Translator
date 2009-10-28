use strict;
use warnings;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME                => 'JCVI::Translator',
    AUTHOR              => 'Kevin Galinsky <kgalinsk@jcvi.org>',
    VERSION_FROM        => 'lib/JCVI/Translator.pm',
    ABSTRACT_FROM       => 'lib/JCVI/Translator.pm',
    PL_FILES            => {},
    PREREQ_PM => {
        'Test::More' => 0,
    },
    dist                => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
    clean               => { FILES => 'JCVI-Translator-*' },
);
