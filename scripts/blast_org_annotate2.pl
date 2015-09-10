#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use Getopt::Long 'GetOptions';

my $accessions = '';
GetOptions(
    'accessions=s' => \$accessions
);

unless ($accessions) {
    die "Missing -a (accesssions file) argument\n";
}

unless (-s $accessions) {
    die "Bad accesssions file ($accessions)\n";
}

open my $fh, '<', $accessions;
my %rh_map = map { chomp; split(/\t/, $_) } <$fh>;
close $fh;

while (my $line = <>) {
    chomp($line);
    if ($line =~ /\#/) {
        print "$line\n";
    }
    else {
        my @flds = split(/\t/, $line);
        my $acc_ver = $flds[-1]; # last column should always be acc.ver

        if (my $val = $rh_map{ $acc_ver }) {
            print join("\t", $line, $val), "\n";
        }
        else {
            die "ERROR: Accession '$acc_ver' not found!\n";
        }
    }
}
