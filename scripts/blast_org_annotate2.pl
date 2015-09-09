#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

open my $fh, '<', 'accessions.txt';
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
