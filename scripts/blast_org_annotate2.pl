#!/usr/bin/perl

use strict;
use warnings;
use FileHandle;

our $MAP_FILE = 'accessions.txt';

MAIN_CODE: {
    my $rh_map = load_map_into_hash();
    while(<>) {
	my $i = $_;
	chomp $i;
	if ($i =~ /\#/) {
	    print $i, "\n"; next; #print comment line and move on.
	}
	my @i = split(/\t/, $i);
	my $acc_ver = $i[-1]; #last column should always be acc.ver
	unless(defined $rh_map->{$acc_ver}) {
	    die "ERROR: accession, $acc_ver, does not appear to be in map file.  Should not see this.\n";
	}
	print $i, "\t", $rh_map->{$acc_ver}, "\n";
    }
}

sub load_map_into_hash {
    my $fh_IN = new FileHandle;
    $fh_IN->open($MAP_FILE) || die "ERROR: Cannot open mapaccession input file.\n";
    my %hash;
    while(<$fh_IN>) {
	chomp $_;
	my $i = $_;
	my @i = split(/\t/, $i);
	$hash{$i[0]} = $i[1];
    }
    $fh_IN->close();
    return \%hash;
}
