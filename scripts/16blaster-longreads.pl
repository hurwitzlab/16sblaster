#!/usr/bin/env perl

use common::sense;
use autodie;
use Cwd qw(cwd);
use Data::Dump qw(dump);
use File::Spec::Functions qw(canonpath catdir catfile);
use File::Path qw(make_path remove_tree);
use File::Basename qw(basename);
use File::Find::Rule;
use File::Which qw(which);
use List::Util qw(max);
use Perl6::Slurp qw(slurp);
use Getopt::Long;
use Pod::Usage;
use Readonly;

my $DEBUG = 0;

main();

# --------------------------------------------------
sub main {
    my $in_dir     = '';
    my $blast_db   = '';
    my $out_dir    = catdir(cwd(), 'RESULTS');
    my $similarity = 0.98;
    my $percent    = 0.01;
    my $length     = 25;
    my $glob       = '*.fasta';
    my ($help, $man_page);

    GetOptions(
        'help'           => \$help,
        'man'            => \$man_page,
        'd|dir=s'        => \$in_dir,
        'o|out:s'        => \$out_dir,
        'c|similarity:f' => \$similarity,
        'p|percent:f'    => \$percent,
        'l|length:i'     => \$length,
        'g|glob:s'       => \$glob,
        'b|blast-db=s'   => \$blast_db,
        'debug'          => \$DEBUG,
    ) or pod2usage(2);

    if ($help || $man_page) {
        pod2usage({
            -exitval => 0,
            -verbose => $man_page ? 2 : 1
        });
    }; 

    if ($in_dir) {
        $in_dir = canonpath($in_dir);
    }
    else {
        pod2usage('Missing "-d" input directory');
    }

    unless (-d $in_dir) {
        pod2usage("Bad directory ($in_dir)");
    }

    unless ($blast_db) {
        pod2usage('Missing --blast-db argument');
    }

    unless (-s $blast_db) {
        pod2usage("Bad --blast-db ($blast_db)");
    }

    my @files = File::Find::Rule->file()->name($glob)->in($in_dir)
        or die "Cannot find any '$glob' files in '$in_dir'\n";

    printf "Found %s files in '$in_dir'\n", scalar @files, $in_dir;

    for my $exe ('cd-hit-est', 'blastn') {
        unless (which($exe)) {
            pod2usage("Cannot find '$exe'");
        }
    }

    process(
        out_dir    => $out_dir, 
        in_dir     => $in_dir, 
        files      => \@files,
        similarity => $similarity,
        percent    => $percent,
        length     => $length,
        blast_db   => $blast_db,
    );
}

# --------------------------------------------------
sub process {
    my %args       = @_;
    my $in_dir     = $args{'in_dir'};
    my $out_dir    = $args{'out_dir'};
    my $files      = $args{'files'};
    my $similarity = $args{'similarity'};
    my $percent    = $args{'percent'};
    my $length     = $args{'length'};
    my $blast_db   = $args{'blast_db'};

    if (-d $out_dir) {
        say "Previous '$out_dir' directory deleted.";
        remove_tree($out_dir, { keep_root => 1 });
    }
    else {
        make_path($out_dir);
    }

    # Perform two rounds of clustering using CD-HIT
    my @arraycounts;
    my @readcounts;
    my $i = 0;

    FILE:
    for my $path (@$files) {
        my $filename = basename($path);
        my $out_file = catfile($out_dir, $filename);
        my $num_seqs = `grep '>' $path | wc -l`;

        unless ($num_seqs > 0) {
            say STDERR "No sequences in '$path'?";
            next FILE;        
        }

        printf "%5d: cd-hit-est '%s' -> '%s'\n", ++$i, $filename, $out_file;

        execute(
            'cd-hit-est',
            '-i', $path, 
            '-o', $out_file, 
            '-c', $similarity, 
            '-l', $length, 
            qw(-n 9 -d 0 -g 1 -T 8 -M 32000 -s 0.98)
        );

        my $take      = int($num_seqs * $percent);
        my $multi_dir = catdir($out_dir, $filename . "MAKEMULTI$take");

        execute(
            'make_multi_seq.pl',
            $path,
            catfile($out_dir, $filename . '.clstr'),
            $multi_dir,
            $take,
        );

        my @multi_files = File::Find::Rule->file()->in($multi_dir);

        unless (@multi_files) {
            say STDERR "No results from 'make_multi_seq.pl', skipping.";
            next FILE;
        }

        my $contigs_file = catfile($out_dir, $filename . '.contigs.fasta');
        execute("cat $multi_dir/* > $contigs_file");

        my $reclustered_file = catfile($out_dir, $filename . '.reclustered');
        execute(
            'cd-hit-est',
            '-i', $contigs_file,
            '-o', $reclustered_file,
            '-c', $similarity,
            '-l', $length,
            qw(-n 9 -d 0 -g 1 -T 8 -M 32000 -s 0.98),
        );

        #
        # While we're in this loop, count the size of the clusters (in
        # reads) so we can report that later on
        #
        my $cluster_file = catfile($out_dir, "$filename.recluster.fasta.clstr");

        unless (-e $cluster_file) {
            say "Cannot find expected cluster file '$cluster_file'";
            next FILE;
        }

        open my $CLSTR_FH, '<', $cluster_file;
        chomp(my @clusters = <$CLSTR_FH>);
        close $CLSTR_FH;
        my $element_count = 0;
        my $count;

        for my $line (@clusters) {
            $element_count++;

            if (
                ($line =~ /^\>/ && $element_count > 1) # not the first
                ||                                     # or
                $element_count == scalar @clusters     # EOF
            ) {
                # look at previous line
                my $number = (split /\t/, $clusters[ $element_count - 1 ])[0];
                push @arraycounts, $number;
                $count++;
            }
        }

        while ($count > 0) {
            my $number = ($element_count - (scalar @arraycounts));
            push @readcounts, $number;
            $count--;
        }

        say 'Finished clustering';
    }

    #
    # After two rounds of clustering we have a fasta file of
    # representative reads, make a .sta (just a fasta with one
    # read) file for each one so they can be blasted one at a
    # time...
    #
    my @reclustered =
        File::Find::Rule->file()->name('*.reclustered.fasta')->in($out_dir);
    my $read_num = 0;
    my @read_ids = ();

    for my $file (@reclustered) {
        open my $fh, '<', $file;
        local $/ = '>';
        while (my $rec = <$fh>) { 
            chomp($rec);
            next unless $rec;
            my ($read_id, @seq) = split(/\n/, $rec);
            push @read_ids, $read_id;
            my $out_name = join('.', $file, $i, $read_id, 'sta');
            open my $out, '>', catfile($out_dir, $out_name);
            print $out, join("\n", ">$out_name", join('', @seq), '');
            close $out;
            $i++;
        }
    }

    my %hash;
    for my $file (
        File::Find::Rule->file()->name('*.sta')->in($out_dir)
    ) {
        my $key        = (split /\./,   $file)[6];
        my $samplename = (split /\/|_/, $file)[2];
        $hash{ $key } = "$file";
    }

    my @sorted;
    for my $key (sort { $a <=> $b } keys %hash) {
        push @sorted, $hash{ $key };
    }

    #
    # Blast each of the .sta files
    #
    print "\nBlasting clusters and parsing the results...\n";
    my $num_alignments = 250;

    for my $file (@sorted) {
        my $cmd = join(' ',
            'blastn',
            '-db', $blast_db,
            '-query', $file, 
            '-num_threads', '9',
            '-outfmt', '"7 qacc sallseqid evalue bitscore pident qstart '
                     . 'qend sstart send saccver"',
            '-num_alignments', $num_alignments,
            'grep -v "#" | perl ./blast_org_annotate2.pl'
        );
        my $hits = `$cmd`;

        debug("hits: $hits");

        my @besthits;
        my @maxhits;
        my @genuslist;
        my @specieslist;
        my $sequence = '';

        if (defined $hits && $hits ne '') {
            my (@pidents, @bitscores);
            my @blastlines = split(/\n/, $hits);
            for my $blastline (@blastlines) {
                my @res = split(/\t/, $blastline);
                push @pidents,   $res[4];
                push @bitscores, $res[3];
            }

            my $max          = max(@bitscores);
            my $bitscoresize = scalar(@bitscores);
            (my $filename    = $file) =~ s{^$out_dir/?}{};
            open my $fh, '<', $file;
            my @contents = <$fh>;
            close $fh;
            $sequence = @contents[1] or die 'No sequence';

            #
            # Add a header indicating the file we are reporting on, 
            # consensus, sequence, number of reads represented.
            #
            $readcounts[0]++ if $readcounts[0] < 1;

            push @besthits,
                "\n\n$filename",
                "\nReadID and sequence:\n>$read_ids[0]\n$sequence",
                "\nWith:\t",
                $arraycounts[0],
                " reads out of $readcounts[0] total = ",
                (sprintf '%.1f', (100 * ($arraycounts[0] / $readcounts[0]))),
                "% of reads";

            shift @read_ids;
            my $y = 0;
            while ($y < $bitscoresize) {
                if ($bitscores[$y] == $max) {
                    my @results = split(/\t/, $blastlines[$y]);
                    push(@maxhits, $results[10]);
                    if ($pidents[$y] >= 98.00) {
                        push(@besthits,
                            "\nFound bitscore ",
                            $results[3],
                            " with ",
                            $results[4],
                            "% match at the species level (>98%) with: ",
                            "\t$results[10]\taccession: $results[9]");
                    }
                    else {
                        push(@besthits,
                            "\nFound bitscore ",
                            $results[3],
                            " with ",
                            $results[4],
                            "% match at the Genus level (i.e. <98%) with: ",
                            "\t$results[10]\taccession: $results[9]");
                    }
                }
                $y++;
                foreach my $result (@maxhits) {
                    my $genus = (split / /, $result)[0];

                    #$genus = $genus.",";
                    my $species = (split / /, $result)[1];
                    my $genusspecies = join " ", $genus, $species;
                    if (not $genus ~~ @genuslist) {
                        push(@genuslist,   $genus);
                        push(@specieslist, $species);
                    }
                    elsif (not $species ~~ @specieslist) {
                        push(@specieslist, $species);
                    }
                    else { }
                }
            }
        }
        else {
            open my $fh, '<', $file;
            my @contents = <$fh>;
            $sequence = $contents[1];
            (my $filename = $file) =~ s{$out_dir/?}{};

            $readcounts[0]++ if $readcounts[0] < 1;

            push @besthits,
                "\n\n$filename",
                "\nReadID and sequence:\n>$read_ids[0]\n$sequence",
                "\nWith:\t",
                $arraycounts[0],
                " reads out of $readcounts[0] total = ",
                (sprintf '%.1f', (100 * ($arraycounts[0] / $readcounts[0]))),
                "% of reads",
                "\nNo result: no useful sequence or no match in database.\n";

            push @genuslist, $file, ":No BLAST result";
            shift @read_ids;
        }

        open my $results_fh, '>', catfile($out_dir, "$file.besthits.txt");
        print $results_fh "@besthits";
        close $results_fh;

        open my $summary_fh, '>', catfile($out_dir, "$file.summary.txt");
        my @array = split /.clustered.fasta.recluster.fasta./, $file;
        my $sample        = $array[0];
        my $clusternumber = (split /\./, $array[1])[0];
        my $readID2       = (split /\./, $array[1])[1];
        print $summary_fh
        "\n$sample\t$clusternumber\t$readID2\t$sequence\t$arraycounts[0] reads";
        shift @arraycounts;
        shift @readcounts;

        if (scalar @genuslist > 1) {
            print $summary_fh "\tFamily or higher hit\t", 
                join(", ", @genuslist),
                "\tNot Applicable";
        }
        elsif (scalar @specieslist > 1) {
            print $summary_fh "\tGenus hit\t", @genuslist, "\t";
            for my $species (@specieslist) {
                print $summary_fh @genuslist, " ", $species, ", ";
            }
        }
        else {
            print $summary_fh "\tSpecies hit\t", @genuslist, "\t", 
                @genuslist, " ",
                @specieslist;
        }
        close $summary_fh;
    }

    my $commandline = join ' ', $0, @ARGV;

    open my $header_fh, ">", catfile($out_dir, 'HEADER.txt');
    print $header_fh join("\n",
        $commandline,
        '',
        "Similarity score for clustering: $similarity",
        "Fraction of total reads that a cluster must contain to be " .
        "included in re-clustering: $percent",
        "Length of read to be included in clustering: $length",
        ''
    );
    close $header_fh;

    open my $summaryheader_fh, '>', catfile($out_dir, 'SUMMARYHEADER.txt');
    print $summaryheader_fh join("\n",
        $commandline,
        '',
        "Similarity score for clustering: $similarity",
        "Fraction of total reads that a cluster must contain to be " .
        "included in re-clustering: $percent",
        "Length of read to be included in clustering: $length",
        '',
        join("\t", 
            'Sample', 
            'Cluster', 
            'Sequence', 
            'Taxonomic level', 
            'genera/genus matched', 
            'species matched(if applicable)',
        ),
    );

    close $summaryheader_fh;
    print join("\n",
        '',
        'Top BLAST results by bitscore for all input .fasta files are in',
        "$out_dir/ALL_RESULTS.txt",
        '',
        "Summary results for all input .fasta files are in ",
        "$out_dir/SUMMARY.txt",
        '',
    );

    for my $cmd (
        "cat $out_dir/HEADER.txt $out_dir/*.besthits.txt > " .
        "$out_dir/ALL_RESULTS.txt",

        "cat $out_dir/SUMMARYHEADER.txt $out_dir/*.summary.txt > " .
        "$out_dir/SUMMARY.txt"
    ) {
        system($cmd) == 0 or die $?;
    }
}

# --------------------------------------------------
sub execute {
    my @cmd = @_ or return;
    debug("\n\n>>>>>>\n\n", join(' ', @cmd), "\n\n<<<<<<\n\n");

    unless (system(@cmd) == 0) {
        die sprintf(
            "FATAL ERROR! Could not execute command:\n%s\n", 
            join(' ', @cmd)
        );
    }
}


# --------------------------------------------------
sub debug {
    say @_ if $DEBUG;
}

__END__

# --------------------------------------------------

=pod

=head1 NAME

16blaster-longreads.pl - short description here

=head1 SYNOPSIS

  16blaster.longreads.pl -d /path/to/fasta --blast-db /path/to/db

Required Arguments:

  -d|--dir          Directory containing FASTA files processing

  --blast-db        Directory containing the BLAST 16s db

Options:

  -o|--out          Where to write the results (default "$PWD/RESULTS")

  -c|--similarity   Value between 0.80 and 1.00 for percent similarity 
                    for cd-hit clustering. Strongly recommend 0.98 
                    (default) to start.

  -p|--percent      Value between 0 and 1.00 for percent of total reads 
                    in order for clusters to be included second round 
                    of clustering.
                    Recommend 0.01 (default) to start.
                    Though running this script with several reasonable
                    iterations may be necessary to optimize results
                    (e.g. try 0.01, 0.03, 0.10) and then check results
                    to see if they vary.
                    Using too high a value will result in no clusters
                    being included in the re-clustering and an error in
                    the cat step where all results are combined (\"No
                    such file or directory)\"
                    Using too low a value will result in many redundant
                    results, a better outcome than using too high a
                    number, so start small and work upward.

  -l|--length       Use this option to exlude short reads from being
                    included in the clustering. Default is 25. A
                    reasonable length restriction would be 75% of
                    expected full length read of the amplicon in
                    question. This prevents short contaminants or a
                    frequent short read from creating clusters that will
                    be reported in results.

  -g|--glob         File pattern glob to find the FASTA files 
                    (default "\*.fasta")

  --debug           Show extra messages

  --help            Show brief help and exit
  --man             Show full documentation

=head1 DESCRIPTION

This script relies on the "cd-hit" programs and "blast_org_annotate.pl"
which must be in your $PATH.  Also required is the file
"accession_organism_map.txt."

This script uses cd-hit to create clusters of representative
bacterial 16s reads and then BLASTs the representative reads to
identify the species that match the read.
              
All FASTA files in the directory provided with the "-d" option 
will be processed.

Indicate FASTA file extension other than ".fasta" with "-g" option.

Remember to remove the primer sequence used to generate the 16s
amplicons from the FASTA files before using this script. 

=head1 AUTHORS

=over 4

=item * George Watts E<lt>gwatts@email.arizona.eduE<gt>

=item * Ken Youens-Clark E<lt>kyclark@email.arizona.eduE<gt>

=back

=head1 COPYRIGHT

Copyright (c) 2015 George Watts

This module is free software; you can redistribute it and/or
modify it under the terms of the GPL (either version 1, or at
your option, any later version) or the Artistic License 2.0.
Refer to LICENSE for the full license text and to DISCLAIMER for
additional warranty disclaimers.

=cut
