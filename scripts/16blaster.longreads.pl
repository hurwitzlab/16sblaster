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
use Getopt::Long;
use Pod::Usage;
use Readonly;

main();

# --------------------------------------------------
sub main {
    my $in_dir     = '';
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

    say "glob ($glob)";
    my @files = File::Find::Rule->file()->name($glob)->in($in_dir)
        or die "Cannot find any '$glob' files in '$in_dir'\n";

    printf "Found %s files in '$in_dir'\n", scalar @files, $in_dir;

    unless (which('cd-hit-est')) {
        pod2usage('Cannot find cd-hit-est');
    }

    process(
        in_dir     => $in_dir, 
        out_dir    => $out_dir, 
        files      => \@files,
        similarity => $similarity,
        percent    => $percent,
        length     => $length,
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

#    my $command = `ps -o args -C perl`;
#    $command =~ s/COMMAND\n//;
#    chomp $command;
#    my $commandline = join "", "command line: \"", $command,
#      "\" was run with options: ", " -c: ", $similarity, " -p: ", $percent,
#      " -l: ", $length;
#

    if (-d $out_dir) {
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

        last;
    }   
}

# --------------------------------------------------
sub execute {
    my @cmd = @_ or return;
    say "\n\n>>>>>>\n\n", join(' ', @cmd), "\n\n<<<<<<\n\n";

    unless (system(@cmd) == 0) {
        die sprintf("Could not execute command:\n%s\n%s\n", 
            join(' ', @cmd), $?
        );
    }
}

#
#    #While we're in this loop, count the size of the clusters (in reads) so we can report that later on
#        open(CLSTR, "<", "$directory\/RESULTS\/$filename.recluster.fasta.clstr")
#          || die
#    "\n\nCan't open: \$directory\/RESULTS\/\$filename.recluster.fasta.clstr:$directory\/RESULTS\/$filename.recluster.fasta.clstr";
#        local $/;
#        my $clstr = <CLSTR>;
#        close(CLSTR);
#        chomp $clstr;
#        my @cluster = split(/\n/, $clstr);
#        my $elementcount = 0;
#        my $number;
#        my $count;
#
#        foreach my $part (@cluster) {
#            if ($part =~ /^\>/ && $elementcount > 0) {
#                $number = (split /\t/, $cluster[ $elementcount - 1 ])[0];
#                push(@arraycounts, $number);
#                $elementcount++;
#                $count++;
#            }
#            else {
#                $elementcount++;
#            }
#        }
#        $number = (split(/\t/, $cluster[ $elementcount - 1 ]))[0];
#        push(@arraycounts, $number);
#        $count++;    #add the last count to the @arraycounts
#        while ($count > 0) {
#            my $number = ($elementcount - (scalar @arraycounts));
#            push(@readcounts, $number);
#            $count = ($count - 1);
#        }
#        print "\tDone Clustering $filepath\n";
#    }
#
#    #After two rounds of clustering we have a fasta file of representative reads, make a .sta (just a fasta with one read) file for each one so they can be blasted one at a time...
#    @filepath = <$directory/RESULTS/*.recluster.fasta>;
#    my $i = 0;
#    my $readID;
#    my @readIDarray;
#    for my $filepath (@filepath) {
#        open(CONTENTS, "<", "$filepath")
#          || die "\nCan't open \$filepath: $filepath $!\n";
#        local $/;
#        my $contents = <CONTENTS>;
#        close(CONTENTS);
#        chomp($contents);
#        my @array = split(/\>/, $contents);
#        shift(@array);
#        foreach my $headerandsequence (@array) {
#            $readID = (split /\n/, $headerandsequence)[0];
#            push(@readIDarray, $readID);
#            open(OUT2, ">", "$filepath.$i.$readID.sta")
#              or die "Can't open $filepath.$i.$readID.sta for writing: $!\n";
#            print(OUT2 ">", $filepath, ".", $i, ".", $readID, ".sta\t",
#                $headerandsequence, "\n");
#            close(OUT2);
#            $i++;
#        }
#    }
#    my @files = <$directory/RESULTS/*.sta>;
#    my %hash;
#    for my $file (@files) {
#        my $key        = (split /\./,   $file)[6];
#        my $samplename = (split /\/|_/, $file)[2];
#        $hash{$key} = "$file";
#    }
#    my @sorted;
#    for my $key (sort { $a <=> $b } keys %hash) {
#        push(@sorted, $hash{$key});
#    }
#
#    #Blast each of the .sta files
#    print "\nBlasting clusters and parsing the results...\n";
#    my $num_alignments = 250;
#    for my $file (@sorted) {
#        my $hits =
#    `ncbi-blast-2.2.29+/bin/blastn -db ncbi-blast-2.2.29+/db/16SMicrobial -query $file -num_threads 9 -outfmt \"7 qacc sallseqid evalue bitscore pident qstart qend sstart send saccver\" -num_alignments $num_alignments | grep -v \"\#\" | perl ./blast_org_annotate2.pl`;
#
#        #print "hits:\n$hits\n";
#        my @besthits;
#        my @maxhits;
#        my @genuslist;
#        my @specieslist;
#        my $sequence;
#        if (defined $hits && $hits ne '') {
#            my @blastlines = split(/\n/, $hits);
#            my @pidents;
#            my @bitscores;
#            my @blastresults;
#            foreach my $blastline (@blastlines) {
#                @blastresults = split(/\t/, $blastline);
#                push(@pidents,   $blastresults[4]);
#                push(@bitscores, $blastresults[3]);
#            }
#            my @largest      = @bitscores;
#            my $max          = $largest[0];
#            my $bitscoresize = scalar(@bitscores);
#            $max = $_ > $max ? $_ : $max foreach (@largest);
#            $filename = $file;
#            $filename =~ s/$directory\/RESULTS\///;
#            open(CONTENTS, "<", "$file") or die "Can't open \%file: $file $!\n";
#            local $/;
#            my $contents = <CONTENTS>;
#            close(CONTENTS);
#            $sequence = (split /\n/, $contents)[1];
#
#    #add a header indicating the file we are reporting on, consensus seuquence, number of reads represented.
#            if ($readcounts[0] < 1) { $readcounts[0]++ }
#            push(@besthits,
#                "\n\n$filename",
#                "\nReadID and sequence:\n>$readIDarray[0]\n$sequence",
#                "\nWith:\t",
#                $arraycounts[0],
#                " reads out of $readcounts[0] total = ",
#                (sprintf '%.1f', (100 * ($arraycounts[0] / $readcounts[0]))),
#                "% of reads");
#            shift(@readIDarray);
#            my $y = 0;
#            while ($y < $bitscoresize) {
#                if ($bitscores[$y] == $max) {
#                    my @results = split(/\t/, $blastlines[$y]);
#                    push(@maxhits, $results[10]);
#                    if ($pidents[$y] >= 98.00) {
#                        push(@besthits,
#                            "\nFound bitscore ",
#                            $results[3],
#                            " with ",
#                            $results[4],
#                            "% match at the species level (>98%) with: ",
#                            "\t$results[10]\taccession: $results[9]");
#                    }
#                    else {
#                        push(@besthits,
#                            "\nFound bitscore ",
#                            $results[3],
#                            " with ",
#                            $results[4],
#                            "% match at the Genus level (i.e. <98%) with: ",
#                            "\t$results[10]\taccession: $results[9]");
#                    }
#                }
#                $y++;
#                foreach my $result (@maxhits) {
#                    my $genus = (split / /, $result)[0];
#
#                    #$genus = $genus.",";
#                    my $species = (split / /, $result)[1];
#                    my $genusspecies = join " ", $genus, $species;
#                    if (not $genus ~~ @genuslist) {
#                        push(@genuslist,   $genus);
#                        push(@specieslist, $species);
#                    }
#                    elsif (not $species ~~ @specieslist) {
#                        push(@specieslist, $species);
#                    }
#                    else { }
#                }
#            }
#        }
#        else {
#            open(CONTENTS, "<", "$file") or die "Can't open \%file: $file $!\n";
#            local $/;
#            my $contents = <CONTENTS>;
#            close(CONTENTS);
#            $sequence = (split /\n/, $contents)[1];
#            $filename = $file;
#            $filename =~ s/$directory\/RESULTS\///;
#
#            #print "\nIn else loop readcounts = $readcounts[0].\n";
#            if ($readcounts[0] < 1) { $readcounts[0]++ }
#            push(@besthits,
#                "\n\n$filename",
#                "\nReadID and sequence:\n>$readIDarray[0]\n$sequence",
#                "\nWith:\t",
#                $arraycounts[0],
#                " reads out of $readcounts[0] total = ",
#                (sprintf '%.1f', (100 * ($arraycounts[0] / $readcounts[0]))),
#                "% of reads",
#                "\nNo result: no useful sequence or no match in database.\n");
#            push(@genuslist, $file, ":No BLAST result");
#            shift(@readIDarray);
#        }
#        open(RESULTS, ">", "$file.besthits.txt")
#          || die
#    "\n\nCan't open \$file.besthits.txt: $file.besthits.txt for writing $!\n\n";
#        print RESULTS "@besthits";
#        close RESULTS;
#        open(SUMMARY, ">", "$file.summary.txt")
#          || die
#    "\n\nCan't open \$file.besthits.txt: $file.besthits.txt for writing $!\n\n";
#        my @array         = split /.clustered.fasta.recluster.fasta./, $filename;
#        my $sample        = $array[0];
#        my $clusternumber = (split /\./, $array[1])[0];
#        my $readID2       = (split /\./, $array[1])[1];
#        print SUMMARY
#          "\n$sample\t$clusternumber\t$readID2\t$sequence\t$arraycounts[0] reads";
#        shift(@arraycounts);
#        shift(@readcounts);
#
#        if (scalar @genuslist > 1) {
#            print SUMMARY "\tFamily or higher hit\t", join(", ", @genuslist),
#              "\tNot Applicable";
#        }
#        elsif (scalar @specieslist > 1) {
#            print SUMMARY "\tGenus hit\t", @genuslist, "\t";
#            foreach my $species (@specieslist) {
#                print SUMMARY @genuslist, " ", $species, ", ";
#            }
#        }
#        else {
#            print SUMMARY "\tSpecies hit\t", @genuslist, "\t", @genuslist, " ",
#              @specieslist;
#        }
#        close SUMMARY;
#    }
#
#    open(HEADER, ">", "$directory\/RESULTS\/HEADER.txt");
#    print HEADER
#    "$commandline\n\nSimilarity score for clustering: $similarity\nFraction of total reads that a cluster must contain to be included in re-clustering: $percent\nLength of read to be included in clustering: $length\n\n";
#    close HEADER;
#    open(SUMMARYHEADER, ">", "$directory\/RESULTS\/SUMMARYHEADER.txt");
#    print SUMMARYHEADER
#    "$commandline\n\nSimilarity score for clustering: $similarity\nFraction of total reads that a cluster must contain to be included in re-clustering: $percent\nLength of read to be included in clustering: $length\n\nSample\tCluster\tSequence\tTaxonomic level\tgenera/genus matched\tspecies matched(if applicable)";
#    close SUMMARYHEADER;
#    print
#    "\nTop BLAST results by bitscore for all input .fasta files are in $directory/RESULTS/ALL_RESULTS.txt\n";
#    print
#    "\nSummary results for all input .fasta files are in $directory/RESULTS/SUMMARY.txt\n\n";
#    system(
#    "cat $directory/RESULTS/HEADER.txt $directory/RESULTS/*.besthits.txt > $directory/RESULTS/ALL_RESULTS.txt"
#    );
#    system(
#    "cat $directory/RESULTS/SUMMARYHEADER.txt $directory/RESULTS/*.summary.txt > $directory/RESULTS/SUMMARY.txt"
#    );

__END__

# --------------------------------------------------

=pod

=head1 NAME

16blaster.longreads.pl - a script

=head1 SYNOPSIS

  16blaster.longreads.pl -i /path/to/fasta-files

Required Arguments:

  -d|dir            Directory containing FASTA files processing

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

=head1 AUTHOR

George Watts E<lt>gwatts@email.arizona.eduE<gt>,
Ken Youens-Clark E<lt>kyclark@email.arizona.eduE<gt>.

=head1 COPYRIGHT

Copyright (c) 2015 George Watts

This module is free software; you can redistribute it and/or
modify it under the terms of the GPL (either version 1, or at
your option, any later version) or the Artistic License 2.0.
Refer to LICENSE for the full license text and to DISCLAIMER for
additional warranty disclaimers.

=cut
