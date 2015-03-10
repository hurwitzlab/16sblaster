#!usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $help = '';
my $directory = '';
my $similarity = '';
my $percent = '';
my $length = '';

GetOptions ("help|h|man|?" => \$help, "i=s" => \$directory, "c=f" => \$similarity, "p=f" => \$percent, "l=i" => \$length) or die "Incorrect usage. Help on using this script is available using the option -h or --help \n";

my $help_text = "

PURPOSE:        This script relies on other programs and scripts: cd-hit, blast_org_annotate.pl and the file: accession_organism_map.txt.
                This script utilizes cd-hit to create clusters of representative bacterial 16s reads and then blasts the representative reads to identify the species that match the read.
              
Input:          required:
                -i <directory> Directory containing .fasta files processing.
		
		optional
		-c <float> Value between 0.80 and 1.00 for percent similarity for cd-hit clustering. Strongly recommend 0.98 (default) to start.
		-p <float> Value between 0 and 1.00 for percent of total reads in order for clusters to be included second round of clustering. Recommend 0.01 (default) to start.
		Though running this script with several reasonable iterations may be necessary to optimize results (e.g. try 0.01, 0.03, 0.10) and then check results to see if they vary.
		Using too high a value will result in no clusters being included in the re-clustering and an error in the cat step where all results are combined (\"No such file or directory)\"
		Using too low a value will result in many redundant results, a better outcome than using too high a number, so start small and work upward.
		-l <integer> Use this option to exlude short reads from being included in the clustering. Default is 25. A reasonable length restriction would be 75% of expected full length read
		of the amplicon in question. This prevents short contaminants or a frequent short read from creating clusters that will be reported in results.

		All .fasta files in the directory provided with the -i option will be processed.
                Data filenames must end with \".fasta\"

                For example, if Ion Torrent IonXpress barcode 17 was used with a sample to generate
                reads from 16s rRNA hypervariable region V1, a good filename would be: \"17.V1.fasta\"

		Remember to remove the primer sequence used to generate the 16s amplicons from the .fasta file(s) before using this script. 

Output:         A directory called: \"RESULTS\".

Usage:          \"perl blaster.pl -i <directory>\" where directory contains .fasta files for blasting.

";

if ( $help ) {
	print $help_text; exit;
}

chomp ( $directory );
$directory =~ s/\/$//; #remove any / the user may have put at the end of the directory given with -i

#check that the user provided directory doesn't contain any periods - this script has a bug that will cause errors if the directory name contains a period.
if ( $directory =~ /\./ ) {
	print "Error: this script cannot use directory names with periods in them, please rename your input directory.\n"; exit;
}

#check that the user provided a directory to work with.
if ( $directory eq '' ) {
        print "\nError: Input directory must be indicated with \"-i <directory>\".\nUsage: \"perl test8.pl -i <directory>\"\nSee more help with the -h option.\nExiting.\n"; exit;
}

#check that the directory provided exists and contains FASTA files using .fasta extension
my @filepath;
if ( -d $directory && -e $directory ) {
	@filepath = <$directory/*.fasta>;
	if ( scalar @filepath == 0 ) {
		print "Directory \"$directory\" does not contain any FASTA (.fasta) files. Exiting...\n"; exit;
	}
} else {
	print "Directory \"$directory\" does not exist or is not a directory. Exiting...\n"; exit;
}

if ( $similarity eq '' ) {
	$similarity = 0.98;
}

if ( $percent eq '' ) {
	$percent = 0.01;
}

if ( $length eq '' ) {
	$length = 25;
}

my $command = `ps -o args -C perl`;
$command =~ s/COMMAND\n//;
chomp $command;
my $commandline = join "", "command line: \"", $command, "\" was run with options: ", " -c: ", $similarity, " -p: ", $percent, " -l: ", $length;

#Delete previous results
system ( "rm -rf $directory\/RESULTS" );
print "\n\tPrevious $directory/RESULTS directory deleted.\n";

#Make a directory for output
system ( "mkdir $directory\/RESULTS" );

#Perform two rounds of clustering using CD-HIT
my $filename;
my @arraycounts;
my @readcounts;

for my $filepath ( @filepath ) {
	$filename = $filepath; $filename =~ s/$directory\///; $filename =~ s/.fasta/.clustered.fasta/;
	print "\n\tClustering:\t\t$filepath in cd-hit-est.\n\tSending output to:\t$directory\/RESULTS\/$filename\n";
	`cd-hit-v4.6.1-2012-08-27/cd-hit-est -i $filepath -o $directory\/RESULTS\/$filename -c $similarity -n 9 -d 0 -g 1 -T 8 -l $length -M 32000 -s 0.98`;
	my $numlines = int(( split (/ /, `wc -l $filepath`))[0] / ( 2*( 1/$percent ))); #use the number of reads in the fasta file to limit the clusters analyzed to those containing a certain percentage of total reads
	`cd-hit-v4.6.1-2012-08-27/make_multi_seq.pl $filepath $directory\/RESULTS\/$filename.clstr $directory\/RESULTS\/$filename.MAKEMULTI$numlines $numlines`;
	system ( "cat $directory\/RESULTS\/$filename.MAKEMULTI$numlines\/* > $directory\/RESULTS\/$filename.contigs.fasta" );
	`cd-hit-v4.6.1-2012-08-27/cd-hit-est -i $directory/RESULTS/$filename.contigs.fasta -o $directory\/RESULTS\/$filename.recluster.fasta -c $similarity -n 9 -d 0 -g 1 -T 8 -l $length -M 32000 -s 0.98`;
	#While we're in this loop, count the size of the clusters (in reads) so we can report that later on
	open( CLSTR, "<", "$directory\/RESULTS\/$filename.recluster.fasta.clstr" ) || die "\n\nCan't open: \$directory\/RESULTS\/\$filename.recluster.fasta.clstr:$directory\/RESULTS\/$filename.recluster.fasta.clstr";
	local $/;
	my $clstr = <CLSTR>;
	close (CLSTR);
	chomp $clstr;
	my @cluster = split (/\n/, $clstr);
	my $elementcount = 0;
	my $number;
	my $count;
	foreach my $part ( @cluster ) {
		if ( $part =~ /^\>/ && $elementcount > 0) {
			$number = (split /\t/, $cluster[$elementcount-1])[0];
			push ( @arraycounts, $number );
			$elementcount++;
			$count++
		}
		else {
			$elementcount++;
		}
	}
	$number = ( split (/\t/, $cluster[$elementcount-1]))[0];
	push ( @arraycounts, $number ); $count++; #add the last count to the @arraycounts			
	while ( $count > 0 ) {
		my $number = ( $elementcount - ( scalar @arraycounts ));
		push ( @readcounts, $number );
		$count = ( $count - 1 );
	}
	print "\tDone Clustering $filepath\n";
}


#After two rounds of clustering we have a fasta file of representative reads, make a .sta (just a fasta with one read) file for each one so they can be blasted one at a time...
@filepath = <$directory/RESULTS/*.recluster.fasta>;
my $i=0; my $readID; my @readIDarray;
for my $filepath ( @filepath ) {
	open( CONTENTS, "<", "$filepath" ) || die "\nCan't open \$filepath: $filepath $!\n";
	local $/;
	my $contents = <CONTENTS>;
	close ( CONTENTS );
	chomp ( $contents );
	my @array = split ( /\>/, $contents );
	shift(@array);
	foreach my $headerandsequence ( @array ) {
		$readID = (split /\n/, $headerandsequence)[0];
		push ( @readIDarray, $readID );
		open ( OUT2, ">", "$filepath.$i.$readID.sta" ) or die "Can't open $filepath.$i.$readID.sta for writing: $!\n";
		print ( OUT2 ">", $filepath, ".", $i, ".", $readID, ".sta\t", $headerandsequence, "\n" );
		close ( OUT2 );
		$i++;
	}
}
my @files = <$directory/RESULTS/*.sta>;
my %hash;
for my $file ( @files ) {
	my $key = (split /\./, $file)[6];
	my $samplename = (split /\/|_/, $file)[2];
	$hash{$key} = "$file";
}
my @sorted;
for my $key ( sort {$a<=>$b} keys %hash ) {
	push ( @sorted, $hash{$key} );
}


#Blast each of the .sta files
print "\nBlasting clusters and parsing the results...\n";
my $num_alignments = 250;
for my $file ( @sorted ) {
	my $hits = `ncbi-blast-2.2.29+/bin/blastn -db ncbi-blast-2.2.29+/db/16SMicrobial -query $file -num_threads 9 -outfmt \"7 qacc sallseqid evalue bitscore pident qstart qend sstart send saccver\" -num_alignments $num_alignments | grep -v \"\#\" | perl ./blast_org_annotate2.pl`;
	#print "hits:\n$hits\n";
	my @besthits;
	my @maxhits;
	my @genuslist;
	my @specieslist;
	my $sequence;
	if ( defined $hits && $hits ne '' ) {
		my @blastlines = split (/\n/, $hits);
		my @pidents;
		my @bitscores;
		my @blastresults;
			foreach my $blastline ( @blastlines ) {
			@blastresults = split (/\t/, $blastline);
			push ( @pidents, $blastresults[4] );
			push ( @bitscores, $blastresults[3] );
		}           
		my @largest = @bitscores;
		my $max = $largest[0];
		my $bitscoresize = scalar ( @bitscores );
		$max= $_>$max ? $_ : $max foreach ( @largest );
		$filename = $file;
		$filename =~ s/$directory\/RESULTS\///;
		open ( CONTENTS, "<", "$file" ) or die "Can't open \%file: $file $!\n";
		local $/;
		my $contents = <CONTENTS>;
		close (CONTENTS);
		$sequence = (split /\n/, $contents)[1];
		#add a header indicating the file we are reporting on, consensus seuquence, number of reads represented.
                if ( $readcounts[0] < 1 ) { $readcounts[0]++ };
		push ( @besthits,
			"\n\n$filename", 
			"\nReadID and sequence:\n>$readIDarray[0]\n$sequence", "\nWith:\t", $arraycounts[0], " reads out of $readcounts[0] total = ", ( sprintf '%.1f', ( 100*( $arraycounts[0]/$readcounts[0] ))), "% of reads" );
		shift ( @readIDarray );
		my $y = 0;
		while ( $y<$bitscoresize ) {
			if ( $bitscores[$y] == $max ) {
				my @results = split (/\t/, $blastlines[$y]);
				push ( @maxhits, $results[10] );
               	                if ( $pidents[$y] >= 98.00 ) {
					push ( @besthits, "\nFound bitscore ", $results[3], " with ", $results[4], "% match at the species level (>98%) with: ", "\t$results[10]\taccession: $results[9]" );
                       	        }
                       	        else {
					push ( @besthits, "\nFound bitscore ", $results[3], " with ", $results[4], "% match at the Genus level (i.e. <98%) with: ", "\t$results[10]\taccession: $results[9]" );
				}
			}
			$y++;
			foreach my $result (@maxhits) {
				my $genus = (split / /, $result)[0];
				#$genus = $genus.","; 
				my $species = (split / /, $result)[1];
				my $genusspecies = join " ", $genus, $species;
				if ( not $genus ~~ @genuslist) {
					push ( @genuslist, $genus);
					push ( @specieslist, $species);
				} elsif ( not $species ~~ @specieslist) {
					push ( @specieslist, $species);
				} else {}		
			}
		}
	} else {
                open ( CONTENTS, "<", "$file" ) or die "Can't open \%file: $file $!\n";
                local $/;
                my $contents = <CONTENTS>;
                close (CONTENTS);
                $sequence = (split /\n/, $contents)[1];
		$filename = $file;
		$filename =~ s/$directory\/RESULTS\///;
		#print "\nIn else loop readcounts = $readcounts[0].\n";
		if ( $readcounts[0] < 1 ) { $readcounts[0]++ };
		push ( @besthits,
			"\n\n$filename",
			"\nReadID and sequence:\n>$readIDarray[0]\n$sequence", "\nWith:\t", $arraycounts[0], " reads out of $readcounts[0] total = ",
			 ( sprintf '%.1f', ( 100*( $arraycounts[0]/$readcounts[0] ))), "% of reads", "\nNo result: no useful sequence or no match in database.\n" );
		push ( @genuslist, $file, ":No BLAST result");
		shift ( @readIDarray );
	}
	open ( RESULTS, ">", "$file.besthits.txt" ) || die "\n\nCan't open \$file.besthits.txt: $file.besthits.txt for writing $!\n\n";
	print RESULTS "@besthits";
	close RESULTS;
	open ( SUMMARY, ">", "$file.summary.txt" ) || die "\n\nCan't open \$file.besthits.txt: $file.besthits.txt for writing $!\n\n";
	my @array = split /.clustered.fasta.recluster.fasta./, $filename;
	my $sample = $array[0];
	my $clusternumber = (split /\./, $array[1])[0];
	my $readID2 = (split /\./, $array[1])[1];
	print SUMMARY "\n$sample\t$clusternumber\t$readID2\t$sequence\t$arraycounts[0] reads";
	shift ( @arraycounts ); shift ( @readcounts );
	if ( scalar @genuslist > 1 ) {
		print SUMMARY "\tFamily or higher hit\t", join(", ", @genuslist), "\tNot Applicable";
	} elsif ( scalar @specieslist > 1) {
		print SUMMARY "\tGenus hit\t", @genuslist, "\t";
		foreach my $species (@specieslist) {
			print SUMMARY @genuslist, " ", $species, ", ";
			}
	} else {
		print SUMMARY "\tSpecies hit\t", @genuslist, "\t", @genuslist, " ", @specieslist;
	}
	close SUMMARY;
}

open ( HEADER, ">", "$directory\/RESULTS\/HEADER.txt" );
print HEADER "$commandline\n\nSimilarity score for clustering: $similarity\nFraction of total reads that a cluster must contain to be included in re-clustering: $percent\nLength of read to be included in clustering: $length\n\n";
close HEADER;
open ( SUMMARYHEADER, ">", "$directory\/RESULTS\/SUMMARYHEADER.txt" );
print SUMMARYHEADER "$commandline\n\nSimilarity score for clustering: $similarity\nFraction of total reads that a cluster must contain to be included in re-clustering: $percent\nLength of read to be included in clustering: $length\n\nSample\tCluster\tSequence\tTaxonomic level\tgenera/genus matched\tspecies matched(if applicable)";
close SUMMARYHEADER;
print "\nTop BLAST results by bitscore for all input .fasta files are in $directory/RESULTS/ALL_RESULTS.txt\n";
print "\nSummary results for all input .fasta files are in $directory/RESULTS/SUMMARY.txt\n\n";
system ( "cat $directory/RESULTS/HEADER.txt $directory/RESULTS/*.besthits.txt > $directory/RESULTS/ALL_RESULTS.txt" );
system ( "cat $directory/RESULTS/SUMMARYHEADER.txt $directory/RESULTS/*.summary.txt > $directory/RESULTS/SUMMARY.txt" );

