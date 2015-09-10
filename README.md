# 16S BLASTER

Perl code, etc., for George's paper.

# Docker

To build Docker image, first build these programs into the "bin" directory:

* BLAST+ 
* cutadapt (https://pypi.python.org/pypi/cutadapt/)
* Samtools
* cd-hit-est, make_multi_seq.pl

To build:

    docker build -t 16sblaster .

To run:

    docker run --rm -v $(pwd):/work -w /work 16sblaster \
    --blast-db /data/blast -d fasta -o out \
    --accessions /data/accession_organism_map.txt
