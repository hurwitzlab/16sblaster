FROM perl:latest

MAINTAINER Ken Youens-Clark <kyclark@email.arizona.edu>

RUN cpanm File::Find::Rule

RUN cpanm File::Which

COPY scripts/16blaster-longreads.pl /usr/local/bin/

COPY scripts/blast_org_annotate2.pl /usr/local/bin/

COPY scripts/16blaster.longreads.pl.orig /usr/local/bin/

COPY bin /usr/local/bin/

COPY blast /data/blast/

COPY data/accession_organism_map.txt /data/

ENTRYPOINT ["16blaster-longreads.pl"]

CMD ["-h"]
