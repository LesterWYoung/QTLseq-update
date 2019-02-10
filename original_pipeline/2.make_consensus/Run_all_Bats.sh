#! /bin/sh

./Bat_make_symbolic_link_of_qualified_fastq.sh
./Bat_bwa2bam.sh
./Bat_run_coval-refine-bam.pl.sh
./Bat_run_coval-call-pileup.pl.sh
./Bat_RYKMSWBDHV_to_ACGT.pl.sh
./Bat_make_consensus.sh
