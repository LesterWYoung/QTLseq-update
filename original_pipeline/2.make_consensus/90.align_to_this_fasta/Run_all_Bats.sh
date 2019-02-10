#! /bin/sh

./Bat_make_symbolic_link_of_your_consensus.fa.sh
./Bat_make_symbolic_link_of_qualified_fastq_again.sh
./Bat_bwa2bam.sh
./Bat_run_coval-refine-bam.pl.sh
./Bat_run_coval-call-pileup.pl.sh
