#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Pod::Usage;
use File::Tee qw(tee);

my $reference_file = '';         # -r
my $read_type = 'PE';		 # -t
my $limited_mismatch_num = 2;	 # -n
my $limited_mismatch_num_pairs = 0;	# -f
my $limited_mismatch_rate = 0;	 # -m
my $minimum_mapq = 0;            # -mq
my $base_call_type = 'auto';     # -e
my $min_mismatch_qual = 41;	 # -k
my $minimum_ave_alt_qual = 10;   # -xa
my $minimum_allele_freq = 0;     # -xf
my $ave_insert_size_SD = 'auto'; # -a
my $maximum_insert_SD_fold = 5;	 # -i
my $multi_sample_mode = 'none';  # -ms
my $filter_pairs = 0;	 	 # -g
my $dis_realign = 0;	         # -d
my $read_correct = 0;            # -x
my $targeted_align_mode = 0;     # -tm
my $error_correct_mode = 0;      # -xm
my $out_sam = 0;		 # -s
my $reference_type = 'COMPLETE'; # -rt
my $include_disc_read = 0;       # -c
my $include_unmap_read = 0;	 # -u
my $long_insert_read = 0;	 # -l
my $out_prefix = 'coval-filter'; # -p
my $help;			 # -h
my $soap_aligner = 0;            # -b

my $multiple_allele_sample = 0;  # -ma

my $read_length = 0;
my $ave_read_length = 0;

my $max_mismatch_num = 0;
my $MDtag_sam = 0;

GetOptions(
    'ref|r=s' => \$reference_file,
    'num|n=i' => \$limited_mismatch_num,
    'fnum|f=i' => \$limited_mismatch_num_pairs,
    'mrate|m=f' => \$limited_mismatch_rate,
    'mapq|mq=i' => \$minimum_mapq,
    'mallel|ma' => \$multiple_allele_sample,
    'qcal|e=s' => \$base_call_type,
    'minq|k=i' => \$min_mismatch_qual,
    'qual_ave|xa=i' => \$minimum_ave_alt_qual,
    'mfreq|xf=f' => \$minimum_allele_freq,
    'avins|a=s' => \$ave_insert_size_SD,
    'insd|i=i' => \$maximum_insert_SD_fold,
    'type|t=s' => \$read_type,
    'fpair|g' => \$filter_pairs,
    'msamp|ms=s' => \$multi_sample_mode,
    'dis_align|d' => \$dis_realign,
    'err_cor|x' => \$read_correct,
    'ec_md|xm' => \$error_correct_mode,
    'talign_md|tm' => \$targeted_align_mode,
    'sam|s' => \$out_sam,
    'reftype|rt=s' => \$reference_type,
    'cdisc|c' => \$include_disc_read,
    'unmap|u' => \$include_unmap_read,
    'lins|l' => \$long_insert_read,
    'pref|p=s' => \$out_prefix,
    'soap|b' => \$soap_aligner,
    'help|h' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

$limited_mismatch_num_pairs = $limited_mismatch_num * 1.7 if ($limited_mismatch_num_pairs == 0);

if ($error_correct_mode == 1){
    $read_correct = 1;
    $filter_pairs = 1;
}
if ($limited_mismatch_num_pairs < $limited_mismatch_num){
    die "The number specified with --fnum should be equal or larger than the value specified with --num.\n$!";
}

if ($targeted_align_mode == 1){
    if ($read_correct == 1){
	print STDERR "\n########## Warning! ##########\nThe job is changed from the error correction mode to the basic mode when the --talign option is set\n";
	$read_correct = 0;
    }
    $filter_pairs = 1;
    $limited_mismatch_num_pairs = $limited_mismatch_num;
}

if ($multi_sample_mode eq 'homo'){
    $read_correct = 1;
    $multi_sample_mode = 1;
    $minimum_allele_freq = 0.8 if ($minimum_allele_freq == 0);
}
elsif ($multi_sample_mode eq 'hetero'){
    $read_correct = 1;
    $multi_sample_mode = 1;
    $minimum_allele_freq = 0.3 if ($minimum_allele_freq == 0);
}
elsif ($multi_sample_mode eq 'none'){
    $multi_sample_mode = 0;
}

$include_disc_read = 1 if ($reference_type eq 'DRAFT') or ($reference_type eq 'draft');

unless (@ARGV){
print "Your sam/bam file name should be added to the argument.\n";
exit;
}
my $input_file = shift (@ARGV);
if ($input_file eq '-'){
}
elsif (!-f $input_file){
    die "inputfile does not exist: $!\n";
}

die "Reference file is not found: $!\n" unless (-f $reference_file);
die "Prefix of output is not specified: $!\n" if ($out_prefix eq '');
die "Read type is not properly specified: $!\n" unless (uc $read_type eq 'PE') or (uc $read_type eq 'SE');
    
open(my $log_fh, '>', "$out_prefix.log");
tee(STDERR, $log_fh);


=head1 SYNOPSIS

 coval refine [options] <input_sorted_bam/sam_file> (Read names and reference sequence names in sam file should not contain characters '||' or '='.)
  Output:
   out_prefix.bam/sam

  Options:
   --ref or -r <STR>	   reference fasta file used for the alignment
   --pref or -p <STR>      prefix of output file
   --num or -n <INT>       maximum number of mismatches contained in a read [default: 2] (incompatible with --mrate)
   --mrate or -m <FLOAT>   maximum rate of mismatches contained in a read [0..1.0] (incompatible with --num)
   --fnum or -f <INT>	   maximum number of total mismatches contained in two paired reads [default: 1.7 * <INT> specified with --num] (incompatible with --mrate)
   --fpair or -g	   filter out the other pair of a filtered paired-end read [default: false]
   --mapq or -mq <INT>     minimum mapping quality not filtered [default: 0]
   --mflq or -fq <FLOAT>   maximum allele frequency for removing reads with a low mapping quality (specified with -mq) [default: 0.8]
   --qcall or -e <STR>	   base call quality format of reads represented in sam/bam file; illumina (Phred+64) or sanger (Phred+33) [default: auto]
   --minq or -k <INT>	   minimum base-call quality of mismatch base (counted as 'mismatch' only mismatched bases with the base quality lower than the specified value) [default: 41]
   --avins or -a <STR>     average distance of paired-end reads (between mapped 5' ends) and its SD, separated with comma [default: auto]
   --insd or -i <INT>      value M for setting a maximum read distance to be filtered: mean read distance + SD * M [default: 5]
   --type or -t <STR>      read type, paired-end (PE) or single-end (SE) [default: PE]
   --lins or -l		   input reads contain mate-pair or paired-end reads with a long insert size [default: false]
   --sam or -s		   output is a sam file [default: false]
   --reftype or -rt <COMPLETE | DRAFT> force --cdisc option when specifying DRAFT; DRAFT: draft-level of genomes or cDNAs/ESTs, COMPLETE: finished genome sequences [default: COMPLETE]
   --cdisc or -c	   include (not filter) discordant paired-end reads (either mate is improperly aligned) [default: false]
   --unmap or -u           include (not filter) unmapped reads [default: false]
   --soap or -b            when alignment was generated with SOAP aligner, use this option [default: false]
   --dis_align or -d       disable realignmnet [default: false]
   --err_cor or -x         correct potential sequencing errors in reads [default: false]
   --help or -h            output help message
   
   << only when -x is set >>
   --qave or -xa <INT>     minimum mean base-call quality of mismatch bases covered at the non-reference position in error correction mode [default: 10]
   --mfreq or -xf <FLOAT>  minimum allele frequency at a potential non-reference allele in error correction mode [0..1.0] [default: 0]
			   (Mismatch bases with lower quality and lower frequency than the values specified with -xb and -xf are corrected.)
   --msamp or -ms <STR>    multi-sample mode; calculate allele frequency for each sample based on RG tag. Read bases lower than the allele frequency specified with -xf are corrected. < none | homo | hetero >
                          (homo: homozygous sample, hetero: heterozygous sample, -x -xf 0.8 and -x -xf 0.3 are automatically set, respectively. The -xf value can be overwritten by setting -xf separately.) [default: none]
   --mallel or -ma         allow multiple non-reference alleles for each sample in multi-sample mode (e.g., alleles A, C, and T for samples-1, -2, and -3) [default: false]
   
   ### Preset options ###
   --ec_md or -xm	   set options suitable for 'error correction mode' (equal to '--err_cor --fpair') [default: false]
   --talign_md or -tm      set options suitable for 'targeted' alignment (equal to '--fpair --fnum <INT specified with --num>') [default: false]

=cut

my %chr_seq;
open (FILE, $reference_file) or die "$!";
    my $ref_chr_name;
    my $count_chr = 0;
    my $chrseq = '';
    while (my $line = <FILE>){
        chomp $line;
	if ($line =~ /^>(\S*)/){
	    if ($count_chr > 0){
		$chr_seq{$ref_chr_name} = $chrseq;
		$chrseq = '';
	    }
	    $ref_chr_name = $1;
	    if (($ref_chr_name =~ /==/) or ($ref_chr_name =~ /\|\|/)){
		die "Reference chromosome/contig names contain unacceptable characters (== or \|\|): $!";
	    }
	    $count_chr ++;
	}
	else{
	    $chrseq .= uc $line;
	}
    }
    $chr_seq{$ref_chr_name} = $chrseq;
    $chrseq = '';
close (FILE);

my $terminal_mismatch = 0;
my $read_number = 0;
my $read_number_2 = 0;
my $discordant_read_num = 0;
my $filtered_read_num = 0;
my $count_realigned_indel = 0;
my $count_realigned_read = 0;
my $count_indel = 0;
my $sum_read_len = 0;
my $chrom;
my $pre_chr = '-';
my %mismatch_read;
my %filtered_mismatch_read;
my %match_read;
my %indel_read;
my %read_pair;
my %change_pos_read_F;
my %change_pos_read_R;
my %filtered_pairs;
my %num_mismatch_pair;
my @unmapped;

my ($m1, $m2, $m3, $m4, $m5, $m6, $m7, $m8, $m9, $m10, $m10o) = (0,0,0,0,0,0,0,0,0,0,0);
my ($n1, $n2, $n3, $n4, $n5, $n6, $n7, $n8, $n9, $n10, $n10o) = (0,0,0,0,0,0,0,0,0,0,0);
my $ave_dist = 0;
my $SD_dist = 0;
if ($ave_insert_size_SD =~ /^(\d+),(\d+)/){
    $ave_dist = $1;
    $SD_dist = $2;
    print STDERR 'Average read distance = ', $ave_dist, ' +/- ', $SD_dist, "\n";
}

my $qual_seq;
my $sanger_type;
my $solexa_type;
my $RGtag_sam = 0;
open (FILE, $input_file) or die "$input_file: $!";
    while (<FILE>){
	next if ($_ =~ /^\@/);
        my @line = split(/\s+/, $_);
        next if ($line[5] =~ /\*/);
        $qual_seq .= $line[10];
	$MDtag_sam = 1 if ($_ =~ /MD:Z:[\dACGT]/);
	$RGtag_sam = 1 if ($_ =~ /RG:Z:\S+/);
        last if (length $qual_seq >= 1000);
    }
    $sanger_type += $qual_seq =~ /[!Ó#%&'\(\)\*,\.\/0123456789:;<=>]/g;
    $solexa_type += $qual_seq =~ /JKLMNOPQRSTUVWXYZ\[\\\]\^_`abcdefgh/g;
close (FILE);
die "RG tag must be added to your input alignmen file when specifying -ms option:\n$!" if ($multi_sample_mode == 1) and ($RGtag_sam == 0);

if ($sanger_type > $solexa_type){
    if ($base_call_type eq 'illumina'){
        print STDERR "\n########## Warning! ##########\nYour base quality format appears a sanger type, but the job will be processed with illumina type.\n To change it use -qcall option.\n";
    }
    else{
	$base_call_type = 'sanger';
	print STDERR "Base call quality format: sanger type\n";
    }
}
else{
    if ($base_call_type eq 'sanger'){
        print STDERR "\n########## Warning! ##########\nYour base quality format appears a illumina type, but the job will be processed with sanger type.\n To change use -qcall option.\n";
    }
    else{
	$base_call_type = 'illumina';
	print STDERR "Base call quality format: illumina type\n";
    }
}

my $realignment = 'ON' if ($dis_realign == 0);
$realignment = 'OFF' if ($dis_realign == 1);

open (FILE, $input_file) or die "$input_file: $!";
    while (my $line = <FILE>){
	if ($line =~ /^\@/){
	    print $line;
	    next;
	}
        chomp $line;
        my @line = split (/\s+/, $line);
	next if (@line < 11);
	if ($line[5] =~ /\*/){
	    push (@unmapped, $line) if ($include_unmap_read == 1);
	    next;
	}
	$read_number++;
	if (($line[5] =~ /H/) or ($line[5] =~ /^\d+S.+S$/) or ($line[5] =~ /^\d+S\d+M\d+[DI]/) or ($line[5] =~ /\d+[DI]\d+M\d+S$/) or ($line[5] =~ /\d+[DI]\d+M\d+[DI]\d+M\d+[DI]/)){
	    if ($filter_pairs == 1){
		my $read_name = $line[0];
		$filtered_pairs{$read_name} = 1;
	    }
	    next;
	}
	if ($minimum_mapq > 0){
	    next if $line[4] < $minimum_mapq;
	}
        $line[2] =~ /^(\S*)/;
        my $chrom = $1;
	$pre_chr = $chrom if ($pre_chr eq '-');
        my $pos = $line[3];
        my $pos09d = sprintf ("%09d", $pos);
        my $read_id = $line[0];
	my $sam_flag = $line[1];
	my $direct = 1;
	$read_id = $1 if ($read_id =~ /(\S+)\/[12]$/);
        my $chr_pos = $chrom . '==' . $pos09d;
	my $mismatch_num = 0;
	$read_length = length $line[9];
	$max_mismatch_num = int ($read_length * 0.1);
	$line[9] = uc $line[9];
	if ($limited_mismatch_rate != 0){
	    $limited_mismatch_num = int ($read_length * $limited_mismatch_rate + 0.5);
	    if ($targeted_align_mode == 0){
		$limited_mismatch_num_pairs = $limited_mismatch_num * 1.7;
	    }
	    else{
		$limited_mismatch_num_pairs = $limited_mismatch_num;
	    }
	}
	$read_number_2++;
	if ($read_number_2 == 1){
	    if (($read_id =~ /==/) or ($read_id =~ /\|\|/)){
		die "Read names contain unacceptable characters (== or \|\|): $!";
	    }
	}
	$sum_read_len += $read_length if ($read_number_2 <= 100);
	if ($read_number_2 == 100){
	    $ave_read_length = int ($sum_read_len / 100 + 0.5);
	}
        
        # remove paired reads with an extraordinary insert size and orphand reads
	if ($read_type eq 'PE'){
	    $discordant_read_num ++ if (($soap_aligner == 0) and ($long_insert_read == 0) and (($sam_flag == 97) or ($sam_flag == 161) or ($sam_flag == 81) or ($sam_flag == 145) or ($sam_flag == 129) or ($sam_flag == 65) or ($sam_flag == 113) or ($sam_flag == 177)
					or ($line[8] == 0)));
	    $discordant_read_num ++ if ((($soap_aligner == 1) or ($long_insert_read == 1)) and (($sam_flag == 129) or ($sam_flag == 65) or ($sam_flag == 113) or ($sam_flag == 177) or ($line[8] == 0)));
	    $discordant_read_num ++ if ((($line[8] < -($ave_dist + $maximum_insert_SD_fold * $SD_dist)) or ($line[8] > $ave_dist + $maximum_insert_SD_fold * $SD_dist)) and ($ave_dist > 0));
	    $direct = 2 if ($sam_flag == 83) or ($sam_flag == 83) or ($sam_flag == 147) or ($sam_flag == 81) or ($sam_flag == 145) or ($sam_flag == 113) or ($sam_flag == 177);
	    
	    if ($include_disc_read == 0){
		if ($line[8] == 0){
		    if ($filter_pairs == 1){
			my $read_name = $line[0];
			$filtered_pairs{$read_name} = 1;
		    }
		    next;
		}
		if ($sam_flag >= 1000){
		    if ($filter_pairs == 1){
			my $read_name = $line[0];
			$filtered_pairs{$read_name} = 1;
		    }
		    next;
		}
		if ($maximum_insert_SD_fold <= 5){
		    if (($soap_aligner == 0) and ($long_insert_read == 0) and (($sam_flag == 97) or ($sam_flag == 161) or ($sam_flag == 81) or ($sam_flag == 145) or ($sam_flag == 129) or ($sam_flag == 65) or ($sam_flag == 113) or ($sam_flag == 177))){
			if ($filter_pairs == 1){
			    my $read_name = $line[0];
			    $filtered_pairs{$read_name} = 1;
			}
			next;
		    }
		    elsif ((($soap_aligner == 1) or ($long_insert_read == 1)) and (($sam_flag == 129) or ($sam_flag == 65) or ($sam_flag == 113) or ($sam_flag == 177))){
			if ($filter_pairs == 1){
			    my $read_name = $line[0];
			    $filtered_pairs{$read_name} = 1;
			}
			next;
		    }
		}
		if ((($line[8] < -($ave_dist + $maximum_insert_SD_fold * $SD_dist)) or ($line[8] > $ave_dist + $maximum_insert_SD_fold * $SD_dist)) and ($ave_dist > 0)){
		    if ($filter_pairs == 1){
			my $read_name = $line[0];
			$filtered_pairs{$read_name} = 1;
		    }
		    next;
		}
	    }
	}
        my $chr_pos_id = $chrom . '==' . $pos09d . '==' . $read_id . '==' . $direct;
	
        if ($pre_chr ne $chrom){
	    my %mismatch_pos;
	    my %mismatch_pos_ms;
	    my %mismatch_pos_id;
	    my %mismatch_pos_id_ms;
	    my %read_depth;
	    my %read_depth_ms;
	    if ($realignment eq 'ON'){
		print STDERR "Realigning reads around indels in $pre_chr ... \n";
		&realign_indel ();
		if ($read_correct == 0){
		    foreach my $key (keys %filtered_mismatch_read){
			unless (exists $match_read{$key}){
			    $match_read{$key} = $filtered_mismatch_read{$key};
			}
		    }
		}
	    }
	    if ($read_correct == 1){				# for read error correction and/or remove alignments with low mapping quality at heterogygous allele
		if ($read_correct == 1){
		    print STDERR "Correcting potential errors in reads ... \n" if ($dis_realign == 0);
		    print STDERR "Correcting potential errors in reads in $pre_chr ... \n" if ($dis_realign == 1);
		}
		if ($minimum_allele_freq > 0){
		    foreach my $key (keys %match_read){			# record read number of mismatch-free reads at each chromosome position
			my ($ch, $ps, $id) = split (/==/, $key);
			$ps =~ s/^0*//;
			my $MDtag = $1 if ($match_read{$key} =~ /MD:Z:(\S+)/);
			my $RGtag = $1 if ($match_read{$key} =~ /RG:Z:(\S+)/) and ($multi_sample_mode == 1);
			my $start_pos = 0;
			while ($MDtag =~ /(\d+)([A-Z]*)(\^[A-Z]+)*/g){
			    for (my $i = $ps + $start_pos; $i < $ps + $start_pos + $1; $i++){
				$read_depth{$i} ++;
				${$read_depth_ms{$i}}{$RGtag} ++ if ($multi_sample_mode == 1);
			    }
			    $start_pos += $1;
			    my $del = $3;	
			    if ($2 ne ''){
				my $mismatch_pos = $ps + $start_pos;
				$read_depth{$mismatch_pos} ++;
				${$read_depth_ms{$mismatch_pos}}{$RGtag} ++ if ($multi_sample_mode == 1);
				$start_pos ++;
			    }
			    if ((defined $del) and ($del =~ /\^([A-Z]+)/)){
				$start_pos += length ($1);
			    }
			}
		    }
		}
		foreach my $key (keys %filtered_mismatch_read){		# record reference position and base quality of mismatch base of mismatch-containing reads and the number of the same mismatch base covered at the site
		    my ($ch, $ps, $id, $dir) = split (/==/, $key);		# record read number at each chromosome position
		    my @mismatch_line = split (/\s+/, $filtered_mismatch_read{$key});
		    my $pos = $mismatch_line[3];
		    my $pos_09d = sprintf ("%09d", $pos);
		    my $MDtag = $1 if ($filtered_mismatch_read{$key} =~ /MD:Z:(\S+)/);
		    my $RGtag = $1 if ($filtered_mismatch_read{$key} =~ /RG:Z:(\S+)/) and ($multi_sample_mode == 1);
		    my $CIAGR = $mismatch_line[5];
		    my $start_pos = 0;
		    my $read_start_pos = 0;
		    while ($MDtag =~ /(\d+)([A-Z]*)(\^[A-Z]+)*/g){
			if ($minimum_allele_freq > 0){
			    for (my $i = $pos + $start_pos; $i < $pos + $start_pos + $1; $i++){
				$read_depth{$i} ++;
				${$read_depth_ms{$i}}{$RGtag} ++ if ($multi_sample_mode == 1) and ($minimum_allele_freq > 0);
			    }
			}
			$start_pos += $1;
			$read_start_pos += $1;
			my $del = $3;	
			if ($2 ne ''){
			    my $mismatch_pos = $pos + $start_pos;
			    $start_pos ++;
			    $read_start_pos ++;		    
			    my ($mismatch_base, $mismatch_qual) = &correct_read ($mismatch_line[9], $CIAGR, $read_start_pos, $2, 'rbase', $mismatch_line[10]);		
			    ${$mismatch_pos{$mismatch_pos}}{$mismatch_base} .= $mismatch_qual;
			    ${${$mismatch_pos_ms{$mismatch_pos}}{$RGtag}}{$mismatch_base} .= $mismatch_qual if ($multi_sample_mode == 1);
			    $read_depth{$mismatch_pos} ++ if ($minimum_allele_freq > 0);
			    ${$read_depth_ms{$mismatch_pos}}{$RGtag} ++ if ($multi_sample_mode == 1) and ($minimum_allele_freq > 0);
			    if (!exists ${$mismatch_pos_id{$mismatch_pos}}{$mismatch_base}){
				${$mismatch_pos_id{$mismatch_pos}}{$mismatch_base} = $id . '==' . $pos_09d . '==' . $read_start_pos . '==' . $dir;
			    }
			    else{
				${$mismatch_pos_id{$mismatch_pos}}{$mismatch_base} .= '||' . $id . '==' . $pos_09d . '==' . $read_start_pos . '==' . $dir;
			    }
			    if ($multi_sample_mode == 1){
				if (!exists ${${$mismatch_pos_id_ms{$mismatch_pos}}{$RGtag}}{$mismatch_base}){
				    ${${$mismatch_pos_id_ms{$mismatch_pos}}{$RGtag}}{$mismatch_base} = $id . '==' . $pos_09d . '==' . $read_start_pos . '==' . $dir;
				}
				else{
				    ${${$mismatch_pos_id_ms{$mismatch_pos}}{$RGtag}}{$mismatch_base} .= '||' . $id . '==' . $pos_09d . '==' . $read_start_pos . '==' . $dir;
				}
			    }
			}
			if ((defined $del) and ($del =~ /\^([A-Z]+)/)){
			    $start_pos += length ($1);
			}
		    }
		    if ($ps ne $pos_09d){
			my $new_key = $ch . '==' . $pos_09d . '==' . $id . '==' . $dir;
			$filtered_mismatch_read{$new_key} = $filtered_mismatch_read{$key};
			delete $filtered_mismatch_read{$key};
		    }
		}
		my %top_allele_pos;
		if (($multi_sample_mode == 0) or (($multi_sample_mode == 1) and ($multiple_allele_sample == 0))){
		    foreach my $mismatch_pos (keys %mismatch_pos){
			my $allele_1 = '';
			my $allele_2 = '';
			my $top_allele = '';
			my $count = 0;
			foreach my $allele (sort {length ${$mismatch_pos{$mismatch_pos}}{$b} <=> length ${$mismatch_pos{$mismatch_pos}}{$a}} keys %{$mismatch_pos{$mismatch_pos}}){
			    $count++;
			    $allele_1 = $allele if ($count == 1);
			    $allele_2 = $allele if ($count == 2);
			}
			if (($allele_2 eq '') or (($allele_2 ne '') and (length ${$mismatch_pos{$mismatch_pos}}{$allele_1} > length ${$mismatch_pos{$mismatch_pos}}{$allele_2}))){
			    $top_allele = $allele_1;
			    my $sum_alt_qual_1 = 0;
			    my $ave_alt_qual = 0;
			    my @alt_qual_1 = split (//, ${$mismatch_pos{$mismatch_pos}}{$allele_1});
			    foreach (@alt_qual_1){
				my $base_qual = ord ($_) - 64 if ($base_call_type eq 'illumina');
				$base_qual = ord ($_) - 33 if ($base_call_type eq 'sanger');
				$sum_alt_qual_1 += $base_qual;
			    }
			    $ave_alt_qual = $sum_alt_qual_1 / @alt_qual_1 if (@alt_qual_1 > 0);
			    if ($ave_alt_qual < $minimum_ave_alt_qual){
				$top_allele = '';
			    }
			}
			elsif (($allele_2 ne '') and (length ${$mismatch_pos{$mismatch_pos}}{$allele_1} == length ${$mismatch_pos{$mismatch_pos}}{$allele_2})){
			    my $sum_alt_qual_1 = 0;
			    my $sum_alt_qual_2 = 0;
			    my $ave_alt_qual = 0;
			    my @alt_qual_1 = split (//, ${$mismatch_pos{$mismatch_pos}}{$allele_1});
			    my @alt_qual_2 = split (//, ${$mismatch_pos{$mismatch_pos}}{$allele_2});
			    foreach (@alt_qual_1){
				my $base_qual = ord ($_) - 64 if ($base_call_type eq 'illumina');
				$base_qual = ord ($_) - 33 if ($base_call_type eq 'sanger');
				$sum_alt_qual_1 += $base_qual;
			    }
			    foreach (@alt_qual_2){
				my $base_qual = ord ($_) - 64 if ($base_call_type eq 'illumina');
				$base_qual = ord ($_) - 33 if ($base_call_type eq 'sanger');
				$sum_alt_qual_2 += $base_qual;
			    }
			    if ($sum_alt_qual_1 >= $sum_alt_qual_2){
				$top_allele = $allele_1;
				$ave_alt_qual = $sum_alt_qual_1 / @alt_qual_1 if (@alt_qual_1 > 0);
			    }
			    else{
				$top_allele = $allele_2;
				$ave_alt_qual = $sum_alt_qual_2 / @alt_qual_2 if (@alt_qual_2 > 0);
			    }
			    if ($ave_alt_qual < $minimum_ave_alt_qual){
				$top_allele = '';
			    }
			}
			$top_allele_pos{$mismatch_pos} = $top_allele;
			if ($multi_sample_mode == 0){
			    foreach my $allele (keys %{$mismatch_pos{$mismatch_pos}}){
				my $allele_freq = length (${$mismatch_pos{$mismatch_pos}}{$allele}) / $read_depth{$mismatch_pos} if ($minimum_allele_freq > 0);
				if (($allele ne $top_allele) or (($allele eq $top_allele) and ($minimum_allele_freq > 0) and ($allele_freq < $minimum_allele_freq))){
				    my @id_pos = split (/\|\|/, ${$mismatch_pos_id{$mismatch_pos}}{$allele});
				    foreach (@id_pos){
					my ($id, $ref_pos, $read_pos, $dir) = split (/==/, $_);
					my $read_key = $pre_chr . '==' . $ref_pos . '==' . $id . '==' . $dir;
					next if (!exists $filtered_mismatch_read{$read_key});
					my @mismatch_line = split (/\s+/, $filtered_mismatch_read{$read_key});
					my $ref_base = substr ($chr_seq{$pre_chr}, $mismatch_pos - 1, 1);
					my $corrected_read = &correct_read ($mismatch_line[9], $mismatch_line[5], $read_pos, $ref_base, 'correct');
					$mismatch_line[9] = $corrected_read;
					$filtered_mismatch_read{$read_key} = join ("\t", @mismatch_line);
					my $MDtag = $1 if ($filtered_mismatch_read{$read_key} =~ /MD:Z:(\S+)/);
					my $start_pos = 0;
					my $pre_match_len = 0;
					my $pre_del = '';
					my $MDtag_cor = '';
					while ($MDtag =~ /(\d+)([A-Z]*)(\^[A-Z]+)*/g){
					    $start_pos += $1;
					    if ($pre_match_len > 0){
						if ($pre_del eq ''){
						    $MDtag_cor .= ($pre_match_len + $1)
						}
						if ($2 ne ''){
						    $start_pos ++;
						    $MDtag_cor .= $2 if (!defined $3);
						    $MDtag_cor .= $2 . $3 if (defined $3);
						}
						else{
						    $MDtag_cor .= $3 if (defined $3);
						}
						$pre_match_len = 0;
					    }
					    else{
						if ($2 ne ''){
						    $start_pos ++;
						    if ($start_pos == $read_pos){
							if (defined $3){
							    $pre_del = $3;
							    $MDtag_cor .= ($1 + 1) . $3;
							}
							else{
							    $pre_match_len = $1 + 1;
							}
						    }
						    else{
							$MDtag_cor .= $1 . $2 if (!defined $3);
							$MDtag_cor .= $1 . $2 . $3 if (defined $3);
						    }
						}
						else{
						    $MDtag_cor .= $1 if (!defined $3);
						    $MDtag_cor .= $1 . $3 if (defined $3);
						}
					    }
					}
					$filtered_mismatch_read{$read_key} =~ s/MD:Z:\S+/MD:Z:$MDtag_cor/;				
				    }
				}
			    }
			}
		    }
		}
		if ($multi_sample_mode == 1){
		    foreach my $mismatch_pos (keys %mismatch_pos_ms){
			foreach my $RG (keys %{$mismatch_pos_ms{$mismatch_pos}}){
			    my $allele_1 = '';
			    my $allele_2 = '';
			    my $top_allele = '';
			    if ($multiple_allele_sample == 1){
				my $count = 0;
				foreach my $allele (sort {length ${${$mismatch_pos_ms{$mismatch_pos}}{$RG}}{$b} <=> length ${${$mismatch_pos_ms{$mismatch_pos}}{$RG}}{$a}} keys %{${$mismatch_pos_ms{$mismatch_pos}}{$RG}}){
				    $count++;
				    $allele_1 = $allele if ($count == 1);
				    $allele_2 = $allele if ($count == 2);
				}
				if (($allele_2 eq '') or (($allele_2 ne '') and (length ${${$mismatch_pos_ms{$mismatch_pos}}{$RG}}{$allele_1} > length ${${$mismatch_pos_ms{$mismatch_pos}}{$RG}}{$allele_2}))){
				    $top_allele = $allele_1;
				    my $sum_alt_qual_1 = 0;
				    my $ave_alt_qual = 0;
				    my @alt_qual_1 = split (//, ${${$mismatch_pos_ms{$mismatch_pos}}{$RG}}{$allele_1});
				    foreach (@alt_qual_1){
					my $base_qual = ord ($_) - 64 if ($base_call_type eq 'illumina');
					$base_qual = ord ($_) - 33 if ($base_call_type eq 'sanger');
					$sum_alt_qual_1 += $base_qual;
				    }
				    $ave_alt_qual = $sum_alt_qual_1 / @alt_qual_1 if (@alt_qual_1 > 0);
				    if ($ave_alt_qual < $minimum_ave_alt_qual){
					$top_allele = '';
				    }
				}
				elsif (($allele_2 ne '') and (length ${${$mismatch_pos_ms{$mismatch_pos}}{$RG}}{$allele_1} == length ${${$mismatch_pos_ms{$mismatch_pos}}{$RG}}{$allele_2})){
				    my $sum_alt_qual_1 = 0;
				    my $sum_alt_qual_2 = 0;
				    my $ave_alt_qual = 0;
				    my @alt_qual_1 = split (//, ${${$mismatch_pos_ms{$mismatch_pos}}{$RG}}{$allele_1});
				    my @alt_qual_2 = split (//, ${${$mismatch_pos_ms{$mismatch_pos}}{$RG}}{$allele_2});
				    foreach (@alt_qual_1){
					my $base_qual = ord ($_) - 64 if ($base_call_type eq 'illumina');
					$base_qual = ord ($_) - 33 if ($base_call_type eq 'sanger');
					$sum_alt_qual_1 += $base_qual;
				    }
				    foreach (@alt_qual_2){
					my $base_qual = ord ($_) - 64 if ($base_call_type eq 'illumina');
					$base_qual = ord ($_) - 33 if ($base_call_type eq 'sanger');
					$sum_alt_qual_2 += $base_qual;
				    }
				    if ($sum_alt_qual_1 >= $sum_alt_qual_2){
					$top_allele = $allele_1;
					$ave_alt_qual = $sum_alt_qual_1 / @alt_qual_1 if (@alt_qual_1 > 0);
				    }
				    else{
					$top_allele = $allele_2;
					$ave_alt_qual = $sum_alt_qual_2 / @alt_qual_2 if (@alt_qual_2 > 0);
				    }
				    if ($ave_alt_qual < $minimum_ave_alt_qual){
					$top_allele = '';
				    }
				}
			    }
			    foreach my $allele (keys %{${$mismatch_pos_ms{$mismatch_pos}}{$RG}}){
				$top_allele = $top_allele_pos{$mismatch_pos} if ($multiple_allele_sample == 0);
				my $allele_freq = length (${${$mismatch_pos_ms{$mismatch_pos}}{$RG}}{$allele}) / ${$read_depth_ms{$mismatch_pos}}{$RG} if ($minimum_allele_freq > 0);
#print STDERR "$mismatch_pos\t$allele\t$top_allele\t$allele_freq\t$pre_chr\t$RG\n" if ($pre_chr =~ /chr01/);
				if (($allele ne $top_allele) or (($allele eq $top_allele) and ($minimum_allele_freq > 0) and ($allele_freq < $minimum_allele_freq))){
				    my @id_pos = split (/\|\|/, ${${$mismatch_pos_id_ms{$mismatch_pos}}{$RG}}{$allele});
				    foreach (@id_pos){
					my ($id, $ref_pos, $read_pos, $dir) = split (/==/, $_);
					my $read_key = $pre_chr . '==' . $ref_pos . '==' . $id . '==' . $dir;
					next if (!exists $filtered_mismatch_read{$read_key});
					my @mismatch_line = split (/\s+/, $filtered_mismatch_read{$read_key});
					my $ref_base = substr ($chr_seq{$pre_chr}, $mismatch_pos - 1, 1);
					my $corrected_read = &correct_read ($mismatch_line[9], $mismatch_line[5], $read_pos, $ref_base, 'correct');
					$mismatch_line[9] = $corrected_read;
					$filtered_mismatch_read{$read_key} = join ("\t", @mismatch_line);
					my $MDtag = $1 if ($filtered_mismatch_read{$read_key} =~ /MD:Z:(\S+)/);
					my $start_pos = 0;
					my $pre_match_len = 0;
					my $pre_del = '';
					my $MDtag_cor = '';
					while ($MDtag =~ /(\d+)([A-Z]*)(\^[A-Z]+)*/g){
					    $start_pos += $1;
					    if ($pre_match_len > 0){
						if ($pre_del eq ''){
						    $MDtag_cor .= ($pre_match_len + $1)
						}
						if ($2 ne ''){
						    $start_pos ++;
						    $MDtag_cor .= $2 if (!defined $3);
						    $MDtag_cor .= $2 . $3 if (defined $3);
						}
						else{
						    $MDtag_cor .= $3 if (defined $3);
						}
						$pre_match_len = 0;
					    }
					    else{
						if ($2 ne ''){
						    $start_pos ++;
						    if ($start_pos == $read_pos){
							if (defined $3){
							    $pre_del = $3;
							    $MDtag_cor .= ($1 + 1) . $3;
							}
							else{
							    $pre_match_len = $1 + 1;
							}
						    }
						    else{
							$MDtag_cor .= $1 . $2 if (!defined $3);
							$MDtag_cor .= $1 . $2 . $3 if (defined $3);
						    }
						}
						else{
						    $MDtag_cor .= $1 if (!defined $3);
						    $MDtag_cor .= $1 . $3 if (defined $3);
						}
					    }
					}
					$filtered_mismatch_read{$read_key} =~ s/MD:Z:\S+/MD:Z:$MDtag_cor/;				
				    }
				}
			    }
			}
		    }
		}
		if ($read_correct == 1){
		    foreach my $key (keys %filtered_mismatch_read){
			my ($ch, $ps, $id) = split (/==/, $key);
			my @mismatch_line = split (/\s+/, $filtered_mismatch_read{$key});
			my $CIAGR = $mismatch_line[5];
			my $indel_num = 0;
			$indel_num = $CIAGR =~ s/[IDS]/A/g;
			my $MDtag = $1 if ($filtered_mismatch_read{$key} =~ /MD:Z:(\S+)/);
			my $alt_num = $MDtag =~ s/\d+[A-Z]//g;
			$alt_num = 0 if ($alt_num eq '');
			my $alt_indel_num = $alt_num + $indel_num;
			if ($alt_indel_num <= $limited_mismatch_num){
			    $filtered_mismatch_read{$key} =~ s/mismatch:\d+/mismatch:$alt_num/;
			    $match_read{$key} = $filtered_mismatch_read{$key};
			}
			else{
			    $filtered_pairs{$id} = 1 if ($filter_pairs == 1);
			}
		    }
		}
	    }
	    
	    foreach my $key (keys %match_read){
		my ($chr, $ps, $id, $dir) = split (/==/, $key);
		my @line = split (/\s+/, $match_read{$key});
		my $pos09d_2 = sprintf ("%09d", $line[3]);
		if ($filter_pairs == 1){
		    my @key_items = split (/==/, $key);
		    if (exists $filtered_pairs{$key_items[2]}){
			delete $match_read{$key};
			next;
		    }
		}
		if (($read_type eq 'PE') and ($limited_mismatch_num_pairs < $limited_mismatch_num * 2)){
		    $num_mismatch_pair{$id} = 0;
		    if ($match_read{$key} =~ /mismatch:(\d+)/){
			$num_mismatch_pair{$id} += $1;
		    }
		    if (($match_read{$key} =~ /\t\d+M\d+[IDN]\d+M\t/) or ($match_read{$key} =~ /\t\d+S\d+M\t/) or ($match_read{$key} =~ /\t\d+M\d+S\t/)){
			$num_mismatch_pair{$id} += 1;
		    }
		    if ($match_read{$key} =~ /\t\d+M\d+[IDN]\d+M\d+[ID]\d+M\t/){
			$num_mismatch_pair{$id} += 2;
		    }
		    if (exists $read_pair{$id}){
			if ($num_mismatch_pair{$id} > $limited_mismatch_num_pairs){
			    delete $match_read{$key};
			    next;
			}
		    }
		    else{
			$read_pair{$id} = 1;
		    }
		}
		if ($match_read{$key} =~ /mismatch:(\d+)/){
		    my $mismatch_num2 = $1;                
		    $n1++ if ($mismatch_num2 == 1);
		    $n2++ if ($mismatch_num2 == 2);
		    $n3++ if ($mismatch_num2 == 3);
		    $n4++ if ($mismatch_num2 == 4);
		    $n5++ if ($mismatch_num2 == 5);
		    $n6++ if ($mismatch_num2 == 6);
		    $n7++ if ($mismatch_num2 == 7);
		    $n8++ if ($mismatch_num2 == 8);
		    $n9++ if ($mismatch_num2 == 9);
		    $n10++ if ($mismatch_num2 == 10);
		    $n10o++ if ($mismatch_num2 > 10);
		}
		if ($match_read{$key} =~ /\trealigned:-*(\d+)/){
		    my $change_pos = $1;
		    if ($change_pos != 0){
			if ($read_type eq 'PE'){
			    $change_pos_read_F{$id} = $change_pos if (($line[1] == 99) or ($line[1] == 163) or ($line[1] == 97) or ($line[1] == 161) or ($line[1] == 65) or ($line[1] == 129));
			    $change_pos_read_R{$id} = $change_pos if (($line[1] == 83) or ($line[1] == 147) or ($line[1] == 81) or ($line[1] == 145) or ($line[1] == 113) or ($line[1] == 177));
			}
		    }
		}
		$match_read{$key} =~ s/\tmismatch:\d*//;
		$match_read{$key} =~ s/\tindel:\d*//;
		$match_read{$key} =~ s/\tsoftclip:\d*//;
		$match_read{$key} =~ s/\trealigned:-*\d*//;
		$match_read{$key} =~ s/[\t\s]*$//;
		if ($ps ne $pos09d_2){
		    my $new_key = $chr . '==' . $pos09d_2 . '==' . $id . '==' . $dir;
		    $match_read{$new_key} = $match_read{$key};
		    delete $match_read{$key};
		}
	    }
	    foreach my $key (sort keys %match_read){   
		$filtered_read_num ++;
		if ($read_type eq 'PE'){
		    my ($chrom, $posit, $readid) = split (/==/, $key);
		    my @line = split (/\s+/, $match_read{$key});
		    my $count = 0;
		    if ((exists $change_pos_read_F{$readid}) and (exists $change_pos_read_R{$readid})){
			@line = split (/\s+/, $match_read{$key});
			if (($line[1] == 99) or ($line[1] == 163) or ($line[1] == 97) or ($line[1] == 161) or ($line[1] == 65) or ($line[1] == 129)){
			    $line[8] = $line[8] - $change_pos_read_F{$readid} + $change_pos_read_R{$readid};
			    $line[7] += $change_pos_read_R{$readid};
			}
			elsif (($line[1] == 83) or ($line[1] == 147) or ($line[1] == 81) or ($line[1] == 145) or ($line[1] == 113) or ($line[1] == 177)){
			    $line[8] = $line[8] - $change_pos_read_R{$readid} + $change_pos_read_F{$readid};
			    $line[7] += $change_pos_read_F{$readid};
			}
			print join ("\t", @line), "\n";
			next;
		    }
		    elsif (exists $change_pos_read_F{$readid}){
			@line = split (/\s+/, $match_read{$key});
			if (($line[1] == 99) or ($line[1] == 163) or ($line[1] == 97) or ($line[1] == 161) or ($line[1] == 65) or ($line[1] == 129)){
			    $line[8] = $line[8] - $change_pos_read_F{$readid};
			}
			elsif (($line[1] == 83) or ($line[1] == 147) or ($line[1] == 81) or ($line[1] == 145) or ($line[1] == 113) or ($line[1] == 177)){
			    $line[8] = $line[8] + $change_pos_read_F{$readid};
			    $line[7] += $change_pos_read_F{$readid};
			}
			print join ("\t", @line), "\n";
			next;
		    }
		    elsif (exists $change_pos_read_R{$readid}){
			@line = split (/\s+/, $match_read{$key});
			if (($line[1] == 99) or ($line[1] == 163) or ($line[1] == 97) or ($line[1] == 161) or ($line[1] == 65) or ($line[1] == 129)){
			    $line[8] = $line[8] + $change_pos_read_R{$readid};
			    $line[7] += $change_pos_read_R{$readid};
			}
			elsif (($line[1] == 83) or ($line[1] == 147) or ($line[1] == 81) or ($line[1] == 145) or ($line[1] == 113) or ($line[1] == 177)){
			    $line[8] = $line[8] - $change_pos_read_R{$readid};
			}
			print join ("\t", @line), "\n";
			next;
		    }
		}
		print $match_read{$key}, "\n";
	    }
	    if (@unmapped > 0){
		foreach (@unmapped){
		    print $_, "\n";
		}
	    }
	    $count_indel += scalar (keys (%indel_read));
	    %indel_read = ();
	    %filtered_mismatch_read = ();
	    %match_read = ();
	    %mismatch_read = ();
	    %filtered_pairs = ();
	    %num_mismatch_pair = ();
	    %read_pair = ();
	    %change_pos_read_F = ();
	    %change_pos_read_R = ();
	    @unmapped = ();
        }
	$pre_chr = $chrom;
	my $chr_len = length $chr_seq{$chrom};
	if ($line[5] eq $read_length . 'M'){
	    my $ref_seq = substr ($chr_seq{$chrom}, $line[3] - 1, $read_length);
	    if (length ($ref_seq) < $read_length){
		next;
	    }
	    if ($line[9] eq $ref_seq){
		my $new_line = $line . "\tmismatch:0" if ($MDtag_sam == 1);
		$new_line = $line . "\tMD:Z:$read_length\tmismatch:0" if ($MDtag_sam == 0);
		$match_read{$chr_pos_id} = $new_line;
	    }
	    else{
		my ($total_mismatch, $lowq_mismatch, $terminal_mismatch_N, $terminal_mismatch_C, $MD) = &count_mismatch_T ($line[9], $ref_seq, $line[10], $chrom, $line[3]);
		my $new_line = $line . "\tmismatch:$total_mismatch" if ($MDtag_sam == 1);
		$new_line = $line . "\tMD:Z:$MD\tmismatch:$total_mismatch" if ($MDtag_sam == 0);
		if ($read_correct == 1){
		    if ($lowq_mismatch > $max_mismatch_num){
			$filtered_pairs{$read_id} = 1 if ($filter_pairs == 1);
		    }
		    else{
			if ($realignment eq 'ON'){
			    unless (exists $mismatch_read{$chr_pos}){
				$mismatch_read{$chr_pos} = $new_line;
			    }
			    elsif (exists $mismatch_read{$chr_pos}){
				$mismatch_read{$chr_pos} = $mismatch_read{$chr_pos} . '||' . $new_line;
			    }
			    $filtered_mismatch_read{$chr_pos_id} = $new_line;
			}
			elsif ($realignment eq 'OFF'){
			    $filtered_mismatch_read{$chr_pos_id} = $new_line;
			}
		    }
		}
		elsif ($read_correct == 0){
		    if ($realignment eq 'ON'){
			unless (exists $mismatch_read{$chr_pos}){
			    $mismatch_read{$chr_pos} = $new_line;
			}
			elsif (exists $mismatch_read{$chr_pos}){
			    $mismatch_read{$chr_pos} = $mismatch_read{$chr_pos} . '||' . $new_line;
			}
			if ((($terminal_mismatch_N > 0) and ($terminal_mismatch_C > 0)) or ($lowq_mismatch > $limited_mismatch_num)){
			    $filtered_pairs{$read_id} = 1 if ($filter_pairs == 1);
			}
			else{
			    $filtered_mismatch_read{$chr_pos_id} = $new_line;
			}
		    }
		    elsif ($realignment eq 'OFF'){
			if ((($terminal_mismatch_N > 0) and ($terminal_mismatch_C > 0)) or ($lowq_mismatch > $limited_mismatch_num)){
			    $filtered_pairs{$read_id} = 1 if ($filter_pairs == 1);
			}
			else{
			    $match_read{$chr_pos_id} = $new_line;
			}
		    }
		}
	    }
	}
	elsif ($line[5] =~ /^(\d+)M(\d+)I(\d+)M$/){
	    my $ref1 = substr ($chr_seq{$chrom}, $line[3] - 1, $1);
	    my $read1 = substr ($line[9], 0, $1);
	    my $qual1 = substr ($line[10], 0, $1);
	    my $ref12 = '';
	    my $read12 = '';
	    my $qual12 = '';
	    my $indel = '';
	    my $clip_tag = 0;
	    my $SC = '';
	    if ($chr_len - $line[3] - $1 + 1 - $3 >= 0){
		my $ref2 = substr ($chr_seq{$chrom}, $line[3] + $1 - 1, $3);
		$ref12 = $ref1 . $ref2;
		my $read2 = substr ($line[9], $1 + $2, $3);
		$read12 = $read1 . $read2;
		my $qual2 = substr ($line[10], $1 + $2, $3);
		$qual12 = $qual1 . $qual2;
		my $ins_pos = $pos + $1 - 1;
		my $ins_base = substr ($line[9], $1, $2);
		my $ins_info = $chrom . '==' . $ins_pos . '==' . $2 . '==I==' . $ins_base;
		my $M1_size = $1;
		my $M2_size = $3;
		my $indel = $2;
		$indel_read{$ins_info} = 1 if ($realignment eq 'ON');
	    }
	    else{
		$ref12 = $ref1;
		$read12 = $read1;
		$qual12 = $qual1;
		$SC = $read_length - $1;
		$line[5] = $1 . 'M' . $SC . 'S';
		$clip_tag = 1;
	    }
	    if ($read12 eq $ref12){
		my $new_line = $line . "\tmismatch:0\tindel:$indel" if ($MDtag_sam == 1) and ($clip_tag == 0);
		$new_line = join ("\t", @line) . "\tmismatch:0\tsoftclip:$SC" if ($MDtag_sam == 1) and ($clip_tag == 1);
		$new_line = $line . "\tMD:Z:" . ($1 + $3) . "\tmismatch:0\tindel:$indel" if ($MDtag_sam == 0) and ($clip_tag == 0);
		$new_line = join ("\t", @line) . "\tMD:Z:" . $1 . "\tmismatch:0\tsoftclip:$SC" if ($MDtag_sam == 0) and ($clip_tag == 1);
		$match_read{$chr_pos_id} = $new_line;
		$mismatch_read{$chr_pos} = $new_line if ($realignment eq 'ON');
	    }
	    else{
		my ($total_mismatch, $lowq_mismatch, $terminal_mismatch_N, $terminal_mismatch_C, $MD) = &count_mismatch_T ($read12, $ref12, $qual12, $chrom, $line[3]);
		my $new_line = $line . "\tmismatch:$total_mismatch\tindel:$indel" if ($MDtag_sam == 1) and ($clip_tag == 0);
		$new_line = join ("\t", @line) . "\tmismatch:$total_mismatch\tsoftclip:$SC" if ($MDtag_sam == 1) and ($clip_tag == 1);
		$new_line = $line . "\tMD:Z:$MD\tmismatch:$total_mismatch\tindel:$indel" if ($MDtag_sam == 0) and ($clip_tag == 0);
		$new_line = join ("\t", @line) . "\tMD:Z:$MD\tmismatch:$total_mismatch\tsoftclip:$SC" if ($MDtag_sam == 0) and ($clip_tag == 1);
		if ($read_correct == 1){
		    if ($lowq_mismatch > $max_mismatch_num){
			$filtered_pairs{$read_id} = 1 if ($filter_pairs == 1);
		    }
		    else{
			if ($realignment eq 'ON'){
			    unless (exists $mismatch_read{$chr_pos}){
				$mismatch_read{$chr_pos} = $new_line;
			    }
			    elsif (exists $mismatch_read{$chr_pos}){
				$mismatch_read{$chr_pos} = $mismatch_read{$chr_pos} . '||' . $new_line;
			    }
			    $filtered_mismatch_read{$chr_pos_id} = $new_line;
			}
			elsif ($realignment eq 'OFF'){
			    $filtered_mismatch_read{$chr_pos_id} = $new_line;
			}
		    }
		}
		elsif ($read_correct == 0){
		    if ($realignment eq 'ON'){
			unless (exists $mismatch_read{$chr_pos}){
			    $mismatch_read{$chr_pos} = $new_line;
			}
			elsif (exists $mismatch_read{$chr_pos}){
			    $mismatch_read{$chr_pos} = $mismatch_read{$chr_pos} . '||' . $new_line;
			}
			if ((($terminal_mismatch_N > 0) and ($terminal_mismatch_C > 0)) or ($lowq_mismatch > $limited_mismatch_num - 1)){
			    $filtered_pairs{$read_id} = 1 if ($filter_pairs == 1);
			}
			else{
			    $filtered_mismatch_read{$chr_pos_id} = $new_line;
			}
		    }
		    elsif ($realignment eq 'OFF'){
			if ((($terminal_mismatch_N > 0) and ($terminal_mismatch_C > 0)) or ($lowq_mismatch > $limited_mismatch_num - 1)){
			    $filtered_pairs{$read_id} = 1 if ($filter_pairs == 1);
			}
			else{
			    $match_read{$chr_pos_id} = $new_line;
			}
		    }
		}
	    }
	}
	elsif ($line[5] =~ /^(\d+)M(\d+)D(\d+)M$/){
	    my $ref1 = substr ($chr_seq{$chrom}, $line[3] - 1, $1);
	    my $read1 = substr ($line[9], 0, $1);
	    my $qual1 = substr ($line[10], 0, $1);
	    my $ref2 = '';
	    my $read2 = '';
	    my $qual2 = '';
	    my $ref12 = '';
	    my $read12 = '';
	    my $qual12 = '';
	    my $indel = '';
	    my $ref_del = '';
	    my $clip_tag = 0;
	    my $SC = '';
	    if ($chr_len - $line[3] - $1 - $2 + 1 - $3 >= 0){
		$ref2 = substr ($chr_seq{$chrom}, $line[3] + $1 + $2 - 1, $3);
		$ref12 = $ref1 . $ref2;
		$ref_del = substr ($chr_seq{$chrom}, $line[3] + $1 - 1, $2) if ($MDtag_sam == 0);
		$read2 = substr ($line[9], $1, $3);
		$read12 = $read1 . $read2;
		$qual2 = substr ($line[10], $1, $3);
		$qual12 = $qual1 . $qual2;
		my $del_pos = $pos + $1 - 1;
		my $del_info = $chrom . '==' . $del_pos . '==' . $2 . '==D==xxx';
		my $M1_size = $1;
		my $M2_size = $3;
		$indel = $2;
		$indel_read{$del_info} = 1 if ($realignment eq 'ON');
	    }
	    else{
		$ref12 = $ref1;
		$read12 = $read1;
		$qual12 = $qual1;
		$SC = $read_length - $1;
		$line[5] = $1 . 'M' . $SC . 'S';
		$clip_tag = 1;
	    }
	    if ($read12 eq $ref12){
		my $new_line = $line . "\tmismatch:0\tindel:$indel" if ($MDtag_sam == 1) and ($clip_tag == 0);
		$new_line = join ("\t", @line) . "\tmismatch:0\tsoftclip:$SC" if ($MDtag_sam == 1) and ($clip_tag == 1);
		$new_line = $line . "\tMD:Z:$1\^$ref_del$3\tmismatch:0\tindel:$indel" if ($MDtag_sam == 0) and ($clip_tag == 0);
		$new_line = join ("\t", @line) . "\tMD:Z:$1\^$ref_del$3\tmismatch:0\tsoftclip:$SC" if ($MDtag_sam == 0) and ($clip_tag == 1);
		$match_read{$chr_pos_id} = $new_line;
		$mismatch_read{$chr_pos} = $new_line if ($realignment eq 'ON');
	    }
	    else{
		my ($total_mismatch, $lowq_mismatch, $terminal_mismatch_N, $terminal_mismatch_C, $MD) = &count_mismatch_T ($read12, $ref12, $qual12) if ($MDtag_sam == 1);
		if ($MDtag_sam == 0){
		    my ($total_mismatch1, $lowq_mismatch1, $terminal_mismatch_N1, $terminal_mismatch_C1, $MD1) = &count_mismatch_T ($read1, $ref1, $qual1, $chrom, $line[3]);
		    my ($total_mismatch2, $lowq_mismatch2, $terminal_mismatch_N2, $terminal_mismatch_C2, $MD2) = (0, 0, 0, 0, '');
		    if ($clip_tag == 0){
			my $ref_start2 = $line[3] + $1 + $2;
			($total_mismatch2, $lowq_mismatch2, $terminal_mismatch_N2, $terminal_mismatch_C2, $MD2) = &count_mismatch_T ($read2, $ref2, $qual2, $chrom, $ref_start2);
		    }
		    $total_mismatch = $total_mismatch1 + $total_mismatch2;
		    $lowq_mismatch = $lowq_mismatch1 + $lowq_mismatch2;
		    $terminal_mismatch_N = $terminal_mismatch_N1 + $terminal_mismatch_N2;
		    $terminal_mismatch_C = $terminal_mismatch_C1 + $terminal_mismatch_C2;
		    $MD = $MD1 . '^' . $ref_del . $MD2 if ($clip_tag == 0);
		    $MD = $MD1 if ($clip_tag == 1);
		}
		my $new_line = $line . "\tmismatch:$total_mismatch\tindel:$indel" if ($MDtag_sam == 1) and ($clip_tag == 0);
		$new_line = join ("\t", @line) . "\tmismatch:$total_mismatch\tsoftclip:$SC" if ($MDtag_sam == 1) and ($clip_tag == 1);
		$new_line = $line . "\tMD:Z:$MD\tmismatch:$total_mismatch\tindel:$indel" if ($MDtag_sam == 0) and ($clip_tag == 0);
		$new_line = join ("\t", @line) . "\tMD:Z:$MD\tmismatch:$total_mismatch\tsoftclip:$SC" if ($MDtag_sam == 0) and ($clip_tag == 1);
		if ($read_correct == 1){
		    if ($lowq_mismatch > $max_mismatch_num){
			$filtered_pairs{$read_id} = 1 if ($filter_pairs == 1);
		    }
		    else{
			if ($realignment eq 'ON'){
			    unless (exists $mismatch_read{$chr_pos}){
				$mismatch_read{$chr_pos} = $new_line;
			    }
			    elsif (exists $mismatch_read{$chr_pos}){
				$mismatch_read{$chr_pos} = $mismatch_read{$chr_pos} . '||' . $new_line;
			    }
			    $filtered_mismatch_read{$chr_pos_id} = $new_line;
			}
			elsif ($realignment eq 'OFF'){
			    $filtered_mismatch_read{$chr_pos_id} = $new_line;
			}
		    }
		}
		elsif ($read_correct == 0){
		    if ($realignment eq 'ON'){
			unless (exists $mismatch_read{$chr_pos}){
			    $mismatch_read{$chr_pos} = $new_line;
			}
			elsif (exists $mismatch_read{$chr_pos}){
			    $mismatch_read{$chr_pos} = $mismatch_read{$chr_pos} . '||' . $new_line;
			}
			if ((($terminal_mismatch_N > 0) and ($terminal_mismatch_C > 0)) or ($lowq_mismatch > $limited_mismatch_num - 1)){
			    $filtered_pairs{$read_id} = 1 if ($filter_pairs == 1);
			}
			else{
			    $filtered_mismatch_read{$chr_pos_id} = $new_line;
			}
		    }
		    elsif ($realignment eq 'OFF'){
			if ((($terminal_mismatch_N > 0) and ($terminal_mismatch_C > 0)) or ($lowq_mismatch > $limited_mismatch_num - 1)){
			    $filtered_pairs{$read_id} = 1 if ($filter_pairs == 1);
			}
			else{
			    $match_read{$chr_pos_id} = $new_line;
			}
		    }
		}
	    }
	}
	elsif ($line[5] =~ /(\d+)M(\d+)([ID])(\d+)M(\d+)([ID])(\d+)M/){
	    my $ref1 = substr ($chr_seq{$chrom}, $line[3] - 1, $1);
	    my $read1 = substr ($line[9], 0, $1);
	    my $qual1 = substr ($line[10], 0, $1);
	    my $ref2 = '';
	    my $ref3 = '';
	    my $read2 = '';
	    my $read3 = '';
	    my $qual2 = '';
	    my $qual3 = '';
	    my $ref12 = '';
	    my $read12 = '';
	    my $qual12 = '';
	    my $indel = '';
	    my $ref_del = '';
	    my $clip_tag = 0;
	    my $SC = '';
	    if ((($3 eq 'I') and ($chr_len - $line[3] - $1 + 1 - $4 >= 0)) or (($3 eq 'D') and ($chr_len - $line[3] - $1 - $2 + 1 - $4 >= 0))){
		if ((($3 eq 'I') and ($6 eq 'I') and ($chr_len - $line[3] - $1 + 1 - $4 - $7 >= 0)) or (($3 eq 'I') and ($6 eq 'D') and ($chr_len - $line[3] - $1 + 1 - $4 - $5 - $7 >= 0)) or (($3 eq 'D') and ($6 eq 'I') and ($chr_len - $line[3] - $1 - $2 + 1 - $4 - $7 >= 0)) or (($3 eq 'D') and ($6 eq 'D') and ($chr_len - $line[3] - $1 - $2 + 1 - $4 - $5 - $7 >= 0))){
		    $ref2 = substr ($chr_seq{$chrom}, $line[3] + $1 - 1, $4) if ($3 eq 'I');
		    $ref2 = substr ($chr_seq{$chrom}, $line[3] + $1 + $2 - 1, $4) if ($3 eq 'D');	    
		    $ref3 = substr ($chr_seq{$chrom}, $line[3] + $1 + $4 - 1, $7) if (($3 eq 'I') and ($6 eq 'I'));
		    $ref3 = substr ($chr_seq{$chrom}, $line[3] + $1 + $4 + $5 - 1, $7) if (($3 eq 'I') and ($6 eq 'D'));
		    $ref3 = substr ($chr_seq{$chrom}, $line[3] + $1 + $2 + $4 - 1, $7) if (($3 eq 'D') and ($6 eq 'I'));
		    $ref3 = substr ($chr_seq{$chrom}, $line[3] + $1 + $2 + $4 + $5 - 1, $7) if (($3 eq 'D') and ($6 eq 'D'));
		    $ref12 = $ref1 . $ref2 . $ref3;
		    $read2 = substr ($line[9], $1 + $2, $4) if ($3 eq 'I');
		    $read2 = substr ($line[9], $1, $4) if ($3 eq 'D');
		    next if (length $ref12 != $1 + $4 + $7);
		    $read3 = substr ($line[9], $1 + $2 + $4 + $5, $7) if (($3 eq 'I') and ($6 eq 'I'));
		    $read3 = substr ($line[9], $1 + $2 + $4, $7) if (($3 eq 'I') and ($6 eq 'D'));
		    $read3 = substr ($line[9], $1 + $4 + $5, $7) if (($3 eq 'D') and ($6 eq 'I'));
		    $read3 = substr ($line[9], $1 + $4, $7) if (($3 eq 'D') and ($6 eq 'D'));
		    $read12 = $read1 . $read2 . $read3;
		    
		    $qual2 = substr ($line[10], $1 + $2, $4) if ($3 eq 'I');
		    $qual2 = substr ($line[10], $1, $4) if ($3 eq 'D');
		    $qual3 = substr ($line[10], $1 + $2 + $4 + $5, $7) if (($3 eq 'I') and ($6 eq 'I'));
		    $qual3 = substr ($line[10], $1 + $2 + $4, $7) if (($3 eq 'I') and ($6 eq 'D'));
		    $qual3 = substr ($line[10], $1 + $4 + $5, $7) if (($3 eq 'D') and ($6 eq 'I'));
		    $qual3 = substr ($line[10], $1 + $4, $7) if (($3 eq 'D') and ($6 eq 'D'));
		    $qual12 = $qual1 . $qual2 . $qual3;
		    $indel = $2 + $5;
		}
		else{
		    $ref12 = $ref1;
		    $read12 = $read1;
		    $qual12 = $qual1;
		    $SC = $read_length - $1;
		    $line[5] = $1 . 'M' . $SC . 'S';
		    $clip_tag = 1;
		}
	    }
	    else{
		$ref12 = $ref1;
		$read12 = $read1;
		$qual12 = $qual1;
		$SC = $read_length - $1;
		$line[5] = $1 . 'M' . $SC . 'S';
		$clip_tag = 1;
	    }
	    my $MDtag = '';	
	    if ($read12 eq $ref12){
		if ($MDtag_sam == 1){
		    my $new_line = $line . "\tmismatch:0\tindel:$indel" if ($clip_tag == 0);
		    $new_line = $line . "\tmismatch:0\tsoftclip:$SC" if ($clip_tag == 1);
		    $match_read{$chr_pos_id} = $new_line;
		}
		elsif ($MDtag_sam == 0){
		    if ($3 eq 'D'){
			my $ref_del1 = substr ($chr_seq{$chrom}, $line[3] + $1 - 1, $2);
			if ($6 eq 'D'){
			    my $ref_del2 = substr ($chr_seq{$chrom}, $line[3] + $1 + $2 + $4 - 1, $5);
			    $MDtag = $1 . '^' . $ref_del1 . $4 . '^' . $ref_del2 . $7;
			}
			else{
			    $MDtag = $1 . '^' . $ref_del1 . ($4 + $7);
			}
		    }
		    else{
			if ($6 eq 'D'){
			    my $ref_del2 = substr ($chr_seq{$chrom}, $line[3] + $1 + $4 - 1, $5);
			    $MDtag = ($1 + $4) . '^' . $ref_del2 . $7;
			}
			else{
			    $MDtag = $1 + $4 + $7;
			}
		    }
		    my $new_line = $line . "\tMD:Z:$MDtag\tmismatch:0\tindel:$indel" if ($clip_tag == 0);
		    $new_line = join ("\t", @line) . "\tMD:Z:$1\tmismatch:0\tsoftclip:$SC" if ($clip_tag == 1);
		    $match_read{$chr_pos_id} = $new_line;
		}
	    }
	    else{
		my ($total_mismatch, $lowq_mismatch, $terminal_mismatch_N, $terminal_mismatch_C, $MD) = &count_mismatch_T ($read12, $ref12, $qual12) if ($MDtag_sam == 1);
		if (($MDtag_sam == 0) and ($clip_tag == 0)){
		    if ($3 eq 'D'){
			my $ref_del1 = substr ($chr_seq{$chrom}, $line[3] + $1 - 1, $2);
			if ($6 eq 'D'){
			    my $ref_del2 = substr ($chr_seq{$chrom}, $line[3] + $1 + $2 + $4 - 1, $5);
			    my $ref_start2 = $line[3] + $1 + $2;
			    my $ref_start3 = $line[3] + $1 + $2 + $4 + $5;
			    my ($total_mismatch1, $lowq_mismatch1, $terminal_mismatch_N1, $terminal_mismatch_C1, $MD1) = &count_mismatch_T ($read1, $ref1, $qual1, $chrom, $line[3]);
			    my ($total_mismatch2, $lowq_mismatch2, $terminal_mismatch_N2, $terminal_mismatch_C2, $MD2) = &count_mismatch_T ($read2, $ref2, $qual2, $chrom, $ref_start2);
			    my ($total_mismatch3, $lowq_mismatch3, $terminal_mismatch_N3, $terminal_mismatch_C3, $MD3) = &count_mismatch_T ($read3, $ref3, $qual3, $chrom, $ref_start3);
			    $total_mismatch = $total_mismatch1 + $total_mismatch2 + $total_mismatch3;
			    $lowq_mismatch = $lowq_mismatch1 + $lowq_mismatch2 + $lowq_mismatch3;
			    $terminal_mismatch_N = $terminal_mismatch_N1 + $terminal_mismatch_N2 + $terminal_mismatch_N3;
			    $terminal_mismatch_C = $terminal_mismatch_C1 + $terminal_mismatch_C2 + $terminal_mismatch_C3;
			    $MD = $MD1 . '^' . $ref_del1 . $MD2 . '^' . $ref_del2 . $MD3;
			}
			else{
			    my $read2a = $read2 . $read3;
			    my $ref2a = $ref2 . $ref3;
			    my $qual2a = $qual2 . $qual3;
			    my $ref_start2 = $line[3] + $1 + $2;
			    my ($total_mismatch1, $lowq_mismatch1, $terminal_mismatch_N1, $terminal_mismatch_C1, $MD1) = &count_mismatch_T ($read1, $ref1, $qual1, $chrom, $line[3]);
			    my ($total_mismatch2, $lowq_mismatch2, $terminal_mismatch_N2, $terminal_mismatch_C2, $MD2) = &count_mismatch_T ($read2a, $ref2a, $qual2a, $chrom, $ref_start2);
			    $total_mismatch = $total_mismatch1 + $total_mismatch2;
			    $lowq_mismatch = $lowq_mismatch1 + $lowq_mismatch2;
			    $terminal_mismatch_N = $terminal_mismatch_N1 + $terminal_mismatch_N2;
			    $terminal_mismatch_C = $terminal_mismatch_C1 + $terminal_mismatch_C2;
			    $MD = $MD1 . '^' . $ref_del1 . $MD2;
			}
		    }
		    else{
			if ($6 eq 'D'){
			    my $ref_del2 = substr ($chr_seq{$chrom}, $line[3] + $1 + $4 - 1, $5);
			    my $read1a = $read1 . $read2;
			    my $ref1a = $ref1 . $ref2;
			    my $qual1a = $qual1 . $qual2;
			    my $ref_start3 = $line[3] + $1 + $2 + $4 + $5;
			    my ($total_mismatch1, $lowq_mismatch1, $terminal_mismatch_N1, $terminal_mismatch_C1, $MD1) = &count_mismatch_T ($read1a, $ref1a, $qual1a, $chrom, $line[3]);
			    my ($total_mismatch2, $lowq_mismatch2, $terminal_mismatch_N2, $terminal_mismatch_C2, $MD2) = &count_mismatch_T ($read3, $ref3, $qual3, $chrom, $ref_start3);
			    $total_mismatch = $total_mismatch1 + $total_mismatch2;
			    $lowq_mismatch = $lowq_mismatch1 + $lowq_mismatch2;
			    $terminal_mismatch_N = $terminal_mismatch_N1 + $terminal_mismatch_N2;
			    $terminal_mismatch_C = $terminal_mismatch_C1 + $terminal_mismatch_C2;
			    $MD = $MD1 . '^' . $ref_del2 . $MD2;
			}
			else{
			    ($total_mismatch, $lowq_mismatch, $terminal_mismatch_N, $terminal_mismatch_C, $MD) = &count_mismatch_T ($read12, $ref12, $qual12, $chrom, $line[3]);
			}
		    }
		}
		my $new_line = $line . "\tmismatch:$total_mismatch\tindel:$indel" if ($MDtag_sam == 1) and ($clip_tag == 0);
		$new_line = join ("\t", @line) . "\tmismatch:$total_mismatch\tsoftclip:$SC" if ($MDtag_sam == 1) and ($clip_tag == 1);
		$new_line = $line . "\tMD:Z:$MD\tmismatch:$total_mismatch\tindel:$indel" if ($MDtag_sam == 0) and ($clip_tag == 0);
		$new_line = join ("\t", @line) . "\tMD:Z:$MD\tmismatch:$total_mismatch\tsoftclip:$SC" if ($MDtag_sam == 0) and ($clip_tag == 1);
		if ($read_correct == 1){
		    if ($lowq_mismatch > $max_mismatch_num){
			$filtered_pairs{$read_id} = 1 if ($filter_pairs == 1);
		    }
		    else{
			$filtered_mismatch_read{$chr_pos_id} = $new_line;
		    }
		}
		elsif ($read_correct == 0){
		    if ($realignment eq 'ON'){
			if ((($terminal_mismatch_N > 0) and ($terminal_mismatch_C > 0)) or ($lowq_mismatch > $limited_mismatch_num - 2)){
			    $filtered_pairs{$read_id} = 1 if ($filter_pairs == 1);
			}
			else{
			    $filtered_mismatch_read{$chr_pos_id} = $new_line;
			}
		    }
		    elsif ($realignment eq 'OFF'){
			if ((($terminal_mismatch_N > 0) and ($terminal_mismatch_C > 0)) or ($lowq_mismatch > $limited_mismatch_num - 2)){
			    $filtered_pairs{$read_id} = 1 if ($filter_pairs == 1);
			}
			else{
			    $match_read{$chr_pos_id} = $new_line;
			}
		    }
		}
	    }
	}
	elsif ($line[5] =~ /^(\d+)S(\d+)M$/){
	    $pos -= $1;
	    next if ($pos <= 0);
	    my $pos09d = sprintf ("%09d", $pos);
	    my $read_id = $line[0];
	    my $chr_pos = $chrom . '==' . $pos09d;
	    my $chr_pos_id = $chrom . '==' . $pos09d . '==' . $read_id . '==' . $direct;
	    my $ref = substr ($chr_seq{$chrom}, $line[3] - 1, $2);
	    my $read = substr ($line[9], $1, $2);
	    my $qual = substr ($line[10], $1, $2);
	    my $SC = $1;
	    if ($read eq $ref){
		my $new_line = $line . "\tmismatch:0\tsoftclip:$SC" if ($MDtag_sam == 1);
		$new_line = $line . "\tMD:Z:$2\tmismatch:0\tsoftclip:$SC" if ($MDtag_sam == 0);
		if ($realignment eq 'ON'){
		    unless (exists $mismatch_read{$chr_pos}){
			$mismatch_read{$chr_pos} = $new_line;
		    }
		    elsif (exists $mismatch_read{$chr_pos}){
			$mismatch_read{$chr_pos} = $mismatch_read{$chr_pos} . '||' . $new_line;
		    }
		}
		elsif ($realignment eq 'OFF'){
		    $match_read{$chr_pos_id} = $new_line;
		}
	    }
	    else{
		my ($total_mismatch, $lowq_mismatch, $terminal_mismatch_N, $terminal_mismatch_C, $MD) = &count_mismatch_T ($read, $ref, $qual, $chrom, $line[3]);
		my $new_line = $line . "\tmismatch:$total_mismatch\tsoftclip:$SC" if ($MDtag_sam == 1);
		$new_line = $line . "\tMD:Z:$MD\tmismatch:$total_mismatch\tsoftclip:$SC" if ($MDtag_sam == 0);
		if ($read_correct == 1){
		    if ($lowq_mismatch > $max_mismatch_num){
			$filtered_pairs{$read_id} = 1 if ($filter_pairs == 1);
		    }
		    else{
			if ($realignment eq 'ON'){
			    unless (exists $mismatch_read{$chr_pos}){
				$mismatch_read{$chr_pos} = $new_line;
			    }
			    elsif (exists $mismatch_read{$chr_pos}){
				$mismatch_read{$chr_pos} = $mismatch_read{$chr_pos} . '||' . $new_line;
			    }
			    $filtered_mismatch_read{$chr_pos_id} = $new_line;
			}
			elsif ($realignment eq 'OFF'){
			    $filtered_mismatch_read{$chr_pos_id} = $new_line;
			}
		    }
		}
		elsif ($read_correct == 0){
		    if ($realignment eq 'ON'){
			unless (exists $mismatch_read{$chr_pos}){
			    $mismatch_read{$chr_pos} = $new_line;
			}
			elsif (exists $mismatch_read{$chr_pos}){
			    $mismatch_read{$chr_pos} = $mismatch_read{$chr_pos} . '||' . $new_line;
			}
			if ((($terminal_mismatch_N > 0) and ($terminal_mismatch_C > 0)) or ($lowq_mismatch > $limited_mismatch_num - 1)){
			    $filtered_pairs{$read_id} = 1 if ($filter_pairs == 1);
			}
			else{
			    $filtered_mismatch_read{$chr_pos_id} = $new_line;
			}
		    }
		    elsif ($realignment eq 'OFF'){
			if ((($terminal_mismatch_N > 0) and ($terminal_mismatch_C > 0)) or ($lowq_mismatch > $limited_mismatch_num - 1)){
			    $filtered_pairs{$read_id} = 1 if ($filter_pairs == 1);
			}
			else{
			    $match_read{$chr_pos_id} = $new_line;
			}
		    }
		}
	    }
	}
	elsif ($line[5] =~ /^(\d+)M(\d+)S$/){
	    my $ref = substr ($chr_seq{$chrom}, $line[3] - 1, $1);
	    my $read = substr ($line[9], 0, $1);
	    my $qual = substr ($line[10], 0, $1);
	    my $SC = $2;
	    if ($read eq $ref){
		my $new_line = $line . "\tmismatch:0\tsoftclip:$SC" if ($MDtag_sam == 1);
		$new_line = $line . "\tMD:Z:$1\tmismatch:0\tsoftclip:$SC" if ($MDtag_sam == 0);
		if ($realignment eq 'ON'){
		    unless (exists $mismatch_read{$chr_pos}){
			$mismatch_read{$chr_pos} = $new_line;
		    }
		    elsif (exists $mismatch_read{$chr_pos}){
			$mismatch_read{$chr_pos} = $mismatch_read{$chr_pos} . '||' . $new_line;
		    }
		}
		elsif ($realignment eq 'OFF'){
		    $match_read{$chr_pos_id} = $new_line;
		}
	    }
	    else{
		my ($total_mismatch, $lowq_mismatch, $terminal_mismatch_N, $terminal_mismatch_C, $MD) = &count_mismatch_T ($read, $ref, $qual, $chrom, $line[3]);
		my $new_line = $line . "\tmismatch:$total_mismatch\tsoftclip:$SC" if ($MDtag_sam == 1);
		$new_line = $line . "\tMD:Z:$MD\tmismatch:$total_mismatch\tsoftclip:$SC" if ($MDtag_sam == 0);
		if ($read_correct == 1){
		    if ($lowq_mismatch > $max_mismatch_num){
			$filtered_pairs{$read_id} = 1 if ($filter_pairs == 1);
		    }
		    else{
			if ($realignment eq 'ON'){
			    unless (exists $mismatch_read{$chr_pos}){
				$mismatch_read{$chr_pos} = $new_line;
			    }
			    elsif (exists $mismatch_read{$chr_pos}){
				$mismatch_read{$chr_pos} = $mismatch_read{$chr_pos} . '||' . $new_line;
			    }
			    $filtered_mismatch_read{$chr_pos_id} = $new_line;
			}
			elsif ($realignment eq 'OFF'){
			    $filtered_mismatch_read{$chr_pos_id} = $new_line;
			}
		    }
		}
		elsif ($read_correct == 0){
		    if ($realignment eq 'ON'){
			unless (exists $mismatch_read{$chr_pos}){
			    $mismatch_read{$chr_pos} = $new_line;
			}
			elsif (exists $mismatch_read{$chr_pos}){
			    $mismatch_read{$chr_pos} = $mismatch_read{$chr_pos} . '||' . $new_line;
			}
			if ((($terminal_mismatch_N > 0) and ($terminal_mismatch_C > 0)) or ($lowq_mismatch > $limited_mismatch_num - 1)){
			    $filtered_pairs{$read_id} = 1 if ($filter_pairs == 1);
			}
			else{
			    $filtered_mismatch_read{$chr_pos_id} = $new_line;
			}
		    }
		    elsif ($realignment eq 'OFF'){
			if ((($terminal_mismatch_N > 0) and ($terminal_mismatch_C > 0)) or ($lowq_mismatch > $limited_mismatch_num - 1)){
			    $filtered_pairs{$read_id} = 1 if ($filter_pairs == 1);
			}
			else{
			    $match_read{$chr_pos_id} = $new_line;
			}
		    }
		}
	    }
	}
	elsif ($line[5] =~ /N/){
	    my $ref12 = '';
	    my $read12 = $line[9];
	    my $qual12 = $line[10];
	    if ($line[5] =~ /^(\d+)M(\d+)N(\d+)M$/){
		my $ref1 = substr ($chr_seq{$chrom}, $line[3] - 1, $1);
		my $ref2 = substr ($chr_seq{$chrom}, $line[3] + $1 + $2 - 1, $3);
	    }
	    elsif ($line[5] =~ /^(\d+)M(\d+)N(\d+)M(\d+)N(\d+)M$/){
		my $ref1 = substr ($chr_seq{$chrom}, $line[3] - 1, $1);
		my $ref2 = substr ($chr_seq{$chrom}, $line[3] + $1 + $2 - 1, $3);
		my $ref3 = substr ($chr_seq{$chrom}, $line[3] + $1 + $2 + $3 + $4 - 1, $5);
		my $ref12 = $ref1 . $ref2 . $ref3;
	    }
	    if ($read12 eq $ref12){
		my $new_line = $line . "\tmismatch:0" if ($MDtag_sam == 1);
		$new_line = $line . "\tMD:Z:$read_length\tmismatch:0" if ($MDtag_sam == 0);
		$match_read{$chr_pos_id} = $new_line;
	    }
	    else{
		my ($total_mismatch, $lowq_mismatch, $terminal_mismatch_N, $terminal_mismatch_C, $MD) = &count_mismatch_T ($read12, $ref12, $qual12) if ($MDtag_sam == 1);
		my $new_line = $line . "\tmismatch:$total_mismatch" if ($MDtag_sam == 1);
		$new_line = $line . "\tMD:Z:$MD\tmismatch:$total_mismatch" if ($MDtag_sam == 0);
		if ((($terminal_mismatch_N > 0) and ($terminal_mismatch_C > 0)) or ($lowq_mismatch > $limited_mismatch_num - 1)){
		    $filtered_pairs{$read_id} = 1 if ($filter_pairs == 1);
		}
		else{
		    $match_read{$chr_pos_id} = $new_line;
		}
	    }
	}
    }
close (FILE);

my %mismatch_pos;
my %mismatch_pos_ms;
my %mismatch_pos_id;
my %mismatch_pos_id_ms;
my %read_depth;
my %read_depth_ms;
if ($realignment eq 'ON'){
    print STDERR "Realigning reads around indels in $pre_chr ... \n";
    &realign_indel ();
    if ($read_correct == 0){
	foreach my $key (keys %filtered_mismatch_read){
	    unless (exists $match_read{$key}){
		$match_read{$key} = $filtered_mismatch_read{$key};
	    }
	}
    }
}
if ($read_correct == 1){				# for read error correction and/or remove alignments with low mapping quality at heterogygous allele
    if ($read_correct == 1){
	print STDERR "Correcting potential errors in reads ... \n" if ($dis_realign == 0);
	print STDERR "Correcting potential errors in reads in $pre_chr ... \n" if ($dis_realign == 1);
    }
    if ($minimum_allele_freq > 0){
	foreach my $key (keys %match_read){			# record read number of mismatch-free reads at each chromosome position
	    my ($ch, $ps, $id) = split (/==/, $key);
	    $ps =~ s/^0*//;
	    my $MDtag = $1 if ($match_read{$key} =~ /MD:Z:(\S+)/);
	    my $RGtag = $1 if ($match_read{$key} =~ /RG:Z:(\S+)/) and ($multi_sample_mode == 1);
	    my $start_pos = 0;
	    while ($MDtag =~ /(\d+)([A-Z]*)(\^[A-Z]+)*/g){
		for (my $i = $ps + $start_pos; $i < $ps + $start_pos + $1; $i++){
		    $read_depth{$i} ++;
		    ${$read_depth_ms{$i}}{$RGtag} ++ if ($multi_sample_mode == 1);
		}
		$start_pos += $1;
		my $del = $3;	
		if ($2 ne ''){
		    my $mismatch_pos = $ps + $start_pos;
		    $read_depth{$mismatch_pos} ++;
		    ${$read_depth_ms{$mismatch_pos}}{$RGtag} ++ if ($multi_sample_mode == 1);
		    $start_pos ++;
		}
		if ((defined $del) and ($del =~ /\^([A-Z]+)/)){
		    $start_pos += length ($1);
		}
	    }
	}
    }
    foreach my $key (keys %filtered_mismatch_read){		# record reference position and base quality of mismatch base of mismatch-containing reads and the number of the same mismatch base covered at the site
	my ($ch, $ps, $id, $dir) = split (/==/, $key);		# record read number at each chromosome position
	my @mismatch_line = split (/\s+/, $filtered_mismatch_read{$key});
	my $pos = $mismatch_line[3];
	my $pos_09d = sprintf ("%09d", $pos);
	my $MDtag = $1 if ($filtered_mismatch_read{$key} =~ /MD:Z:(\S+)/);
	my $RGtag = $1 if ($filtered_mismatch_read{$key} =~ /RG:Z:(\S+)/) and ($multi_sample_mode == 1);
	my $CIAGR = $mismatch_line[5];
	my $start_pos = 0;
	my $read_start_pos = 0;
	while ($MDtag =~ /(\d+)([A-Z]*)(\^[A-Z]+)*/g){
	    if ($minimum_allele_freq > 0){
		for (my $i = $pos + $start_pos; $i < $pos + $start_pos + $1; $i++){
		    $read_depth{$i} ++;
		    ${$read_depth_ms{$i}}{$RGtag} ++ if ($multi_sample_mode == 1) and ($minimum_allele_freq > 0);
		}
	    }
	    $start_pos += $1;
	    $read_start_pos += $1;
	    my $del = $3;	
	    if ($2 ne ''){
		my $mismatch_pos = $pos + $start_pos;
		$start_pos ++;
		$read_start_pos ++;		    
		my ($mismatch_base, $mismatch_qual) = &correct_read ($mismatch_line[9], $CIAGR, $read_start_pos, $2, 'rbase', $mismatch_line[10]);		
		${$mismatch_pos{$mismatch_pos}}{$mismatch_base} .= $mismatch_qual;
		${${$mismatch_pos_ms{$mismatch_pos}}{$RGtag}}{$mismatch_base} .= $mismatch_qual if ($multi_sample_mode == 1);
		$read_depth{$mismatch_pos} ++ if ($minimum_allele_freq > 0);
		${$read_depth_ms{$mismatch_pos}}{$RGtag} ++ if ($multi_sample_mode == 1) and ($minimum_allele_freq > 0);
		if (!exists ${$mismatch_pos_id{$mismatch_pos}}{$mismatch_base}){
		    ${$mismatch_pos_id{$mismatch_pos}}{$mismatch_base} = $id . '==' . $pos_09d . '==' . $read_start_pos . '==' . $dir;
		}
		else{
		    ${$mismatch_pos_id{$mismatch_pos}}{$mismatch_base} .= '||' . $id . '==' . $pos_09d . '==' . $read_start_pos . '==' . $dir;
		}
		if ($multi_sample_mode == 1){
		    if (!exists ${${$mismatch_pos_id_ms{$mismatch_pos}}{$RGtag}}{$mismatch_base}){
			${${$mismatch_pos_id_ms{$mismatch_pos}}{$RGtag}}{$mismatch_base} = $id . '==' . $pos_09d . '==' . $read_start_pos . '==' . $dir;
		    }
		    else{
			${${$mismatch_pos_id_ms{$mismatch_pos}}{$RGtag}}{$mismatch_base} .= '||' . $id . '==' . $pos_09d . '==' . $read_start_pos . '==' . $dir;
		    }
		}
	    }
	    if ((defined $del) and ($del =~ /\^([A-Z]+)/)){
		$start_pos += length ($1);
	    }
	}
	if ($ps ne $pos_09d){
	    my $new_key = $ch . '==' . $pos_09d . '==' . $id . '==' . $dir;
	    $filtered_mismatch_read{$new_key} = $filtered_mismatch_read{$key};
	    delete $filtered_mismatch_read{$key};
	}
    }
    my %top_allele_pos;
    if (($multi_sample_mode == 0) or (($multi_sample_mode == 1) and ($multiple_allele_sample == 0))){
	foreach my $mismatch_pos (keys %mismatch_pos){
	    my $allele_1 = '';
	    my $allele_2 = '';
	    my $top_allele = '';
	    my $count = 0;
	    foreach my $allele (sort {length ${$mismatch_pos{$mismatch_pos}}{$b} <=> length ${$mismatch_pos{$mismatch_pos}}{$a}} keys %{$mismatch_pos{$mismatch_pos}}){
		$count++;
		$allele_1 = $allele if ($count == 1);
		$allele_2 = $allele if ($count == 2);
	    }
	    if (($allele_2 eq '') or (($allele_2 ne '') and (length ${$mismatch_pos{$mismatch_pos}}{$allele_1} > length ${$mismatch_pos{$mismatch_pos}}{$allele_2}))){
		$top_allele = $allele_1;
		my $sum_alt_qual_1 = 0;
		my $ave_alt_qual = 0;
		my @alt_qual_1 = split (//, ${$mismatch_pos{$mismatch_pos}}{$allele_1});
		foreach (@alt_qual_1){
		    my $base_qual = ord ($_) - 64 if ($base_call_type eq 'illumina');
		    $base_qual = ord ($_) - 33 if ($base_call_type eq 'sanger');
		    $sum_alt_qual_1 += $base_qual;
		}
		$ave_alt_qual = $sum_alt_qual_1 / @alt_qual_1 if (@alt_qual_1 > 0);
		if ($ave_alt_qual < $minimum_ave_alt_qual){
		    $top_allele = '';
		}
	    }
	    elsif (($allele_2 ne '') and (length ${$mismatch_pos{$mismatch_pos}}{$allele_1} == length ${$mismatch_pos{$mismatch_pos}}{$allele_2})){
		my $sum_alt_qual_1 = 0;
		my $sum_alt_qual_2 = 0;
		my $ave_alt_qual = 0;
		my @alt_qual_1 = split (//, ${$mismatch_pos{$mismatch_pos}}{$allele_1});
		my @alt_qual_2 = split (//, ${$mismatch_pos{$mismatch_pos}}{$allele_2});
		foreach (@alt_qual_1){
		    my $base_qual = ord ($_) - 64 if ($base_call_type eq 'illumina');
		    $base_qual = ord ($_) - 33 if ($base_call_type eq 'sanger');
		    $sum_alt_qual_1 += $base_qual;
		}
		foreach (@alt_qual_2){
		    my $base_qual = ord ($_) - 64 if ($base_call_type eq 'illumina');
		    $base_qual = ord ($_) - 33 if ($base_call_type eq 'sanger');
		    $sum_alt_qual_2 += $base_qual;
		}
		if ($sum_alt_qual_1 >= $sum_alt_qual_2){
		    $top_allele = $allele_1;
		    $ave_alt_qual = $sum_alt_qual_1 / @alt_qual_1 if (@alt_qual_1 > 0);
		}
		else{
		    $top_allele = $allele_2;
		    $ave_alt_qual = $sum_alt_qual_2 / @alt_qual_2 if (@alt_qual_2 > 0);
		}
		if ($ave_alt_qual < $minimum_ave_alt_qual){
		    $top_allele = '';
		}
	    }
	    $top_allele_pos{$mismatch_pos} = $top_allele;
	    if ($multi_sample_mode == 0){
		foreach my $allele (keys %{$mismatch_pos{$mismatch_pos}}){
		    my $allele_freq = length (${$mismatch_pos{$mismatch_pos}}{$allele}) / $read_depth{$mismatch_pos} if ($minimum_allele_freq > 0);
		    if (($allele ne $top_allele) or (($allele eq $top_allele) and ($minimum_allele_freq > 0) and ($allele_freq < $minimum_allele_freq))){
			my @id_pos = split (/\|\|/, ${$mismatch_pos_id{$mismatch_pos}}{$allele});
			foreach (@id_pos){
			    my ($id, $ref_pos, $read_pos, $dir) = split (/==/, $_);
			    my $read_key = $pre_chr . '==' . $ref_pos . '==' . $id . '==' . $dir;
			    next if (!exists $filtered_mismatch_read{$read_key});
			    my @mismatch_line = split (/\s+/, $filtered_mismatch_read{$read_key});
			    my $ref_base = substr ($chr_seq{$pre_chr}, $mismatch_pos - 1, 1);
			    my $corrected_read = &correct_read ($mismatch_line[9], $mismatch_line[5], $read_pos, $ref_base, 'correct');
			    $mismatch_line[9] = $corrected_read;
			    $filtered_mismatch_read{$read_key} = join ("\t", @mismatch_line);
			    my $MDtag = $1 if ($filtered_mismatch_read{$read_key} =~ /MD:Z:(\S+)/);
			    my $start_pos = 0;
			    my $pre_match_len = 0;
			    my $pre_del = '';
			    my $MDtag_cor = '';
			    while ($MDtag =~ /(\d+)([A-Z]*)(\^[A-Z]+)*/g){
				$start_pos += $1;
				if ($pre_match_len > 0){
				    if ($pre_del eq ''){
					$MDtag_cor .= ($pre_match_len + $1)
				    }
				    if ($2 ne ''){
					$start_pos ++;
					$MDtag_cor .= $2 if (!defined $3);
					$MDtag_cor .= $2 . $3 if (defined $3);
				    }
				    else{
					$MDtag_cor .= $3 if (defined $3);
				    }
				    $pre_match_len = 0;
				}
				else{
				    if ($2 ne ''){
					$start_pos ++;
					if ($start_pos == $read_pos){
					    if (defined $3){
						$pre_del = $3;
						$MDtag_cor .= ($1 + 1) . $3;
					    }
					    else{
						$pre_match_len = $1 + 1;
					    }
					}
					else{
					    $MDtag_cor .= $1 . $2 if (!defined $3);
					    $MDtag_cor .= $1 . $2 . $3 if (defined $3);
					}
				    }
				    else{
					$MDtag_cor .= $1 if (!defined $3);
					$MDtag_cor .= $1 . $3 if (defined $3);
				    }
				}
			    }
			    $filtered_mismatch_read{$read_key} =~ s/MD:Z:\S+/MD:Z:$MDtag_cor/;				
			}
		    }
		}
	    }
	}
    }
    if ($multi_sample_mode == 1){
	foreach my $mismatch_pos (keys %mismatch_pos_ms){
	    foreach my $RG (keys %{$mismatch_pos_ms{$mismatch_pos}}){
		my $allele_1 = '';
		my $allele_2 = '';
		my $top_allele = '';
		if ($multiple_allele_sample == 1){
		    my $count = 0;
		    foreach my $allele (sort {length ${${$mismatch_pos_ms{$mismatch_pos}}{$RG}}{$b} <=> length ${${$mismatch_pos_ms{$mismatch_pos}}{$RG}}{$a}} keys %{${$mismatch_pos_ms{$mismatch_pos}}{$RG}}){
			$count++;
			$allele_1 = $allele if ($count == 1);
			$allele_2 = $allele if ($count == 2);
		    }
		    if (($allele_2 eq '') or (($allele_2 ne '') and (length ${${$mismatch_pos_ms{$mismatch_pos}}{$RG}}{$allele_1} > length ${${$mismatch_pos_ms{$mismatch_pos}}{$RG}}{$allele_2}))){
			$top_allele = $allele_1;
			my $sum_alt_qual_1 = 0;
			my $ave_alt_qual = 0;
			my @alt_qual_1 = split (//, ${${$mismatch_pos_ms{$mismatch_pos}}{$RG}}{$allele_1});
			foreach (@alt_qual_1){
			    my $base_qual = ord ($_) - 64 if ($base_call_type eq 'illumina');
			    $base_qual = ord ($_) - 33 if ($base_call_type eq 'sanger');
			    $sum_alt_qual_1 += $base_qual;
			}
			$ave_alt_qual = $sum_alt_qual_1 / @alt_qual_1 if (@alt_qual_1 > 0);
			if ($ave_alt_qual < $minimum_ave_alt_qual){
			    $top_allele = '';
			}
		    }
		    elsif (($allele_2 ne '') and (length ${${$mismatch_pos_ms{$mismatch_pos}}{$RG}}{$allele_1} == length ${${$mismatch_pos_ms{$mismatch_pos}}{$RG}}{$allele_2})){
			my $sum_alt_qual_1 = 0;
			my $sum_alt_qual_2 = 0;
			my $ave_alt_qual = 0;
			my @alt_qual_1 = split (//, ${${$mismatch_pos_ms{$mismatch_pos}}{$RG}}{$allele_1});
			my @alt_qual_2 = split (//, ${${$mismatch_pos_ms{$mismatch_pos}}{$RG}}{$allele_2});
			foreach (@alt_qual_1){
			    my $base_qual = ord ($_) - 64 if ($base_call_type eq 'illumina');
			    $base_qual = ord ($_) - 33 if ($base_call_type eq 'sanger');
			    $sum_alt_qual_1 += $base_qual;
			}
			foreach (@alt_qual_2){
			    my $base_qual = ord ($_) - 64 if ($base_call_type eq 'illumina');
			    $base_qual = ord ($_) - 33 if ($base_call_type eq 'sanger');
			    $sum_alt_qual_2 += $base_qual;
			}
			if ($sum_alt_qual_1 >= $sum_alt_qual_2){
			    $top_allele = $allele_1;
			    $ave_alt_qual = $sum_alt_qual_1 / @alt_qual_1 if (@alt_qual_1 > 0);
			}
			else{
			    $top_allele = $allele_2;
			    $ave_alt_qual = $sum_alt_qual_2 / @alt_qual_2 if (@alt_qual_2 > 0);
			}
			if ($ave_alt_qual < $minimum_ave_alt_qual){
			    $top_allele = '';
			}
		    }
		}
		foreach my $allele (keys %{${$mismatch_pos_ms{$mismatch_pos}}{$RG}}){
		    $top_allele = $top_allele_pos{$mismatch_pos} if ($multiple_allele_sample == 0);
		    my $allele_freq = length (${${$mismatch_pos_ms{$mismatch_pos}}{$RG}}{$allele}) / ${$read_depth_ms{$mismatch_pos}}{$RG} if ($minimum_allele_freq > 0);
		    if (($allele ne $top_allele) or (($allele eq $top_allele) and ($minimum_allele_freq > 0) and ($allele_freq < $minimum_allele_freq))){
			my @id_pos = split (/\|\|/, ${${$mismatch_pos_id_ms{$mismatch_pos}}{$RG}}{$allele});
			foreach (@id_pos){
			    my ($id, $ref_pos, $read_pos, $dir) = split (/==/, $_);
			    my $read_key = $pre_chr . '==' . $ref_pos . '==' . $id . '==' . $dir;
			    next if (!exists $filtered_mismatch_read{$read_key});
			    my @mismatch_line = split (/\s+/, $filtered_mismatch_read{$read_key});
			    my $ref_base = substr ($chr_seq{$pre_chr}, $mismatch_pos - 1, 1);
			    my $corrected_read = &correct_read ($mismatch_line[9], $mismatch_line[5], $read_pos, $ref_base, 'correct');
			    $mismatch_line[9] = $corrected_read;
			    $filtered_mismatch_read{$read_key} = join ("\t", @mismatch_line);
			    my $MDtag = $1 if ($filtered_mismatch_read{$read_key} =~ /MD:Z:(\S+)/);
			    my $start_pos = 0;
			    my $pre_match_len = 0;
			    my $pre_del = '';
			    my $MDtag_cor = '';
			    while ($MDtag =~ /(\d+)([A-Z]*)(\^[A-Z]+)*/g){
				$start_pos += $1;
				if ($pre_match_len > 0){
				    if ($pre_del eq ''){
					$MDtag_cor .= ($pre_match_len + $1)
				    }
				    if ($2 ne ''){
					$start_pos ++;
					$MDtag_cor .= $2 if (!defined $3);
					$MDtag_cor .= $2 . $3 if (defined $3);
				    }
				    else{
					$MDtag_cor .= $3 if (defined $3);
				    }
				    $pre_match_len = 0;
				}
				else{
				    if ($2 ne ''){
					$start_pos ++;
					if ($start_pos == $read_pos){
					    if (defined $3){
						$pre_del = $3;
						$MDtag_cor .= ($1 + 1) . $3;
					    }
					    else{
						$pre_match_len = $1 + 1;
					    }
					}
					else{
					    $MDtag_cor .= $1 . $2 if (!defined $3);
					    $MDtag_cor .= $1 . $2 . $3 if (defined $3);
					}
				    }
				    else{
					$MDtag_cor .= $1 if (!defined $3);
					$MDtag_cor .= $1 . $3 if (defined $3);
				    }
				}
			    }
			    $filtered_mismatch_read{$read_key} =~ s/MD:Z:\S+/MD:Z:$MDtag_cor/;				
			}
		    }
		}
	    }
	}
    }
    if ($read_correct == 1){
	foreach my $key (keys %filtered_mismatch_read){
	    my ($ch, $ps, $id) = split (/==/, $key);
	    my @mismatch_line = split (/\s+/, $filtered_mismatch_read{$key});
	    my $CIAGR = $mismatch_line[5];
	    my $indel_num = 0;
	    $indel_num = $CIAGR =~ s/[IDS]/A/g;
	    my $MDtag = $1 if ($filtered_mismatch_read{$key} =~ /MD:Z:(\S+)/);
	    my $alt_num = $MDtag =~ s/\d+[A-Z]//g;
	    $alt_num = 0 if ($alt_num eq '');
	    my $alt_indel_num = $alt_num + $indel_num;
	    if ($alt_indel_num <= $limited_mismatch_num){
		$filtered_mismatch_read{$key} =~ s/mismatch:\d+/mismatch:$alt_num/;
		$match_read{$key} = $filtered_mismatch_read{$key};
	    }
	    else{
		$filtered_pairs{$id} = 1 if ($filter_pairs == 1);
	    }
	}
    }
}

foreach my $key (keys %match_read){
    my ($chr, $ps, $id, $dir) = split (/==/, $key);
    my @line = split (/\s+/, $match_read{$key});
    my $pos09d = sprintf ("%09d", $line[3]);
    if ($filter_pairs == 1){
	my @key_items = split (/==/, $key);
	if (exists $filtered_pairs{$key_items[2]}){
	    delete $match_read{$key};
	    next;
	}
    }
    if ($line[5] =~ /^(\d+)M$/){
    }
    elsif ($line[5] =~ /^(\d+)M(\d+)I(\d+)M$/){
    }
    elsif ($line[5] =~ /^(\d+)M(\d+)D(\d+)M$/){
    }
    elsif ($line[5] =~ /^(\d+)S(\d+)M$/){
    }
    elsif ($line[5] =~ /^(\d+)M(\d+)S$/){
    }
    
    
    if (($read_type eq 'PE') and ($limited_mismatch_num_pairs < $limited_mismatch_num * 2)){
	$num_mismatch_pair{$id} = 0;
	if ($match_read{$key} =~ /mismatch:(\d+)/){
	    $num_mismatch_pair{$id} += $1;
	}
	if (($match_read{$key} =~ /\t\d+M\d+[IDN]\d+M\t/) or ($match_read{$key} =~ /\t\d+S\d+M\t/) or ($match_read{$key} =~ /\t\d+M\d+S\t/)){
	    $num_mismatch_pair{$id} += 1;
	}
	if ($match_read{$key} =~ /\t\d+M\d+[IDN]\d+M\d+[ID]\d+M\t/){
	    $num_mismatch_pair{$id} += 2;
	}
	if (exists $read_pair{$id}){
	    if ($num_mismatch_pair{$id} > $limited_mismatch_num_pairs){
		delete $match_read{$key};
		next;
	    }
	}
	else{
	    $read_pair{$id} = 1;
	}
    }
    if ($match_read{$key} =~ /mismatch:(\d+)/){
	my $mismatch_num2 = $1;                
	$n1++ if ($mismatch_num2 == 1);
	$n2++ if ($mismatch_num2 == 2);
	$n3++ if ($mismatch_num2 == 3);
	$n4++ if ($mismatch_num2 == 4);
	$n5++ if ($mismatch_num2 == 5);
	$n6++ if ($mismatch_num2 == 6);
	$n7++ if ($mismatch_num2 == 7);
	$n8++ if ($mismatch_num2 == 8);
	$n9++ if ($mismatch_num2 == 9);
	$n10++ if ($mismatch_num2 == 10);
	$n10o++ if ($mismatch_num2 > 10);
    }
    if ($match_read{$key} =~ /\trealigned:-*(\d+)/){
	my $change_pos = $1;
	if ($change_pos != 0){
	    if ($read_type eq 'PE'){
		$change_pos_read_F{$id} = $change_pos if (($line[1] == 99) or ($line[1] == 163) or ($line[1] == 97) or ($line[1] == 161) or ($line[1] == 65) or ($line[1] == 129));
		$change_pos_read_R{$id} = $change_pos if (($line[1] == 83) or ($line[1] == 147) or ($line[1] == 81) or ($line[1] == 145) or ($line[1] == 113) or ($line[1] == 177));
	    }
	}
    }
    $match_read{$key} =~ s/\tmismatch:\d*//;
    $match_read{$key} =~ s/\tindel:\d*//;
    $match_read{$key} =~ s/\tsoftclip:\d*//;
    $match_read{$key} =~ s/\trealigned:-*\d*//;
    $match_read{$key} =~ s/[\t\s]*$//;
    if ($ps ne $pos09d){
	my $new_key = $chr . '==' . $pos09d . '==' . $id . '==' . $dir;
	$match_read{$new_key} = $match_read{$key};
	delete $match_read{$key};
    }
}
foreach my $key (sort keys %match_read){
    $filtered_read_num ++;
    if ($read_type eq 'PE'){
	my ($chrom, $posit, $readid) = split (/==/, $key);
	my @line = ();
	if ((exists $change_pos_read_F{$readid}) and (exists $change_pos_read_R{$readid})){
	    @line = split (/\s+/, $match_read{$key});
	    if (($line[1] == 99) or ($line[1] == 163) or ($line[1] == 97) or ($line[1] == 161) or ($line[1] == 65) or ($line[1] == 129)){
		$line[8] = $line[8] - $change_pos_read_F{$readid} + $change_pos_read_R{$readid};
		$line[7] += $change_pos_read_R{$readid};
	    }
	    elsif (($line[1] == 83) or ($line[1] == 147) or ($line[1] == 81) or ($line[1] == 145) or ($line[1] == 113) or ($line[1] == 177)){
		$line[8] = $line[8] - $change_pos_read_R{$readid} + $change_pos_read_F{$readid};
		$line[7] += $change_pos_read_F{$readid};
	    }
	    print join ("\t", @line), "\n";
	    next;
	}
	elsif (exists $change_pos_read_F{$readid}){
	    @line = split (/\s+/, $match_read{$key});
	    if (($line[1] == 99) or ($line[1] == 163) or ($line[1] == 97) or ($line[1] == 161) or ($line[1] == 65) or ($line[1] == 129)){
		$line[8] = $line[8] - $change_pos_read_F{$readid};
	    }
	    elsif (($line[1] == 83) or ($line[1] == 147) or ($line[1] == 81) or ($line[1] == 145) or ($line[1] == 113) or ($line[1] == 177)){
		$line[8] = $line[8] + $change_pos_read_F{$readid};
		$line[7] += $change_pos_read_F{$readid};
	    }
	    print join ("\t", @line), "\n";
	    next;
	}
	elsif (exists $change_pos_read_R{$readid}){
	    @line = split (/\s+/, $match_read{$key});
	    if (($line[1] == 99) or ($line[1] == 163) or ($line[1] == 97) or ($line[1] == 161) or ($line[1] == 65) or ($line[1] == 129)){
		$line[8] = $line[8] + $change_pos_read_R{$readid};
		$line[7] += $change_pos_read_R{$readid};
	    }
	    elsif (($line[1] == 83) or ($line[1] == 147) or ($line[1] == 81) or ($line[1] == 145) or ($line[1] == 113) or ($line[1] == 177)){
		$line[8] = $line[8] - $change_pos_read_R{$readid};
	    }
	    print join ("\t", @line), "\n";
	    next;
	}
    }
    print $match_read{$key}, "\n";
}
$count_indel += scalar (keys (%indel_read));
if (@unmapped > 0){
    foreach (@unmapped){
	print $_, "\n";
    }
}



print STDERR "\n<< Number of mismatch reads before filtering >>\n";
print STDERR "1-mismatch reads = ", $m1, "\n";
print STDERR "2-mismatch reads = ", $m2, "\n";
print STDERR "3-mismatch reads = ", $m3, "\n";
print STDERR "4-mismatch reads = ", $m4, "\n";
print STDERR "5-mismatch reads = ", $m5, "\n";
print STDERR "6-mismatch reads = ", $m6, "\n";
print STDERR "7-mismatch reads = ", $m7, "\n";
print STDERR "8-mismatch reads = ", $m8, "\n";
print STDERR "9-mismatch reads = ", $m9, "\n";
print STDERR "10-mismatch reads = ", $m10, "\n";
print STDERR ">10-mismatch reads = ", $m10o, "\n";

#print STDERR "terminal mismatch reads = $terminal_mismatch\n";
print STDERR "discordant paired-end reads = $discordant_read_num\n";
print STDERR "total aligned read number = $read_number\n\n";


print STDERR "<< Number of mismatch reads after filtering >>\n";
print STDERR "1-mismatch reads = ", $n1, "\n";
print STDERR "2-mismatch reads = ", $n2, "\n";
print STDERR "3-mismatch reads = ", $n3, "\n";
print STDERR "4-mismatch reads = ", $n4, "\n";
print STDERR "5-mismatch reads = ", $n5, "\n";
print STDERR "6-mismatch reads = ", $n6, "\n";
print STDERR "7-mismatch reads = ", $n7, "\n";
print STDERR "8-mismatch reads = ", $n8, "\n";
print STDERR "9-mismatch reads = ", $n9, "\n";
print STDERR "10-mismatch reads = ", $n10, "\n";
print STDERR ">10-mismatch reads = ", $n10o, "\n";

print STDERR "total indels = $count_indel\n" if ($realignment eq 'ON');
print STDERR "realigned reads = $count_realigned_read\n" if ($realignment eq 'ON');
print STDERR "read number after filtering = $filtered_read_num\n";


###############################################################################

sub realign_indel{      
    foreach my $key (sort keys %indel_read){
	my ($chr, $indelpos, $indelsize, $indel, $insbase) = split (/==/, $key);
	$count_realigned_indel ++;
	my $del_size = $indelsize if ($indel eq 'D');
	$del_size = 0 if ($indel eq 'I');
	for (my $posi = $indelpos - $ave_read_length + 2; $posi <= $indelpos + $del_size; $posi++){
	    my $pos09d = sprintf ("%09d", $posi);
	    my $chr_pos = $chr . '==' . $pos09d;
	    if (exists $mismatch_read{$chr_pos}){
		my @mismatch_read_lines = split (/\|\|/, $mismatch_read{$chr_pos});
		my $count_mismatch_line = 0;
		my @corrected_read_pos;
		foreach my $mismatch_line (@mismatch_read_lines){
		    $count_mismatch_line ++;
		    my $correct_mismatch_flag = 0;
		    my @line = split (/\s+/, $mismatch_line);
		    $read_length = length $line[9];
		    next if ($indelpos - $posi + 1 > $read_length - 1);
		    $mismatch_line =~ /mismatch:(\d+)/;
		    my $mismatch_num = $1;
		    $mismatch_line =~ /softclip:(\d+)/;
		    my $SC_num = $1;
		    my $read_id = $line[0];
		    my $flag = $line[1];
		    my $dir = 1;
		    $dir = 2 if ($flag == 83) or ($flag == 83) or ($flag == 147) or ($flag == 81) or ($flag == 145) or ($flag == 113) or ($flag == 177);
		    my $chr_pos_id = $chr_pos . '==' . $read_id . '==' . $dir;
		    my $ref12 = '';
		    my $read12 = '';
		    my $read1 = '';
		    my $read2 = '';
		    my $read3 = '';
		    my $ref1 = '';
		    my $ref2 = '';
		    my $ref3 = '';
		    my $mismatch_count = 0;		
		    if ($line[5] eq $read_length . 'M'){
			if ($indel eq 'I'){
			    if (($indelpos - $posi + 1) < $read_length / 2){
				if (($indelpos - $posi + 1) > $indelsize){
				    my $ins_read = substr ($line[9], $indelpos - $posi + 1 - $indelsize, $indelsize);
				    next if ($ins_read ne $insbase);
				    $ref12 = substr ($chr_seq{$chr}, $posi - 1 + $indelsize, $read_length - $indelsize);
				    $read1 = substr ($line[9], 0, $indelpos - $posi + 1 - $indelsize);
				    next if ($indelpos - $posi + 1 > $read_length - 1);
				    $read2 = substr ($line[9], $indelpos - $posi + 1);
				    $read12 = $read1 . $read2;
				    $line[3] += $indelsize;
				    next if (length $read12 != length $ref12);
				    if ($read12 eq $ref12){
					$line[5] = ($indelpos - $posi + 1 - $indelsize) . 'M' . $indelsize . 'I' . ($read_length - $indelpos + $posi - 1) . 'M';
					my $new_line = join ("\t", @line);
					$new_line =~ s/mismatch:\d+/mismatch:0\trealigned:$indelsize/;
					my $MDtag = $read_length - $indelsize;
					$new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
					$match_read{$chr_pos_id} = $new_line;
					$count_realigned_read ++;
					$correct_mismatch_flag = 1;
					delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
					delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
				    }
				    else{
					$mismatch_count = &count_mismatch_I ($read12, $ref12);
					if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
					    $line[5] = ($indelpos - $posi + 1 - $indelsize) . 'M' . $indelsize . 'I' . ($read_length - $indelpos + $posi - 1) . 'M';
					    my $new_line = join ("\t", @line);
					    $new_line =~ s/mismatch:\d+/mismatch:$mismatch_count\trealigned:$indelsize/;
					    my $MDtag = &get_mismatch_info_I ($read12, $ref12);
					    $new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
					    $filtered_mismatch_read{$chr_pos_id} = $new_line;
					    $count_realigned_read ++;
					    $correct_mismatch_flag = 1;
					    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					}
				    }
				}
				else{
				    my $ins_read = substr ($line[9], 0, $indelpos - $posi + 1);
				    next if ($ins_read ne substr ($insbase, length ($insbase) - length ($ins_read)));
				    $ref12 = substr ($chr_seq{$chr}, $indelpos, $read_length - $indelpos + $posi - 1);
				    $read12 = substr ($line[9], $indelpos - $posi + 1);
				    my $indelsize_term = $indelpos - $posi + 1;
				    $line[3] = $indelpos + 1;
				    next if (length $read12 != length $ref12);
				    if ($read12 eq $ref12){
					$line[5] = ($indelpos - $posi + 1) . 'I' . ($read_length - $indelpos + $posi - 1) . 'M';
					my $new_line = join ("\t", @line);
					$new_line =~ s/mismatch:\d+/mismatch:0\trealigned:$indelsize_term/;
					my $MDtag = $read_length - $indelsize;
					$new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
					$match_read{$chr_pos_id} = $new_line;
					$count_realigned_read ++;
					$correct_mismatch_flag = 1;
					delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
					delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
				    }
				    else{
					$mismatch_count = &count_mismatch_I ($read12, $ref12);
					if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
					    $line[5] = ($indelpos - $posi + 1) . 'I' . ($read_length - $indelpos + $posi - 1) . 'M';
					    my $new_line = join ("\t", @line);
					    $new_line =~ s/mismatch:\d+/mismatch:$mismatch_count\trealigned:$indelsize_term/;
					    my $MDtag = &get_mismatch_info_I ($read12, $ref12);
					    $new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
					    $filtered_mismatch_read{$chr_pos_id} = $new_line;
					    $count_realigned_read ++;
					    $correct_mismatch_flag = 1;
					    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					}
				    }
				}
			    }
			    else{
				if ($read_length - $indelpos + $posi - 1 > $indelsize){
				    my $ins_read = substr ($line[9], $indelpos - $posi + 1, $indelsize);
				    next if ($ins_read ne $insbase);
				    $ref12 = substr ($chr_seq{$chr}, $posi - 1, $read_length - $indelsize);
				    $read1 = substr ($line[9], 0, $indelpos - $posi + 1);
				    next if ($indelpos - $posi + 1 + $indelsize > $read_length - 1);
				    $read2 = substr ($line[9], $indelpos - $posi + 1 + $indelsize);
				    $read12 = $read1 . $read2;
				    next if (length $read12 != length $ref12);
				    if ($read12 eq $ref12){
					$line[5] = ($indelpos - $posi + 1) . 'M' . $indelsize . 'I' . ($read_length - $indelpos + $posi - 1 - $indelsize) . 'M';
					my $new_line = join ("\t", @line);
					$new_line =~ s/mismatch:\d+/mismatch:0\trealigned:0/;
					my $MDtag = $read_length - $indelsize;
					$new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
					$match_read{$chr_pos_id} = $new_line;
					$count_realigned_read ++;
					$correct_mismatch_flag = 1;
					delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
					delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
				    }
				    else{
					$mismatch_count = &count_mismatch_I ($read12, $ref12);
					if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
					    $line[5] = ($indelpos - $posi + 1) . 'M' . $indelsize . 'I' . ($read_length - $indelpos + $posi - 1 - $indelsize) . 'M';
					    my $new_line = join ("\t", @line);
					    $new_line =~ s/mismatch:\d+/mismatch:$mismatch_count\trealigned:0/;
					    my $MDtag = &get_mismatch_info_I ($read12, $ref12);
					    $new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
					    $filtered_mismatch_read{$chr_pos_id} = $new_line;
					    $count_realigned_read ++;
					    $correct_mismatch_flag = 1;
					    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					}
				    }
				}
				else{
				    my $ins_read = substr ($line[9], $indelpos - $posi + 1, $read_length - $indelpos + $posi - 1);
				    next if ($ins_read ne substr ($insbase, 0, length ($ins_read)));
				    $ref12 = substr ($chr_seq{$chr}, $posi - 1, $indelpos - $posi + 1);
				    $read12 = substr ($line[9], 0, $indelpos - $posi + 1);
				    next if (length $read12 != length $ref12);
				    if ($read12 eq $ref12){
					$line[5] = ($indelpos - $posi + 1) . 'M' . ($read_length - $indelpos + $posi - 1) . 'I';
					my $new_line = join ("\t", @line);
					$new_line =~ s/mismatch:\d+/mismatch:0\trealigned:0/;
					my $MDtag = $read_length - $indelsize;
					$new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
					$match_read{$chr_pos_id} = $new_line;
					$count_realigned_read ++;
					$correct_mismatch_flag = 1;
					delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
					delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
				    }
				    else{
					$mismatch_count = &count_mismatch_I ($read12, $ref12);
					if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
					    $line[5] = ($indelpos - $posi + 1) . 'M' . ($read_length - $indelpos + $posi - 1) . 'I';
					    my $new_line = join ("\t", @line);
					    $new_line =~ s/mismatch:\d+/mismatch:$mismatch_count\trealigned:0/;
					    my $MDtag = &get_mismatch_info_I ($read12, $ref12);
					    $new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
					    $filtered_mismatch_read{$chr_pos_id} = $new_line;
					    $count_realigned_read ++;
					    $correct_mismatch_flag = 1;
					    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					}
				    }
				}
			    }
			}
			elsif ($indel eq 'D'){
			    if (($indelpos <= $posi) or ($indelpos - $posi + 1 <= $read_length - $indelsize - $indelpos + $posi - 1)){
				$ref1 = substr ($chr_seq{$chr}, $posi - 1 - $indelsize, $indelpos - $posi + 1 + $indelsize);
				$ref2 = substr ($chr_seq{$chr}, $indelpos + $indelsize, $read_length - length $ref1);
				$ref12 = $ref1 . $ref2;
				$read12 = $line[9];
				$line[3] -= $indelsize;
				next if (length $read12 != length $ref12);
				if ($read12 eq $ref12){
				    $line[5] = ($indelpos - $posi + 1 + $indelsize) . 'M' . $indelsize . 'D' . ($read_length - $indelpos + $posi - 1 - $indelsize) . 'M';
				    my $new_line = join ("\t", @line);
				    $new_line =~ s/mismatch:\d+/mismatch:0\trealigned:-$indelsize/;
				    my $delseq = substr ($chr_seq{$chr}, $indelpos, $indelsize);
				    my $MDtag = length $read1 . '^' . $delseq . length $read2;
				    $new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
				    $match_read{$chr_pos_id} = $new_line;
				    $count_realigned_read ++;
				    $correct_mismatch_flag = 1;
				    delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
				    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
				}
				else{
				    $mismatch_count = &count_mismatch_I ($read12, $ref12);
				    if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
					$line[5] = ($indelpos - $posi + 1 + $indelsize) . 'M' . $indelsize . 'D' . ($read_length - $indelpos + $posi - 1 - $indelsize) . 'M';
					my $new_line = join ("\t", @line);
					$new_line =~ s/mismatch:\d+/mismatch:$mismatch_count\trealigned:-$indelsize/;
					my $delseq = substr ($chr_seq{$chr}, $indelpos, $indelsize);		
					$read1 = substr ($read12, 0, $indelpos - $posi + 1 + $indelsize);
					$read2 = substr ($read12, $indelpos - $posi + 1 + $indelsize);
					next if (length $read1 != length $ref1) or (length $read2 != length $ref2);
					my $MDtag = &get_mismatch_info_D ($read1, $ref1, $read2, $ref2, $delseq);
					$new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
					$filtered_mismatch_read{$chr_pos_id} = $new_line;
					$count_realigned_read ++;
					$correct_mismatch_flag = 1;
					delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
				    }
				}
			    }
			    else{
				$ref1 = substr ($chr_seq{$chr}, $posi - 1, $indelpos - $posi + 1);
				$ref2 = substr ($chr_seq{$chr}, $indelpos + $indelsize, $read_length - length $ref1);
				$ref12 = $ref1 . $ref2;
				$read12 = $line[9];
				next if (length $read12 != length $ref12);
				if ($read12 eq $ref12){
				    $line[5] = ($indelpos - $posi + 1) . 'M' . $indelsize . 'D' . ($read_length - $indelpos + $posi - 1) . 'M';
				    my $new_line = join ("\t", @line);
				    $new_line =~ s/mismatch:\d+/mismatch:0\trealigned:0/;
				    my $delseq = substr ($chr_seq{$chr}, $indelpos, $indelsize);
				    my $MDtag = length $read1 . '^' . $delseq . length $read2;
				    $new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
				    $match_read{$chr_pos_id} = $new_line;
				    $count_realigned_read ++;
				    $correct_mismatch_flag = 1;
				    delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
				    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
				}
				else{
				    $mismatch_count = &count_mismatch_I ($read12, $ref12);
				    if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
					$line[5] = ($indelpos - $posi + 1) . 'M' . $indelsize . 'D' . ($read_length - $indelpos + $posi - 1) . 'M';
					my $new_line = join ("\t", @line);
					$new_line =~ s/mismatch:\d+/mismatch:$mismatch_count\trealigned:0/;
					my $delseq = substr ($chr_seq{$chr}, $indelpos, $indelsize);					
					$read1 = substr ($read12, 0, $indelpos - $posi + 1);
					$read2 = substr ($read12, $indelpos - $posi + 1);
					next if (length $read1 != length $ref1) or (length $read2 != length $ref2);
					my $MDtag = &get_mismatch_info_D ($read1, $ref1, $read2, $ref2, $delseq);
					$new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
					$filtered_mismatch_read{$chr_pos_id} = $new_line;
					$count_realigned_read ++;
					$correct_mismatch_flag = 1;
					delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
				    }
				}
			    }
			}
		    }
		    elsif ($line[5] =~ /^(\d+)M(\d+)I(\d+)M$/){
			next if (($posi + $1 - 1 == $indelpos) and ($indel eq 'I') and ($indelsize == $2));
			if ($indel eq 'I'){
			    if (($indelpos - $posi + 1) < $read_length / 2){				
				if ($posi + $1 - 1 > $indelpos){
				    if (($indelpos - $posi + 1) > $indelsize){
					my $ins_read = substr ($line[9], $indelpos - $posi + 1 - $indelsize, $indelsize);
					next if ($ins_read ne $insbase);
					$ref12 = substr ($chr_seq{$chr}, $posi - 1 + $indelsize, $read_length - $indelsize - $2);
					$read1 = substr ($line[9], 0, $indelpos - $posi + 1 - $indelsize);
					next if ($indelpos - $posi + 1 + $indelsize > $read_length - 1);
					$read2 = substr ($line[9], $indelpos - $posi + 1, $1 - $indelsize - length $read1);
					$read3 = substr ($line[9], $1 + $2, $3);
					$read12 = $read1 . $read2 . $read3;
					$line[3] += $indelsize;
					next if (length $read12 != length $ref12);
					if ($read12 eq $ref12){
					    $line[5] = ($indelpos - $posi + 1 - $indelsize) . 'M' . $indelsize . 'I' . ($1 - $indelpos + $posi - 1) . 'M' . $2 . 'I' . $3 . 'M';
					    my $new_line = join ("\t", @line);
					    $new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:0\trealigned:$indelsize/;
					    $match_read{$chr_pos_id} = $new_line;
					    $count_realigned_read ++;
					    $correct_mismatch_flag = 1;
					    delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
					    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					}
					else{
					    $mismatch_count = &count_mismatch_I ($read12, $ref12);
					    if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
						$line[5] = ($indelpos - $posi + 1 - $indelsize) . 'M' . $indelsize . 'I' . ($1 - $indelpos + $posi - 1) . 'M' . $2 . 'I' . $3 . 'M';
						my $new_line = join ("\t", @line);
						$new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:$mismatch_count\trealigned:$indelsize/;
						my $MDtag = &get_mismatch_info_I ($read12, $ref12);
						$new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
						$filtered_mismatch_read{$chr_pos_id} = $new_line;
						$count_realigned_read ++;
						$correct_mismatch_flag = 1;
						delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					    }
					}
				    }
				    else{
					my $ins_read = substr ($line[9], 0, $indelpos - $posi + 1);
					next if ($ins_read ne substr ($insbase, length ($insbase) - length ($ins_read)));
					$read1 = substr ($line[9], $indelpos - $posi + 1, $1 - $indelpos + $posi - 1);
					next if ($indelpos - $posi + 1 + $indelsize > $read_length - 1);
					$read2 = substr ($line[9], $1 + $2, $3);
					$read12 = $read1 . $read2;
					$ref12 = substr ($chr_seq{$chr}, $indelpos, length $read12);
					my $indelsize_term = $indelpos - $posi + 1;
					$line[3] = $indelpos + 1;
					next if (length $read12 != length $ref12);
					if ($read12 eq $ref12){
					    $line[5] = ($indelpos - $posi + 1) . 'I' . ($1 - $indelpos + $posi - 1) . 'M' . $2 . 'I' . $3 . 'M';
					    my $new_line = join ("\t", @line);
					    $new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:0\trealigned:$indelsize_term/;
					    $match_read{$chr_pos_id} = $new_line;
					    $count_realigned_read ++;
					    $correct_mismatch_flag = 1;
					    delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
					    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					}
					else{
					    $mismatch_count = &count_mismatch_I ($read12, $ref12);
					    if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
						$line[5] = ($indelpos - $posi + 1) . 'I' . ($1 - $indelpos + $posi - 1) . 'M' . $2 . 'I' . $3 . 'M';
						my $new_line = join ("\t", @line);
						$new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:$mismatch_count\trealigned:$indelsize_term/;
						my $MDtag = &get_mismatch_info_I ($read12, $ref12);
						$new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
						$filtered_mismatch_read{$chr_pos_id} = $new_line;
						$count_realigned_read ++;
						$correct_mismatch_flag = 1;
						delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					    }
					}
				    }
				}
				else{				# insersion in read is neglected and the read is treaed as a normal read
				    if (($indelpos - $posi + 1 + $2) > $indelsize){
					my $ins_read = substr ($line[9], $indelpos - $posi + 1 - $indelsize + $2, $indelsize);
					next if ($ins_read ne $insbase);
					$ref12 = substr ($chr_seq{$chr}, $posi - 1 - $indelsize - $2, $read_length - $indelsize);
					$read1 = substr ($line[9], 0, $indelpos - $posi + 1 - $indelsize + $2);
					next if ($indelpos - $posi + 1 + $indelsize > $read_length - 1);
					$read2 = substr ($line[9], $indelpos - $posi + 1 + $2);
					$read12 = $read1 . $read2;
					my $indelsize_term = $indelsize - $2;
					$line[3] += $indelsize_term;
					next if (length $read12 != length $ref12);
					if ($read12 eq $ref12){
					    $line[5] = ($indelpos - $posi + 1 - $indelsize + $2) . 'M' . $indelsize . 'I' . ($3 - $indelpos + $posi - 1 + $1) . 'M';
					    my $new_line = join ("\t", @line);
					    $new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:0\trealigned:$indelsize_term/;
					    $match_read{$chr_pos_id} = $new_line;
					    $count_realigned_read ++;
					    $correct_mismatch_flag = 1;
					    delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
					    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					}
					else{
					    $mismatch_count = &count_mismatch_I ($read12, $ref12);
					    if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
						$line[5] = ($indelpos - $posi + 1 - $indelsize + $2) . 'M' . $indelsize . 'I' . ($3 - $indelpos + $posi - 1 + $1) . 'M';
						my $new_line = join ("\t", @line);
						$new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:$mismatch_count\trealigned:$indelsize_term/;
						my $MDtag = &get_mismatch_info_I ($read12, $ref12);
						$new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
						$filtered_mismatch_read{$chr_pos_id} = $new_line;
						$count_realigned_read ++;
						$correct_mismatch_flag = 1;
						delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					    }
					}
				    }
				    else{
					my $ins_read = substr ($line[9], 0, $indelpos - $posi + 1 + $2);
					next if ($ins_read ne substr ($insbase, length ($insbase) - length ($ins_read)));
					$read12 = substr ($line[9], $indelpos - $posi + 1 + $2);
					$ref12 = substr ($chr_seq{$chr}, $indelpos, length $read12);
					my $indelsize_term = $indelpos - $posi + 1 + $2;
					$line[3] = $indelpos + 1;
					next if (length $read12 != length $ref12);
					if ($read12 eq $ref12){
					    $line[5] = ($indelpos - $posi + 1 + $2) . 'I' . ($3 - $indelpos + $posi - 1 + $1) . 'M';
					    my $new_line = join ("\t", @line);
					    $new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:0\trealigned:$indelsize_term/;
					    $match_read{$chr_pos_id} = $new_line;
					    $count_realigned_read ++;
					    $correct_mismatch_flag = 1;
					    delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
					    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					}
					else{
					    $mismatch_count = &count_mismatch_I ($read12, $ref12);
					    if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
						$line[5] = ($indelpos - $posi + 1 + $2) . 'I' . ($3 - $indelpos + $posi - 1 + $1) . 'M';
						my $new_line = join ("\t", @line);
						$new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:$mismatch_count\trealigned:$indelsize_term/;
						my $MDtag = &get_mismatch_info_I ($read12, $ref12);
						$new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
						$filtered_mismatch_read{$chr_pos_id} = $new_line;
						$count_realigned_read ++;
						$correct_mismatch_flag = 1;
						delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					    }
					}
				    }
				}
			    }
			    else{
				if ($posi + $1 - 1 >= $indelpos){	# insersion in read is neglected and the read is treaed as a normal read
				    if ($read_length - $indelpos + $posi - 1 > $indelsize){
					my $ins_read = substr ($line[9], $indelpos - $posi + 1, $indelsize);
					next if ($ins_read ne $insbase);
					$ref12 = substr ($chr_seq{$chr}, $posi - 1, $read_length - $indelsize);
					$read1 = substr ($line[9], 0, $indelpos - $posi + 1);
					next if ($indelpos - $posi + 1 + $indelsize > $read_length - 1);
					$read2 = substr ($line[9], $indelpos - $posi + 1 + $indelsize);
					$read12 = $read1 . $read2;
					next if (length $read12 != length $ref12);
					if ($read12 eq $ref12){
					    $line[5] = ($indelpos - $posi + 1) . 'M' . $indelsize . 'I' . ($1 + $2 + $3 - $indelsize - $indelpos + $posi - 1) . 'M';
					    my $new_line = join ("\t", @line);
					    $new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:0\trealigned:0/;
					    $match_read{$chr_pos_id} = $new_line;
					    $count_realigned_read ++;
					    $correct_mismatch_flag = 1;
					    delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
					    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					}
					else{
					    $mismatch_count = &count_mismatch_I ($read12, $ref12);
					    if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
						$line[5] = ($indelpos - $posi + 1) . 'M' . $indelsize . 'I' . ($1 + $2 + $3 - $indelsize - $indelpos + $posi - 1) . 'M';
						my $new_line = join ("\t", @line);
						$new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:$mismatch_count\trealigned:0/;
						my $MDtag = &get_mismatch_info_I ($read12, $ref12);
						$new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
						$filtered_mismatch_read{$chr_pos_id} = $new_line;
						$count_realigned_read ++;
						$correct_mismatch_flag = 1;
						delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					    }
					}
				    }
				    else{
					my $ins_read = substr ($line[9], $indelpos - $posi + 1, $read_length - $indelpos + $posi - 1);
					next if ($ins_read ne substr ($insbase, 0, length ($ins_read)));
					$read12 = substr ($line[9], 0, $indelpos - $posi + 1);
					$ref12 = substr ($chr_seq{$chr}, $posi - 1, length $read12);
					my $indelsize_term = $read_length - $indelpos + $posi - 1;
					next if (length $read12 != length $ref12);
					if ($read12 eq $ref12){
					    $line[5] = ($indelpos - $posi + 1) . 'M' . ($read_length - $indelpos + $posi - 1) . 'I';
					    my $new_line = join ("\t", @line);
					    $new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:0\trealigned:0/;
					    $match_read{$chr_pos_id} = $new_line;
					    $count_realigned_read ++;
					    $correct_mismatch_flag = 1;
					    delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
					    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					}
					else{
					    $mismatch_count = &count_mismatch_I ($read12, $ref12);
					    if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
						$line[5] = ($indelpos - $posi + 1) . 'M' . ($read_length - $indelpos + $posi - 1) . 'I';
						my $new_line = join ("\t", @line);
						$new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:$mismatch_count\trealigned:0/;
						my $MDtag = &get_mismatch_info_I ($read12, $ref12);
						$new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
						$filtered_mismatch_read{$chr_pos_id} = $new_line;
						$count_realigned_read ++;
						$correct_mismatch_flag = 1;
						delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					    }
					}
				    }
				}
				else{
				    if ($read_length - $indelpos + $posi - 1 - $2 > $indelsize){
					my $ins_read = substr ($line[9], $2 + $indelpos - $posi + 1, $indelsize);
					next if ($ins_read ne $insbase);
					$ref12 = substr ($chr_seq{$chr}, $posi - 1, $read_length - $indelsize - $2);
					$read1 = substr ($line[9], 0, $1);
					$read2 = substr ($line[9], $1 + $2, $indelpos - $posi + 1 - $1);
					next if (length ($read1) + length ($read2) + $2 + $indelsize > $read_length - 1);
					$read3 = substr ($line[9], length ($read1) + length ($read2) + $2 + $indelsize);
					$read12 = $read1 . $read2 . $read3;
					next if (length $read12 != length $ref12);
					if ($read12 eq $ref12){
					    $line[5] = $1 . 'M' . $2 . 'I' . ($indelpos - $posi + 1 - $1) . 'M' . $indelsize . 'I' . ($1 + $3 - $indelpos + $posi - 1 - $indelsize) . 'M';
					    my $new_line = join ("\t", @line);
					    $new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:0\trealigned:0/;
					    $match_read{$chr_pos_id} = $new_line;
					    $count_realigned_read ++;
					    $correct_mismatch_flag = 1;
					    delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
					    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					}
					else{
					    $mismatch_count = &count_mismatch_I ($read12, $ref12);
					    if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
						$line[5] = $1 . 'M' . $2 . 'I' . ($indelpos - $posi + 1 - $1) . 'M' . $indelsize . 'I' . ($1 + $3 - $indelpos + $posi - 1 - $indelsize) . 'M';
						my $new_line = join ("\t", @line);
						$new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:$mismatch_count\trealigned:0/;
						my $MDtag = &get_mismatch_info_I ($read12, $ref12);
						$new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
						$filtered_mismatch_read{$chr_pos_id} = $new_line;
						$count_realigned_read ++;
						$correct_mismatch_flag = 1;
						delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					    }
					}
				    }
				    elsif ($read_length - $indelpos + $posi - 1 - $2 > 0){
					my $ins_read = substr ($line[9], $2 + $indelpos - $posi + 1, $read_length - $indelpos + $posi - 1 - $2);
					next if ($ins_read ne substr ($insbase, 0, length ($ins_read)));
					$read1 = substr ($line[9], 0, $1);
					$read2 = substr ($line[9], $1 + $2, $indelpos - $posi + 1 - $1);
					$read12 = $read1 . $read2;
					$ref12 = substr ($chr_seq{$chr}, $posi - 1, length $read12);
					next if (length $read12 != length $ref12);
					if ($read12 eq $ref12){
					    $line[5] = $1 . 'M' . $2 . 'I' . ($indelpos - $posi + 1 - $1) . 'M' . ($read_length - $indelpos + $posi - 1 - $2) . 'I';
					    my $new_line = join ("\t", @line);
					    $new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:0\trealigned:0/;
					    $match_read{$chr_pos_id} = $new_line;
					    $count_realigned_read ++;
					    $correct_mismatch_flag = 1;
					    delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
					    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					}
					else{
					    $mismatch_count = &count_mismatch_I ($read12, $ref12);
					    if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
						$line[5] = $1 . 'M' . $2 . 'I' . ($indelpos - $posi + 1 - $1) . 'M' . ($read_length - $indelpos + $posi - 1 - $2) . 'I';
						my $new_line = join ("\t", @line);
						$new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:$mismatch_count\trealigned:0/;
						my $MDtag = &get_mismatch_info_I ($read12, $ref12);
						$new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
						$filtered_mismatch_read{$chr_pos_id} = $new_line;
						$count_realigned_read ++;
						$correct_mismatch_flag = 1;
						delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					    }
					}
				    }
				}
			    }
			}
			elsif ($indel eq 'D'){
			    if (($indelpos <= $posi) or ($indelpos - $posi + 1 <= $read_length - $indelsize - $indelpos + $posi - 1)){
				if ($posi + $1 - 1 > $indelpos + $indelsize){
				    $ref1 = substr ($chr_seq{$chr}, $posi - 1 - $indelsize, $indelpos - $posi + 1 + $indelsize);
				    $ref2 = substr ($chr_seq{$chr}, $indelpos + $indelsize, $read_length - length ($ref1) - $2);
				    $ref12 = $ref1 . $ref2;
				    $read1 = substr ($line[9], 0, $1);
				    $read2 = substr ($line[9], $1 + $2, $3);
				    $read12 = $read1 . $read2;
				    $line[3] -= $indelsize;
				    next if (length $read12 != length $ref12);
				    if ($read12 eq $ref12){
					$line[5] = ($indelpos - $posi + 1 + $indelsize) . 'M' . $indelsize . 'D' . ($1 - $indelsize - $indelpos + $posi - 1) . 'M' . $2 . 'I' . $3 . 'M';
					my $new_line = join ("\t", @line);
					$new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:0\trealigned:-$indelsize/;
					$match_read{$chr_pos_id} = $new_line;
					$count_realigned_read ++;
					$correct_mismatch_flag = 1;
					delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
					delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
				    }
				    else{
					$mismatch_count = &count_mismatch_I ($read12, $ref12);
					if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
					    $line[5] = ($indelpos - $posi + 1 + $indelsize) . 'M' . $indelsize . 'D' . ($1 - $indelsize - $indelpos + $posi - 1) . 'M' . $2 . 'I' . $3 . 'M';
					    my $new_line = join ("\t", @line);
					    $new_line =~ s/mismatch:\d+/mismatch:$mismatch_count\trealigned:-$indelsize/;
					    my $delseq = substr ($chr_seq{$chr}, $indelpos, $indelsize);
					    $read1 = substr ($read12, 0, $indelpos - $posi + 1 + $indelsize);
					    $read2 = substr ($read12, $indelpos - $posi + 1 + $indelsize);
					    next if (length $read1 != length $ref1) or (length $read2 != length $ref2);
					    my $MDtag = &get_mismatch_info_D ($read1, $ref1, $read2, $ref2, $delseq);
					    $new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
					    $filtered_mismatch_read{$chr_pos_id} = $new_line;
					    $count_realigned_read ++;
					    $correct_mismatch_flag = 1;
					    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					}
				    }
				}
				else{				# insersion in read is neglected and the read is treaed as a normal read
				    $ref1 = substr ($chr_seq{$chr}, $posi - 1 - $indelsize - $2, $indelpos - $posi + 1 + $indelsize + $2);
				    $ref2 = substr ($chr_seq{$chr}, $indelpos + $indelsize, $read_length - length $ref1);
				    $ref12 = $ref1 . $ref2;
				    $read12 = $line[9];
				    my $shift_size = $indelsize + $2;
				    $line[3] -= $shift_size;
				    next if (length $read12 != length $ref12);
				    if ($read12 eq $ref12){
					$line[5] = ($indelpos - $posi + 1 + $2 + $indelsize) . 'M' . $indelsize . 'D' . ($3 - $indelsize - $indelpos + $posi - 1 + $1) . 'M';
					my $new_line = join ("\t", @line);
					$new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:0\trealigned:-$shift_size/;
					$match_read{$chr_pos_id} = $new_line;
					$count_realigned_read ++;
					$correct_mismatch_flag = 1;
					delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
					delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
				    }
				    else{
					$mismatch_count = &count_mismatch_I ($read12, $ref12);
					if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
					    $line[5] = ($indelpos - $posi + 1 + $2 + $indelsize) . 'M' . $indelsize . 'D' . ($3 - $indelsize - $indelpos + $posi - 1 + $1) . 'M';
					    my $delseq = substr ($chr_seq{$chr}, $indelpos, $indelsize);
					    $read1 = substr ($read12, 0, $indelpos - $posi + 1 + $indelsize + $2);
					    $read2 = substr ($read12, $indelpos - $posi + 1 + $indelsize + $2);
					    next if (length $read1 != length $ref1) or (length $read2 != length $ref2);
					    my $MDtag = &get_mismatch_info_D ($read1, $ref1, $read2, $ref2, $delseq);
					    my $new_line = join ("\t", @line);
					    $new_line =~ s/mismatch:\d+/mismatch:$mismatch_count\trealigned:-$shift_size/;
					    $new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
					    $filtered_mismatch_read{$chr_pos_id} = $new_line;
					    $count_realigned_read ++;
					    $correct_mismatch_flag = 1;
					    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					}
				    }
				}
			    }
			    else{
				if ($posi + $1 - 1 >= $indelpos){	# insersion in read is neglected and the read is treaed as a normal read
				    $ref1 = substr ($chr_seq{$chr}, $posi - 1, $indelpos - $posi + 1);
				    $ref2 = substr ($chr_seq{$chr}, $indelpos + $indelsize, $read_length - length $ref1);
				    $ref12 = $ref1 . $ref2;
				    $read12 = $line[9];
				    next if (length $read12 != length $ref12);
				    if ($read12 eq $ref12){
					$line[5] = ($indelpos - $posi + 1) . 'M' . $indelsize . 'D' . ($1 + $2 + $3 - $indelpos + $posi - 1) . 'M';
					my $new_line = join ("\t", @line);
					$new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:0\trealigned:0/;
					$match_read{$chr_pos_id} = $new_line;
					$count_realigned_read ++;
					$correct_mismatch_flag = 1;
					delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
					delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
				    }
				    else{
					$mismatch_count = &count_mismatch_I ($read12, $ref12);
					if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
					    $line[5] = ($indelpos - $posi + 1) . 'M' . $indelsize . 'D' . ($1 + $2 + $3 - $indelpos + $posi - 1) . 'M';
					    my $new_line = join ("\t", @line);
					    $new_line =~ s/mismatch:\d+/mismatch:$mismatch_count\trealigned:0/;
					    my $delseq = substr ($chr_seq{$chr}, $indelpos, $indelsize);
					    $read1 = substr ($read12, 0, $indelpos - $posi + 1);
					    $read2 = substr ($read12, $indelpos - $posi + 1);
					    next if (length $read1 != length $ref1) or (length $read2 != length $ref2);
					    my $MDtag = &get_mismatch_info_D ($read1, $ref1, $read2, $ref2, $delseq);
					    $new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
					    $filtered_mismatch_read{$chr_pos_id} = $new_line;
					    $count_realigned_read ++;
					    $correct_mismatch_flag = 1;
					    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					}
				    }
				}
				elsif ($indelpos - $posi + 1 + $2 < $read_length){
				    $ref1 = substr ($chr_seq{$chr}, $posi - 1, $indelpos - $posi + 1);
				    $ref2 = substr ($chr_seq{$chr}, $indelpos + $indelsize, $read_length - $2 - length $ref1);
				    $ref12 = $ref1 . $ref2;
				    $read1 = substr ($line[9], 0, $1);
				    $read2 = substr ($line[9], $1 + $2, $3);
				    $read12 = $read1 . $read2;
				    next if (length $read12 != length $ref12);
				    if ($read12 eq $ref12){
					$line[5] = $1 . 'M' . $2 . 'I' . ($indelpos - $posi + 1 - $1) . 'M' . $indelsize . 'D' . ($read_length - $indelpos + $posi - 1 - $2) . 'M';
					my $new_line = join ("\t", @line);
					$new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:0\trealigned:0/;
					$match_read{$chr_pos_id} = $new_line;
					$count_realigned_read ++;
					$correct_mismatch_flag = 1;
					delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
					delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
				    }
				    else{
					$mismatch_count = &count_mismatch_I ($read12, $ref12);
					if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
					    $line[5] = $1 . 'M' . $2 . 'I' . ($indelpos - $posi + 1 - $1) . 'M' . $indelsize . 'D' . ($read_length - $indelpos + $posi - 1 - $2) . 'M';
					    my $new_line = join ("\t", @line);
					    $new_line =~ s/mismatch:\d+/mismatch:$mismatch_count\trealigned:0/;
					    my $delseq = substr ($chr_seq{$chr}, $indelpos, $indelsize);			    
					    $read1 = substr ($read12, 0, $indelpos - $posi + 1);
					    $read2 = substr ($read12, $indelpos - $posi + 1);
					    next if (length $read1 != length $ref1) or (length $read2 != length $ref2);
					    my $MDtag = &get_mismatch_info_D ($read1, $ref1, $read2, $ref2, $delseq);
					    $new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
					    $filtered_mismatch_read{$chr_pos_id} = $new_line;
					    $count_realigned_read ++;
					    $correct_mismatch_flag = 1;
					    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					}
				    }
				}
			    }
			}
		    }
		    elsif ($line[5] =~ /^(\d+)M(\d+)D(\d+)M$/){
			next if (($posi + $1 - 1 == $indelpos) and ($indel eq 'D') and ($indelsize == $2));
			if ($indel eq 'I'){
			    if (($indelpos - $posi + 1) < $read_length / 2){
				if ($posi + $1 - 1 > $indelpos){
				    if (($indelpos - $posi + 1) > $indelsize){
					my $ins_read = substr ($line[9], $indelpos - $posi + 1 - $indelsize, $indelsize);
					next if ($ins_read ne $insbase);
					$ref1 = substr ($chr_seq{$chr}, $posi - 1 + $indelsize, $1 - $indelsize);
					$ref2 = substr ($chr_seq{$chr}, $posi - 1 + $1 + $2, $3);
					$ref12 = $ref1 . $ref2;
					$read1 = substr ($line[9], 0, $indelpos - $posi + 1 - $indelsize);
					next if ($indelpos - $posi + 1 + $indelsize > $read_length - 1);
					$read2 = substr ($line[9], $indelpos - $posi + 1, $1 - $indelsize - length $read1);
					$read3 = substr ($line[9], $1, $3);
					$read12 = $read1 . $read2 . $read3;
					$line[3] += $indelsize;
					next if (length $read12 != length $ref12);
					if ($read12 eq $ref12){
					    $line[5] = ($indelpos - $posi + 1 - $indelsize) . 'M' . $indelsize . 'I' . ($1 - $indelpos + $posi - 1) . 'M' . $2 . 'D' . $3 . 'M';
					    my $new_line = join ("\t", @line);
					    $new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:0\trealigned:$indelsize/;
					    $match_read{$chr_pos_id} = $new_line;
					    $count_realigned_read ++;
					    $correct_mismatch_flag = 1;
					    delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
					    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					}
					else{
					    $mismatch_count = &count_mismatch_I ($read12, $ref12);
					    if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
						$line[5] = ($indelpos - $posi + 1 - $indelsize) . 'M' . $indelsize . 'I' . ($1 - $indelpos + $posi - 1) . 'M' . $2 . 'D' . $3 . 'M';
						my $new_line = join ("\t", @line);
						$new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:$mismatch_count\trealigned:$indelsize/;
						my $MDtag = &get_mismatch_info_I ($read12, $ref12);
						$new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
						$filtered_mismatch_read{$chr_pos_id} = $new_line;
						$count_realigned_read ++;
						$correct_mismatch_flag = 1;
						delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					    }
					}
				    }
				    else{
					my $ins_read = substr ($line[9], 0, $indelpos - $posi + 1);
					next if ($ins_read ne substr ($insbase, length ($insbase) - length ($ins_read)));
					$ref1 = substr ($chr_seq{$chr}, $indelpos, $1 - $indelpos + $posi - 1);
					$ref2 = substr ($chr_seq{$chr}, $posi - 1 + $1 + $2, $3);
					$ref12 = $ref1 . $ref2;
					$read1 = substr ($line[9], $indelpos - $posi + 1, $1 - $indelpos + $posi - 1);
					$read2 = substr ($line[9], $1, $3);
					$read12 = $read1 . $read2;
					next if (length $read12 != length $ref12);
					my $indelsize_term = $indelpos - $posi + 1;
					$line[3] = $indelpos + 1;
					if ($read12 eq $ref12){
					    $line[5] = ($indelpos - $posi + 1) . 'I' . ($1 - $indelpos + $posi - 1) . 'M' . $2 . 'D' . $3 . 'M';
					    my $new_line = join ("\t", @line);
					    $new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:0\trealigned:$indelsize_term/;
					    $match_read{$chr_pos_id} = $new_line;
					    $count_realigned_read ++;
					    $correct_mismatch_flag = 1;
					    delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
					    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					}
					else{
					    $mismatch_count = &count_mismatch_I ($read12, $ref12);
					    if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
						$line[5] = ($indelpos - $posi + 1) . 'I' . ($1 - $indelpos + $posi - 1) . 'M' . $2 . 'D' . $3 . 'M';
						my $new_line = join ("\t", @line);
						$new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:$mismatch_count\trealigned:$indelsize_term/;
						my $MDtag = &get_mismatch_info_I ($read12, $ref12);
						$new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
						$filtered_mismatch_read{$chr_pos_id} = $new_line;
						$count_realigned_read ++;
						$correct_mismatch_flag = 1;
						delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					    }
					}
				    }
				}
				else{				# deletion in read is neglected and the read is treaed as a normal read
				    if ($indelpos - $posi + 1 - $2 > $indelsize){
					my $ins_read = substr ($line[9], $indelpos - $posi + 1 - $indelsize - $2, $indelsize);
					next if ($ins_read ne $insbase);
					$ref12 = substr ($chr_seq{$chr}, $posi - 1 + $indelsize + $2, $read_length - $indelsize);
					$read1 = substr ($line[9], 0, $indelpos - $posi + 1 - $indelsize - $2);
					next if ($indelpos - $posi + 1 - $2 > $read_length - 1);
					$read2 = substr ($line[9], $indelpos - $posi + 1 - $2);
					$read12 = $read1 . $read2;
					my $indelsize_2 = $indelsize + $2;
					$line[3] += $indelsize_2;
					next if (length $read12 != length $ref12);
					if ($read12 eq $ref12){
					    $line[5] = ($indelpos - $posi + 1 - $indelsize - $2) . 'M' . $indelsize . 'I' . ($3 - $indelpos + $posi - 1 + $1 + $2) . 'M';
					    my $new_line = join ("\t", @line);
					    $new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:0\trealigned:$indelsize_2/;
					    $match_read{$chr_pos_id} = $new_line;
					    $count_realigned_read ++;
					    $correct_mismatch_flag = 1;
					    delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
					    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					}
					else{
					    $mismatch_count = &count_mismatch_I ($read12, $ref12);
					    if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
						$line[5] = ($indelpos - $posi + 1 - $indelsize - $2) . 'M' . $indelsize . 'I' . ($3 - $indelpos + $posi - 1 + $1 + $2) . 'M';
						my $new_line = join ("\t", @line);
						$new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:$mismatch_count\trealigned:$indelsize_2/;
						my $MDtag = &get_mismatch_info_I ($read12, $ref12);
						$new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
						$filtered_mismatch_read{$chr_pos_id} = $new_line;
						$count_realigned_read ++;
						$correct_mismatch_flag = 1;
						delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					    }
					}
				    }
				    elsif ($indelpos - $posi + 1 - $2 > 0){
					my $ins_read = substr ($line[9], 0, $indelpos - $posi + 1 - $2);
					next if ($ins_read ne substr ($insbase, length ($insbase) - length ($ins_read)));
					$read12 = substr ($line[9], $indelpos - $posi + 1 - $2);
					$ref12 = substr ($chr_seq{$chr}, $indelpos, length $read12);
					my $indelsize_term = $indelpos - $posi + 1;
					$line[3] = $indelpos + 1;
					next if (length $read12 != length $ref12);
					if ($read12 eq $ref12){
					    $line[5] = ($indelpos - $posi + 1 - $2) . 'I' . ($read_length - $indelpos + $posi - 1 + $2) . 'M';
					    my $new_line = join ("\t", @line);
					    $new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:0\trealigned:$indelsize_term/;
					    $match_read{$chr_pos_id} = $new_line;
					    $count_realigned_read ++;
					    $correct_mismatch_flag = 1;
					    delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
					    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					}
					else{
					    $mismatch_count = &count_mismatch_I ($read12, $ref12);
					    if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
						$line[5] = ($indelpos - $posi + 1 - $2) . 'I' . ($read_length - $indelpos + $posi - 1 + $2) . 'M';
						my $new_line = join ("\t", @line);
						$new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:$mismatch_count\trealigned:$indelsize_term/;
						my $MDtag = &get_mismatch_info_I ($read12, $ref12);
						$new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
						$filtered_mismatch_read{$chr_pos_id} = $new_line;
						$count_realigned_read ++;
						$correct_mismatch_flag = 1;
						delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					    }
					}
				    }
				}
			    }
			    else{
				if ($posi + $1 - 1 + $2 >= $indelpos){	# deletion in read is neglected and the read is treaed as a normal read
				    if ($read_length - $indelpos + $posi - 1 > $indelsize){
					my $ins_read = substr ($line[9], $indelpos - $posi + 1, $indelsize);
					next if ($ins_read ne $insbase);
					$ref1 = substr ($chr_seq{$chr}, $posi - 1, $1 - $indelsize);
					$ref2 = substr ($chr_seq{$chr}, $posi - 1 + $1 - $indelsize, $3);
					$ref12 = $ref1 . $ref2;
					$read1 = substr ($line[9], 0, $indelpos - $posi + 1);
					next if ($indelpos - $posi + 1 + $indelsize > $read_length - 1);
					$read2 = substr ($line[9], $indelpos - $posi + 1 + $indelsize);
					$read12 = $read1 . $read2;
					next if (length $read12 != length $ref12);
					if ($read12 eq $ref12){
					    $line[5] = ($indelpos - $posi + 1) . 'M' . $indelsize . 'I' . ($1 + $3 - $indelsize - $indelpos + $posi - 1) . 'M';
					    my $new_line = join ("\t", @line);
					    $new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:0\trealigned:0/;
					    $match_read{$chr_pos_id} = $new_line;
					    $count_realigned_read ++;
					    $correct_mismatch_flag = 1;
					    delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
					    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					}
					else{
					    $mismatch_count = &count_mismatch_I ($read12, $ref12);
					    if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
						$line[5] = ($indelpos - $posi + 1) . 'M' . $indelsize . 'I' . ($1 + $3 - $indelsize - $indelpos + $posi - 1) . 'M';
						my $new_line = join ("\t", @line);
						$new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:$mismatch_count\trealigned:0/;
						my $MDtag = &get_mismatch_info_I ($read12, $ref12);
						$new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
						$filtered_mismatch_read{$chr_pos_id} = $new_line;
						$count_realigned_read ++;
						$correct_mismatch_flag = 1;
						delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					    }
					}
				    }
				    else{
					my $ins_read = substr ($line[9], $indelpos - $posi + 1, $read_length - $indelpos + $posi - 1);
					next if ($ins_read ne substr ($insbase, 0, length ($ins_read)));
					$read12 = substr ($line[9], 0, $indelpos - $posi + 1);
					$ref12 = substr ($chr_seq{$chr}, $posi - 1, length $read12);
					next if (length $read12 != length $ref12);
					if ($read12 eq $ref12){
					    $line[5] = ($indelpos - $posi + 1) . 'M' . ($read_length - $indelpos + $posi - 1) . 'I';
					    my $new_line = join ("\t", @line);
					    $new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:0\trealigned:0/;
					    $match_read{$chr_pos_id} = $new_line;
					    $count_realigned_read ++;
					    $correct_mismatch_flag = 1;
					    delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
					    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					}
					else{
					    $mismatch_count = &count_mismatch_I ($read12, $ref12);
					    if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
						$line[5] = ($indelpos - $posi + 1) . 'M' . ($read_length - $indelpos + $posi - 1) . 'I';
						my $new_line = join ("\t", @line);
						$new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:$mismatch_count\trealigned:0/;
						my $MDtag = &get_mismatch_info_I ($read12, $ref12);
						$new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
						$filtered_mismatch_read{$chr_pos_id} = $new_line;
						$count_realigned_read ++;
						$correct_mismatch_flag = 1;
						delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					    }
					}
				    }
				}
				else{
				    if ($read_length - $indelpos + $posi - 1 + $2 > $indelsize){
					my $ins_read = substr ($line[9], $indelpos - $posi + 1 - $2, $indelsize);
					next if ($ins_read ne $insbase);
					$ref1 = substr ($chr_seq{$chr}, $posi - 1, $1);
					$ref2 = substr ($chr_seq{$chr}, $posi - 1 + $1 + $2, $3 - $indelsize);
					$ref12 = $ref1 . $ref2;
					$read1 = substr ($line[9], 0, $indelpos - $posi + 1 - $2);
					next if ($indelpos - $posi + 1 - $2 + $indelsize > $read_length - 1);
					$read2 = substr ($line[9], $indelpos - $posi + 1 - $2 + $indelsize);
					$read12 = $read1 . $read2;
					next if (length $read12 != length $ref12);
					if ($read12 eq $ref12){
					    $line[5] = $1 . 'M' . $2 . 'D' . ($indelpos - $posi + 1 - $1 - $2) . 'M' . $indelsize . 'I' . ($1 + $3 + $2 - $indelpos + $posi - 1 - $indelsize) . 'M';
					    my $new_line = join ("\t", @line);
					    $new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:0\trealigned:0/;
					    $match_read{$chr_pos_id} = $new_line;
					    $count_realigned_read ++;
					    $correct_mismatch_flag = 1;
					    delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
					    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					}
					else{
					    $mismatch_count = &count_mismatch_I ($read12, $ref12);
					    if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
						$line[5] = $1 . 'M' . $2 . 'D' . ($indelpos - $posi + 1 - $1 - $2) . 'M' . $indelsize . 'I' . ($1 + $3 + $2 - $indelpos + $posi - 1 - $indelsize) . 'M';
						my $new_line = join ("\t", @line);
						$new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:$mismatch_count\trealigned:0/;
						my $MDtag = &get_mismatch_info_I ($read12, $ref12);
						$new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
						$filtered_mismatch_read{$chr_pos_id} = $new_line;
						$count_realigned_read ++;
						$correct_mismatch_flag = 1;
						delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					    }
					}
				    }
				    else{
					my $ins_read = substr ($line[9], $indelpos - $posi + 1 - $2, $read_length - $indelpos + $posi - 1  + $2);
					next if ($ins_read ne substr ($insbase, 0, length ($ins_read)));
					$ref1 = substr ($chr_seq{$chr}, $posi - 1, $1);
					$ref2 = substr ($chr_seq{$chr}, $posi - 1 + $1 + $2, $indelpos - $posi + 1 - $2 - $1);
					$ref12 = $ref1 . $ref2;
					$read12 = substr ($line[9], 0, $indelpos - $posi + 1 - $2);
					next if (length $read12 != length $ref12);
					if ($read12 eq $ref12){
					    $line[5] = $1 . 'M' . $2 . 'D' . ($indelpos - $posi + 1 - $1 - $2) . 'M' . ($read_length - $indelpos + $posi - 1  + $2) . 'I';
					    my $new_line = join ("\t", @line);
					    $new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:0\trealigned:0/;
					    $match_read{$chr_pos_id} = $new_line;
					    $count_realigned_read ++;
					    $correct_mismatch_flag = 1;
					    delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
					    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					}
					else{
					    $mismatch_count = &count_mismatch_I ($read12, $ref12);
					    if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
						$line[5] = $1 . 'M' . $2 . 'D' . ($indelpos - $posi + 1 - $1 - $2) . 'M' . ($read_length - $indelpos + $posi - 1  + $2) . 'I';
						my $new_line = join ("\t", @line);
						$new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:$mismatch_count\trealigned:0/;
						my $MDtag = &get_mismatch_info_I ($read12, $ref12);
						$new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
						$filtered_mismatch_read{$chr_pos_id} = $new_line;
						$count_realigned_read ++;
						$correct_mismatch_flag = 1;
						delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					    }
					}
				    }
				}
			    }
			}
			elsif ($indel eq 'D'){
			    if (($indelpos <= $posi) or ($indelpos - $posi + 1 <= $read_length - $indelsize - $indelpos + $posi - 1)){
				if ($posi + $1 - 1 > $indelpos + $indelsize){
				    $ref1 = substr ($chr_seq{$chr}, $posi - 1 - $indelsize, $indelpos - $posi + 1 + $indelsize);
				    $ref2 = substr ($chr_seq{$chr}, $indelpos + $indelsize, $1 - length $ref1);
				    $ref3 = substr ($chr_seq{$chr}, $posi + $1 - 1 + $2, $3);
				    my $ref2_2 = $ref2 . $ref3;
				    $ref12 = $ref1 . $ref2 . $ref3;
				    $read12 = $line[9];
				    $line[3] -= $indelsize;
				    next if (length $read12 != length $ref12);
				    if ($read12 eq $ref12){
					$line[5] = ($indelpos - $posi + 1 + $indelsize) . 'M' . $indelsize . 'D' . ($1 - $indelsize - $indelpos + $posi - 1) . 'M' . $2 . 'D' . $3 . 'M';
					my $new_line = join ("\t", @line);
					$new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:0\trealigned:-$indelsize/;
					$match_read{$chr_pos_id} = $new_line;
					$count_realigned_read ++;
					$correct_mismatch_flag = 1;
					delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
					delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
				    }
				    else{
					$mismatch_count = &count_mismatch_I ($read12, $ref12);
					if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
					    $line[5] = ($indelpos - $posi + 1 + $indelsize) . 'M' . $indelsize . 'D' . ($1 - $indelsize - $indelpos + $posi - 1) . 'M' . $2 . 'D' . $3 . 'M';
					    my $new_line = join ("\t", @line);
					    $new_line =~ s/mismatch:\d+/mismatch:$mismatch_count\trealigned:-$indelsize/;
					    my $delseq = substr ($chr_seq{$chr}, $indelpos, $indelsize);					    
					    $read1 = substr ($read12, 0, $indelpos - $posi + 1 + $indelsize);
					    $read2 = substr ($read12, $indelpos - $posi + 1 + $indelsize);
					    next if (length $read1 != length $ref1) or (length $read2 != length $ref2);
					    my $MDtag = &get_mismatch_info_D ($read1, $ref1, $read2, $ref2_2, $delseq);
					    $new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
					    $filtered_mismatch_read{$chr_pos_id} = $new_line;
					    $count_realigned_read ++;
					    $correct_mismatch_flag = 1;
					    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					}
				    }
				}
				else{				# deletion in read is neglected and the read is treaed as a normal read
				    next if ($indelpos - $posi + 1 + $indelsize - $2 <= 0);
				    $ref1 = substr ($chr_seq{$chr}, $posi - 1 - $indelsize + $2, $indelpos - $posi + 1 + $indelsize - $2);
				    $ref2 = substr ($chr_seq{$chr}, $indelpos + $indelsize, $read_length - length $ref1);				
				    $ref12 = $ref1 . $ref2;
				    $read12 = $line[9];
				    my $shift_size = $indelsize - $2;
				    $line[3] -= $shift_size;
				    next if (length $read12 != length $ref12);
				    if ($read12 eq $ref12){
					$line[5] = ($indelpos - $posi + 1 + $indelsize - $2) . 'M' . $indelsize . 'D' . ($1 + $3 - $indelpos + $posi - 1 - $indelsize + $2) . 'M';
					my $new_line = join ("\t", @line);
					$new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:0\trealigned:-$shift_size/ if ($shift_size >= 0);
					$new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:0\trealigned:$shift_size/ if ($shift_size < 0);
					$match_read{$chr_pos_id} = $new_line;
					$count_realigned_read ++;
					$correct_mismatch_flag = 1;
					delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
					delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
				    }
				    else{
					$mismatch_count = &count_mismatch_I ($read12, $ref12);
					if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
					    $line[5] = ($indelpos - $posi + 1 + $indelsize - $2) . 'M' . $indelsize . 'D' . ($1 + $3 - $indelpos + $posi - 1 - $indelsize + $2) . 'M';
					    my $delseq = substr ($chr_seq{$chr}, $indelpos, $indelsize);					    
					    $read1 = substr ($read12, 0, $indelpos - $posi + 1 + $indelsize - $2);
					    $read2 = substr ($read12, $indelpos - $posi + 1 + $indelsize - $2);
					    next if (length $read1 != length $ref1) or (length $read2 != length $ref2);
					    my $MDtag = &get_mismatch_info_D ($read1, $ref1, $read2, $ref2, $delseq);
					    my $new_line = join ("\t", @line);
					    $new_line =~ s/mismatch:\d+/mismatch:$mismatch_count\trealigned:-$shift_size/ if ($shift_size >= 0);
					    $new_line =~ s/mismatch:\d+/mismatch:$mismatch_count\trealigned:$shift_size/ if ($shift_size < 0);
					    $new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
					    $filtered_mismatch_read{$chr_pos_id} = $new_line;
					    $count_realigned_read ++;
					    $correct_mismatch_flag = 1;
					    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					}
				    }
				}
			    }
			    else{
				if ($posi + $1 - 1 + $2 >= $indelpos){	# deletion in read is neglected and the read is treaed as a normal read
				    $ref1 = substr ($chr_seq{$chr}, $posi - 1, $indelpos - $posi + 1);
				    $ref2 = substr ($chr_seq{$chr}, $indelpos + $indelsize, $read_length - length $ref1);
				    $ref12 = $ref1 . $ref2;
				    $read12 = $line[9];
				    next if (length $read12 != length $ref12);
				    if ($read12 eq $ref12){
					$line[5] = ($indelpos - $posi + 1) . 'M' . $indelsize . 'D' . ($1 + $3 - $indelpos + $posi - 1) . 'M';
					my $new_line = join ("\t", @line);
					$new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:0\trealigned:0/;
					$match_read{$chr_pos_id} = $new_line;
					$count_realigned_read ++;
					$correct_mismatch_flag = 1;
					delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
					delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
				    }
				    else{
					$mismatch_count = &count_mismatch_I ($read12, $ref12);
					if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
					    $line[5] = ($indelpos - $posi + 1) . 'M' . $indelsize . 'D' . ($1 + $3 - $indelpos + $posi - 1) . 'M';
					    my $new_line = join ("\t", @line);
					    $new_line =~ s/mismatch:\d+/mismatch:$mismatch_count\trealigned:0/;
					    my $delseq = substr ($chr_seq{$chr}, $indelpos, $indelsize);					    
					    $read1 = substr ($read12, 0, $indelpos - $posi + 1);
					    $read2 = substr ($read12, $indelpos - $posi + 1);
					    next if (length $read1 != length $ref1) or (length $read2 != length $ref2);
					    my $MDtag = &get_mismatch_info_D ($read1, $ref1, $read2, $ref2, $delseq);
					    $new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
					    $filtered_mismatch_read{$chr_pos_id} = $new_line;
					    $count_realigned_read ++;
					    $correct_mismatch_flag = 1;
					    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					}
				    }
				}
				else{
				    $ref1 = substr ($chr_seq{$chr}, $posi - 1, $1);
				    next if ($indelpos - $posi + 1 - $1 - $2 <= 0);
				    $ref2 = substr ($chr_seq{$chr}, $posi - 1 + $1 + $2, $indelpos - $posi + 1 - $1 - $2);
				    $ref3 = substr ($chr_seq{$chr}, $indelpos + $indelsize, $read_length - $1 - length $ref2);
				    my $ref1_2 = $ref1 . $ref2;
				    $ref12 = $ref1 . $ref2 . $ref3;
				    $read12 = $line[9];
				    next if (length $read12 != length $ref12);
				    if ($read12 eq $ref12){
					$line[5] = $1 . 'M' . $2 . 'D' . ($indelpos - $posi + 1 - $1) . 'M' . $indelsize . 'D' . ($1 + $3 - $indelpos + $posi - 1) . 'M';
					my $new_line = join ("\t", @line);
					$new_line =~ s/mismatch:\d+\tindel:\d+/mismatch:0\trealigned:0/;
					$match_read{$chr_pos_id} = $new_line;
					$count_realigned_read ++;
					$correct_mismatch_flag = 1;
					delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
					delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
				    }
				    else{
					$mismatch_count = &count_mismatch_I ($read12, $ref12);
					if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
					    $line[5] = $1 . 'M' . $2 . 'D' . ($indelpos - $posi + 1 - $1) . 'M' . $indelsize . 'D' . ($1 + $3 - $indelpos + $posi - 1) . 'M';
					    my $delseq = substr ($chr_seq{$chr}, $indelpos, $indelsize);					    
					    $read1 = substr ($read12, 0, $indelpos - $posi + 1 - $2);
					    $read2 = substr ($read12, $indelpos - $posi + 1 - $2);
					    $ref1 = substr ($chr_seq{$chr}, $posi - 1, $1);
					    next if (length $read1 != length $ref1) or (length $read2 != length $ref2);
					    my $MDtag = &get_mismatch_info_D ($read1, $ref1_2, $read2, $ref3, $delseq);
					    my $new_line = join ("\t", @line);
					    $new_line =~ s/mismatch:\d+/mismatch:$mismatch_count\trealigned:0/;
					    $new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
					    $filtered_mismatch_read{$chr_pos_id} = $new_line;
					    $count_realigned_read ++;
					    $correct_mismatch_flag = 1;
					    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					}
				    }
				}
			    }
			}
		    }
		    elsif (($line[5] =~ /^(\d+)S(\d+)M$/) or ($line[5] =~ /^(\d+)M(\d+)S$/)){
			my $total_mismatch = $mismatch_num + $SC_num;
			my $shift_pos;
			if ($indel eq 'I'){
			    if (($indelpos - $posi + 1) < $read_length / 2){
				if (($indelpos - $posi + 1) > $indelsize){
				    my $ins_read = substr ($line[9], $indelpos - $posi + 1 - $indelsize, $indelsize);
				    next if ($ins_read ne $insbase);
				    $ref12 = substr ($chr_seq{$chr}, $posi - 1 + $indelsize, $read_length - $indelsize);
				    $read1 = substr ($line[9], 0, $indelpos - $posi + 1 - $indelsize);
				    next if ($indelpos - $posi + 1 + $indelsize > $read_length - 1);
				    $read2 = substr ($line[9], $indelpos - $posi + 1);
				    $read12 = $read1 . $read2;
				    $line[3] = $posi + $indelsize;
				    $shift_pos = $indelsize - $SC_num if ($line[5] =~ /^\d+S\d+M$/);
				    $shift_pos = $indelsize if ($line[5] =~ /^\d+M\d+S$/);
				    next if (length $read12 != length $ref12);
				    if ($read12 eq $ref12){
					$line[5] = ($indelpos - $posi + 1 - $indelsize) . 'M' . $indelsize . 'I' . ($read_length - $indelpos + $posi - 1) . 'M';
					my $new_line = join ("\t", @line);
					$new_line =~ s/mismatch:\d+\tsoftclip:\d+/mismatch:0\trealigned:$shift_pos/;
					my $MDtag = $read_length - $indelsize;
					$new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
					$match_read{$chr_pos_id} = $new_line;
					$count_realigned_read ++;
					$correct_mismatch_flag = 1;
					delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
					delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
				    }
				    else{
					$mismatch_count = &count_mismatch_I ($read12, $ref12);
					if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
					    $line[5] = ($indelpos - $posi + 1 - $indelsize) . 'M' . $indelsize . 'I' . ($read_length - $indelpos + $posi - 1) . 'M';
					    my $new_line = join ("\t", @line);
					    $new_line =~ s/mismatch:\d+\tsoftclip:\d+/mismatch:$mismatch_count\trealigned:$shift_pos/;
					    my $MDtag = &get_mismatch_info_I ($read12, $ref12);
					    $new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
					    $filtered_mismatch_read{$chr_pos_id} = $new_line;
					    $count_realigned_read ++;
					    $correct_mismatch_flag = 1;
					    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					}
				    }
				}
				else{
				    my $ins_read = substr ($line[9], 0, $indelpos - $posi + 1);
				    next if ($ins_read ne substr ($insbase, length ($insbase) - length ($ins_read)));
				    $ref12 = substr ($chr_seq{$chr}, $indelpos, $read_length - $indelpos + $posi - 1);
				    $read12 = substr ($line[9], $indelpos - $posi + 1);
				    $line[3] = $indelpos + 1;
				    $shift_pos = $indelpos - $posi + 1 - $SC_num if ($line[5] =~ /^\d+S\d+M$/);
				    $shift_pos = $indelpos - $posi + 1 if ($line[5] =~ /^\d+M\d+S$/);
				    next if (length $read12 != length $ref12);
				    if ($read12 eq $ref12){
					$line[5] = ($indelpos - $posi + 1) . 'I' . ($read_length - $indelpos + $posi - 1) . 'M';
					my $new_line = join ("\t", @line);
					$new_line =~ s/mismatch:\d+/mismatch:0\trealigned:$shift_pos/;
					my $MDtag = $read_length - $indelsize;
					$new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
					$match_read{$chr_pos_id} = $new_line;
					$count_realigned_read ++;
					$correct_mismatch_flag = 1;
					delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
					delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
				    }
				    else{
					$mismatch_count = &count_mismatch_I ($read12, $ref12);
					if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
					    $line[5] = ($indelpos - $posi + 1) . 'I' . ($read_length - $indelpos + $posi - 1) . 'M';
					    my $new_line = join ("\t", @line);
					    $new_line =~ s/mismatch:\d+/mismatch:$mismatch_count\trealigned:$shift_pos/;
					    my $MDtag = &get_mismatch_info_I ($read12, $ref12);
					    $new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
					    $filtered_mismatch_read{$chr_pos_id} = $new_line;
					    $count_realigned_read ++;
					    $correct_mismatch_flag = 1;
					    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					}
				    }
				}
			    }
			    else{
				if ($read_length - $indelpos + $posi - 1 > $indelsize){
				    my $ins_read = substr ($line[9], $indelpos - $posi + 1, $indelsize);
				    next if ($ins_read ne $insbase);
				    $ref12 = substr ($chr_seq{$chr}, $posi - 1, $read_length - $indelsize);
				    $read1 = substr ($line[9], 0, $indelpos - $posi + 1);
				    next if ($indelpos - $posi + 1 + $indelsize > $read_length - 1);
				    $read2 = substr ($line[9], $indelpos - $posi + 1 + $indelsize);
				    $read12 = $read1 . $read2;
				    $line[3] = $posi;
				    $shift_pos = $SC_num if ($line[5] =~ /^(\d+)S(\d+)M$/);
				    $shift_pos = 0 if ($line[5] =~ /^(\d+)M(\d+)S$/);
				    next if (length $read12 != length $ref12);
				    if ($read12 eq $ref12){
					$line[5] = ($indelpos - $posi + 1) . 'M' . $indelsize . 'I' . ($read_length - $indelpos + $posi- 1 - $indelsize) . 'M';
					my $new_line = join ("\t", @line);
					$new_line =~ s/mismatch:\d+\tsoftclip:\d+/mismatch:0\trealigned:$shift_pos/;
					my $MDtag = $read_length - $indelsize;
					$new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
					$match_read{$chr_pos_id} = $new_line;
					$count_realigned_read ++;
					$correct_mismatch_flag = 1;
					delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
					delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
				    }
				    else{
					$mismatch_count = &count_mismatch_I ($read12, $ref12);
					if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
					    $line[5] = ($indelpos - $posi + 1) . 'M' . $indelsize . 'I' . ($read_length - $indelpos + $posi - 1 - $indelsize) . 'M';
					    my $new_line = join ("\t", @line);
					    $new_line =~ s/mismatch:\d+\tsoftclip:\d+/mismatch:$mismatch_count\trealigned:$shift_pos/;
					    my $MDtag = &get_mismatch_info_I ($read12, $ref12);
					    $new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
					    $filtered_mismatch_read{$chr_pos_id} = $new_line;
					    $count_realigned_read ++;
					    $correct_mismatch_flag = 1;
					    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					}
				    }
				}
				else{					    
				    my $ins_read = substr ($line[9], $indelpos - $posi + 1);
				    next if ($ins_read ne substr ($insbase, 0, length ($ins_read)));
				    $ref12 = substr ($chr_seq{$chr}, $posi - 1, $indelpos - $posi + 1);
				    $read12 = substr ($line[9], 0, $indelpos - $posi + 1);
				    $line[3] = $posi;
				    $shift_pos = -$SC_num if ($line[5] =~ /^\d+S\d+M$/);
				    $shift_pos = 0 if ($line[5] =~ /^\d+M\d+S$/);
				    next if (length $read12 != length $ref12);
				    if ($read12 eq $ref12){
					$line[5] = ($indelpos - $posi + 1) . 'M' . ($read_length - $indelpos + $posi - 1) . 'I';
					my $new_line = join ("\t", @line);
					$new_line =~ s/mismatch:\d+/mismatch:0\trealigned:$shift_pos/;
					my $MDtag = $read_length - $indelsize;
					$new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
					$match_read{$chr_pos_id} = $new_line;
					$count_realigned_read ++;
					$correct_mismatch_flag = 1;
					delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
					delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
				    }
				    else{
					$mismatch_count = &count_mismatch_I ($read12, $ref12);
					if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
					    $line[5] = ($indelpos - $posi + 1) . 'M' . ($read_length - $indelpos + $posi - 1) . 'I';
					    my $new_line = join ("\t", @line);
					    $new_line =~ s/mismatch:\d+/mismatch:$mismatch_count\trealigned:$shift_pos/;
					    my $MDtag = &get_mismatch_info_I ($read12, $ref12);
					    $new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
					    $filtered_mismatch_read{$chr_pos_id} = $new_line;
					    $count_realigned_read ++;
					    $correct_mismatch_flag = 1;
					    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
					}
				    }
				}
			    }
			}
			elsif ($indel eq 'D'){
			    if (($indelpos <= $posi) or ($indelpos - $posi + 1 <= $read_length - $indelsize - $indelpos + $posi - 1)){
				$ref1 = substr ($chr_seq{$chr}, $posi - 1 - $indelsize, $indelpos - $posi + 1 + $indelsize);
				$ref2 = substr ($chr_seq{$chr}, $indelpos + $indelsize, $read_length - length $ref1);
				$ref12 = $ref1 . $ref2;
				$read12 = $line[9];
				$line[3] = $posi - $indelsize;
				$shift_pos = -$indelsize - $SC_num if ($line[5] =~ /^\d+S\d+M$/);
				$shift_pos = -$indelsize if ($line[5] =~ /^\d+M\d+S$/);
				next if (length $read12 != length $ref12);
				if ($read12 eq $ref12){
				    $line[5] = ($indelpos - $posi + 1 + $indelsize) . 'M' . $indelsize . 'D' . ($read_length - $indelpos + $posi - 1 - $indelsize) . 'M';
				    my $new_line = join ("\t", @line);
				    $new_line =~ s/mismatch:\d+\tsoftclip:\d+/mismatch:0\trealigned:$shift_pos/;
				    my $delseq = substr ($chr_seq{$chr}, $indelpos, $indelsize);
				    my $MDtag = length $read1 . '^' . $delseq . length $read2;
				    $new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
				    $match_read{$chr_pos_id} = $new_line;
				    $count_realigned_read ++;
				    $correct_mismatch_flag = 1;
				    delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
				    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
				}
				else{
				    $mismatch_count = &count_mismatch_I ($read12, $ref12);
				    if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
					$line[5] = ($indelpos - $posi + 1 + $indelsize) . 'M' . $indelsize . 'D' . ($read_length - $indelpos + $posi - 1 - $indelsize) . 'M';
					my $new_line = join ("\t", @line);
					$new_line =~ s/mismatch:\d+\tsoftclip:\d+/mismatch:$mismatch_count\trealigned:$shift_pos/;
					my $delseq = substr ($chr_seq{$chr}, $indelpos, $indelsize);					
					$read1 = substr ($read12, 0, $indelpos - $posi + 1 + $indelsize);
					$read2 = substr ($read12, $indelpos - $posi + 1 + $indelsize);
					next if (length $read1 != length $ref1) or (length $read2 != length $ref2);
					my $MDtag = &get_mismatch_info_D ($read1, $ref1, $read2, $ref2, $delseq);
					$new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
					$filtered_mismatch_read{$chr_pos_id} = $new_line;
					$count_realigned_read ++;
					$correct_mismatch_flag = 1;
					delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
				    }
				}
			    }
			    else{
				$ref1 = substr ($chr_seq{$chr}, $posi - 1, $indelpos - $posi + 1);
				$ref2 = substr ($chr_seq{$chr}, $indelpos + $indelsize, $read_length - length $ref1);
				$ref12 = $ref1 . $ref2;
				$read12 = $line[9];
				$line[3] = $posi;
				$shift_pos = -$SC_num if ($line[5] =~ /^\d+S\d+M$/);
				$shift_pos = 0 if ($line[5] =~ /^\d+M\d+S$/);
				next if (length $read12 != length $ref12);
				if ($read12 eq $ref12){
				    $line[5] = ($indelpos - $posi + 1) . 'M' . $indelsize . 'D' . ($read_length - $indelpos + $posi - 1) . 'M';
				    my $new_line = join ("\t", @line);
				    $new_line =~ s/mismatch:\d+\tsoftclip:\d+/mismatch:0\trealigned:$shift_pos/;
				    my $delseq = substr ($chr_seq{$chr}, $indelpos, $indelsize);
				    my $MDtag = length $read1 . '^' . $delseq . length $read2;
				    $new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
				    $match_read{$chr_pos_id} = $new_line;
				    $count_realigned_read ++;
				    $correct_mismatch_flag = 1;
				    delete $filtered_mismatch_read{$chr_pos_id} if ($read_correct == 1);
				    delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
				}
				else{
				    $mismatch_count = &count_mismatch_I ($read12, $ref12);
				    if (($mismatch_num >= $mismatch_count) and (($read_correct == 1) or (($read_correct == 0) and ($mismatch_count <= $limited_mismatch_num - 1)))){
					$line[5] = ($indelpos - $posi + 1) . 'M' . $indelsize . 'D' . ($read_length - $indelpos + $posi - 1) . 'M';
					my $new_line = join ("\t", @line);
					$new_line =~ s/mismatch:\d+\tsoftclip:\d+/mismatch:$mismatch_count\trealigned:$shift_pos/;
					my $delseq = substr ($chr_seq{$chr}, $indelpos, $indelsize);
					$read1 = substr ($read12, 0, $indelpos - $posi + 1);
					$read2 = substr ($read12, $indelpos - $posi + 1);
					next if (length $read1 != length $ref1) or (length $read2 != length $ref2);
					my $MDtag = &get_mismatch_info_D ($read1, $ref1, $read2, $ref2, $delseq);
					$new_line =~ s/MD:Z:\S+/MD:Z:$MDtag/;
					$filtered_mismatch_read{$chr_pos_id} = $new_line;
					$count_realigned_read ++;
					$correct_mismatch_flag = 1;
					delete $filtered_pairs{$read_id} if (($filter_pairs == 1) and (exists $filtered_pairs{$read_id}));
				    }
				}
			    }
			}
		    }
		    push (@corrected_read_pos, $count_mismatch_line) if ($correct_mismatch_flag == 1);
		}
		if (@corrected_read_pos > 0){
		    if (@corrected_read_pos == @mismatch_read_lines){
			delete $mismatch_read{$chr_pos};
		    }
		    else{
			$mismatch_read{$chr_pos} = '-';
			for (my $i = 0; $i < @mismatch_read_lines; $i++){
			    my $match_pos = 0;
			    foreach my $del_pos (@corrected_read_pos){
				$match_pos = 1 if ($i == $del_pos - 1);
			    }
			    if ($match_pos == 0){
				if ($mismatch_read{$chr_pos} eq '-'){
				    $mismatch_read{$chr_pos} = $mismatch_read_lines[$i];
				}
				else{
				    $mismatch_read{$chr_pos} .= '||' . $mismatch_read_lines[$i];
				}
			    }
			}
		    }
		}
	    }
	}
    }
}

sub count_mismatch_T {
    my ($read, $ref, $qual, $chr, $ref_start) = @_;
    my $count = 0;
    my $count_lq = 0;
    my $read_pos = 0;
    my $terminal_mismatch_N = 0;
    my $terminal_mismatch_C = 0;
    my @qual = split (//, $qual);
    my $MD = '';
    return ($limited_mismatch_num + 1, 0, 0, $MD) if (length $read != length $ref);
    my $mismatch_pos = 0;
    for (my $i = 0; $i < length $read; $i++){
	$read_pos ++;
	my $read_i = substr ($read, $i, 1);
	my $ref_i = substr ($ref, $i, 1);
	if ($read_i ne $ref_i){
	    my $mismatch_qual = ord($qual[$i]) - 64 if ($base_call_type eq 'illumina');
	    $mismatch_qual = ord($qual[$i]) - 33 if ($base_call_type eq 'sanger');
	    $count++;
	    $count_lq++ if ($mismatch_qual < $min_mismatch_qual);
	    $terminal_mismatch_N ++ if (($read_pos == 1) or ($read_pos == 2));
	    $terminal_mismatch_C ++ if (($read_pos == length $read) or ($read_pos == length ($read) - 1));
	    if ($MDtag_sam == 0){
		$MD .= ($i - $mismatch_pos) . $ref_i;
		$mismatch_pos = $i + 1;
	    }
	}
    }
    $MD .= length ($read) - $mismatch_pos;
    $m1++ if ($count == 1);
    $m2++ if ($count == 2);
    $m3++ if ($count == 3);
    $m4++ if ($count == 4);
    $m5++ if ($count == 5);
    $m6++ if ($count == 6);
    $m7++ if ($count == 7);
    $m8++ if ($count == 8);
    $m9++ if ($count == 9);
    $m10++ if ($count == 10);
    $m10o++ if ($count > 10);
    return ($count, $count_lq, $terminal_mismatch_N, $terminal_mismatch_C, $MD);
}
            
sub count_mismatch_I {
    my ($read, $ref) = @_;
    my $count = 0;
    return $limited_mismatch_num + 1 if (length $read != length $ref);
    for (my $i = 0; $i < length $read; $i++){
	my $read_i = substr ($read, $i, 1);
	my $ref_i = substr ($ref, $i, 1);
	if ($read_i ne $ref_i){
	    $count++;
	}
    }
    return $count;
}

sub get_mismatch_info_I {
    my ($read_12, $ref_12) = @_;
    my $MD = '';
    my $mismatch_pos = 0;
    for (my $i = 0; $i < length $read_12; $i++){
	my $read_i = substr ($read_12, $i, 1);
	my $ref_i = substr ($ref_12, $i, 1);
	if ($read_i ne $ref_i){
	    $MD = $MD . ($i - $mismatch_pos) . $ref_i;
	    $mismatch_pos = $i + 1;
	}
    }
    $MD = $MD . (length ($read_12) - $mismatch_pos);
    return $MD;
}

sub get_mismatch_info_D {
    my ($read_1, $ref_1, $read_2, $ref_2, $del_seq) = @_;
    my $MD = '';
    my $mismatch_pos = 0;
    for (my $i = 0; $i < length $read_1; $i++){
	my $read_i = substr ($read_1, $i, 1);
	my $ref_i = substr ($ref_1, $i, 1);
	if ($read_i ne $ref_i){
	    $MD = $MD . ($i - $mismatch_pos) . $ref_i;
	    $mismatch_pos = $i + 1;
	}
    }
    $MD = $MD . (length ($read_1) - $mismatch_pos) . '^' . $del_seq;
    $mismatch_pos = 0;
    for (my $i = 0; $i < length $read_2; $i++){
	my $read_i = substr ($read_2, $i, 1);
	my $ref_i = substr ($ref_2, $i, 1);
	if ($read_i ne $ref_i){
	    $MD = $MD . ($i - $mismatch_pos) . $ref_i;
	    $mismatch_pos = $i + 1;
	}
    }
    $MD = $MD . (length ($read_2) - $mismatch_pos);
    return $MD;
}

sub correct_read{
    my ($read, $ciagr, $start, $ref_base, $tag, $qual) = @_;
    my $mod_start = -1;
    if (($ciagr =~ /^\d+M$/) or ($ciagr =~ /^\d+M\d+S$/) or ($ciagr =~ /^\d+M\d+I$/) or ($ciagr =~ /^\d+M\d+D\d+M$/) or ($ciagr =~ /^\d+M\d+D\d+M\d+S$/) or ($ciagr =~ /^(\d+)M(\d+)D(\d+)M\d+I$/) or ($ciagr =~ /^\d+M\d+D\d+M\d+D\d+M$/)){
	$mod_start = $start;
    }
    elsif (($ciagr =~ /^(\d+)S\d+M$/) or ($ciagr =~ /^(\d+)S\d+M\d+D\d+M$/) or ($ciagr =~ /^(\d+)S\d+M\d+I$/)){
	$mod_start = $1 + $start;
    }
    elsif (($ciagr =~ /^(\d+)M(\d+)I\d+M$/) or ($ciagr =~ /^(\d+)M(\d+)I\d+M\d+D\d+M$/) or ($ciagr =~ /^(\d+)M(\d+)I\d+M\d+I$/)){
	if ($start > $1){
	    $mod_start = $2 + $start;
	}
	else{
	    $mod_start = $start;
	}
    }
    elsif (($ciagr =~ /^(\d+)I\d+M$/) or ($ciagr =~ /^(\d+)I\d+M\d+S$/) or ($ciagr =~ /^(\d+)I\d+M\d+I$/) or ($ciagr =~ /^(\d+)I\d+M\d+D\d+M$/)){
	$mod_start = $1 + $start;
    }
    elsif ($ciagr =~ /^(\d+)S(\d+)M(\d+)I\d+M$/){
	if ($start > $2){
	    $mod_start = $1 + $3 + $start;
	}
	else{
	    $mod_start = $1 + $start;
	}
    }
    elsif ($ciagr =~ /^(\d+)I(\d+)M(\d+)I\d+M$/){
	if ($start > $2){
	    $mod_start = $1 + $3 + $start;
	}
	else{
	    $mod_start = $1 + $start;
	}
    }
    elsif ($ciagr =~ /^(\d+)M(\d+)I(\d+)M(\d+)I\d+M$/){
	if ($start > $1){
	    $mod_start = $2 + $start;
	}
	elsif ($start > $1 + $3){
	    $mod_start = $2 + $4 + $start;
	}
	else{
	    $mod_start = $start;
	}
    }
    elsif ($ciagr =~ /^(\d+)M\d+D(\d+)M(\d+)I\d+M$/){
	if ($start > $1 + $2){
	    $mod_start = $3 + $start;
	}
	else{
	    $mod_start = $start;
	}
    }
    my $read_base = '';
    my $read_qual = '';
    return ($read) if (length $read < $mod_start) and ($tag eq 'correct');
    return ($ref_base, substr ($qual, 0, 1)) if (length $read < $mod_start) and ($tag eq 'rbase');
    $read_base = substr ($read, $mod_start - 1, 1, $ref_base) if ($mod_start >= 0);
    $read_qual = substr ($qual, $mod_start - 1, 1) if ($mod_start >= 0) and ($tag eq 'rbase');
    if ($tag eq 'rbase'){
	return ($read_base, $read_qual);
    }
    elsif ($tag eq 'correct'){
	return ($read);
    }
}
