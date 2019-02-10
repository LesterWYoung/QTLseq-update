#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;

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
my $out_prefix = 'coval-refine'; # -p
my $help;			 # -h
my $soap_aligner = 0;            # -b

my $multiple_allele_sample = 0;  # -ma
my $com_path = '';

my $arg = join (' ', @ARGV);

if ($0 =~ /^(.*)\/coval/){
    $com_path = $1;
}

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

die "Reference file is not found: $!\n" unless (-f $reference_file);
die "Prefix of output is not specified: $!\n" if ($out_prefix eq '');
die "Read type is not properly specified: $!\n" unless (uc $read_type eq 'PE') or (uc $read_type eq 'SE');
die "Option --reftype is not properly specified: $!\n" if ($reference_type ne 'COMPLETE') and ($reference_type ne 'DRAFT') and ($reference_type ne 'draft');

$limited_mismatch_num_pairs = $limited_mismatch_num * 1.7 if ($limited_mismatch_num_pairs == 0);
if ($error_correct_mode == 1){
    $read_correct = 1;
    $filter_pairs = 1;
}

if ($targeted_align_mode == 1){
    if ($read_correct == 1){
	print STDERR "\n########## Warning! ##########\nThe job is changed from the error correction mode to the basic mode when the --talign option is set\n";
	$read_correct = 0;
    }
    $filter_pairs = 1;
    $limited_mismatch_num_pairs = $limited_mismatch_num;
}

unless (@ARGV){
print "Your input file should be added to the argument.\n";
exit;
}
my $input_file = shift (@ARGV);
if ($input_file eq '-'){
}
elsif (!-f $input_file){
    die "inputfile does not exist: $!\n";
}

my $multi_sample = 'OFF' if ($multi_sample_mode eq 'none');
$multi_sample = 'ON: homo' if ($multi_sample_mode eq 'homo');
$multi_sample = 'ON: hetero' if ($multi_sample_mode eq 'hetero');
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

my $current_mode = 'Error correction mode';
$current_mode = 'Error correction (multi-sample) mode' if ($multi_sample_mode == 1);
$current_mode = 'Basic mode' if ($read_correct == 0) and ($targeted_align_mode == 0);
$current_mode = 'Basic mode (Targeted align mode)' if ($read_correct == 0) and ($targeted_align_mode == 1);
my $realign_mode = 'ON';
$realign_mode = 'OFF' if ($dis_realign == 1);
my $filter_pair_func = 'ON';
$filter_pair_func = 'OFF' if ($filter_pairs == 0);
my $filt_discordant_read = 'ON';
$filt_discordant_read = 'OFF' if ($include_disc_read == 1);
my $filter_unmapped_read = 'ON';
$filter_unmapped_read = 'OFF' if ($include_unmap_read == 1);


if ($limited_mismatch_rate == 0){
    if ($read_correct == 0){
	print STDERR "########## coval refine options ##########\n<$current_mode>\nlocal_realignment=$realign_mode filter_pairs=$filter_pair_func filter_discordant-read=$filt_discordant_read filter_unmapped_read=$filter_unmapped_read\n";
        print STDERR "num=$limited_mismatch_num fnum=$limited_mismatch_num_pairs qmap=$minimum_mapq avins=$ave_insert_size_SD insd=$maximum_insert_SD_fold type=$read_type minq=$min_mismatch_qual out_sam=$out_sam\n";
	
    }
    else{
        print STDERR "########## coval refine options ##########\n<$current_mode>\tmulti-sample=$multi_sample\tqual_ave=$minimum_ave_alt_qual allele_freq=$minimum_allele_freq\nlocal_realignment=$realign_mode filter_pairs=$filter_pair_func filter_discordant_read=$filt_discordant_read filter_unmapped_read=$filter_unmapped_read\n";
        print STDERR "num=$limited_mismatch_num fnum=$limited_mismatch_num_pairs qmap=$minimum_mapq avins=$ave_insert_size_SD insd=$maximum_insert_SD_fold type=$read_type minq=$min_mismatch_qual out_sam=$out_sam\n";
    }
    print STDERR "ref=$reference_file\ninput=$input_file\nprefix_out=$out_prefix\n\n";
}
else{
    if ($read_correct == 0){
	print STDERR "########## coval refine options ##########\n<$current_mode>\nlocal_realignment=$realign_mode filter_pairs=$filter_pair_func filter_discordant-read=$filt_discordant_read filter_unmapped_read=$filter_unmapped_read\n";
        print STDERR "mrate=$limited_mismatch_rate qmap=$minimum_mapq avins=$ave_insert_size_SD insd=$maximum_insert_SD_fold type=$read_type minq=$min_mismatch_qual out_sam=$out_sam\n";
    }
    else{
	print STDERR "########## coval refine options ##########\n<$current_mode>\tmulti-sample=$multi_sample\tqual_ave=$minimum_ave_alt_qual allele_freq=$minimum_allele_freq\nlocal_realignment=$realign_mode filter_pairs=$filter_pair_func filter_discordant_read=$filt_discordant_read filter_unmapped_read=$filter_unmapped_read\n";
        print STDERR "mrate=$limited_mismatch_rate qmap=$minimum_mapq avins=$ave_insert_size_SD insd=$maximum_insert_SD_fold type=$read_type minq=$min_mismatch_qual out_sam=$out_sam\n";
    }
    print STDERR "ref=$reference_file\ninput=$input_file\nprefix_out=$out_prefix\n\n";
}
=head1 SYNOPSIS

 coval refine [options] <input_sorted_bam/sam_file> (Read names and reference sequence names in sam file should not contain characters '||' or '='.)
  Output:
   out_prefix.bam/sam

  Options:
   --ref or -r <STR>	   reference fasta file used for the alignment
   --pref or -p <STR>      prefix of output file
   --num or -n <INT>       maximum number of mismatches contained in a read [default: 2] (incompatible with --mrate)
   --mrate or -m <FLOAT>   maximum rate of mismatches contained in a read [0..1.0] (incompatible with --num)
   --fnum or -f <INT>	   maximum number of total mismatches contained in two paired reads [default: 2 * <INT> specified with --num] (incompatible with --mrate)
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
                          (-ms homo: -x -xf 0.8, -ms hetero: -x -xf 0.3 is automatically set. The -xf value can be overwritten by setting -xf separately.) [default: none]
   --mallel or -ma         allow multiple non-reference alleles for each sample in multi-sample mode (e.g., alleles A, C, and T for samples-1, -2, and -3) [default: false]
   
   ### Preset options ###
   --ec_md or -xm	   set options suitable for 'error correction mode' (equal to '--err_cor --fpair') [default: false]
   --talign_md or -tm      set options suitable for 'targeted' alignment (equal to '--fpair --fnum <INT specified with --num>') [default: false]

=cut

$arg =~ s/$input_file//;
$arg =~ s/[\s]{2}/ /;

my $stderr1 = `file $input_file 2>&1`;

if ($stderr1 =~ /text/){
    print STDERR "input file is SAM\n";
    if (($dis_realign == 0) or ($read_correct == 1)){
        my $order = &judge_sort($input_file);
        if ($order == 1){
            die "Input sam file is not sorted.\n";
        }
    }
    if ($out_sam == 1){
        if ($read_type eq 'PE'){
            my $new_arg = &stat_ins_size('sam');
            system ("$com_path/coval-refine-sam.pl $new_arg $input_file > $out_prefix.sam") if ($out_prefix ne '');
            system ("$com_path/coval-refine-sam.pl $new_arg $input_file") if ($out_prefix eq '');
        }
        else{
            system ("$com_path/coval-refine-sam.pl $arg $input_file > $out_prefix.sam") if ($out_prefix ne '');
            system ("$com_path/coval-refine-sam.pl $arg $input_file") if ($out_prefix eq '');
        }
    }
    else{
        if ((!-e "$reference_file.fai") or (-z "$reference_file.fai")){
            print STDERR "making an index file of the reference file --\n";
            system ("samtools faidx $reference_file");
	    my $file_size = -s "$reference_file.fai";
	    if ($file_size <= 1){
		print STDERR "Adjusting line length in reference file\n";
		system ("rm $reference_file.fai");
		system ("mv $reference_file $reference_file.2");
		system ("$com_path/adjust_line_fasta.pl $reference_file.2 > $reference_file");
		system ("samtools faidx $reference_file");
		system "mv $reference_file $reference_file.adj";
		system "mv $reference_file.2 $reference_file";
	    }
        }
        if ($read_type eq 'PE'){
            my $new_arg = &stat_ins_size('sam');
#            system ("$com_path/coval-refine-sam.pl $new_arg $input_file | samtools view -uSt $reference_file.fai - 2>/dev/null | samtools sort - $out_prefix 2>/dev/null");
	    system ("$com_path/coval-refine-sam.pl $new_arg $input_file | samtools view -bSt $reference_file.fai - > $out_prefix.bam");
        }
        else{
            system ("$com_path/coval-refine-sam.pl $arg $input_file | samtools view -bSt $reference_file.fai -  > $out_prefix.bam");
        }
    }
}
else{
    if (($dis_realign == 0) or ($read_correct == 1)){
        my $stderr2 = `samtools index $input_file 2>&1`;
        if (length $stderr2 > 0){
            die "Input bam file is not sorted.\n";
        }
    }
    if ($out_sam == 1){
        if ($read_type eq 'PE'){
            my $new_arg = &stat_ins_size('bam');
            system ("samtools view $input_file | $com_path/coval-refine-sam.pl $new_arg - > $out_prefix.sam") if ($out_prefix ne '');
            system ("samtools view $input_file | $com_path/coval-refine-sam.pl $new_arg -") if ($out_prefix eq '');
        }
        else{
            system ("samtools view $input_file | $com_path/coval-refine-sam.pl $arg - > $out_prefix.sam") if ($out_prefix ne '');
            system ("samtools view $input_file | $com_path/coval-refine-sam.pl $arg -") if ($out_prefix eq '');
        }
    }
    else{
        if ((!-e "$reference_file.fai") or (-z "$reference_file.fai")){
            print STDERR "making an index file of the reference file --\n";
            system ("samtools faidx $reference_file");
	    my $file_size = -s "$reference_file.fai";
	    if ($file_size <= 1){
		print STDERR "Adjusting line length in reference file\n";
		system ("rm $reference_file.fai");
		system ("mv $reference_file $reference_file.2");
		system ("$com_path/adjust_line_fasta.pl $reference_file.2 > $reference_file");
		system ("samtools faidx $reference_file");
		system "mv $reference_file $reference_file.adj";
		system "mv $reference_file.2 $reference_file";
	    }
        }
        if ($read_type eq 'PE'){
            my $new_arg = &stat_ins_size('bam');
#            system ("samtools view $input_file | $com_path/coval-refine-sam.pl $new_arg - | samtools view -uSt $reference_file.fai - 2>/dev/null | samtools sort - $out_prefix 2>/dev/null");
	    system ("samtools view $input_file | $com_path/coval-refine-sam.pl $new_arg - | samtools view -bSt $reference_file.fai - > $out_prefix.bam");
        }
        else{
            system ("samtools view $input_file | $com_path/coval-refine-sam.pl $arg - | samtools view -bSt $reference_file.fai - > $out_prefix.bam");
        }
    }
}

sub judge_sort {
    my ($file) = @_;
    open (FILE, $file) or die "$!";
    my $pre_chr = '-';
    my $pre_pos = 1;
    my $count_line = 0;
    my $diff_chr = 0;
    while (<FILE>){
        my @line = split (/\s+/, $_);
        my $chr = $line[2];
        my $pos = $line[3];
        $pre_chr = $chr if ($pre_chr eq '-');
        $diff_chr++ if ($pre_chr ne $chr);
        return '1' if (($pre_chr eq $chr) and ($pre_pos > $pos));
        return '1' if ($diff_chr > 500);
        return '0' if ($count_line > 1000);
        $pre_chr = $chr;
        $pre_pos = $pos;
        $count_line++;
    }
    close (FILE);
    return '0';
}

sub stat_ins_size {
    my ($format) = @_;
    my $stat_ins = `cat $input_file | $com_path/stat-insert.pl $arg -` if ($format eq 'sam');
    $stat_ins = `samtools view $input_file | $com_path/stat-insert.pl $arg -` if ($format eq 'bam');
    $ave_insert_size_SD = $stat_ins if (($ave_insert_size_SD eq 'auto') or ($ave_insert_size_SD !~ /^\d+,\d+$/));
    $arg =~ s/\s*-+a\w*\s\w*\s*//;
    my $arg_2 = "$arg -a $ave_insert_size_SD";
    return $arg_2;
}
