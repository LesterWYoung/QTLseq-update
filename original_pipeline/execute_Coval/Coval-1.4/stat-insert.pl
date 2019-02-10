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
my $out_prefix = 'coval-filter'; # -p
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
    $limited_mismatch_num_pairs = $limited_mismatch_num if ($limited_mismatch_num_pairs == $limited_mismatch_num * 2);
}

unless (@ARGV){
print "Your sam/bam file name should be added to the argument.\n";
exit;
}
my $input_file = shift (@ARGV);
if ($input_file eq '-'){
}
elsif (!-f $input_file){
    die "Input file does not exist: $!\n";
}

die "Reference file is not found: $!\n" unless (-f $reference_file);
die "Prefix of output is not specified: $!\n" if ($out_prefix eq '');
die "Read type is not properly specified: $!\n" unless (uc $read_type eq 'PE') or (uc $read_type eq 'SE');

$include_disc_read = 1 if ($reference_type eq 'DRAFT') or ($reference_type eq 'draft');

my $read_number = 0;
my $read_length;
my @read_dist;
my $ave_dist = 0;
my $SD_dist = 0;
my $concordant_read = 0;
my $discordant_read = 0;

if ($read_type eq 'PE'){
    open (FILE, $input_file) or die "$input_file: $!";
    while (my $line = <FILE>){
	next if ($line =~ /^\@/);
        chomp $line;
        my @line = split (/\s+/, $line);
	next if (@line < 11);
	$read_length = length $line[9];
	next if ($line[5] =~ /\*/);
	$read_number++;
	if (($line[8] > 0) and (($line[1] == 73) or ($line[1] == 99) or ($line[1] == 137) or ($line[1] == 153))){
	    $concordant_read ++;
	}
	elsif (($line[1] == 81) or ($line[1] == 97) or ($line[1] == 145) or ($line[1] == 161)){
	    $discordant_read ++;
	}
	if ($read_number >= 800){
	    if (($soap_aligner == 0) and ($long_insert_read == 0)){
		if ((($line[1] == 99) or ($line[1] == 163)) and ($line[4] >= 20) and ($line[5] eq $read_length . 'M')){
		    push (@read_dist, $line[8]);
		}
	    }
	    else{
		if ((($line[1] == 97) or ($line[1] == 161) or ($line[1] == 99) or ($line[1] == 163)) and ($line[4] >= 20) and ($line[5] eq $read_length . 'M')){
		    push (@read_dist, $line[8]);
		}
	    }
	}
	last if ($concordant_read == 5000);
    }
    close (FILE);
}


if (@read_dist > 1){
    my $ave_dist_SD = &ave_SD (@read_dist);
    ($ave_dist, $SD_dist) = split (/=/, $ave_dist_SD);
}
#	    print STDERR 'Average read distance = ', $ave_dist, ' +/- ', $SD_dist, "\n";

my $ins_SD = $ave_dist . ',' . $SD_dist;
print $ins_SD;


sub ave_SD {
    my @insert_size = @_;
    my $sum_insert = 0;
    foreach (@insert_size){
	$sum_insert += $_;
    }
    my $average = $sum_insert / @insert_size;
    my $sq_total = 0;
    foreach (@insert_size){
	$sq_total += ($average - $_) ** 2;
    }
    my $SD = int (($sq_total / @insert_size) ** 0.5 + 0.5);
    my $average2 = int ($average + 0.5);
    my $average_SD = $average2 . '=' . $SD;
    return $average_SD;
}
