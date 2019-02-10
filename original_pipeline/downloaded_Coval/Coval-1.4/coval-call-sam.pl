#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;

# Call SNPs and short indels from a sorted sam file

my $reference_file;		# -r
my $minimum_alt_num = 2;        # -n
my $minimum_alt_freq = 0.8;     # -f
my $minimum_ave_alt_qual = 20;  # -qa
my $minimun_alt_qual = 3;       # -qb
my $maximum_read_depth = 36;    # -m
my $output_prefix = 'out';      # -p
my $minimum_het_alt_num = 2;    # -t
my $base_call_type = 'auto';    # -c
my $help;                       # -h

$reference_file = '/home/kosugis/reference/fk02IDrr_nipp-snpindel.fa';
#$reference_file = '/home/kosugis/reference/IRGSPb5.fa';
#$reference_file = '/home/kosugis/reference/simulation/fk02IDrr_elegans_snpindel.fa';
#$reference_file = '/home/kosugis/reference/arabi.fa';

GetOptions(
    'ref=s' => \$reference_file,
    'num=i' => \$minimum_alt_num,
    'freq=f' => \$minimum_alt_freq,
    'qual_ave|qa=i' => \$minimum_ave_alt_qual,
    'qual_base|qb=i' => \$minimun_alt_qual,
    'maxr=i' => \$maximum_read_depth,
    'tnum=i' => \$minimum_het_alt_num,
    'calltype=s' => \$base_call_type,
    'pref=s' => \$output_prefix,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;
pod2usage(-verbose => 0) unless (-e $reference_file);
    

unless (@ARGV){
print "Your sam file name should be added to the argument.\n";
exit;
}
my $input_file = shift (@ARGV);
if ($input_file eq '-'){
}
elsif (!-f $input_file){
    die "inputfile does not exist.\n";
}
#$output_prefix = shift (@ARGV);

print STDERR "########## covar call-sam options ##########\nnum=$minimum_alt_num freq=$minimum_alt_freq qual_ave=$minimum_ave_alt_qual qual_base=$minimun_alt_qual maxr=$maximum_read_depth calltype=$base_call_type\npref=$output_prefix input=$input_file ref=$reference_file\n";

=head1 SYNOPSIS

coval call-sam [options] <input_sorted_sam_file>
  Output:
   prefix-snp.txt

  Options:
   --ref or -r <STR>	    reference fasta file used for the alignment
   --num or -n <INT>        minimum number of reads supporting non-reference allele [default: 2]
   --freq or -f <FLOAT>     minimum frequency of non-reference allele [default: 0.8]
   --qual_ave or -qa <INT>  minimum averaged base-call quality at a SNP-called site [default: 20]
   --qual_base or -qb <INT> minimum base-call quality of a non-reference base [default: 3]
   --maxr or -m <INT>       maximum read number covering non-reference allele [default: 10000]
   --tnum or -t <INT>       minimum number of reads supprting heterozygous non-reference allele [default: 2]
   --calltype or -c <STR>   quality format of the fastq file used; illumina (Phred+64) or sanger (Phred+33) [default: auto]
   --pref or -p <STR>       prefix of output files [default: out]
   --help or -h             output help message

=cut


my $snp_outfile = "$output_prefix-snp-C.txt";
my $indel_outfile = "$output_prefix-indel-C.txt";

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
	    $count_chr ++;
	}
	else{
	    $chrseq .= uc $line;
	}
    }
    $chr_seq{$ref_chr_name} = $chrseq;
    $chrseq = '';
close (FILE);

my @line;
my $read_number = 0;
my $count_aligned_pos = 0;
my $chrom;
my $pre_chr = '-';
my %snp_pos;
my %non_indel_pos;
my @extracted_var;
my @snp_lines;
my @indel_lines;
my $count_Q = 0;
my $count_Q_filt = 0;
my $count_line_filt = 0;
my $calc_line = 10000;

my $read_length = 0;
my $average_read_depth = 0;
my $coverage = 0;
my $genome_size = 0;

keys %snp_pos = 25000;
keys %non_indel_pos = 25000;

my $base_qual;
my $sanger_type;
my $solexa_type;
open (FILE, $input_file) or die "$input_file: $!";
    while (<FILE>){
	next if ($_ =~ /^\@/);
        my @line = split(/\s+/, $_);
        $base_qual = $line[10];
        last if (length $base_qual >= 1000);
    }
    $sanger_type += $base_qual =~ /[!Ó#%&'\(\)\*,\.\/0123456789:;<=>]/g;
    $solexa_type += $base_qual =~ /LMNOPQRSTUVWXYZ\[\\\]\^_`abcdefgh/g;
close (FILE);
if ($sanger_type > $solexa_type){
    if ($base_call_type eq 'illumina'){
        print STDERR "\n########## Warning! ##########\nYour base quality format appears a sanger type, but the job will be processed with illumina type\n. To change it use -calltype option.\n";
    }
    else{
	$base_call_type = 'sanger';
	print STDERR "Base call quality format: sanger type\n";
    }
}
else{
    if ($base_call_type eq 'sanger'){
        print STDERR "\n########## Warning! ##########\nYour base quality format appears a illumina type, but the job will be processed with sanger type\n. To change it use -calltype option.\n";
    }
    else{
	$base_call_type = 'illumina';
	print STDERR "Base call quality format: illumina type\n";
    }
}

my $count_line = 0;

open (FILE, $input_file) or die "$input_file: $!";
    while (my $line = <FILE>){
	if ($line =~ /^\@/){
	    print $line;
	    next;
	}
        chomp $line;
        @line = split (/\s+/, $line);
	next if ($line[5] =~ /\*/);
	$read_number++;
print STDERR 'Filtered reads: ', $read_number, "\t", 'pos = ', $line[3], "\n" if ($read_number % 1000000 == 0);	
        next if ($line[0] =~ /^\@/);
        next if ($line[5] =~ /H/);
	next if ($line[5] =~ /^\d+S.+S$/);
	next if ($line[5] =~ /^\d+S\d+M\d+[DI]/);
	next if ($line[5] =~ /\d+[DI]\d+M\d+S$/);
        next if ($line[5] =~ /\d+[DI]\d+M\d+[DI]\d+M\d+[DI]/);
        next if ($line[1] >= 1000);          #remove PCR-duplicated reads
                 
        $line[2] =~ /^(\S*)/;
        $chrom = $1;
        $pre_chr = $chrom if ($pre_chr eq '-');
        $read_length = length $line[9];

	if ($line[3] >= $calc_line + 10000){
	    &extract_var ();
	    $calc_line += 10000;
	}
	elsif ($chrom ne $pre_chr){
	    &extract_var ();
	    $calc_line = 10000;
	    %snp_pos = ();
	    %non_indel_pos = ();
	}

	if ($line[5] eq $read_length . 'M'){
	    my $ref = substr ($chr_seq{$chrom}, $line[3] - 1, $read_length);
	    if (length $ref != $read_length){
		next;
	    }
	    if ($line[9] eq $ref){
		&pile_match_pos ($line[3], $read_length);
	    }
	    else{
		&pile_mismatch_pos ($line[9], $ref, $line[3]);
	    }
	}
	elsif ($line[5] =~ /^(\d+)M(\d+)I(\d+)M$/){
	    my $ref1 = substr ($chr_seq{$chrom}, $line[3] - 1, $1);
	    my $ref2 = substr ($chr_seq{$chrom}, $line[3] + $1 - 1, $3);
	    my $read1 = substr ($line[9], 0, $1);
	    my $read2 = substr ($line[9], $1 + $2, $3);
	    next if (length $ref1 != length $read1);
	    if ($read1 eq $ref1){
		&pile_match_pos ($line[3], length $read1);
	    }
	    else{
		&pile_mismatch_pos ($read1, $ref1, $line[3]);
	    }
	    next if (length $ref2 != length $read2);
	    if ($read2 eq $ref2){
		&pile_match_pos ($line[3] + $1, length $read2);
	    }
	    else{
		&pile_mismatch_pos ($read2, $ref2, $line[3] + $1);
	    }
	    my $ins_pos = $line[3] + $1 - 1;
	    my $ins_seq = substr ($line[9], $1, $2);
	    &pile_indel_pos ($ins_pos, $ins_seq, 'I');
	}
	elsif ($line[5] =~ /^(\d+)M(\d+)D(\d+)M$/){
	    my $ref1 = substr ($chr_seq{$chrom}, $line[3] - 1, $1);
	    my $ref2 = substr ($chr_seq{$chrom}, $line[3] + $1 + $2 - 1, $3);
	    my $read1 = substr ($line[9], 0, $1);
	    my $read2 = substr ($line[9], $1, $3);
	    next if (length $ref1 != length $read1);
	    if ($read1 eq $ref1){
		&pile_match_pos ($line[3], length $read1);
	    }
	    else{
		&pile_mismatch_pos ($read1, $ref1, $line[3]);
	    }
	    next if (length $ref2 != length $read2);
	    if ($read2 eq $ref2){
		&pile_match_pos ($line[3] + $1 + $2, length $read2);
	    }
	    else{
		&pile_mismatch_pos ($read2, $ref2, $line[3] + $1 + $2);
	    }
	    my $del_pos = $line[3] + $1 - 1;
	    my $del_seq = substr ($chr_seq{$chrom}, $line[3] + $1 - 1, $2);
	    &pile_indel_pos ($del_pos, $del_seq, 'D');
	}
	elsif ($line[5] =~ /^(\d+)M(\d+)I(\d+)M(\d+)I(\d+)M$/){
	    my $ref1 = substr ($chr_seq{$chrom}, $line[3] - 1, $1);
	    my $ref2 = substr ($chr_seq{$chrom}, $line[3] + $1 - 1, $3);
	    my $ref3 = substr ($chr_seq{$chrom}, $line[3] + $1 + $3 - 1, $5);
	    my $read1 = substr ($line[9], 0, $1);
	    my $read2 = substr ($line[9], $1 + $2, $3);
	    my $read3 = substr ($line[9], $1 + $2 + $3 + $4, $5);
	    next if (length $ref1 != length $read1);
	    if ($read1 eq $ref1){
		&pile_match_pos ($line[3], length $read1);
	    }
	    else{
		&pile_mismatch_pos ($read1, $ref1, $line[3]);
	    }
	    next if (length $ref2 != length $read2);
	    if ($read2 eq $ref2){
		&pile_match_pos ($line[3] + $1, length $read2);
	    }
	    else{
		&pile_mismatch_pos ($read2, $ref2, $line[3] + $1);
	    }
	    next if (length $ref3 != length $read3);
	    if ($read3 eq $ref3){
		&pile_match_pos ($line[3] + $1 + $3, length $read3);
	    }
	    else{
		&pile_mismatch_pos ($read3, $ref3, $line[3] + $1 + $3);
	    }
	    my $ins_pos_1 = $line[3] + $1 - 1;
	    my $ins_seq_1 = substr ($line[9], $1, $2);
	    my $ins_pos_2 = $line[3] + $1 + $3 - 1;
	    my $ins_seq_2 = substr ($line[9], $1 + $2 + $3, $4);
	    &pile_indel_pos ($ins_pos_1, $ins_seq_1, 'I');
	    &pile_indel_pos ($ins_pos_2, $ins_seq_2, 'I');
	}
	elsif ($line[5] =~ /^(\d+)M(\d+)I(\d+)M(\d+)D(\d+)M$/){
	    my $ref1 = substr ($chr_seq{$chrom}, $line[3] - 1, $1);
	    my $ref2 = substr ($chr_seq{$chrom}, $line[3] + $1 - 1, $3);
	    my $ref3 = substr ($chr_seq{$chrom}, $line[3] + $1 + $3 + $4 - 1, $5);
	    my $read1 = substr ($line[9], 0, $1);
	    my $read2 = substr ($line[9], $1 + $2, $3);
	    my $read3 = substr ($line[9], $1 + $2 + $3, $5);
	    next if (length $ref1 != length $read1);
	    if ($read1 eq $ref1){
		&pile_match_pos ($line[3], length $read1);
	    }
	    else{
		&pile_mismatch_pos ($read1, $ref1, $line[3]);
	    }
	    next if (length $ref2 != length $read2);
	    if ($read2 eq $ref2){
		&pile_match_pos ($line[3] + $1, length $read2);
	    }
	    else{
		&pile_mismatch_pos ($read2, $ref2, $line[3] + $1);
	    }
	    next if (length $ref3 != length $read3);
	    if ($read3 eq $ref3){
		&pile_match_pos ($line[3] + $1 + $3 + $4, length $read3);
	    }
	    else{
		&pile_mismatch_pos ($read3, $ref3, $line[3] + $1 + $3 + $4);
	    }
	    my $ins_pos = $line[3] + $1 - 1;
	    my $ins_seq = substr ($line[9], $1, $2);
	    my $del_pos = $line[3] + $1 + $3 - 1;
	    my $del_seq = substr ($chr_seq{$chrom}, $line[3] + $1 + $3 - 1, $4);
	    &pile_indel_pos ($ins_pos, $ins_seq, 'I');
	    &pile_indel_pos ($del_pos, $del_seq, 'D');
	}
	elsif ($line[5] =~ /^(\d+)M(\d+)D(\d+)M(\d+)I(\d+)M$/){
	    my $ref1 = substr ($chr_seq{$chrom}, $line[3] - 1, $1);
	    my $ref2 = substr ($chr_seq{$chrom}, $line[3] + $1 + $2 - 1, $3);
	    my $ref3 = substr ($chr_seq{$chrom}, $line[3] + $1 + $2 + $3 - 1, $5);
	    my $read1 = substr ($line[9], 0, $1);
	    my $read2 = substr ($line[9], $1, $3);
	    my $read3 = substr ($line[9], $1 + $3 + $4, $5);
	    next if (length $ref1 != length $read1);
	    if ($read1 eq $ref1){
		&pile_match_pos ($line[3], length $read1);
	    }
	    else{
		&pile_mismatch_pos ($read1, $ref1, $line[3]);
	    }
	    next if (length $ref2 != length $read2);
	    if ($read2 eq $ref2){
		&pile_match_pos ($line[3] + $1 + $2, length $read2);
	    }
	    else{
		&pile_mismatch_pos ($read2, $ref2, $line[3] + $1 + $2);
	    }
	    next if (length $ref3 != length $read3);
	    if ($read3 eq $ref3){
		&pile_match_pos ($line[3] + $1 + $2 + $3, length $read3);
	    }
	    else{
		&pile_mismatch_pos ($read3, $ref3, $line[3] + $1 + $2 + $3);
	    }
	    my $ins_pos = $line[3] + $1 + $2 + $3 - 1;
	    my $ins_seq = substr ($line[9], $1 + $2, $4);
	    my $del_pos = $line[3] + $1 - 1;
	    my $del_seq = substr ($chr_seq{$chrom}, $line[3] + $1 - 1, $2);
	    &pile_indel_pos ($ins_pos, $ins_seq, 'I');
	    &pile_indel_pos ($del_pos, $del_seq, 'D');
	}
	elsif ($line[5] =~ /^(\d+)M(\d+)D(\d+)M(\d+)D(\d+)M$/){
	    my $ref1 = substr ($chr_seq{$chrom}, $line[3] - 1, $1);
	    my $ref2 = substr ($chr_seq{$chrom}, $line[3] + $1 + $2 - 1, $3);
	    my $ref3 = substr ($chr_seq{$chrom}, $line[3] + $1 + $2 + $3 + $4 - 1, $5);
	    my $read1 = substr ($line[9], 0, $1);
	    my $read2 = substr ($line[9], $1, $3);
	    my $read3 = substr ($line[9], $1 + $3, $5);
	    next if (length $ref1 != length $read1);
	    if ($read1 eq $ref1){
		&pile_match_pos ($line[3], length $read1);
	    }
	    else{
		&pile_mismatch_pos ($read1, $ref1, $line[3]);
	    }
	    next if (length $ref2 != length $read2);
	    if ($read2 eq $ref2){
		&pile_match_pos ($line[3] + $1 + $2, length $read2);
	    }
	    else{
		&pile_mismatch_pos ($read2, $ref2, $line[3] + $1 + $2);
	    }
	    next if (length $ref3 != length $read3);
	    if ($read3 eq $ref3){
		&pile_match_pos ($line[3] + $1 + $2 + $3 + $4, length $read3);
	    }
	    else{
		&pile_mismatch_pos ($read3, $ref3, $line[3] + $1 + $2 + $3 + $4);
	    }
	    my $del_pos_1 = $line[3] + $1 - 1;
	    my $del_seq_1 = substr ($chr_seq{$chrom}, $line[3] + $1 - 1, $2);
	    my $del_pos_2 = $line[3] + $1 + $2 + $3 - 1;
	    my $del_seq_2 = substr ($chr_seq{$chrom}, $line[3] + $1 + $2 + $3 - 1, $4);
	    &pile_indel_pos ($del_pos_1, $del_seq_1, 'D');
	    &pile_indel_pos ($del_pos_2, $del_seq_2, 'D');
	}
	elsif ($line[5] =~ /^(\d+)S(\d+)M$/){
	    my $ref = substr ($chr_seq{$chrom}, $line[3] - 1, $2);
	    my $read = substr ($line[9], $1, $2);
	    if (length $ref != length $read){
		next;
	    }
	    if ($read eq $ref){
		&pile_match_pos ($line[3], length $read);
	    }
	    else{
		&pile_mismatch_pos ($read, $ref, $line[3]);
	    }
	}
	elsif ($line[5] =~ /^(\d+)M(\d+)S$/){
	    my $ref = substr ($chr_seq{$chrom}, $line[3] - 1, $1);
	    my $read = substr ($line[9], 0, $1);
	    if (length $ref != length $read){
		next;
	    }
	    if ($read eq $ref){
		&pile_match_pos ($line[3], length $read);
	    }
	    else{
		&pile_mismatch_pos ($read, $ref, $line[3]);
	    }
	}
	elsif ($line[5] =~ /^(\d+)I(\d+)M$/){
	    my $ref = substr ($chr_seq{$chrom}, $line[3] - 1, $1);
	    my $read = substr ($line[9], $1, $2);
	    next if (length $ref != length $read);
	    if ($read eq $ref){
		&pile_match_pos ($line[3], length $read);
	    }
	    else{
		&pile_mismatch_pos ($read, $ref, $line[3]);
	    }
	    my $ins_pos = $line[3] - $1 - 1;
	    my $ins_seq = substr ($line[9], 0, $1);
	    &pile_ins_pos_T ($ins_pos, $ins_seq, 'I');
	}
	elsif ($line[5] =~ /^(\d+)I(\d+)M(\d+)I(\d+)M$/){
	    my $ref1 = substr ($chr_seq{$chrom}, $line[3] - 1, $2);
	    my $ref2 = substr ($chr_seq{$chrom}, $line[3] + $2 - 1, $4);
	    my $read1 = substr ($line[9], $1, $2);
	    my $read2 = substr ($line[9], $1 + $2 + $3, $4);
	    next if (length $ref1 != length $read1);
	    if ($read1 eq $ref1){
		&pile_match_pos ($line[3], length $read1);
	    }
	    else{
		&pile_mismatch_pos ($read1, $ref1, $line[3]);
	    }
	    next if (length $ref2 != length $read2);
	    if ($read2 eq $ref2){
		&pile_match_pos ($line[3] + $2, length $read2);
	    }
	    else{
		&pile_mismatch_pos ($read2, $ref2, $line[3] + $2);
	    }
	    my $ins_pos_1 = $line[3] - $1 - 1;
	    my $ins_seq_1 = substr ($line[9], 0, $1);
	    my $ins_pos_2 = $line[3] + $2 - 1;
	    my $ins_seq_2 = substr ($line[9], $1 + $2, $3);
	    &pile_ins_pos_T ($ins_pos_1, $ins_seq_1, 'I');
	    &pile_indel_pos ($ins_pos_2, $ins_seq_2, 'I');
	}
	elsif ($line[5] =~ /^(\d+)M(\d+)I(\d+)M(\d+)I$/){
	    my $ref1 = substr ($chr_seq{$chrom}, $line[3] - 1, $1);
	    my $ref2 = substr ($chr_seq{$chrom}, $line[3] + $1 - 1, $3);
	    my $read1 = substr ($line[9], 0, $1);
	    my $read2 = substr ($line[9], $1 + $2, $3);
	    next if (length $ref1 != length $read1);
	    if ($read1 eq $ref1){
		&pile_match_pos ($line[3], length $read1);
	    }
	    else{
		&pile_mismatch_pos ($read1, $ref1, $line[3]);
	    }
	    next if (length $ref2 != length $read2);
	    if ($read2 eq $ref2){
		&pile_match_pos ($line[3] + $1, length $read2);
	    }
	    else{
		&pile_mismatch_pos ($read2, $ref2, $line[3] + $1);
	    }
	    my $ins_pos_1 = $line[3] + $1 - 1;
	    my $ins_seq_1 = substr ($line[9], $1, $2);
	    my $ins_pos_2 = $line[3] + $1 + $3 - 1;
	    my $ins_seq_2 = substr ($line[9], $1 + $2 + $3, $4);
	    &pile_indel_pos ($ins_pos_1, $ins_seq_1, 'I');
	    &pile_indel_pos ($ins_pos_2, $ins_seq_2, 'I');
	}
	elsif ($line[5] =~ /^(\d+)I(\d+)M(\d+)D(\d+)M$/){
	    my $ref1 = substr ($chr_seq{$chrom}, $line[3] - 1, $2);
	    my $ref2 = substr ($chr_seq{$chrom}, $line[3] + $2 + $3 - 1, $4);
	    my $read1 = substr ($line[9], $1, $2);
	    my $read2 = substr ($line[9], $1 + $2, $4);
	    next if (length $ref1 != length $read1);
	    if ($read1 eq $ref1){
		&pile_match_pos ($line[3], length $read1);
	    }
	    else{
		&pile_mismatch_pos ($read1, $ref1, $line[3]);
	    }
	    next if (length $ref2 != length $read2);
	    if ($read2 eq $ref2){
		&pile_match_pos ($line[3] + $2 + $3, length $read2);
	    }
	    else{
		&pile_mismatch_pos ($read2, $ref2, $line[3] + $2 + $3);
	    }
	    my $ins_pos = $line[3] - $1 - 1;
	    my $ins_seq = substr ($line[9], 0, $1);
	    my $del_pos = $line[3] + $2 - 1;
	    my $del_seq = substr ($chr_seq{$chrom}, $line[3] + $2 - 1, $3);
	    &pile_ins_pos_T ($ins_pos, $ins_seq, 'I');
	    &pile_indel_pos ($del_pos, $del_seq, 'D');
	}
	elsif ($line[5] =~ /^(\d+)M(\d+)D(\d+)M(\d+)I$/){
	    my $ref1 = substr ($chr_seq{$chrom}, $line[3] - 1, $1);
	    my $ref2 = substr ($chr_seq{$chrom}, $line[3] + $1 + $2 - 1, $3);
	    my $read1 = substr ($line[9], 0, $1);
	    my $read2 = substr ($line[9], $1, $3);
	    next if (length $ref1 != length $read1);
	    if ($read1 eq $ref1){
		&pile_match_pos ($line[3], length $read1);
	    }
	    else{
		&pile_mismatch_pos ($read1, $ref1, $line[3]);
	    }
	    next if (length $ref2 != length $read2);
	    if ($read2 eq $ref2){
		&pile_match_pos ($line[3] + $1 + $2, length $read2);
	    }
	    else{
		&pile_mismatch_pos ($read2, $ref2, $line[3] + $1 + $2);
	    }
	    my $ins_pos = $line[3] + $1 + $2 + $3 - 1;
	    my $ins_seq = substr ($line[9], $1 + $3, $4);
	    my $del_pos = $line[3] + $1 - 1;
	    my $del_seq = substr ($chr_seq{$chrom}, $line[3] + $1 - 1, $2);
	    &pile_indel_pos ($ins_pos, $ins_seq, 'I');
	    &pile_indel_pos ($del_pos, $del_seq, 'D');
	}
        $pre_chr = $chrom;
	@line = ();
    }
    &extract_var ();
close (FILE);

foreach my $key (keys %chr_seq){
    $genome_size += length $chr_seq{$key};
}

$average_read_depth = int ($read_number * $read_length / $count_aligned_pos * 10) / 10;
$coverage = int ($count_aligned_pos / $genome_size * 1000) / 10;


open (OUTFILE, ">$snp_outfile");
foreach my $line (@snp_lines){
    print OUTFILE $line, "\n";
}
close (OUTFILE);

open (OUTFILE, ">$indel_outfile");
foreach my $line (@indel_lines){
    print OUTFILE $line, "\n";
}

print STDERR "read number = $read_number\n";
print STDERR "genome size = $genome_size\n";
print STDERR "total genome positions covered by reads = $count_aligned_pos\n";
print STDERR "genome coverage = $coverage %\n";
print STDERR "average read depth = $average_read_depth", ' X', "\n";


print STDERR "average Q for all SNPs/indels with the minimum number/frequency of SNPs/indels = ", $count_Q / $count_line, "\n";
print STDERR "average Q for selected SNPs/indels = ", $count_Q_filt / $count_line_filt, "\n";
print STDERR "total selected SNPs/indels = ", $count_line_filt, "\n";


####################################################################

sub pile_mismatch_pos{
    my ($read, $ref, $pos) = @_;
    for (my $i = 0; $i < length $read; $i++){
	my $read_i = substr ($read, $i, 1);
	my $ref_i = substr ($ref, $i, 1);
	next if ($ref_i eq 'N');
	my $snp_base;
	my $position = $pos + $i;
        my $mapQ = $line[4];
        my $pos09d = sprintf ("%09d", $position);
	if ($read_i ne $ref_i){
	    $snp_base = $read_i;
	    my $qualities = substr ($line[10], $i, 1);
	    if (!exists $snp_pos{$pos09d}){ 
		my $new_line = '|m' . $mapQ . '|a' . $snp_base . '|q' . $qualities;
		$snp_pos{$pos09d} = $new_line;
	    }
	    else{
		$snp_pos{$pos09d} .= '|m' . $mapQ . '|a' . $snp_base . '|q' . $qualities;
	    }
	}
        else{
            $snp_base = ':';
	    if (!exists $snp_pos{$pos09d}){ 
		my $new_line = '|m' . $mapQ . '|a' . $snp_base;
		$snp_pos{$pos09d} = $new_line;
	    }
	    else{
		$snp_pos{$pos09d} .= '|m' . $mapQ . '|a' . $snp_base;
	    }
        }
	$non_indel_pos{$position} ++ if ($i < length ($read) - 1);
    }
}

sub pile_match_pos{
    my ($pos, $read_length) = @_;
    for (my $i = 0; $i < $read_length; $i++){
	my $snp_base = ':';
	my $position = $pos + $i;
        my $mapQ = $line[4];
        my $pos09d = sprintf ("%09d", $position);
	if (!exists $snp_pos{$pos09d}){ 
	    my $new_line = '|m' . $mapQ . '|a' . $snp_base;
	    $snp_pos{$pos09d} = $new_line;
	}
	else{
	    $snp_pos{$pos09d} .= '|m' . $mapQ . '|a' . $snp_base;
	}
	$non_indel_pos{$position} ++ if ($i < $read_length - 1);
    }
}

sub pile_indel_pos{
    my ($pos, $seq, $indel) = @_;
    my $pos09d = sprintf ("%09d", $pos);
    my $indel_pos = $pos09d . '=' . $indel . '=1';
    my $qualities = substr ($line[10], $pos - $line[3], 1);
    my $mapQ = $line[4];
    if (!exists $snp_pos{$indel_pos}){
	my $new_line = 'i' . $seq . '|m' . $mapQ . '|a1' . '|q' . $qualities;
        $snp_pos{$indel_pos} = $new_line;
    }
    else{
	$snp_pos{$indel_pos} =~ /^i(\w+)\|/;
	if ($1 =~ /$seq/){
            $snp_pos{$indel_pos} .= '|m' . $mapQ . '|a1' . '|q' . $qualities;
        }
        else{
            my $indel_pos_2 = $pos09d . '=' . $indel . '=2';
            if (!exists $snp_pos{$indel_pos_2}){
                my $new_line = 'i' . $seq . '|m' . $mapQ . '|a1' . '|q' . $qualities;
		$snp_pos{$indel_pos_2} = $new_line;
            }
            else{
                $snp_pos{$indel_pos_2} .= '|m' . $mapQ . '|a1' . '|q' . $qualities;
            }
        }
    }
}

sub pile_ins_pos_T{
    my ($pos, $seq, $indel) = @_;
    for (my $i = $pos - 10; $i <= $pos; $i++){
	my $pos09d = sprintf ("%09d", $i);
        my $ins_pos = $pos09d . '=' . $indel . '=1';
        if (exists $snp_pos{$ins_pos}){
	    my $mapQ = $line[4];
	    my $qualities = substr ($line[10], $i - $line[3], 1);
	    $snp_pos{$ins_pos} .= '|m' . $mapQ . '|a1' . '|q' . $qualities;
	}
    }
}

sub calc_qual{
    my ($qual, $vartype) = @_;
    my @snp_qual = split(//, $qual);
    my $Q_base = 0;
    my $count_lowQ = 0;
    my $del_base_pos = '-';
    if ($base_call_type eq 'illumina'){
	for (my $i = 0; $i < @snp_qual; $i++){
	    if ((ord($snp_qual[$i]) - 64 < $minimun_alt_qual) and ($vartype eq 'SNP')){
		if ($del_base_pos eq '-'){
		    $del_base_pos = $i;
		}
		else{
		    $del_base_pos .= '-' . $i;
		}
		next;
	    }
	    if ($snp_qual[$i] =~ /[abcdefgh`]/){
		$Q_base += 40;
	    }
	    elsif ($snp_qual[$i] =~ /B/){
		$Q_base += 0;
	    }
	    else {
		$Q_base += ord($snp_qual[$i]) - 64;
	    }
	}
    }
    elsif ($base_call_type eq 'sanger'){
	for (my $i = 0; $i < @snp_qual; $i++){
	    if ((ord($snp_qual[$i]) - 33 < $minimun_alt_qual) and ($vartype eq 'SNP')){
		if ($del_base_pos eq '-'){
		    $del_base_pos = $i;
		}
		else{
		    $del_base_pos .= '-' . $i;
		}
		next;
	    }
	    if ($snp_qual[$i] =~ /[ABCDEFGHI>\?@]/){
		$Q_base += 40;
	    }
	    else {
		$Q_base += ord($snp_qual[$i]) - 33;
	    }
	}
    }
    return ($Q_base, $del_base_pos);
}

sub extract_var{
#print STDERR scalar (keys %snp_pos), "\n";
    return if (scalar keys %snp_pos < 1);
    foreach my $key (sort keys %snp_pos){
	if ($key =~ /=/){
	    my ($pos09d, $indel, $indel_2) = split (/=/, $key);
	    my $pos = $1 if ($pos09d =~ /^0*(\d+)/);
	    last if (($pos >= $calc_line) and ($chrom eq $pre_chr));
	    my @indel_items = split (/\|/, $snp_pos{$key});
	    delete $snp_pos{$key};
	    my $indel_seq = '';
	    my $indel_num = 0;
	    my $sum_mapQ = 0;
	    my $qual = '';
	    foreach (@indel_items){
		$indel_seq = $1 if ($_ =~ /^i(.+)/);
		$indel_num ++ if ($_ eq 'a1');
		$sum_mapQ += $1 if ($_ =~ /^m(.+)/);
		$qual .= $1 if ($_ =~ /^q(.+)/);
	    }
	    next if ($indel_num < $minimum_alt_num);
	    my $non_indel_read_num = 0;
	    $non_indel_read_num = $non_indel_pos{$pos} if (exists $non_indel_pos{$pos});
	    my $read_num = $indel_num + $non_indel_read_num;
	    my $allele_freq = int ($indel_num / $read_num * 100) / 100;
	    next if ($allele_freq < $minimum_alt_freq);
	    my ($Q_sum, $del_pos) = &calc_qual ($qual, 'indel');
	    my $Q_ave = int ($Q_sum / $indel_num * 10) / 10;
	    my $indel_mapQ = int ($sum_mapQ / $read_num * 10) / 10;
	    my $indel_line = $chrom . "\t" . $pos . "\t" . '*' . "\t" . $indel_seq . "\t" . $indel_num . "\t" . $read_num . "\t" . $allele_freq. "\t" . $indel_mapQ . "\t" . $Q_ave if ($indel eq 'I');
	    $indel_line = $chrom . "\t" . $pos . "\t" . $indel_seq . "\t" . '*' . "\t" . $indel_num . "\t" . $read_num . "\t" . $allele_freq. "\t" . $indel_mapQ . "\t" . $Q_ave if ($indel eq 'D');
	    push (@extracted_var, $indel_line);
	}
	else{
	    last if (($key >= $calc_line) and ($chrom eq $pre_chr));
	    $count_aligned_pos ++;
	    if ($snp_pos{$key} !~ /\|a[A-Z]/){
		delete $snp_pos{$key};
		next;
	    }
	    my @snp_items = split (/\|/, $snp_pos{$key});
	    my $base = '';
	    my $sum_mapQ = 0;
	    my $qual = '';
	    my $pos = $1 if ($key =~ /^0*(\d+)/);
	    delete $snp_pos{$key};
	    
	    foreach (@snp_items){
		$base .= $1 if ($_ =~ /^a(.)/);
		$sum_mapQ += $1 if ($_ =~ /^m(.+)/);
		$qual .= $1 if ($_ =~ /^q(.)/);
	    }
	    next if (length $base < $minimum_alt_num);
	    my ($Q_sum, $del_pos) = &calc_qual ($qual, 'SNP');
	    if ($del_pos ne '-'){
		my @del_pos = split (/-/, $del_pos);
		foreach my $delpos (sort {$b <=> $a} @del_pos){
		    substr ($base, $delpos, 1, '');
		}
	    }
	    my $ref_snp = $base;
	    my $snp_C = 0;
	    my $snp_A = 0;
	    my $snp_G = 0;
	    my $snp_T = 0;
	    my $snp_N = 0;
	    $snp_C = $ref_snp =~ s/C/C/g;
	    $snp_A = $ref_snp =~ s/A/A/g;
	    $snp_G = $ref_snp =~ s/G/G/g;
	    $snp_T = $ref_snp =~ s/T/T/g;
	    $snp_N = $ref_snp =~ s/N/N/g;
	    my %snp = ('C' => $snp_C, 'A' => $snp_A, 'G' => $snp_G, 'T' => $snp_T);
	    my @snp_base;
	    my @snp_num;
	    foreach (sort {$snp{$a} <=> $snp{$b}} keys %snp){
		push (@snp_base, $_);
		push (@snp_num, $snp{$_});
	    }
	    my $read_num = length ($base);
	    next if ($read_num < 1);
	    my $snp_base = $snp_base[3];
	    my $snp_num = $snp_num[3];
	    my $allele_freq = int ($snp_num[3] / $read_num * 100) / 100;
	    next if ($snp_num[3] < $minimum_alt_num);
	    next if ($allele_freq < $minimum_alt_freq);
	    
	    my $Q_ave = int ($Q_sum / $snp_num * 10) / 10;
	    $count_Q += $Q_ave;
	    my $ave_mapQ = int ($sum_mapQ / $read_num * 10) / 10;
	    my $ref = substr ($chr_seq{$chrom}, $pos - 1, 1);
	    my $snp_line = $chrom . "\t" . $pos . "\t" . $ref . "\t" . $snp_base . "\t" . $snp_num . "\t" . $read_num . "\t" . $allele_freq . "\t" . $ave_mapQ . "\t" . $Q_ave;
	    push (@extracted_var, $snp_line);
	}
    }
    &filter_var ();
    foreach my $key_2 (sort {$a <=> $b} keys %non_indel_pos){
	delete $non_indel_pos{$key_2} if ($key_2 <= $calc_line);
	last if ($key_2 >= $calc_line);
    }
}

sub filter_var{
#print STDERR 'extracted_var = ', scalar @extracted_var, "\n";
    my @last_var_line = ();
    my $indel_flag = 0;
    my $pre_indel_pos = 0;
    my $pre_indel_alt_num = 0;
    my $pre1_pos = 0;
    my $pre2_pos = 0;
    my $pre3_pos = 0;
    my $IorD = "";
    my $del_length = 0;
    foreach my $line (@extracted_var){  
	my ($chr, $pos, $ref, $alt, $alt_num, $read_num, $freq, $mapQ, $baseQ) = split(/\s+/, $line);
	next if ($alt_num < 1);
	$count_line ++;
	if (($ref eq '*') or ($alt eq '*')){
	    if (($alt_num >= $minimum_alt_num) and ($freq >=  $minimum_alt_freq) and ($read_num <= $maximum_read_depth)){
		$IorD = 'Ins' if ($ref eq '*');
		$IorD = 'Del' if ($alt eq '*');
		$del_length = length $ref;
		$indel_flag = 1;
		if ($pos - $pre_indel_pos <= 2){
		    if ($alt >= $pre_indel_alt_num){
			pop (@indel_lines);
		    }
		    else{
			next;
		    }
		}
		if ($pos - $pre3_pos <= 3){		# remove snp_lines if their positions are 1~3 bp around indel positions
		    pop (@snp_lines);
		    pop (@snp_lines);
		    pop (@snp_lines);
		}
		elsif ($pos - $pre2_pos <= 2){
		    pop (@snp_lines);
		    pop (@snp_lines);
		}
		elsif ($pos - $pre1_pos <= 1){
		    pop (@snp_lines);
		}
		my $new_line = $chr . "\t" . $pos . "\t" . $ref . "\t" . $alt . "\t" . $alt_num . "\t" . $read_num . "\t" . $freq . "\t" . $mapQ . "\t" . $baseQ;
		if (($pos <= $calc_line) and ($pos >= $calc_line - 2)){
		    push (@last_var_line, $new_line);
		}
		else{
		    push (@indel_lines, $new_line);
		}
		$pre_indel_pos = $pos;
		$pre_indel_alt_num = $alt_num;
	    }
	    next;
	}
	elsif ($indel_flag == 1){								# remove snp_lines if their positions are 1~3 bp around indel positions
	    if ($IorD eq 'Ins'){
		if ($pos - $pre_indel_pos <= 3){			# when the current position is within 3 bp downstream of a reliable insersion position     
		    next;
		}
		else{
		    $indel_flag = 0;
		}
	    }
	    elsif ($IorD eq 'Del'){
		if ($pos - $pre_indel_pos <= $del_length){		# when the current position is within a reliable deletion    
		    next;
		}
		elsif ($pos <= $pre_indel_pos + $del_length + 3){	# when the current position is within 3 bp downstream of a reliable deletion position
		    next;
		}
		else{
		    $indel_flag = 0;
		}
	    }
	}
     
	if ($read_num <= $maximum_read_depth){
	    if (($alt_num >=  $minimum_alt_num) and ($freq >= $minimum_alt_freq) and ($baseQ >= $minimum_ave_alt_qual)){	
		my $new_line = $chr . "\t" . $pos . "\t" . $ref . "\t" . $alt . "\t" . $alt_num . "\t" . $read_num . "\t" . $freq . "\t" . $mapQ . "\t" . $baseQ;
		if ($minimum_alt_freq >= 0.8){
		    if (($pos <= $calc_line) and ($pos >= $calc_line - 3)){
			push (@last_var_line, $new_line);
		    }
		    else{
			push (@snp_lines, $new_line);
		    }
		}
		else{
		    if ($freq >= 0.9){
			if (($pos <= $calc_line) and ($pos >= $calc_line - 3)){
			    push (@last_var_line, $new_line);
			}
			else{
			    push (@snp_lines, $new_line);
			}
		    }
		    elsif (($freq < 0.9) and ($alt_num >= $minimum_het_alt_num)){
			if (($pos <= $calc_line) and ($pos >= $calc_line - 3)){
			    push (@last_var_line, $new_line);
			}
			else{
			    push (@snp_lines, $new_line);
			}
		    }
		}
		
		$pre3_pos = $pre2_pos;
		$pre2_pos = $pre1_pos;
		$pre1_pos = $pos;                   
		$count_Q_filt += $baseQ;
		$count_line_filt++;
	    } 
	}
    }
    @extracted_var = ();
    @extracted_var = @last_var_line;
    close (OUTFILE);
}
