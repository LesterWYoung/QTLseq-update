#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;

# Making a simulated reference genome where a user-defined rate of aritificial SNPs and 1- to 6-bp indels are randomely introduced.
# The SNP bases and the indel frequency depending on the length are automatically selected, based on a naturally occuring SNP compoition and indel frequency.

# @ARGV: arguments: 	-r: a fasta file containing chromosome sequences,
#			-p: a prefix for output files (1. a fasta file containing indel/snp-substituted genome, 2. an indel file containing the incorporated snps/indels)
#			-s: mutation rate for SNP [0.001]
#			-i: mutation rate for indel [0.0001]

my $reference;			# -r
my $output_prefix = 'out';	# -p
my $snp_rate = 0.001;		# -s
my $indel_rate = 0.0001;	# -i
my $help;			# -h

GetOptions(
    'ref=s' => \$reference,
    'pref=s' => \$output_prefix,
    'snp=f' => \$snp_rate,
    'indel=f' => \$indel_rate,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;
pod2usage(-verbose => 0) unless (-e $reference);


print STDERR "########## covar simulate options ##########\nsnp=$snp_rate indel=$indel_rate pref=$output_prefix ref=$reference\n";

=head1 SYNOPSIS

coval simulate [options] -r <reference file> -p <prefix of output files>
  Output:
   prefix.fa		simulated reference fasta file
   prefix.snp		text file showing the reference positions and bases of introduced SNPs and indels

  Options:
   --ref or -r <STR>        reference fasta file
   --pref or -p <STR>       prefix of outputfiles [default: out]
   --snp or -s <FLOAT>      genomic mutation rate for SNPs [default: 0.01]
   --indel or -i <FLOAT>    genomic mutation rate for indels [default: 0.001]
   --help or -h             output help message

=cut


my $out_fasta = $output_prefix . '.' . 'fa';
my $out_indel = $output_prefix . '.' . 'snp';

my $genome_size = 0;
my %chr_seq;
my %chr_length;
my %chr_rate;


open (FILE, $reference) or die "$!";
    my $ref_chr;
    my $pre_chrom = '-';
    my $ref_chr_frac;
    my $count_1M = 1;
    my $chrseq = '';
    my $seq_add = '';
    my $count_line = 0;
    while (my $line = <FILE>){
        chomp $line;
	$count_line ++;
	if ($line =~ /^>(\S*)/){
	    $ref_chr = $1;
	    if ($pre_chrom ne '-'){
		my $count_1M_9d = sprintf ("%09d", $count_1M);
		$ref_chr_frac = $pre_chrom . '=' . $count_1M_9d;
		$chr_seq{$ref_chr_frac} = $chrseq;
		$genome_size += length $chrseq;
		$chrseq = '';
		$count_1M = 1;
	    }
	    $pre_chrom = $ref_chr;
	}
	else{
	    $chrseq .= uc $line;
	}
	if (length $chrseq >= 1000000){
	    $seq_add = substr ($chrseq, 1000000, length ($chrseq) - 1000000, '') if (length $chrseq > 1000000);
	    my $count_1M_9d = sprintf ("%09d", $count_1M);
	    $ref_chr_frac = $ref_chr . '=' . $count_1M_9d;
	    $chr_seq{$ref_chr_frac} = $chrseq;
	    $genome_size += length $chrseq;
	    $chrseq = $seq_add;
	    $seq_add = '';
	    $count_1M += 1000000;
	}
    }
    my $count_1M_9d = sprintf ("%09d", $count_1M);
    $ref_chr_frac = $ref_chr . '=' . $count_1M_9d;
    $chr_seq{$ref_chr_frac} = $chrseq;
    $genome_size += length $chrseq;
    $chrseq = '';
close (FILE);


foreach my $key (keys %chr_seq){
	my ($chr, $chr_frac) = split (/=/, $key);
	$chr_length{$chr} += length $chr_seq{$key};
	$chr_rate{$chr} = 0;
}

foreach (keys %chr_rate){
	$chr_rate{$_} = int (($chr_length{$_} / $genome_size) * 1000) / 1000;
}



my %indel;
my %indel_allpos;
my %indel_allpos_2;
my @indel_info;
my $count_snp = 0;
my $count_snp_2 = 0;
my $count_var = 0;
my $count_overlength = 0;

my $snp_no = $snp_rate * $genome_size;
my $indel_no = $indel_rate * $genome_size;
# the expected distribution of indel 1-6 bp is derived from the japanese human genome (Nature Genet 42, 931 2010 Suppl.Fig.13)
my $indel_1b = $indel_no * 0.5 * 0.66;
my $indel_2b = $indel_no * 0.5 * 0.17;
my $indel_3b = $indel_no * 0.5 * 0.07;
my $indel_4b = $indel_no * 0.5 * 0.07;
my $indel_5b = $indel_no * 0.5 * 0.02;
my $indel_6b = $indel_no * 0.5 * 0.01;

my @nuc = ('A', 'C', 'G', 'T');
my @snp_nuc_A = ('G', 'G', 'G', 'G', 'C', 'T');
my @snp_nuc_G = ('A', 'A', 'A', 'A', 'C', 'T');
my @snp_nuc_C = ('T', 'T', 'T', 'T', 'A', 'G');
my @snp_nuc_T = ('C', 'C', 'C', 'C', 'A', 'G');

&make_indel ($indel_1b, $indel_2b, $indel_3b, $indel_4b, $indel_5b, $indel_6b);

&make_snp ($snp_no);

&subst_genome ();


print STDERR 'total SNPs = ', $count_snp, "\n";
print STDERR 'susbstituted SNPs = ', $count_snp_2, "\n";
print STDERR 'total SNPs&Indels = ', $count_var, "\n";
print STDERR 'number of overlength = ', $count_overlength, "\n";


my $prechrom = '-';
my $chrom;
my $chr_frac;
my $chr_seq = '';
open (FILE1, "> $out_fasta") or die "$!";
foreach my $key (sort keys %chr_seq){
    ($chrom, $chr_frac) = split (/=/, $key);
    $prechrom = $chrom if ($prechrom eq '-');
    if ($prechrom eq $chrom){
	$chr_seq .= $chr_seq{$key};
    }
    else{
	print FILE1 '>', "$prechrom\n";
	while (length ($chr_seq) > 0){
	    my $seq_50 = substr ($chr_seq, 0, 50, '');
	    print FILE1 $seq_50, "\n";
	}
	$chr_seq = $chr_seq{$key};
    }
    $prechrom = $chrom;
}
print FILE1 '>', "$chrom\n";
while (length ($chr_seq) > 0){
    my $seq_50 = substr ($chr_seq, 0, 50, '');
    print FILE1 $seq_50, "\n";
}
close (FILE1);


open (FILE2, "> $out_indel") or die "$out_indel: $!";
print FILE2 'Chr', "\t", 'Pos', "\t", 'Var', "\t", 'Ref', "\t", 'Alt', "\n";
foreach (@indel_info){
	print FILE2 $_, "\n";	
}
close (FILE2);


###########################################################################################################

sub make_indel {
	my @indel_num = @_;
	foreach my $chrno (sort keys %chr_rate){
		for (my $delno = 1; $delno < 7; $delno++){
			for (my $i = 0; $i < int ($indel_num[$delno - 1] * $chr_rate{$chrno} + 0.5); $i++){
				my $pos = int (rand ($chr_length{$chrno})) + 1;
				my $pos08d = sprintf ("%09d", $pos);
				my $key_pos = $chrno . '=' . $pos08d;
				if (exists $indel_allpos{$key_pos}){
					$i = $i - 1;
					next;
				}
				else{
					$indel{$key_pos} = 'D-' . $delno;
					for (my $j = -1; $j <= $delno; $j++){
						my $pos08d = sprintf ("%09d", $pos + $j);
						my $key_pos = $chrno . '=' . $pos08d; 
						$indel_allpos{$key_pos} = 1;
					}
				}
			}
			for (my $i = 0; $i < int ($indel_num[$delno - 1] * $chr_rate{$chrno} + 0.5); $i++){
				my $pos = int (rand ($chr_length{$chrno})) + 1;
				my $pos08d = sprintf ("%09d", $pos);
				my $key_pos = $chrno . '=' . $pos08d;
				my $ins_base = '';
				for (my $j = 0; $j < $delno; $j++){
					$ins_base .= $nuc[(rand (4))];
				}
				if (exists $indel_allpos{$key_pos}){
					$i = $i - 1;
					next;
				}
				else{
					$indel{$key_pos} = 'I' . '-' . $ins_base;
					for (my $j = -1; $j <= $delno; $j++){
						my $pos08d = sprintf ("%09d", $pos + $j);
						my $key_pos = $chrno . '=' . $pos08d; 
						$indel_allpos{$key_pos} = 1;
					}
				}
			}
		}
	}
}


sub make_snp {
	my ($snp_numb) = @_;
	foreach my $chrno (sort keys %chr_rate){
		for (my $i = 0; $i < int ($snp_numb * $chr_rate{$chrno} + 0.5); $i++){
			my $pos = int (rand ($chr_length{$chrno})) + 1;
			my $pos08d = sprintf ("%09d", $pos);
			my $key_pos = $chrno . '=' . $pos08d;
			if (exists $indel_allpos{$key_pos}){
					$i = $i - 1;
					next;
			}
			else{
				$indel{$key_pos} = 'SNP';
				$count_snp++;
			}
		}
	}
}


sub subst_genome {
	my $pre_chr = '-';
	my $adjust_size = 0;
	my %adjust_size_frac;
	foreach my $key (sort keys %indel){
		my ($chr, $pos) = split (/=/, $key);
		$pre_chr = $chr if ($pre_chr eq '-');
		$adjust_size = 0 if ($chr ne $pre_chr);
		$pos =~ s/^0+//;
		$pos += $adjust_size;
$count_overlength++ if ($chr_length{$chr} < $pos);
		next if ($chr_length{$chr} < $pos);
		my $value = $indel{$key};
		my $new_key = $chr . '=' . $pos;
		next if (exists $indel_allpos_2{$new_key});
		$count_var++;
print STDERR "substituted variants - $count_var\n" if ($count_var % 10000 == 0);
		$pre_chr = $chr;
		my $sum_frac_adjust_size = 0;
		my $pre_ref_chr = '-';
		if ($value eq 'SNP'){
			foreach my $key (sort keys %chr_seq){
				my ($ref_chr, $ref_1st_pos) = split (/=/, $key);
				$ref_1st_pos =~ s/^0*//;
				$pre_ref_chr = $ref_chr if ($pre_ref_chr eq '-');
				$sum_frac_adjust_size = 0 if ($pre_ref_chr ne $ref_chr);
				$pre_ref_chr = $ref_chr;
				my $ref_1st_pos_2 = $ref_1st_pos + $sum_frac_adjust_size;
				$sum_frac_adjust_size += $adjust_size_frac{$key} if (exists $adjust_size_frac{$key});
				if (($ref_1st_pos_2 <= $pos) and ($ref_1st_pos_2 + length ($chr_seq{$key}) - 1 > $pos) and ($chr eq $ref_chr)){
#print STDERR $pos, "\t", $indel, "\t", $indel_base, "\t", $ref_1st_pos, "\t", $ref_1st_pos_2, "\t", $sum_frac_adjust_size, "\n";
					my $pos_2 = $pos - $ref_1st_pos_2 + 1;
					my $ref_base = substr ($chr_seq{$key}, $pos_2 - 1, 1);
					last if ($ref_base =~ /[^ACGT]/);
					my $snp_base = $snp_nuc_A[rand (6)] if ($ref_base eq 'A');
					$snp_base = $snp_nuc_G[rand (6)] if ($ref_base eq 'G');
					$snp_base = $snp_nuc_C[rand (6)] if ($ref_base eq 'C');
					$snp_base = $snp_nuc_T[rand (6)] if ($ref_base eq 'T');
					my $ref_base_2 = substr ($chr_seq{$key}, $pos_2 - 1, 1, $snp_base);
					my $new_line = $chr . "\t" . $pos . "\t" . 'S' . "\t" . $ref_base . "\t" . $snp_base;
					push (@indel_info, $new_line);
					$indel_allpos_2{$new_key} = 1;
					$adjust_size += 0;
					$count_snp_2++;
					last;
				}
			}
		}
		elsif ($value =~ /^D/){
			my ($D, $length) = split (/-/, $value);
			next if ($chr_length{$chr} < $pos + $length);
			foreach my $key (sort keys %chr_seq){
				my ($ref_chr, $ref_1st_pos) = split (/=/, $key);
				$ref_1st_pos =~ s/^0*//;
				$pre_ref_chr = $ref_chr if ($pre_ref_chr eq '-');
				$sum_frac_adjust_size = 0 if ($pre_ref_chr ne $ref_chr);
				$pre_ref_chr = $ref_chr;
				my $ref_1st_pos_2 = $ref_1st_pos + $sum_frac_adjust_size;
				$sum_frac_adjust_size += $adjust_size_frac{$key} if (exists $adjust_size_frac{$key});
				if (($ref_1st_pos_2 <= $pos) and ($ref_1st_pos_2 + length ($chr_seq{$key}) - 1 > $pos) and ($chr eq $ref_chr)){
					my $pos_2 = $pos - $ref_1st_pos_2 + 1;
					my $ref_base;
					if ($ref_1st_pos_2 + length ($chr_seq{$key}) - 1 - $pos >= length $length){
					    last if ((substr ($chr_seq{$key}, $pos_2, $length) =~ /[^ACGT]/) and (substr ($chr_seq{$key}, $pos_2 - 1, 1) !~ /[ACGT]/));
						$ref_base = substr ($chr_seq{$key}, $pos_2, $length, '');
						$adjust_size -= $length;
						$adjust_size_frac{$key} -= $length;
					}
					else{
					    last if (substr ($chr_seq{$key}, $pos_2, $ref_1st_pos_2 + length ($chr_seq{$key}) - 1 - $pos =~ /[^ACGT]/) and (substr ($chr_seq{$key}, $pos_2 - 1, 1) !~ /[ACGT]/));
						$ref_base = substr ($chr_seq{$key}, $pos_2, $ref_1st_pos_2 + length ($chr_seq{$key}) - 1 - $pos, '');
						$adjust_size_frac{$key} -= length ($ref_1st_pos_2 + length ($chr_seq{$key}) - 1 - $pos);
						my $count_1M_9d = sprintf ("%09d", $ref_1st_pos + 1000000);
						my $chr_chrfrac_2 = $chr . '=' . $count_1M_9d;
						my $ref_base_2 = substr ($chr_seq{$chr_chrfrac_2}, 0, $length - $ref_1st_pos_2 - length ($chr_seq{$key}) + 1 + $pos, '');
						$ref_base .= $ref_base_2;
						$adjust_size -= $length;
						$adjust_size_frac{$chr_chrfrac_2} -= length ($length - $ref_1st_pos_2 - length ($chr_seq{$key}) + 1 + $pos);
					}
					my $new_line = $chr . "\t" . $pos . "\t" . 'D' . "\t" . '-' . "\t" . $ref_base;
					push (@indel_info, $new_line);
					for (my $i = 0; $i <= $length; $i++){
						my $new_pos = $pos + $i;
						my $key_pos = $chr . '=' . $new_pos;
						$indel_allpos{$key_pos} = 1;
					}
					last;
				}
			}
		}
		elsif ($value =~ /^I/){
			my ($I, $seq) = split (/-/, $value);
			foreach my $key (sort keys %chr_seq){
				my ($ref_chr, $ref_1st_pos) = split (/=/, $key);
				$ref_1st_pos =~ s/^0*//;
				$pre_ref_chr = $ref_chr if ($pre_ref_chr eq '-');
				$sum_frac_adjust_size = 0 if ($pre_ref_chr ne $ref_chr);
				$pre_ref_chr = $ref_chr;
				my $ref_1st_pos_2 = $ref_1st_pos + $sum_frac_adjust_size;
				$sum_frac_adjust_size += $adjust_size_frac{$key} if (exists $adjust_size_frac{$key});
				if (($ref_1st_pos_2 <= $pos) and ($ref_1st_pos_2 + length ($chr_seq{$key}) - 1 > $pos) and ($chr eq $ref_chr)){
					my $pos_2 = $pos - $ref_1st_pos_2 + 1;
					last unless (substr ($chr_seq{$key}, $pos_2 - 1, 1) =~ /[ACGT]/);
					my $ref_base = substr ($chr_seq{$key}, $pos_2, 0, $seq);
					my $new_line = $chr . "\t" . $pos . "\t" . 'I' . "\t" . '+' . "\t" . $seq;
					$adjust_size += length $seq;
					$adjust_size_frac{$key} += length $seq;
					push (@indel_info, $new_line);
					for (my $i = 0; $i <= length $seq; $i++){
						my $new_pos = $pos + $i;
						my $key_pos = $chr . '=' . $new_pos;
						$indel_allpos{$key_pos} = 1;
					}
					last;
				}
			}		
		}
		
	}
}


