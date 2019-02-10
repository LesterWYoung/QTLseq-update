#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;

my $minimum_alt_num = 2;        # -n
my $minimum_alt_freq = 0.8;     # -f
my $minimum_ave_alt_qual = 20;  # -qa
my $minimun_alt_qual = 3;       # -qb   
my $maximum_read_depth = 35;    # -m
my $output_prefix = 'out';      # -p
my $minimum_het_alt_num = 3;    # -t
my $base_call_type = 'auto';  # -c
my $help;                       # -h

GetOptions(
    'num=i' => \$minimum_alt_num,
    'freq=f' => \$minimum_alt_freq,
    'qual_ave|qa=i' => \$minimum_ave_alt_qual,
    'qual_base|qb=i' => \$minimun_alt_qual,
    'maxr=i' => \$maximum_read_depth,
    'tnum=i' => \$minimum_het_alt_num,
    'pref=s' => \$output_prefix,
    'calltype=s' => \$base_call_type,
    'help' => \$help
) or pod2usage(-verbose => 0);
    pod2usage(-verbose => 0) if $help;
    
print STDERR "-numb option should not be smaller than 1\n" if ($minimum_alt_num < 1);
    

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

print STDERR "########## coval call options ##########\nnum=$minimum_alt_num freq=$minimum_alt_freq qual_ave=$minimum_ave_alt_qual qual_base=$minimun_alt_qual maxr=$maximum_read_depth calltype=$base_call_type\npref=$output_prefix\ninput=$input_file\n";

=head1 SYNOPSIS

covar call [options] <input_pileup_file>
  Output:
   prefix-snp.pileup and prefix-indel.pileup

  Options:
   --num or -n <INT>        minimum number of reads supporting non-reference allele [default: 2]
   --freq or -f <FLOAT>     minimum frequency of non-reference allele [default: 0.8]
   --qual_ave or -qa <INT>  minimum averaged base-call quality at a SNP-called site [default: 20]
   --qual_base or -qb <INT> minimum base-call quality of a non-reference base [default: 3]
   --maxr or -m <INT>       maximum read number covering non-reference allele [default: 10000]
   --tnum or -t <INT>       minimum number of reads supprting heterozygous non-reference allele [default: 3]
   --calltype or -c <STR>   quality format of the fastq file used; illumina (Phred+64) or sanger (Phred+33) [default: auto]
   --pref or -p <STR>       prefix of outputfiles [default: out]
   --help or -h             output help message

=cut

my $stderr = `samtools 2>&1`;
if ($stderr =~ /Version:\s\d\.\d\.(\d+)/){
    my $samtool_ver = $1;
    if ($samtool_ver > 8){
        print "\n########## Warning! ##########\nThe version of samtools in your PATH appears to be higher than ver. 0.1.8.\n";
        print "If you generated the pileup file with the samtools 'pileup' from ver. 0.1.9 or later, expected results will not be obtained. ";
        print "If not, ignore this warning.\n";
    }
}


my $base_qual = '';
my $sanger_type = 0;
my $solexa_type = 0;

open (FILE, $input_file) or die "$input_file: $!";
    while (<FILE>){
        my @line = split(/\s+/, $_);
        next if ($line[2] =~ /\*/);
        $base_qual = $line[9];
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


my @snp_lines;
my @indel_lines;
my %SC_indel;
my $indel_flag = 0;
my $pre_indel_chr = "";
my $pre_indel_pos = 0;
my $pre_indel_alt_num = 0;
my $indel_next = 0;
my $pre1_indel_chr = "";
my $pre2_indel_chr = "";
my $pre3_indel_chr = "";
my $pre1_pos = 0;
my $pre2_pos = 0;
my $pre3_pos = 0;
my $IorD = "";
my $indel_length = 0;
my $count_Q = 0;
my $count_Q_filt = 0;
my $count_line_filt = 0;

my $count_line = 0;
open (FILE, $input_file) or die "$input_file: $!";
    while (my $line = <FILE>){
        my $Q_C = 0;
        my $Q_A = 0;
        my $Q_G = 0;
        my $Q_T = 0;
        my @snp_base;
        my @snp_qual;        
        my @line = split(/\s+/, $line);
        my $chr = $line[0];
        my $pos = $line[1];
        my $ref_base = $line[2];
        my $alt_base = $line[3];
        my $read_depth = $line[7];
        my $base_str = $line[8];
        my $qual_str = $line[9];
#        next if ($ref_base !~ /[ACGTacgt\*]/);
        next if (($ref_base !~ /\*/) and ($alt_base =~ /[^ACGTacgt]/) and ($minimum_alt_freq >= 0.8));
        my $chr_pos = $chr . '=' . $pos;
        if (($base_str =~ /[\+-]\d+/) and ($base_str =~ /\$/)){
            my $indel_bases = $base_str;
            my $SC_read_num = $indel_bases =~ s/\$/S/g;
            $SC_indel{$chr_pos} = $SC_read_num;
        }
        next if ($ref_base eq $alt_base);
        if ($ref_base =~ /\*/){
            if (exists $SC_indel{$chr_pos}){
                $read_depth -= $SC_indel{$chr_pos};
                next if ($read_depth == 0);
            }
            $line[7] = $read_depth;
            if ($alt_base =~ /([+-])(\w+)/){
                $IorD = $1;
                $indel_length = length $2;
            }
            if (($base_str !~ /\*/) and ($line[10] / $read_depth >= $minimum_alt_freq) and ($line[10] >= $minimum_alt_num) and ($read_depth <= $maximum_read_depth)){     # remove snp_lines if their positions are 1~3 bp around indel positions
                if (($pos - $pre_indel_pos <= 2) and ($chr eq $pre_indel_chr)){
                    if ($line[10] >= $pre_indel_alt_num){
                        pop (@indel_lines);
                    }
                    else{
                        next;
                    }
                }
                if (($line[0] eq $pre3_indel_chr) and (($pos - $pre3_pos) <= 3)){
                    pop (@snp_lines);
                    pop (@snp_lines);
                    pop (@snp_lines);
                }
                elsif (($line[0] eq $pre2_indel_chr) and (($pos - $pre2_pos) <= 2)){
                    pop (@snp_lines);
                    pop (@snp_lines);
                }
                elsif (($chr eq $pre1_indel_chr) and (($pos - $pre1_pos) <= 1)){
                    pop (@snp_lines);
                }
                $line[7] = $read_depth;
                $line[11] -= $SC_indel{$chr_pos} if (exists $SC_indel{$chr_pos});
                my $new_line = join ("\t", @line);
                push (@indel_lines, $new_line);
                $pre_indel_chr = $chr;
                $pre_indel_pos = $pos;
                $pre_indel_alt_num = $line[10];
                $indel_flag = 1;
            }
            elsif (($qual_str !~ /\*/) and ($line[11] / $read_depth >= $minimum_alt_freq) and ($line[11] >= $minimum_alt_num) and ($read_depth <= $maximum_read_depth)){
                if (($pos - $pre_indel_pos <= 2) and ($chr eq $pre_indel_chr)){
                    if ($line[11] >= $pre_indel_alt_num){
                        pop (@indel_lines);
                    }
                    else{
                        next;
                    }
                }
                if (($chr eq $pre3_indel_chr) and (($pos - $pre3_pos) <= 3)){
                    pop (@snp_lines);
                    pop (@snp_lines);
                    pop (@snp_lines);
                }
                elsif (($chr eq $pre2_indel_chr) and (($pos - $pre2_pos) <= 2)){
                    pop (@snp_lines);
                    pop (@snp_lines);
                }
                elsif (($chr eq $pre1_indel_chr) and (($pos - $pre1_pos) <= 1)){
                    pop (@snp_lines);
                }
                if ($read_depth <= $maximum_read_depth){
                    $line[7] = $read_depth;
                    $line[10] -= $SC_indel{$chr_pos} if (exists $SC_indel{$chr_pos});
                    my $new_line = join ("\t", @line);
                    push (@indel_lines, $new_line);
                }
                $pre_indel_chr = $chr;
                $pre_indel_pos = $pos;
                $pre_indel_alt_num = $line[11];
                $indel_flag = 1;
            }
            next;
        }
        elsif ($indel_flag == 1){                                                                       # remove snp_lines if their positions are 1~3 bp around indel positions
            if ($IorD eq '+'){
                if (($chr eq $pre_indel_chr) and (($pos - $pre_indel_pos) <= 3)){   # when the current position is within 3 bp downstream of a reliable insersion position
                    next;
                }
                else{
                    $indel_flag = 0;
                }
            }
            elsif ($IorD eq '-'){
                if (($chr eq $pre_indel_chr) and (($pos - $pre_indel_pos) <= $indel_length)){       # when the current position is within a reliable deletion 
                    next;
                }
                elsif (($chr eq $pre_indel_chr) and ($pos <= $pre_indel_pos + $indel_length + 3)){  # when the current position is within 3 bp downstream of a reliable deletion position
                    next;
                }
                else{
                    $indel_flag = 0;
                }
            }
        }
        
        if ($base_str =~ /[+-]\d+[A-Za-z]/){
            while ($base_str =~ /[+-](\d+)[A-Za-z]/g){
                $base_str = uc $base_str;
                $base_str =~ s/([+-]\d+[A-Z]{$1})//;
            }
            $base_str =~ s/\^/\./g;
        }
        $base_str =~ s/[^A-Za-z\.,\*]//g;
        
        @snp_base = split(//, $base_str);
        @snp_qual = split(//, $qual_str);
        my @del_base = ();
        if ($base_call_type eq 'illumina'){
            for (my $i = 0; $i < @snp_qual; $i++){
                if (ord($snp_qual[$i]) - 64 < $minimun_alt_qual){
                    $read_depth -= 1;
                    push @del_base, $i;
                    next;
                }
                if ($snp_qual[$i] =~ /[abcdefgh`]/){
                    $Q_C += 40 if (($snp_base[$i] eq "C") or ($snp_base[$i] eq "c"));
                    $Q_A += 40 if (($snp_base[$i] eq "A") or ($snp_base[$i] eq "a"));
                    $Q_G += 40 if (($snp_base[$i] eq "G") or ($snp_base[$i] eq "g"));
                    $Q_T += 40 if (($snp_base[$i] eq "T") or ($snp_base[$i] eq "t"));
                }
                elsif ($snp_qual[$i] =~ /B/){
                    $Q_C += 0 if (($snp_base[$i] eq "C") or ($snp_base[$i] eq "c"));
                    $Q_A += 0 if (($snp_base[$i] eq "A") or ($snp_base[$i] eq "a"));
                    $Q_G += 0 if (($snp_base[$i] eq "G") or ($snp_base[$i] eq "g"));
                    $Q_T += 0 if (($snp_base[$i] eq "T") or ($snp_base[$i] eq "t"));
                }
                else {
                    $Q_C += ord($snp_qual[$i]) - 64 if (($snp_base[$i] eq "C") or ($snp_base[$i] eq "c"));
                    $Q_A += ord($snp_qual[$i]) - 64 if (($snp_base[$i] eq "A") or ($snp_base[$i] eq "a"));
                    $Q_G += ord($snp_qual[$i]) - 64 if (($snp_base[$i] eq "G") or ($snp_base[$i] eq "g"));
                    $Q_T += ord($snp_qual[$i]) - 64 if (($snp_base[$i] eq "T") or ($snp_base[$i] eq "t"));  
                }
            }
        }
        elsif ($base_call_type eq 'sanger'){
            for (my $i = 0; $i < @snp_qual; $i++){
                if (ord($snp_qual[$i]) - 33 < $minimun_alt_qual){
                    $read_depth -= 1;
                    push @del_base, $i;
                    next;
                }
                if ($snp_qual[$i] =~ /[ABCDEFGHI>\?@]/){
                    $Q_C += 40 if (($snp_base[$i] eq "C") or ($snp_base[$i] eq "c"));
                    $Q_A += 40 if (($snp_base[$i] eq "A") or ($snp_base[$i] eq "a"));
                    $Q_G += 40 if (($snp_base[$i] eq "G") or ($snp_base[$i] eq "g"));
                    $Q_T += 40 if (($snp_base[$i] eq "T") or ($snp_base[$i] eq "t"));
                }
                else {
                    $Q_C += ord($snp_qual[$i]) - 33 if (($snp_base[$i] eq "C") or ($snp_base[$i] eq "c"));
                    $Q_A += ord($snp_qual[$i]) - 33 if (($snp_base[$i] eq "A") or ($snp_base[$i] eq "a"));
                    $Q_G += ord($snp_qual[$i]) - 33 if (($snp_base[$i] eq "G") or ($snp_base[$i] eq "g"));
                    $Q_T += ord($snp_qual[$i]) - 33 if (($snp_base[$i] eq "T") or ($snp_base[$i] eq "t"));  
                }
            }
        }
        foreach my $delpos (sort {$b <=> $a} @del_base){
            substr ($base_str, $delpos, 1, '');
            substr ($qual_str, $delpos, 1, '');
        }
        next if ($read_depth == 0);
                    
        my $snp_C = 0;
        my $snp_A = 0;
        my $snp_G = 0;
        my $snp_T = 0;
        
        $snp_C = $base_str =~ s/[Cc]/C/g;
        $snp_A = $base_str =~ s/[Aa]/A/g;
        $snp_G = $base_str =~ s/[Gg]/G/g;
        $snp_T = $base_str =~ s/[Tt]/T/g;
        
        my %snp = ('A' => $snp_A, 'C' => $snp_C, 'G' => $snp_G, 'T' => $snp_T);
        my $top_snp = '';
        my $top_snp_num = 0;
        foreach my $key (sort {$snp{$a} <=> $snp{$b}} keys %snp){
            $top_snp = $key;
            $top_snp_num = $snp{$key};
        }
        
        my %Q = ('C' => $Q_C, 'A' => $Q_A, 'G' => $Q_G, 'T' => $Q_T);
        my $Q_ave = 0;
        $Q_ave = $Q{$top_snp} / $top_snp_num if ($top_snp_num > 0);         
        my $snp_freq = $top_snp_num / $read_depth;
        
        $count_Q += $Q_ave;
        $count_line++;
        
        if ($read_depth <= $maximum_read_depth){
            if (($top_snp_num >= $minimum_alt_num) and ($snp_freq >= $minimum_alt_freq) and ($Q_ave >= $minimum_ave_alt_qual)){
                $line[7] = $read_depth;
                $line[8] = $base_str;
                $line[9] = $qual_str;
                my $new_line = join ("\t", @line);
                $new_line .= "\t" . int ($snp_freq * 100) / 100;
                if ($minimum_alt_freq >= 0.8){
                    push (@snp_lines, $new_line);
                }
                else{
                    if ($snp_freq >= 0.9){
                        push (@snp_lines, $new_line);
                    }
                    elsif (($snp_freq < 0.9) and ($top_snp_num >= $minimum_het_alt_num)){
                        push (@snp_lines, $new_line);
                    }
                }
                $pre3_indel_chr = $pre2_indel_chr;
                $pre2_indel_chr = $pre1_indel_chr;
                $pre1_indel_chr = $chr;
                $pre3_pos = $pre2_pos;
                $pre2_pos = $pre1_pos;
                $pre1_pos = $pos;                   
                $count_Q_filt += $Q_ave;
                $count_line_filt++;
            } 
        }  
    }
close (FILE);


my $snp_outfile = "$output_prefix-snp.pileup";
my $indel_outfile = "$output_prefix-indel.pileup";

open (OUTFILE, ">$snp_outfile");
foreach my $line (@snp_lines){
    print OUTFILE $line, "\n";
#    my @line = split (/\t/, $line);
#    print OUTFILE "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[7]\t$line[10]\n";
}
close (OUTFILE);

open (OUTFILE, ">$indel_outfile");
foreach my $line (@indel_lines){
    print OUTFILE $line, "\n";
}
close (OUTFILE);

print STDERR "average Q for all SNPs/indels = ", $count_Q / $count_line, "\n";
print STDERR "average Q for selected SNPs/indels = ", $count_Q_filt / $count_line_filt, "\n";
print STDERR "total selected SNPs/indels = ", $count_line_filt, "\n";

####################################################################