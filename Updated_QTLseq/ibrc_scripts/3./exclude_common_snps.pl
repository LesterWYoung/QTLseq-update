#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;

# 2011.09.27 mod by summer: output file to current directory
# summer
use File::Basename;
my $my_fname = "";

# exclude SNP lines in pileup files, corresponding to **-common-pos.txt files containig the common SNP positions between 2-13 pileup files 

# $input_file_1 from @ARGV is a filtered pileup file.
# $input_file_2 from @ARGV is a common-snp-pos.txt file for $input_file_1.

my %chr_pos;
my $k = 0;
my $help;

GetOptions(
    'help' => \$help
) or pod2usage(-verbose => 0);
    pod2usage(-verbose => 0) if $help;
    

unless (@ARGV){
print "Your pileup file and its common-snp-pos.txt file should be added to the argument.\n";
exit;
}

=head1 SYNOPSIS

exclude_common_snps.pl <input_pileup_file> <input_common-snp-pos.txt>
  
  Output:
   (prefix of input_pileup_file)-rmc2snp.pileup

  Options:
   --help or -h            output help message

=cut

my $file_number = scalar @ARGV;

my $input_file_1 = shift (@ARGV);
my $input_file_2 = shift (@ARGV);

unless (-e $input_file_2){
die "'$input_file_1' does not exist.\n";
}
elsif (-f $input_file_2){
open (FILE, $input_file_2) or die "$input_file_2: $!";
	while (my $line = <FILE>){
            chomp $line;
	    $chr_pos{$line} = 1;
        }
close (FILE);
}

unless (-e $input_file_1){
die "'$input_file_1' does not exist.\n";
}
elsif (-f $input_file_1){
#summer
$my_fname = basename($input_file_1); $my_fname =~ s/\.pileup$//;
my $out_file = "${my_fname}-rmc2snp.pileup" if ($input_file_1 =~ /^(.*)\.pileup/);
open (FILE, $input_file_1) or die "$input_file_1: $!";
open (NEWFILE, "> $out_file") or die "$!";
	while (my $line = <FILE>){
            chomp $line;
            my @line = split(/\s+/, $line);
            my $item = $line[0].'='.$line[1];
	    if (exists $chr_pos{$item}){
		next;
	    }
            else{
		print NEWFILE $line, "\n";
	    }
        }
close (NEWFILE);
close (FILE);
}
