#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;

# 2011.09.27 mod by summer: output file to current directory
use File::Basename;
my $my_fname = "";

# extract common SNP positions from 2~13 pileup files

my ($file_1, $file_2, $file_3, $file_4, $file_5, $file_6, $file_7, $file_8, $file_9, $file_10, $file_11, $file_12, $file_13, $match, $help);
my (@chr_pos_1, @chr_pos_2, @chr_pos_3, @chr_pos_4, @chr_pos_5, @chr_pos_6, @chr_pos_7, @chr_pos_8, @chr_pos_9, @chr_pos_10, @chr_pos_11, @chr_pos_12, @chr_pos_13);

# $file_1 from @ARGV is 1st filtered pileup file.
# $file_2 from @ARGV is 2d filtered pileup file.
# $file_2 from @ARGV is 3d filtered pileup file.
# $file_2 from @ARGV is 4th filtered pileup file.
#
#
# $file_13 from @ARGV is 13th filtered pileup file.


GetOptions(
    'help' => \$help
) or pod2usage(-verbose => 0);
    pod2usage(-verbose => 0) if $help;
    

unless (@ARGV){
print "Your pileup files should be added to the argument.\n";
exit;
}

=head1 SYNOPSIS

extract_common_snp_pos.pl <input_pileup_file1> <input_pileup_file2>-----
  
    Output:
    (prefix of input_pileup_file1)-common-pos.txt, (prefix of input_pileup_file2)-common-pos.txt, -----

    Options:
    --help or -h            output help message

=cut

my $file_number = scalar @ARGV;
#print "$file_number\n";

$file_1 = shift (@ARGV);
unless (-e $file_1){
die "'$file_1' does not exist.\n";
}
elsif (-f $file_1){
open (FILE, $file_1) or die "$file_1: $!";
	while (my $line = <FILE>){
            chomp $line;
            my @line = split(/\s+/, $line);
            my $item = $line[0].'='.$line[1];
            push (@chr_pos_1, $item);
        }
close (FILE);
}

$file_2 = shift (@ARGV);
unless (-e $file_2){
die "'$file_2' does not exist.\n";
}
elsif (-f $file_2){
open (FILE, $file_2) or die "$file_2: $!";
	while (my $line = <FILE>){
            chomp $line;
            my @line = split(/\s+/, $line);
            my $item = $line[0].'='.$line[1];
            push (@chr_pos_2, $item);
        }
close (FILE);
}

if (@ARGV > 0){
$file_3 = shift (@ARGV) if (@ARGV > 0);
unless (-e $file_3){
die "'$file_3' does not exist.\n";
}
elsif (-f $file_3){
open (FILE, $file_3) or die "$file_3: $!";
	while (my $line = <FILE>){
            chomp $line;
            my @line = split(/\s+/, $line);
            my $item = $line[0].'='.$line[1];
            push (@chr_pos_3, $item);
        }
close (FILE);
}
}

if (@ARGV > 0){
$file_4 = shift (@ARGV);
unless (-e $file_4){
die "'$file_4' does not exist.\n";
}
elsif (-f $file_4){
open (FILE, $file_4) or die "$file_4: $!";
	while (my $line = <FILE>){
            chomp $line;
            my @line = split(/\s+/, $line);
            my $item = $line[0].'='.$line[1];
            push (@chr_pos_4, $item);
        }
close (FILE);
}
}

if (@ARGV > 0){
$file_5 = shift (@ARGV);
unless (-e $file_5){
die "'$file_5' does not exist.\n";
}
elsif (-f $file_5){
open (FILE, $file_5) or die "$file_5: $!";
	while (my $line = <FILE>){
            chomp $line;
            my @line = split(/\s+/, $line);
            my $item = $line[0].'='.$line[1];
            push (@chr_pos_5, $item);
        }
close (FILE);
}
}

if (@ARGV > 0){
$file_6 = shift (@ARGV);
unless (-e $file_6){
die "'$file_6' does not exist.\n";
}
elsif (-f $file_6){
open (FILE, $file_6) or die "$file_6: $!";
	while (my $line = <FILE>){
            chomp $line;
            my @line = split(/\s+/, $line);
            my $item = $line[0].'='.$line[1];
            push (@chr_pos_6, $item);
        }
close (FILE);
}
}

if (@ARGV > 0){
$file_7 = shift (@ARGV);
unless (-e $file_7){
die "'$file_7' does not exist.\n";
}
elsif (-f $file_7){
open (FILE, $file_7) or die "$file_7: $!";
	while (my $line = <FILE>){
            chomp $line;
            my @line = split(/\s+/, $line);
            my $item = $line[0].'='.$line[1];
            push (@chr_pos_7, $item);
        }
close (FILE);
}
}

if (@ARGV > 0){
$file_8 = shift (@ARGV);
unless (-e $file_8){
die "'$file_8' does not exist.\n";
}
elsif (-f $file_8){
open (FILE, $file_8) or die "$file_8: $!";
	while (my $line = <FILE>){
            chomp $line;
            my @line = split(/\s+/, $line);
            my $item = $line[0].'='.$line[1];
            push (@chr_pos_8, $item);
        }
close (FILE);
}
}

if (@ARGV > 0){
$file_9 = shift (@ARGV);
unless (-e $file_9){
die "'$file_9' does not exist.\n";
}
elsif (-f $file_9){
open (FILE, $file_9) or die "$file_9: $!";
	while (my $line = <FILE>){
            chomp $line;
            my @line = split(/\s+/, $line);
            my $item = $line[0].'='.$line[1];
            push (@chr_pos_9, $item);
        }
close (FILE);
}
}

if (@ARGV > 0){
$file_10 = shift (@ARGV);
unless (-e $file_10){
die "'$file_10' does not exist.\n";
}
elsif (-f $file_10){
open (FILE, $file_10) or die "$file_10: $!";
	while (my $line = <FILE>){
            chomp $line;
            my @line = split(/\s+/, $line);
            my $item = $line[0].'='.$line[1];
            push (@chr_pos_10, $item);
        }
close (FILE);
}
}

if (@ARGV > 0){
$file_11 = shift (@ARGV);
unless (-e $file_11){
die "'$file_11' does not exist.\n";
}
elsif (-f $file_11){
open (FILE, $file_11) or die "$file_11: $!";
	while (my $line = <FILE>){
            chomp $line;
            my @line = split(/\s+/, $line);
            my $item = $line[0].'='.$line[1];
            push (@chr_pos_11, $item);
        }
close (FILE);
}
}

if (@ARGV > 0){
$file_12 = shift (@ARGV);
unless (-e $file_12){
die "'$file_12' does not exist.\n";
}
elsif (-f $file_12){
open (FILE, $file_12) or die "$file_12: $!";
	while (my $line = <FILE>){
            chomp $line;
            my @line = split(/\s+/, $line);
            my $item = $line[0].'='.$line[1];
            push (@chr_pos_12, $item);
        }
close (FILE);
}
}

if (@ARGV > 0){
$file_13 = shift (@ARGV);
unless (-e $file_13){
die "'$file_13' does not exist.\n";
}
elsif (-f $file_13){
open (FILE, $file_13) or die "$file_13: $!";
	while (my $line = <FILE>){
            chomp $line;
            my @line = split(/\s+/, $line);
            my $item = $line[0].'='.$line[1];
            push (@chr_pos_13, $item);
        }
close (FILE);
}
}

# extract the chr-positions common between at least 2 files of 2~13 pileup files

my @chr_pos_ref = (\@chr_pos_1, \@chr_pos_2, \@chr_pos_3, \@chr_pos_4, \@chr_pos_5, \@chr_pos_6, \@chr_pos_7, \@chr_pos_8, \@chr_pos_9, \@chr_pos_10, \@chr_pos_11, \@chr_pos_12, \@chr_pos_13);
my (@common_pos_1, @common_pos_2, @common_pos_3, @common_pos_4, @common_pos_5, @common_pos_6, @common_pos_7, @common_pos_8, @common_pos_9, @common_pos_10, @common_pos_11, @common_pos_12, @common_pos_13);
my @common_pos_ref = (\@common_pos_1, \@common_pos_2, \@common_pos_3, \@common_pos_4, \@common_pos_5, \@common_pos_6, \@common_pos_7, \@common_pos_8, \@common_pos_9, \@common_pos_10, \@common_pos_11, \@common_pos_12, \@common_pos_13);

for (my $i = 0; $i < $file_number; $i++){
    for (my $j = 0; $j < $file_number; $j++){
	my %count;
	@{$common_pos_ref[$i]} = (@{$common_pos_ref[$i]}, grep {$count{$_}++} (@{$chr_pos_ref[$i]}, @{$chr_pos_ref[$j]})) if ($i != $j);
    }
# extract >= 3 chr-positions (To extract >= 2 chr-positions, comment out the following 2 lines. To extract >= 4 chr-positions, change '>= 2' to '>= 3')
#    my %count1;
#    @{$common_pos_ref[$i]} = grep {$count1{$_} >= 2} grep {++$count1{$_} > 1} (@{$common_pos_ref[$i]});

    my %count2;
# remove redundant chr-positions
    @{$common_pos_ref[$i]} = grep {!$count2{$_}++} (@{$common_pos_ref[$i]});
print STDERR scalar @{$common_pos_ref[$i]}, "\n";
}



# summer
$my_fname = basename($file_1); $my_fname =~ s/\.pileup$//;
my $out_file_1 = "${my_fname}-common-pos.txt" if ($file_1 =~ /^(.*)\.pileup/);

open (NEWFILE, "> $out_file_1") or die "$!";
foreach (@common_pos_1){
    print NEWFILE $_, "\n";
}
close (NEWFILE);


# summer
$my_fname = basename($file_2); $my_fname =~ s/\.pileup$//;
my $out_file_2 = "${my_fname}-common-pos.txt" if ($file_2 =~ /^(.*)\.pileup/);
open (NEWFILE, "> $out_file_2") or die "$!";
foreach (@common_pos_2){
    print NEWFILE $_, "\n";
}
close (NEWFILE);

if ($file_number > 2){
# summer
$my_fname = basename($file_3); $my_fname =~ s/\.pileup$//;
my $out_file_3 = "${my_fname}-common-pos.txt" if ($file_3 =~ /^(.*)\.pileup/);
open (NEWFILE, "> $out_file_3") or die "$!";
foreach (@common_pos_3){
    print NEWFILE $_, "\n";
}
close (NEWFILE);
}

if ($file_number > 3){
# summer
$my_fname = basename($file_4); $my_fname =~ s/\.pileup$//;
my $out_file_4 = "${my_fname}-common-pos.txt" if ($file_4 =~ /^(.*)\.pileup/);
open (NEWFILE, "> $out_file_4") or die "$!";
foreach (@common_pos_4){
    print NEWFILE $_, "\n";
}
close (NEWFILE);
}

if ($file_number > 4){
# summer
$my_fname = basename($file_5); $my_fname =~ s/\.pileup$//;
my $out_file_5 = "${my_fname}-common-pos.txt" if ($file_5 =~ /^(.*)\.pileup/);
open (NEWFILE, "> $out_file_5") or die "$!";
foreach (@common_pos_5){
    print NEWFILE $_, "\n";
}
close (NEWFILE);
}

if ($file_number > 5){
# summer
$my_fname = basename($file_6); $my_fname =~ s/\.pileup$//;
my $out_file_6 = "${my_fname}-common-pos.txt" if ($file_6 =~ /^(.*)\.pileup/);
open (NEWFILE, "> $out_file_6") or die "$!";
foreach (@common_pos_6){
    print NEWFILE $_, "\n";
}
close (NEWFILE);
}

if ($file_number > 6){
# summer
$my_fname = basename($file_7); $my_fname =~ s/\.pileup$//;
my $out_file_7 = "${my_fname}-common-pos.txt" if ($file_7 =~ /^(.*)\.pileup/);
open (NEWFILE, "> $out_file_7") or die "$!";
foreach (@common_pos_7){
    print NEWFILE $_, "\n";
}
close (NEWFILE);
}


if ($file_number > 7){
# summer
$my_fname = basename($file_8); $my_fname =~ s/\.pileup$//;
my $out_file_8 = "${my_fname}-common-pos.txt" if ($file_8 =~ /^(.*)\.pileup/);
open (NEWFILE, "> $out_file_8") or die "$!";
foreach (@common_pos_8){
    print NEWFILE $_, "\n";
}
close (NEWFILE);
}

if ($file_number > 8){
# summer
$my_fname = basename($file_9); $my_fname =~ s/\.pileup$//;
my $out_file_9 = "${my_fname}-common-pos.txt" if ($file_9 =~ /^(.*)\.pileup/);
open (NEWFILE, "> $out_file_9") or die "$!";
foreach (@common_pos_9){
    print NEWFILE $_, "\n";
}
close (NEWFILE);
}

if ($file_number > 9){
# summer
$my_fname = basename($file_10); $my_fname =~ s/\.pileup$//;
my $out_file_10 = "${my_fname}-common-pos.txt" if ($file_10 =~ /^(.*)\.pileup/);
open (NEWFILE, "> $out_file_10") or die "$!";
foreach (@common_pos_10){
    print NEWFILE $_, "\n";
}
close (NEWFILE);
}

if ($file_number > 10){
# summer
$my_fname = basename($file_11); $my_fname =~ s/\.pileup$//;
my $out_file_11 = "${my_fname}-common-pos.txt" if ($file_11 =~ /^(.*)\.pileup/);
open (NEWFILE, "> $out_file_11") or die "$!";
foreach (@common_pos_11){
    print NEWFILE $_, "\n";
}
close (NEWFILE);
}

if ($file_number > 11){
# summer
$my_fname = basename($file_12); $my_fname =~ s/\.pileup$//;
my $out_file_12 = "${my_fname}-common-pos.txt" if ($file_12 =~ /^(.*)\.pileup/);
open (NEWFILE, "> $out_file_12") or die "$!";
foreach (@common_pos_12){
    print NEWFILE $_, "\n";
}
close (NEWFILE);
}

if ($file_number > 12){
# summer
$my_fname = basename($file_13); $my_fname =~ s/\.pileup$//;
my $out_file_13 = "${my_fname}-common-pos.txt" if ($file_13 =~ /^(.*)\.pileup/);
open (NEWFILE, "> $out_file_13") or die "$!";
foreach (@common_pos_13){
    print NEWFILE $_, "\n";
}
close (NEWFILE);
}
=cut
