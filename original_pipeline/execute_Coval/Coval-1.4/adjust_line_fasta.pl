#!/usr/bin/perl -w
use strict;

unless (@ARGV){
print "Your sam file should be added to the argument.\n";
exit;
}
my $file_1 = shift (@ARGV);
unless (-e $file_1){
die "'$file_1' does not exist.\n";
}

open (FILE, $file_1) or die "$file_1: $!";
my $pre_chr = '-';
my $chr;
my $chr_seq = '';
while (my $line = <FILE>){
        $line =~ s/(\r\n|\n|\r)//g;
        if ($line =~ /^>(\S+)/){
                $chr = $1;
                $chr = $1 if ($chr =~ /(\S+)\|size\d+/);
                if ($chr_seq ne ''){
                        print ">$pre_chr\n";
                        while (length $chr_seq > 0){
                                if (length $chr_seq >= 50){
                                        my $subseq = substr ($chr_seq, 0, 50, '');
                                        print $subseq, "\n";
                                }
                                else{
                                        print $chr_seq, "\n";
                                        $chr_seq = '';
                                }
                        }
                }
                $chr_seq = '';
        }
        else{
                $chr_seq .= uc $line;
        }
        $pre_chr = $chr;
}
close (FILE);
if ($chr_seq ne ''){
        print ">$pre_chr\n";
        while (length $chr_seq > 0){
                if (length $chr_seq >= 50){
                        my $subseq = substr ($chr_seq, 0, 50, '');
                        print $subseq, "\n";
                }
                else{
                        print $chr_seq, "\n";
                        $chr_seq = '';
                }
        }
}
$chr_seq = '';
