#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;

# generate reference genome that was substituted with SNPs of a filtered pileup file

my $reference = '/home/kosugis/reference/mouse_129S5.fa';
my $help;

GetOptions(
    'ref=s' => \$reference,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;
pod2usage(-verbose => 0) unless (-e $reference);

unless (@ARGV){
print "Your pileup file should be added to the argument.\n";
exit;
}
my $file_1 = shift (@ARGV);
if ($file_1 eq '-'){
    # no operation
}
elsif (!-f $file_1){
    die "'$file_1' does not exist.\n";
}

=head1 SYNOPSIS

make_consensus.pl -r <reference_file> <input_file (filtered pileup file)>
  Output:
   STDOUT		a sunstituted fasta file

  Options:
   --ref or -r <STR>       a reference file to be substituted
   --help or -h            output help message

=cut


my %chr_no;				#2014/01/09 kikuchi
my %chr_seq;
open (FILE, $reference) or die "$!";
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
		
	    $chr_no{$count_chr} = $ref_chr_name;	#2014/01/09 kikuchi
#print "no $count_chr, 	$ref_chr_name"
	}
	else{
	    $chrseq .= uc $line;
	}
    }
    $chr_seq{$ref_chr_name} = $chrseq;
    $chrseq = '';
close (FILE);


open (FILE, $file_1) or die "$file_1: $!";
    while (my $line = <FILE>){
	chomp $line;
	my @line = split (/\s+/, $line);
	$line[0] =~ /(\S*)/;
	my $chr = $1;
	my $ref_base = substr ($chr_seq{$chr}, $line[1] - 1, 1, $line[3]);
	print STDERR 'Ref base is different from that of the pileup file', "\n" if (uc $ref_base ne uc $line[2]);
    }
close (FILE);

#								##########2014/01/09 kikuchi
foreach my $key (sort {$a <=> $b} keys %chr_no){
    print '>', "$chr_no{$key}\n";
    while (length $chr_seq{$chr_no{$key}} > 0){
	my $seq_50 = substr ($chr_seq{$chr_no{$key}}, 0, 50, '');
	print $seq_50, "\n";
    }
}
#								##########2014/01/09 kikuchi


#foreach my $key (sort keys %chr_seq){
#    print '>', "$key\n";
#    while (length $chr_seq{$key} > 0){
#	my $seq_50 = substr ($chr_seq{$key}, 0, 50, '');
#	print $seq_50, "\n";
#    }
#}
	