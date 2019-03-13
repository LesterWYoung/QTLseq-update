#! /usr/bin/perl
use File::Basename;

use strict;
my $Command = basename $0;

if (scalar(@ARGV) < 1) {
	&usage();
}

my $pileup_file1 = shift; # 1st argment	# inputname
my $pileup_file2 = shift; # 2nd argment # outputname

#if (!-e $pileup_file1 || !-e $pileup_file2) {  &usage("invalid file"); }
if (!-e $pileup_file1) {  &usage("invalid file"); }

#----------------------------------------------------------------------------
# options
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
# Main 
#----------------------------------------------------------------------------

&scan_pileup_and_convert($pileup_file1, $pileup_file2);


#----------------------------------------------------------------------------
# Main Routin
#----------------------------------------------------------------------------
sub scan_pileup_and_convert() {
	my($pileup_file1, $pileup_file2) = @_;
	
	open(PILEUP1, $pileup_file1) or die "No $pileup_file1";
	open(PILEUP2, ">$pileup_file2") or die "No $pileup_file2";
	
	print STDERR "start scan $pileup_file1\n";

	while (<PILEUP1>) {
		# my ($chr, $pos, $ref, $cons, $bases) = (split)[0,1,2,3,8];
		
		my $oneline = $_;
		chomp $oneline;
		my @columns = split(/\t/, $oneline);
		my $ncol = @columns;
		# print "n : $ncol\n";
		
        my ($chr, $pos, $ref, $cons, $q1, $q2, $q3, $depth, $bases, $quals, $snpindex) = split(/\t/, $oneline);	
		# print "$chr\n";
		# print "$pos\n";
		# print "$ref\n";
		# print "$cons\n";
		# print "$bases\n";
		
		my $cnvcons = "N";
		
		if ($cons =~ /^[ACGTN]$/) {
			# print         "$_";
			print PILEUP2 "$_";
		}
		else {	
			if ($ref =~ /^[Tt]$/ && $cons =~ /^[Rr]$/) {
				$cnvcons = &get_major($ref, $cons, $bases);
			}
			elsif ($ref =~ /^[Cc]$/ && $cons =~ /^[Rr]$/) {
				$cnvcons = &get_major($ref, $cons, $bases);
			}
			elsif ($ref =~ /^[Aa]$/ && $cons =~ /^[Yy]$/) {
				$cnvcons = &get_major($ref, $cons, $bases);
			}
			elsif ($ref =~ /^[Gg]$/ && $cons =~ /^[Yy]$/) {
				$cnvcons = &get_major($ref, $cons, $bases);
			}
			elsif ($ref =~ /^[Aa]$/ && $cons =~ /^[Kk]$/) {
				$cnvcons = &get_major($ref, $cons, $bases);
			}
			elsif ($ref =~ /^[Cc]$/ && $cons =~ /^[Kk]$/) {
				$cnvcons = &get_major($ref, $cons, $bases);
			}
			elsif ($ref =~ /^[Gg]$/ && $cons =~ /^[Mm]$/) {
				$cnvcons = &get_major($ref, $cons, $bases);
			}
			elsif ($ref =~ /^[Tt]$/ && $cons =~ /^[Mm]$/) {
				$cnvcons = &get_major($ref, $cons, $bases);
			}
			elsif ($ref =~ /^[Aa]$/ && $cons =~ /^[Ss]$/) {
				$cnvcons = &get_major($ref, $cons, $bases);
			}
			elsif ($ref =~ /^[Tt]$/ && $cons =~ /^[Ss]$/) {
				$cnvcons = &get_major($ref, $cons, $bases);
			}
			elsif ($ref =~ /^[Gg]$/ && $cons =~ /^[Ww]$/) {
				$cnvcons = &get_major($ref, $cons, $bases);
			}
			elsif ($ref =~ /^[Cc]$/ && $cons =~ /^[Ww]$/) {
				$cnvcons = &get_major($ref, $cons, $bases);
			}
			# 2013/09/26 bug fixed by Yaegashi
			else{
				$cnvcons = &RYKMSW_to_ACGT($ref, $cons, $bases);	
			}
			
			# # printf STDERR  "%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%.2f\n", $chr, $pos, $ref, $cnvcons, $q1, $q2, $q3, $depth, $bases, $quals, $snpindex;
			# printf PILEUP2 "%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%.2f\n", $chr, $pos, $ref, $cnvcons, $q1, $q2, $q3, $depth, $bases, $quals, $snpindex;

			# printf STDERR  "%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n", $chr, $pos, $ref, $cnvcons, $q1, $q2, $q3, $depth, $bases, $quals, $snpindex;
			printf PILEUP2 "%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n", $chr, $pos, $ref, $cnvcons, $q1, $q2, $q3, $depth, $bases, $quals, $snpindex;
		}
		# next if ($ref =~ /^\*$/);
	}
	
	close(PILEUP1);
	close(PILEUP2);

	print STDERR "scan finished! $pileup_file1\n";

}


#----------------------------------------------------------------------------
# RYKMSW --> ACGT
#----------------------------------------------------------------------------
sub RYKMSW_to_ACGT() {
	my ($ref, $cons, $bases) = @_;

	my $retval = "N";
	if($ref =~ /^[Gg]$/ && $cons =~ /^[Rr]$/){ 
		$retval = "A";
	}
	elsif($ref =~ /^[Aa]$/ && $cons =~ /^[Rr]$/){ 
		$retval = "G";
	}
	elsif($ref =~ /^[Tt]$/ && $cons =~ /^[Yy]$/){ 
		$retval = "C";
	}
	elsif($ref =~ /^[Cc]$/ && $cons =~ /^[Yy]$/){ 
		$retval = "T";
	}
	elsif($ref =~ /^[Tt]$/ && $cons =~ /^[Kk]$/){ 
		$retval = "G";
	}
	elsif($ref =~ /^[Gg]$/ && $cons =~ /^[Kk]$/){ 
		$retval = "T";
	}
	elsif($ref =~ /^[Cc]$/ && $cons =~ /^[Mm]$/){ 
		$retval = "A";
	}
	elsif($ref =~ /^[Aa]$/ && $cons =~ /^[Mm]$/){ 
		$retval = "C";
	}
	elsif($ref =~ /^[Cc]$/ && $cons =~ /^[Ss]$/){ 
		$retval = "G";
	}
	elsif($ref =~ /^[Gg]$/ && $cons =~ /^[Ss]$/){ 
		$retval = "C";
	}
	elsif($ref =~ /^[Tt]$/ && $cons =~ /^[Ww]$/){ 
		$retval = "A";
	}
	elsif($ref =~ /^[Aa]$/ && $cons =~ /^[Ww]$/){ 
		$retval = "T";
	}
	else{
		$retval = "N";
	}

	return $retval;


}

#----------------------------------------------------------------------------
# RYKMSW --> ACGT
#----------------------------------------------------------------------------
sub get_major() {
	my ($ref, $cons, $bases) = @_;
	
	my $len0 = length($bases);
	# print "$len0\n";
	
	my $strxxxx = $bases;
	my $strxCGT = $bases;
	my $strAxGT = $bases;
	my $strACxT = $bases;
	my $strACGx = $bases;

	$strxxxx =~ s/[ACGTacgt]//g;
	$strxCGT =~ s/[Aa]//g;
	$strAxGT =~ s/[Cc]//g;
	$strACxT =~ s/[Gg]//g;
	$strACGx =~ s/[Tt]//g;

	my $lenxxxx = length($strxxxx);
	my $lenxCGT = length($strxCGT);
	my $lenAxGT = length($strAxGT);
	my $lenACxT = length($strACxT);
	my $lenACGx = length($strACGx);

	my $cntA = $len0 - $lenxCGT;
	my $cntC = $len0 - $lenAxGT;
	my $cntG = $len0 - $lenACxT;
	my $cntT = $len0 - $lenACGx;
	
	my $retval = "N";

	if($cons =~ /^[Rr]$/){ 
		if ($cntA > $cntG){
			$retval = "A";
		}
		else{
			$retval = "G";
		}
	}
	elsif($cons =~ /^[Yy]$/){ 
		if ($cntC > $cntT){
			$retval = "C";
		}
		else{
			$retval = "T";
		}
	}
	elsif($cons =~ /^[Kk]$/){ 
		if ($cntG > $cntT){
			$retval = "G";
		}
		else{
			$retval = "T";
		}
	}
	elsif($cons =~ /^[Mm]$/){ 
		if ($cntA > $cntC){
			$retval = "A";
		}
		else{
			$retval = "C";
		}
	}
	elsif($cons =~ /^[Ss]$/){ 
		if ($cntG > $cntC){
			$retval = "G";
		}
		else{
			$retval = "C";
		}
	}
	elsif($cons =~ /^[Ww]$/){ 
		if ($cntG > $cntC){
			$retval = "A";
		}
		else{
			$retval = "T";
		}
	}
	else{
		$retval = "N";
	}

	return $retval;
}

#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
sub usage() {
	my($msg) = @_;
	if($msg ne ""){
		print STDERR $msg . "\n";
	}
	print STDERR "usage: $Command <file1> <file2>\n";  
	exit; 
}


