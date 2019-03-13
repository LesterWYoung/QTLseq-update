#! /usr/bin/perl
use File::Basename;

use strict;
my $Command = basename $0;

if (scalar(@ARGV) < 1) {
	&usage();
}

my $pvalue_file = shift; # 1st argment # input_interval_table
my $pileup_file = shift; # 2nd argment # input_pileup
my $output_file = shift; # 3rd argment # output_pileup

if (!-e $pvalue_file) {  &usage("invalid file"); }
if (!-e $pileup_file) {  &usage("invalid file"); }

#----------------------------------------------------------------------------
# options
#----------------------------------------------------------------------------
my $onesided = 0;
#  1 : one-sided; upper
# -1 : one-sided; lower
#  0 : two-sided

#----------------------------------------------------------------------------
# Main 
#----------------------------------------------------------------------------

&scan_pileup_and_cbind($pvalue_file, $pileup_file, $output_file);


#----------------------------------------------------------------------------
# Main Routin
#----------------------------------------------------------------------------
sub scan_pileup_and_cbind() {
	my($intrvl_file0, $pileup_file1, $pileup_file2) = @_;
	
	open(INTRVL0, $intrvl_file0) or die "No $intrvl_file0";
	open(PILEUP1, $pileup_file1) or die "No $pileup_file1";
	open(PILEUP2, ">$pileup_file2") or die "No $pileup_file2";
	
	my %p_l90 = ();
	my %p_h90 = ();
	my %p_l95 = ();
	my %p_h95 = ();
	my %p_l99 = ();
	my %p_h99 = ();
	
	my $min_dep = 100000;
	my $max_dep = 0;
	
	print STDERR "load confidence interval $intrvl_file0\n";
	while (<INTRVL0>) {
		my $oneline = $_;
		chomp $oneline;
		
		my ($dep, $pl95, $ph95, $pl99, $ph99);
		if($onesided == 1){
			($dep, $ph95, $ph99) = split(/\t/, $oneline);
			$pl95 = 0;
			$pl99 = 0;
		}
		elsif($onesided == -1){
			($dep, $pl95, $pl99) = split(/\t/, $oneline);
			$ph95 = 1;
			$ph99 = 1;
		}
		else{
			($dep, $pl95, $ph95, $pl99, $ph99) = split(/\t/, $oneline);
		}
		if ($dep =~ /^[0-9]+$/){
			$p_l90{$dep} = -1;
			$p_h90{$dep} = 1;
			$p_l95{$dep} = $pl95;
			$p_h95{$dep} = $ph95;
			$p_l99{$dep} = $pl99;
			$p_h99{$dep} = $ph99;
			
			if($dep < $min_dep){
				$min_dep = $dep;
			}
			if($dep > $max_dep){
				$max_dep = $dep;
			}
		}
		
	}
	close(INTRVL0);
	print STDERR "start scan $pileup_file1\n";

	while (<PILEUP1>) {
		
		my $oneline = $_;
		chomp $oneline;
		
        ## my ($chr, $pos, $ref, $cons, $q1, $q2, $q3, $depth, $snpindex, $bases, $quals) = split(/\t/, $oneline);	
        my ($chrA, $posA, $refA, $consA, $q1A, $q2A, $q3A, $depthA, $snpindexA, $basesA, $qualsA,
			$chrB, $posB, $refB, $consB, $q1B, $q2B, $q3B, $depthB, $snpindexB, $basesB, $qualsB) = split(/\t/, $oneline);	
		
		my $depth = $depthA;
		if($depthB < $depthA){ $depth = $depthB; }
		my $snpindex = $snpindexA - $snpindexB;
		
		if($depth < $min_dep){
			$p_l90{$depth} = -1;
			$p_h90{$depth} = 1;
			$p_l95{$depth} = -1;
			$p_h95{$depth} = 1;
			$p_l99{$depth} = -1;
			$p_h99{$depth} = 1;
		}
		if($depth > $max_dep){
			$p_l90{$depth} = $p_l90{$max_dep};
			$p_h90{$depth} = $p_h90{$max_dep};
			$p_l95{$depth} = $p_l95{$max_dep};
			$p_h95{$depth} = $p_h95{$max_dep};
			$p_l99{$depth} = $p_l99{$max_dep};
			$p_h99{$depth} = $p_h99{$max_dep};
		}
		my ($j_l90, $j_m90, $j_h90);
		my ($j_l95, $j_m95, $j_h95);
		my ($j_l99, $j_m99, $j_h99);
		if   ($snpindex < $p_l90{$depth}) {	$j_l90 = -1; $j_m90 = 0; $j_h90 = 0; }
		elsif($snpindex <= $p_h90{$depth}){	$j_l90 = 0;  $j_m90 = 1; $j_h90 = 0; }
		else                              {	$j_l90 = 0;  $j_m90 = 0; $j_h90 = 1; }
		if   ($snpindex < $p_l95{$depth}) {	$j_l95 = -1; $j_m95 = 0; $j_h95 = 0; }
		elsif($snpindex <= $p_h95{$depth}){	$j_l95 = 0;  $j_m95 = 1; $j_h95 = 0; }
		else                              {	$j_l95 = 0;  $j_m95 = 0; $j_h95 = 1; }
		if   ($snpindex < $p_l99{$depth}) {	$j_l99 = -1; $j_m99 = 0; $j_h99 = 0; }
		elsif($snpindex <= $p_h99{$depth}){	$j_l99 = 0;  $j_m99 = 1; $j_h99 = 0; }
		else                              {	$j_l99 = 0;  $j_m99 = 0; $j_h99 = 1; }
		
		printf PILEUP2 "%s\t%d\t%f\t%f\t%f\t%d\t%d\t%d\t%f\t%f\t%d\t%d\t%d\t%f\t%f\t%d\t%d\t%d\n", $oneline, $depth, $snpindex, $p_h90{$depth}, $p_l90{$depth}, $j_l90, $j_m90, $j_h90, $p_h95{$depth}, $p_l95{$depth}, $j_l95, $j_m95, $j_h95, $p_h99{$depth}, $p_l99{$depth}, $j_l99, $j_m99, $j_h99;
		
		# next if ($ref =~ /^\*$/);
	}
	
	close(PILEUP1);
	close(PILEUP2);

	print STDERR "scan finished! $pileup_file1\n";
}


#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
sub usage() {
	my($msg) = @_;
	if($msg ne ""){
		print STDERR $msg . "\n";
	}
	print STDERR "usage: $Command <in.intervaltable.txt> <in.pileup> <out.pileup>\n";  
	exit; 
}



