#! /usr/bin/perl
use File::Basename;

use strict;
my $Command = basename $0;

if (scalar(@ARGV) < 1) {
	&usage();
}

my $pileup_file1 = shift; # 1st argment	# inputname
my $pileup_file2 = shift; # 2nd argment # outputname

if (!-e $pileup_file1) {  &usage("invalid file"); }

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
		
		my $oneline = $_;
		chomp $oneline;
		my @columns = split(/\t/, $oneline);
		my $ncol = @columns;
		
        my ($chr, $pos, $ref, $cons, $q1, $q2, $q3, $depth, $bases, $quals) = split(/\t/, $oneline);	
		
		my $cnvcons = "n";
		my $snpindex = 0.0;
		
		next if ($ref =~ /^\*$/);
		
		my $trim_bases = $bases;
		$trim_bases = &trim_hat_dollar($ref, $cons, $bases);
		$trim_bases = &trim_indels($ref, $cons, $trim_bases);

		$cnvcons = &get_major($ref, $cons, $trim_bases);
		$snpindex = &get_snpindex($ref, $cons, $depth, $trim_bases, $cnvcons);

		#printf PILEUP2 "%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n", $chr, $pos, $ref, $cnvcons, $q1, $q2, $q3, $depth, $bases, $quals, $snpindex;
		printf PILEUP2 "%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n", $chr, $pos, $ref, $cons, $q1, $q2, $q3, $depth, $bases, $quals, $snpindex;

	}
	
	close(PILEUP1);
	close(PILEUP2);

	print STDERR "scan finished! $pileup_file1\n";

}


#----------------------------------------------------------------------------
# 
#----------------------------------------------------------------------------
sub get_major() {
	my ($ref, $cons, $bases) = @_;
	
	my $retval = "N";
	my $act_bases = $bases;
	
	$act_bases =~ s/[^.,\*ACGTNacgtn]//g;
	$act_bases =~ s/[.,]/${ref}/g;
	$act_bases =~ s/\*/${cons}/g;

	my $a_freq = $act_bases;
	my $c_freq = $act_bases;
	my $g_freq = $act_bases;
	my $t_freq = $act_bases;
	
	$a_freq =~ s/[^Aa]//g;
	$c_freq =~ s/[^Cc]//g;
	$g_freq =~ s/[^Gg]//g;
	$t_freq =~ s/[^Tt]//g;
	
	my $a_len = length($a_freq);
	my $c_len = length($c_freq);
	my $g_len = length($g_freq);
	my $t_len = length($t_freq);

	if($a_len >= $c_len){
		if($a_len >= $g_len){
			$retval = ($a_len >= $t_len)? "A" : "T";
		}
		else{
			$retval = ($g_len >= $t_len)? "G" : "T";
		}
	}
	else{
		if($c_len >= $g_len){
			$retval = ($c_len >= $t_len)? "C" : "T";
		}
		else{
			$retval = ($g_len >= $t_len)? "G" : "T";
		}
	}
	# print "$ref\t$cons\t$bases\t$act_bases\tA:$a_len\tC:$c_len\tG:$g_len\tT:$t_len\t$retval\n";
	return $retval;
}

#----------------------------------------------------------------------------
# 
#----------------------------------------------------------------------------
sub get_snpindex() {
	my ($ref, $cons, $depth, $bases, $cnvcons) = @_;
	my $retval = 0.0;
	my $act_bases = $bases;
	
	$act_bases =~ s/[^.,\*ACGTNacgtn]//g;
	$act_bases =~ s/[.,]/${ref}/g;
	$act_bases =~ s/\*/${cons}/g;
	
	my $cnv_freq = $act_bases;
	$cnv_freq =~ s/[^${cnvcons}]//ig;
	
	if($ref =~ /${cnvcons}/i){
		$retval = $depth - length($cnv_freq);
		$retval /= $depth;
	}
	else{
		$retval = length($cnv_freq);
		$retval /= $depth;
	}
	
	return $retval;
}

#----------------------------------------------------------------------------
#
#----------------------------------------------------------------------------
sub trim_hat_dollar(){
	my ($ref, $cons, $bases) = @_;
	my $act_bases = $bases;

	# --------------------------------------------------
	# a symbol '^' + one ASCII or a symol '$'
	# --------------------------------------------------
	if($act_bases =~ /\^.|\$/){
		$act_bases =~ s/\^.|\$//g;
	}
	return $act_bases;
}

#----------------------------------------------------------------------------
# for trim insertions and deletions
#	befor=,,,-4tata,-4tata,-4tata,-4tata,-4tata..,-6tatata.,-4tata
#	after=,,,,,,,..,.,
#----------------------------------------------------------------------------
sub trim_indels(){
	my ($ref, $cons, $bases) = @_;
	my $act_bases = $bases;

	if($act_bases =~ /[+-]\d+[ACGTNacgtn]+/){
		my @ppp = split(/[+-]/, $act_bases);
		$act_bases = "";
		foreach my $qqq (@ppp){
			my $rrr = $qqq;
			if($rrr =~ /(\d+)(.*)/){
				# ----------------------------------------
				# if $rrr = 11CATCATACAAA,g
				#	==> $1 = 11
				#	==> $2 = CATCATACAAA,g
				#	==> substr($2, $1) = ,g
				# ----------------------------------------
			 	if($1 > 0){
			 		# print "$1\t$2\t";
					$act_bases .= substr($2, $1);
			 	}
			}
			else{
				$act_bases .= $qqq;
			}
		}
	}

	return $act_bases;
}

#----------------------------------------------------------------------------
# 
#----------------------------------------------------------------------------
sub usage() {
	my($msg) = @_;
	if($msg ne ""){
		print STDERR $msg . "\n";
	}
	print STDERR "usage: $Command <file1> <file2>\n";  
	exit; 
}


