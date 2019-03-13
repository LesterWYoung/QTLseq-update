#! /usr/bin/perl
use File::Basename;

use strict;
my $Command = basename $0;

if (scalar(@ARGV) < 4) {
	&usage();
}

my $pileup_file1 = shift; 	# 1st argment
my $bam_file2 = shift; 		# 2nd argment
my $ref_file2 = shift; 		# 3rd argment
my $out_path = shift;		# 4th argment 2013/03/22
## my $pileup_out = shift;

# my $pileup_file2; # 2nd argment


if (!-e $pileup_file1 || !-e $bam_file2 || !-e $ref_file2) {  &usage("invalid file"); }


#----------------------------------------------------------------------------
# options
#----------------------------------------------------------------------------
my $sw9 = 0; # show_screen;					must be 0 if &comp_snp_to_pileup()
my $sw0 = 0; # dump to sidebyside_file;		must be 0 if &comp_snp_to_pileup()
my $sw1 = 1; # dump to common_pos_file1;	usually 1
my $sw2 = 1; # dump to common_pos_file2;	usually 1
my $sw3 = 1; # dump to unique_file1;		usually 1
my $sw4 = 0; # dump to unique_file2;		usually 0 if pileup_file2 is 'full-pileup'


# &scan_pileups($pileup_file1, $pileup_file2);
# &comp_snp_to_pileup($pileup_file1, $pileup_file2);
&comp_snp_to_bam($pileup_file1, $bam_file2, $ref_file2);

#----------------------------------------------------------------------------
# comp_snp_to_pileup()の発展形
# 
#----------------------------------------------------------------------------
sub comp_snp_to_bam() {
	my($pileup_file1, $bam_file2, $ref_file2) = @_;

	my $oneline1;
	my $oneline2;
	my @items1 = ();
 	my @items2 = ();

	my %hash1 = ();
	my $chr_pos;

	my $nullcolumn = "\t\t\t\t\t\t\t\t\t\t";


	my $pileup_name1 = `basename ${pileup_file1}`;
	chomp $pileup_name1;
	
	my $pileup_name2 = `basename $bam_file2`;
	chomp $pileup_name2;
	$pileup_name2 = "$pileup_name2.pileup";

	my @output = ();
	$output[0] = "${out_path}/sidebyside_${pileup_name1}_${pileup_name2}";
	$output[1] = "${out_path}/common_pos_${pileup_name1}";
	$output[2] = "${out_path}/common_pos_${pileup_name2}";
	$output[3] = "${out_path}/unique_${pileup_name1}";
	$output[4] = "${out_path}/unique_${pileup_name2}";
	
	my $output_3_ = $output[3];
	## $output[3] = "unsorted_".$output[3];
	$output[3] = "${out_path}/unsorted_unique_${pileup_name1}";

	#----------------------------------------------------------------------------
	# load pileup1 ---> メモリーに格納%hash1
	#----------------------------------------------------------------------------
	open(PILEUP1, $pileup_file1) or die "No $pileup_file1";
	while (<PILEUP1>) {
		$oneline1 = $_;
		chomp $oneline1;
		@items1 = split(/\t/, $oneline1);
		$chr_pos = $items1[0]."_".$items1[1];
		$hash1{$chr_pos} = $oneline1;
	}
	close(PILEUP1);


	open(OUT0, ">$output[0]") or die "can not open $output[0]";
	open(OUT1, ">$output[1]") or die "can not open $output[1]";
	open(OUT2, ">$output[2]") or die "can not open $output[2]";  
	open(OUT3, ">$output[3]") or die "can not open $output[3]";
	open(OUT4, ">$output[4]") or die "can not open $output[4]";


	#----------------------------------------------------------------------------
	# load pileup2 ---> 逐次実行
	#----------------------------------------------------------------------------
	## open(PILEUP2, $pileup_file2) or die "No $pileup_file2";
	open(PILEUP2, "samtools pileup -c -f $ref_file2 $bam_file2 | ") or die "No $bam_file2, $ref_file2";	
	
 	while (<PILEUP2>) {
		$oneline2 = $_;
		chomp $oneline2;
		@items2 = split(/\t/, $oneline2);
		$chr_pos = $items2[0]."_".$items2[1];

		if(exists($hash1{$chr_pos})){
			$oneline1 = $hash1{$chr_pos};
			if($sw9 == 1){ print "1&2\t$oneline1\t$oneline2\n"; }
			if($sw0 == 1){ print OUT0 "$oneline1\t$oneline2\n"; }
			if($sw9 == 1){ print "1\'\t$oneline1\n"; }
			if($sw1 == 1){ print OUT1 "$oneline1\n"; }
			if($sw9 == 1){ print "2\'\t$oneline2\n"; }
			if($sw2 == 1){ print OUT2 "$oneline2\n"; }
			
			delete $hash1{$chr_pos};
		}
		else{
			if($sw9 == 1){ print   "2\t$oneline2\n"; }
			if($sw0 == 1){ print OUT0 "$nullcolumn\t$oneline2\n"; }
			if($sw9 == 1){ print   "2\t$oneline2\n"; }
			if($sw4 == 1){ print OUT4 "$oneline2\n"; }			
		}
	}

	foreach my $key ( keys( %hash1 ) ) {
		$oneline1 = $hash1{$key};
		if($sw9 == 1){ print   "1\t$oneline1\n"; }
		if($sw0 == 1){ print OUT0 "$oneline1\t$nullcolumn\n"; }
		if($sw9 == 1){ print   "1\t$oneline1\n"; }
		if($sw3 == 1){ print OUT3 "$oneline1\n"; }		
	}

	close(OUT0);
	close(OUT1);
	close(OUT2);
	close(OUT3);
	close(OUT4);

	close(PILEUP2);
	
	if($sw3 == 1){
		&my_sort($output[3], $output_3_);
	}	
	
}



#----------------------------------------------------------------------------
# メモリー節約型
# pileup1または2の一方がSNPで、もう一方はフルpileupファイルの場合にはこちら
#----------------------------------------------------------------------------
sub comp_snp_to_pileup() {
	my($pileup_file1, $pileup_file2) = @_;

	my $oneline1;
	my $oneline2;
	my @items1 = ();
 	my @items2 = ();

	my %hash1 = ();
	my $chr_pos;

	my $nullcolumn = "\t\t\t\t\t\t\t\t\t\t";

	my @output = ();
	$output[0] = "sidebyside_${pileup_file1}_${pileup_file2}";
	$output[1] = "common_pos_$pileup_file1";
	$output[2] = "common_pos_$pileup_file2";
	$output[3] = "unique_$pileup_file1";
	$output[4] = "unique_$pileup_file2";
	
	my $output_3_ = $output[3];
	$output[3] = "unsorted_".$output[3];


	#----------------------------------------------------------------------------
	# load pileup1 ---> メモリーに格納%hash1
	#----------------------------------------------------------------------------
	open(PILEUP1, $pileup_file1) or die "No $pileup_file1";
	while (<PILEUP1>) {
		$oneline1 = $_;
		chomp $oneline1;
		@items1 = split(/\t/, $oneline1);
		$chr_pos = $items1[0]."_".$items1[1];
		$hash1{$chr_pos} = $oneline1;
	}
	close(PILEUP1);


	open(OUT0, ">$output[0]") or die "can not open $output[0]";
	open(OUT1, ">$output[1]") or die "can not open $output[1]";
	open(OUT2, ">$output[2]") or die "can not open $output[2]";  
	open(OUT3, ">$output[3]") or die "can not open $output[3]";
	open(OUT4, ">$output[4]") or die "can not open $output[4]";


	#----------------------------------------------------------------------------
	# load pileup2 ---> 逐次実行
	#----------------------------------------------------------------------------
	open(PILEUP2, $pileup_file2) or die "No $pileup_file2";
 	while (<PILEUP2>) {
		$oneline2 = $_;
		chomp $oneline2;
		@items2 = split(/\t/, $oneline2);
		$chr_pos = $items2[0]."_".$items2[1];

		if(exists($hash1{$chr_pos})){
			$oneline1 = $hash1{$chr_pos};
			if($sw9 == 1){ print "1&2\t$oneline1\t$oneline2\n"; }
			if($sw0 == 1){ print OUT0 "$oneline1\t$oneline2\n"; }
			if($sw9 == 1){ print "1\'\t$oneline1\n"; }
			if($sw1 == 1){ print OUT1 "$oneline1\n"; }
			if($sw9 == 1){ print "2\'\t$oneline2\n"; }
			if($sw2 == 1){ print OUT2 "$oneline2\n"; }
			
			delete $hash1{$chr_pos};
		}
		else{
			if($sw9 == 1){ print   "2\t$oneline2\n"; }
			if($sw0 == 1){ print OUT0 "$nullcolumn\t$oneline2\n"; }
			if($sw9 == 1){ print   "2\t$oneline2\n"; }
			if($sw4 == 1){ print OUT4 "$oneline2\n"; }			
		}
	}

	foreach my $key ( keys( %hash1 ) ) {
		$oneline1 = $hash1{$key};
		if($sw9 == 1){ print   "1\t$oneline1\n"; }
		if($sw0 == 1){ print OUT0 "$oneline1\t$nullcolumn\n"; }
		if($sw9 == 1){ print   "1\t$oneline1\n"; }
		if($sw3 == 1){ print OUT3 "$oneline1\n"; }		
	}

	close(OUT0);
	close(OUT1);
	close(OUT2);
	close(OUT3);
	close(OUT4);

	close(PILEUP2);
	
	if($sw3 == 1){
		&my_sort($output[3], $output_3_);
	}	
	
}


#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
sub my_sort(){
	my ($unsort, $sorted) = @_;

	my $command = "sort -k1,1 -k2,2n ".$unsort." > ".$sorted;

	system($command);
	
	system("rm ".$unsort)

}


#----------------------------------------------------------------------------
# 逐次実行型
# pileupファイル両方ともSNPの場合にはこちら
#----------------------------------------------------------------------------
sub scan_pileups() {
	my($pileup_file1, $pileup_file2) = @_;

	open(PILEUP1, $pileup_file1) or die "No $pileup_file1";  my @arr1 = <PILEUP1>;  close(PILEUP1);

	open(PILEUP2, $pileup_file2) or die "No $pileup_file2";  my @arr2 = <PILEUP2>;  close(PILEUP2);

	my $nullcolumn = "\t\t\t\t\t\t\t\t\t\t";

	my @output = ();
	$output[0] = "sidebyside_${pileup_file1}_${pileup_file2}";
	$output[1] = "common_pos_$pileup_file1";
	$output[2] = "common_pos_$pileup_file2";
	$output[3] = "unique_$pileup_file1";
	$output[4] = "unique_$pileup_file2";

	open(OUT0, ">$output[0]") or die "can not open $output[0]";
	open(OUT1, ">$output[1]") or die "can not open $output[1]";
	open(OUT2, ">$output[2]") or die "can not open $output[2]";  
	open(OUT3, ">$output[3]") or die "can not open $output[3]";
	open(OUT4, ">$output[4]") or die "can not open $output[4]";

 	my @items1 = ();
 	my @items2 = ();

 	my $i = 0;
 	my $j = 0;

 	do{
		if($i < @arr1 && $j < @arr2){
			chomp $arr1[$i];
			@items1 = split(/\t/, $arr1[$i]);

			chomp $arr2[$j];
			@items2 = split(/\t/, $arr2[$j]);

			if($items1[0] eq $items2[0]){
				if($items1[1] == $items2[1]){
					if($sw9 == 1){ print "1&2\t$arr1[$i]\t$arr2[$j]\n"; }
					if($sw0 == 1){ print OUT0 "$arr1[$i]\t$arr2[$j]\n"; }
					if($sw9 == 1){ print "1\'\t$arr1[$i]\n"; }
					if($sw1 == 1){ print OUT1 "$arr1[$i]\n"; }
					if($sw9 == 1){ print "2\'\t$arr2[$j]\n"; }
					if($sw2 == 1){ print OUT2 "$arr2[$j]\n"; }
					
					$i++;
					$j++;
				}
				elsif($items1[1] < $items2[1]){
					if($sw9 == 1){ print   "1\t$arr1[$i]\n"; }
					if($sw0 == 1){ print OUT0 "$arr1[$i]\t$nullcolumn\n"; }
					if($sw9 == 1){ print   "1\t$arr1[$i]\n"; }
					if($sw3 == 1){ print OUT3 "$arr1[$i]\n"; }
					
					$i++;
				}
				else{
					if($sw9 == 1){ print   "2\t$arr2[$j]\n"; }
					if($sw0 == 1){ print OUT0 "$nullcolumn\t$arr2[$j]\n"; }
					if($sw9 == 1){ print   "2\t$arr2[$j]\n"; }
					if($sw4 == 1){ print OUT4 "$arr2[$j]\n"; }
					
					$j++;
				}
			}
			elsif($items1[0] lt $items2[0]){
				if($sw9 == 1){ print   "1\t$arr1[$i]\n"; }
				if($sw0 == 1){ print OUT0 "$arr1[$i]\t$nullcolumn\n"; }
				if($sw9 == 1){ print   "1\t$arr1[$i]\n"; }
				if($sw3 == 1){ print OUT3 "$arr1[$i]\n"; }
				
				$i++;
			}
			else{
				if($sw9 == 1){ print   "2\t$arr2[$j]\n"; }
				if($sw0 == 1){ print OUT0 "$nullcolumn\t$arr2[$j]\n"; }
				if($sw9 == 1){ print   "2\t$arr2[$j]\n"; }
				if($sw4 == 1){ print OUT4 "$arr2[$j]\n"; }
				
				$j++;
			}
		}

		elsif($i == @arr1 && $j < @arr2){
			chomp $arr2[$j];
			
			if($sw9 == 1){ print   "2\t$arr2[$j]\n"; }
			if($sw0 == 1){ print OUT0 "$nullcolumn\t$arr2[$j]\n"; }
			if($sw9 == 1){ print   "2\t$arr2[$j]\n"; }
			if($sw4 == 1){ print OUT4 "$arr2[$j]\n"; }
			
			$j++;
		}
		elsif($i < @arr1 && $j == @arr2){
			chomp $arr1[$i];
			
			if($sw9 == 1){ print   "1\t$arr1[$i]\n"; }
			if($sw0 == 1){ print OUT0 "$arr1[$i]\t$nullcolumn\n"; }
			if($sw9 == 1){ print   "1\t$arr1[$i]\n"; }
			if($sw3 == 1){ print OUT3 "$arr1[$i]\n"; }
			
			$i++;
		}

	}while($i < @arr1 || $j < @arr2);

	close(OUT0);
	close(OUT1);
	close(OUT2);
	close(OUT3);
	close(OUT4);

}

#	#----------------------------------------------------------------------------
#	
#	#----------------------------------------------------------------------------
#	sub usage() {
#		my($msg) = @_;
#		if($msg ne ""){
#			print STDERR $msg . "\n";
#		}
#		print STDERR "usage: $Command <file1> <file2>\n";  
#		exit; 
#	}

#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
sub usage() {
	my($msg) = @_;
	if($msg ne ""){
		print STDERR $msg . "\n";
	}
	print STDERR "usage: $Command <snp_pileup1> <bam> <reference> <outpath>\n";  
	exit; 
}

