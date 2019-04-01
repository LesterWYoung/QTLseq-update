#!/usr/bin/env perl

use strict;
#use File::Basename;

my $Command = `basename $0`;
my $arg_cnt = scalar(@ARGV);


if ($arg_cnt < 3) {
	print STDERR "usage: $Command read1 read2 output_dir\n";
	exit;
}

#set_input_file
my $input_file1 = "";
my $input_file2 = "";
my $output_dir = "";

$input_file1 = shift;
$input_file2 = shift;
$output_dir = shift;


#set variable
my $input_basename ="";

my %index_id = ();
my $output_file1 = "";
my $output_file2 = "";
my $output_fileS = "";

#set_ 1st_file ID_array	

my $id_fomat = &set_id_array($input_file1);



#set_output_file
$input_basename = `basename $input_file1 .1.gz`;
chomp $input_basename;

$output_file1 = $output_dir."/".$input_basename."S.1.gz";
$output_file2 = $output_dir."/".$input_basename."S.2.gz";
$output_fileS = $output_dir."/".$input_basename."S.s.gz";

#input_file open
open(INFILE1, "gunzip -q -c $input_file1 | ") or die "No $input_file1";
open(INFILE2, "gunzip -q -c $input_file2 | ") or die "No $input_file2";

#output_file open
open(OUTFILE1, " | gzip -q -f - > $output_file1") or die "No $output_file1";
open(OUTFILE2, " | gzip -q -f - > $output_file2") or die "No $output_file2";
open(OUTFILES, " | gzip -q -f - > $output_fileS") or die "No $output_fileS";


#output file1 file2 fileS

&select_pair($id_fomat) ;


#input_file close
close(INFILE1);
close(INFILE2);

#output_file close
close(OUTFILE1);
close(OUTFILE2);
close(OUTFILES);


#----------------------------------------------------------------------------
# 1st_file ID_array set
#----------------------------------------------------------------------------
sub set_id_array() {
	my($input_file1) = @_;
	my $i = 0;
	my $id_format1 = "";

	#set id_format
	$id_format1 = &check_which_kind_of_id($input_file1);
	
	print "id_format=$id_format1";
	
	open(INFILE1, "gunzip -q -c $input_file1 | ") or die "No $input_file1";
	
	print STDERR "start scan $input_file1\n";

	while (<INFILE1>) {	
		chomp;
		if ( $i % 4 == 0) {	
			my @item = split(/$id_format1/,$_);
			$index_id{$item[0]} = $i;
			
		}
		$i++;
	}	
	close(INFILE1);
	
	print STDERR "id_array finished! $input_file1\n";
	
	return $id_format1;
}

#----------------------------------------------------------------------------
# 1st_file 2nd_file pair select
#----------------------------------------------------------------------------
sub select_pair() {
	
	my $id_format = $_;

	my $file1_line = "";
	my $file2_line = "";
	my @out_file1_line = ();
	my @out_file2_line = ();

	my $id1_cnt = 0;
	my @item = ();
	my $i = 4;

	print STDERR "start scan $input_file2\n";
	
	while ( $file2_line = <INFILE2> ) {

		if ( $i % 4 == 0) {	
			@item = split(/$id_fomat/, $file2_line);
			
			#file2_id = file1_id
			if (exists $index_id{$item[0]}) {
				
				#file1 read
				if (!eof INFILE1) {

					while ($index_id{$item[0]}+3 >= $id1_cnt) {
						
						$file1_line = <INFILE1>;
					
			#output out_fileS  <-- file1 	
						if ($index_id{$item[0]} == $id1_cnt) {
							print OUTFILES @out_file1_line;
							@out_file1_line = ();
						}
					
			#array file1		
						push (@out_file1_line,$file1_line);
						
			#output out_file1   <-- file1 
						if($index_id{$item[0]}+3 == $id1_cnt){
							print OUTFILE1 @out_file1_line;
							@out_file1_line = ();				
						}
						$id1_cnt++;
					}
				}
				
			#output out_fileS  <-- file2 	
				print OUTFILES @out_file2_line;
				@out_file2_line = ();
				$i = 0;
				
			}			
		}
		#array file1		
		push (@out_file2_line,$file2_line);


		#output out_file2 <-- file2
		if($i == 3){
			print OUTFILE2 @out_file2_line;
			@out_file2_line = ();
		}
		$i++;
		
	}

	#output out_files from the rest of out_file2
	print OUTFILES @out_file2_line;

	#output out_files from the rest of out_file1
	while($file1_line = <INFILE1>){
		print OUTFILES $file1_line;
	}

	print STDERR "scan finished! $input_file2\n";

}

# ==================================================
# 
# ==================================================
sub check_which_kind_of_id(){
	my($infastq1) = @_;
	open (my $fh1, "gunzip -q -c $infastq1 |") or die "cannot open $infastq1\n";
	my $oneline1 = <$fh1>;
	chomp($oneline1);
	
	my $id_format2 = &which_kind_of_id($oneline1);
	
	close $fh1;
	
	return $id_format2;
}

# ==================================================
# 
# ==================================================
sub which_kind_of_id(){
	my $retval = "";
	my ($oneline) = @_;
	
	my(@readid1) = split(/\s/,$oneline);
	my $num_array = @readid1;
	
	if($num_array == 1){
		$retval = "/"	# "#";
		
	}else{
		$retval = "\\s";	
	}
	return $retval;
}
