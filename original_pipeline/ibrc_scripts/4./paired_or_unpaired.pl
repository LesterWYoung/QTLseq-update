#! /usr/bin/perl -w

use strict;
#use File::Basename;

my $Command = `basename $0`;
my $arg_cnt = scalar(@ARGV);


if ($arg_cnt < 6) {
	print STDERR "usage: $Command A_SNP B_SNP A_SNP_PAIRED B_SNP_PAIRED A_SNP_UNPAIRED B_SNP_UNPAIRED\n";
	exit;
}

#set_input_file
my $input_file1 = shift;
my $input_file2 = shift;
my $output_file1 = shift;
my $output_file2 = shift;
my $output_file1_un = shift;
my $output_file2_un = shift;



#set variable
#my $input_basename ="";

my %index_id = ();

#set_ 1st_file ID_array	
&set_id_array($input_file1);
#my $id_fomat = &set_id_array($input_file1);


#set_output_file
#$input_basename = `basename $input_file1 .1.gz`;
#chomp $input_basename;


#input_file open
open(INFILE1, "$input_file1") or die "No $input_file1";
open(INFILE2, "$input_file2") or die "No $input_file2";

#output_file open
open(OUTFILE1, " > $output_file1") or die "No $output_file1";
open(OUTFILE2, " > $output_file2") or die "No $output_file2";
open(OUTFILE1_UNPAIR, " > $output_file1_un") or die "No $output_file1_un";
open(OUTFILE2_UNPAIR, "  > $output_file2_un") or die "No $output_file2_un";


#output file1 file2 fileS
#&select_pair($id_fomat) ;

#output file1_pair file2_pair file1_unpair file2_unpair
&select_pair() ;


#input_file close
close(INFILE1);
close(INFILE2);

#output_file close
close(OUTFILE1);
close(OUTFILE2);
close(OUTFILE1_UNPAIR);
close(OUTFILE2_UNPAIR);


#----------------------------------------------------------------------------
# 1st_file ID_array set
#----------------------------------------------------------------------------
sub set_id_array() {
	my($input_file1) = @_;
	my $i = 0;
	my $id_format1 = "";

	#set id_format
#	$id_format1 = &check_which_kind_of_id($input_file1);
#	print "id_format=$id_format1";
	
	open(INFILE1, "$input_file1") or die "No $input_file1";
	
	print STDERR "start scan $input_file1\n";

	while (<INFILE1>) {	
		chomp;
#		if ( $i % 4 == 0) {	
			my @item = split(/\t/,$_);
			$index_id{$item[0]."____".$item[1]} = $i;
			
#		}
		$i++;
	}	
	close(INFILE1);
	
	print STDERR "id_array finished! $input_file1\n";
	
#	return $id_format1;
}

#----------------------------------------------------------------------------
# 1st_file 2nd_file pair select
#----------------------------------------------------------------------------
sub select_pair() {
	
#	my $id_format = $_;

	my $file1_line = "";
	my $file2_line = "";
	my @out_file1_line = ();
#	my @out_file2_line = ();

	my $id1_cnt = 0;
	my @item = ();
	my $key = "";
	
	print STDERR "start scan $input_file2\n";
	
	while ( $file2_line = <INFILE2> ) {

		@item = split(/\t/, $file2_line);
		$key = $item[0]."____".$item[1];
			
		#file2_id = file1_id
		if (exists $index_id{$key}) {
				
				#file1 read
			if (!eof INFILE1) {

				while ($index_id{$key} >= $id1_cnt) {
					
					$file1_line = <INFILE1>;
				
		#output out_file1_unpair  <-- file1 	
					if ($index_id{$key} == $id1_cnt) {
						print OUTFILE1_UNPAIR @out_file1_line;
						@out_file1_line = ();
					}
					
		#array file1		
					push (@out_file1_line,$file1_line);
					
		#output out_file1_pair   <-- file1 
					if($index_id{$key} == $id1_cnt){
						print OUTFILE1 @out_file1_line;
						@out_file1_line = ();				
					}
					$id1_cnt++;
				}
			}
				
		#output out_file2_unpair  <-- file2 	
			print OUTFILE2 $file2_line;
			delete $index_id{$key};
			next;
		}			

		#output out_file2_pair <-- file2
		print OUTFILE2_UNPAIR $file2_line;
		
	}

#	#output out_files from the rest of out_file2
#	print OUTFILE2_UNPAIR $file2_line;

	#output out_files from the rest of out_file1
	while($file1_line = <INFILE1>){
		print OUTFILE1_UNPAIR $file1_line;
	}
	print STDERR "scan finished! $input_file2\n";
}
