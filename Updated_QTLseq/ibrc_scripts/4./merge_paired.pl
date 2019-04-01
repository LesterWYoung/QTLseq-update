#! /usr/bin/perl -w

use strict;
#use File::Basename;

my $Command = `basename $0`;
my $arg_cnt = scalar(@ARGV);


if ($arg_cnt < 5) {
	print STDERR "usage: $Command REF_FAI PAIRED_AB UNPAIRED_A UNPAIRED_B MERGED_PAIRED_AB\n";
	exit;
}

#set_input_file
my $ref_fai = shift;						# REF_FAI
my $input_file = shift;						# PAIRED_AB
my $input_file1 = shift;					# UNPAIRED_A
my $input_file2 = shift;					# UNPAIRED_B
my $output_file = shift;					# MERGED_PAIRED_AB


print STDERR "input_file $input_file\n";
print STDERR "input_file1 $input_file1\n";
print STDERR "input_file2 $input_file2\n";




#input_file open
open(INREF, "$ref_fai") or die "No $ref_fai";
open(INFILE, "$input_file") or die "No $input_file";
open(INFILE1, "$input_file1") or die "No $input_file1";
open(INFILE2, "$input_file2") or die "No $input_file2";

#output_file open
open(OUTFILE, " > $output_file") or die "No $output_file";


# PAIRED_AB UNPAIRED_A UNPAIRED_B -> MERGED_PAIRED_AB
&sortfile() ;


#input_file close
close(INREF);
close(INFILE);
close(INFILE1);
close(INFILE2);

#output_file close
close(OUTFILE);


#----------------------------------------------------------------------------
# sortfile
#----------------------------------------------------------------------------
sub sortfile() {
	
	my %index_id = ();

	my $ref_line = "";
	my $file_line = "";
	my $file1_line = "";
	my $file2_line = "";
	my @out_file1_line = ();
	my @filename = ();
	
	my $id1_cnt = 0;
	my @item0 = ();
	my @item = ();
	my @item1 = ();
	my @item2 = ();
	my $key = "";
	my $max = "";
	
	print STDERR "start scan $ref_line\n";
	

	while ( $ref_line = <INREF>) {

		
		@item0 = split(/\t/, $ref_line);
		$key = $item0[0];
		$max = $item0[1];
		
		print STDERR "key $item0[0]\n";

		until ( eof INFILE && eof INFILE1 && eof INFILE2 && !%index_id){

			# PAIRED_AB read
			if (!exists $index_id{infile}) {
				if (!$item[0] && !eof INFILE){
					$file_line = <INFILE>;				
					@item = split(/\t/, $file_line);
				}

				if ($item[0] && $item[0] eq $key) {
					$index_id{infile} = $item[1];
					@item = ();
				}else{
					$index_id{infile} = $max;
				}
			}

			# UNPAIRED_A read
			if (!exists $index_id{infile1}) {
				if (!$item1[0] && !eof INFILE1){
					$file1_line = <INFILE1>;
					@item1 = split(/\t/, $file1_line);
				}
				if ($item1[0] && $item1[0] eq $key) {
					$index_id{infile1} = $item1[1];
					@item1 = ();
				}else{
					$index_id{infile1} = $max;
				}
			}
			
			# UNPAIRED_B read
			if (!exists $index_id{infile2}) {
				if (!$item2[0] && !eof INFILE2){
					$file2_line = <INFILE2>;
					@item2 = split(/\t/, $file2_line);
				}
				if ($item2[0] && $item2[0] eq $key) {
					$index_id{infile2} = $item2[1];
					@item2 = ();
				}else{
					$index_id{infile2} = $max;
				}
			}
		
			# min_valu sort
			@filename = sort {$index_id{$a} <=> $index_id{$b}} keys%index_id;
	
			# min_valu == $max -> nex_id read
			if ($index_id{$filename[0]} == $max){
				%index_id = ();
				last;
			}	

			# min_valu_file output
			if ($filename[0] eq "infile"){
				print OUTFILE $file_line;
				$file_line = "";
			}elsif($filename[0] eq "infile1"){
				print OUTFILE $file1_line;			
				$file1_line = "";
			}elsif($filename[0] eq "infile2"){			
				print OUTFILE $file2_line;			
				$file2_line = "";	
			}
			
			# min_valu del
			delete($index_id{$filename[0]});
		}
			
	}
	

	print STDERR "scan finished!\n";
}
