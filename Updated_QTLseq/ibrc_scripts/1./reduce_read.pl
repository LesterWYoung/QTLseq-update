#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
my $Command = basename($0);

use Math::Random::MT::Auto 'rand';

# ----------------------------------------
# <参考> 'rand'以外も使う場合
# ----------------------------------------
# use Math::Random::MT::Auto
# 	qw(rand irand shuffle gaussian
# 		exponential erlang poisson binomial
# 		srand get_seed set_seed get_state set_state),
# 		'/dev/urandom' => 256,
# 		'random_org';
# 
# # http://search.cpan.org/~jdhedden/Math-Random-MT-Auto-6.14/lib/Math/Random/MT/Auto.pm
# # http://www.fifthdimension.jp/wiki.cgi?page=Math%3A%3ARandom%3A%3AMT%3A%3AAuto
# # 32bit整数を返すirand()
# # リスト・配列の並び替えを行うshuffle()
# # 正規乱数gaussian()
# # 指数乱数exponential()
# # ポアソン乱数poisson()
# # 二項乱数binomial()
 
# ----------------------------------------
# <参考> 四捨五入
# ----------------------------------------
# # http://d.hatena.ne.jp/end0tknr/20080928/1222581535
# use Math::Round qw/nearest/;



# ==================================================
# Main routin
# ==================================================
if (scalar(@ARGV) == 0) { &usage(); }

my @InputFile = ();
$InputFile[0] = $ARGV[0];
$InputFile[1] = $ARGV[1];
my $ReduceNum = $ARGV[2];			# 減らすべき到達リード数


my $OutFile_Suffix = "";
#-if ($ReduceNum) {			#--2014.05.09 kikuchi
if ($ReduceNum ne "") {
	$ReduceNum =~ s/,//g;	# 1,345,678,901
	$OutFile_Suffix .= "_" if ($OutFile_Suffix ne "");
	$OutFile_Suffix .= "reduced-$ReduceNum";
}
$OutFile_Suffix .= ".gz";

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) count fastq records
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

my $All_Fastq_Cnt = 0;

if ($ReduceNum) {
	print STDERR "Counting input file line ...\n";
	$All_Fastq_Cnt = &count_fastq_record($InputFile[0]);
	print STDERR "All_Fastq_Cnt=$All_Fastq_Cnt\n";
#-	if ($ReduceNum >= $All_Fastq_Cnt) {				#--2014.05.09 kikuchi
	if ($ReduceNum > $All_Fastq_Cnt) {
#-		&usage("ReduceNum >= All_Fastq_Cnt invalid value, $ReduceNum >= $All_Fastq_Cnt");
		&usage("ReduceNum > All_Fastq_Cnt invalid value, $ReduceNum > $All_Fastq_Cnt");
	}
}

print STDERR "All_Fastq_Cnt: $All_Fastq_Cnt\n";
print STDERR "ReduceNum: $ReduceNum\n";
print STDERR "OutFile_Suffix: $OutFile_Suffix\n";


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2) if reduce read, make random reduce line data.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
my %Random_WriteOut_ReadNo = ();	# set random line no
my $Reduce_Mode = "writeout";

if ($ReduceNum) {
	print STDERR "\nStart Random point get .. \n";

	my $pickup_cnt = $ReduceNum;
	my $reduce_rate = $ReduceNum/$All_Fastq_Cnt;

	print STDERR "reduce_rate=ReduceNum/All_Fastq_Cnt= ";
	print STDERR "$ReduceNum/$All_Fastq_Cnt= ";
	print STDERR $ReduceNum/$All_Fastq_Cnt . "\n";

	if ($reduce_rate > 0.5) {
		$Reduce_Mode = "not_writeout";
		$pickup_cnt = $All_Fastq_Cnt - $ReduceNum;
	}
	print STDERR "Reduce_Mode=$Reduce_Mode, pickup_cnt=$pickup_cnt\n";

	&set_ramdom_lineno($All_Fastq_Cnt, $pickup_cnt);

	my $get_cnt = scalar(keys(%Random_WriteOut_ReadNo));
	print STDERR "Done. $get_cnt\n";
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3) write out paired-end read
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
&output_reduced_read(\@InputFile, $OutFile_Suffix);


# ==================================================
# End of Main routin
# ==================================================


# ==================================================
# get line counts --> read counts
# ==================================================
sub count_fastq_record() {
	my ($fastq_file) = @_;
	my $wc_line = "";
	if ($fastq_file =~ /\.gz$/) {
		$wc_line = `gunzip -c $fastq_file | wc -l`;
	} else {
		$wc_line = `wc -l $fastq_file`;
	}
	my ($fastq_linecnt) = split(/\s+/, $wc_line);
	my ($fastq_reccnt) = $fastq_linecnt / 4;
	$fastq_reccnt = int($fastq_reccnt);

	return $fastq_reccnt;
}

# ==================================================
# get random number
# 	$all_fastq_cnt (integer)	fastq's all line cnt
# 	$reduce_cnt (integer)		target cnt to reduce number of rec
# ==================================================
sub set_ramdom_lineno() {
	my ($all_fastq_cnt, $reduce_cnt) = @_;

	# 0からはじめて、reduce_cntを超えたら終わり
	for (my $remain_cnt = 0; $remain_cnt < $reduce_cnt; $remain_cnt++) {

		# greater than or equal to 0 and less than the value of EXPR
		# recno = start at 1
   		my $remain_lno = int(rand($all_fastq_cnt)) + 1;
		if (defined($Random_WriteOut_ReadNo{$remain_lno}) and
			$Random_WriteOut_ReadNo{$remain_lno} == $remain_lno) {
			# more (Duplicate $remain_lno);
			&set_ramdom_lineno($all_fastq_cnt, 1);
		}
		else {
			$Random_WriteOut_ReadNo{$remain_lno} = $remain_lno;
		}
		#print "$remain_cnt<$reduce_cnt: remain_lno=$remain_lno ($all_fastq_cnt)\n";
	}
}

# ==================================================
# write to file
# ==================================================
sub output_reduced_read(){
	my ($p_fastq_files, $outfile_suffix) = @_;

	my $write_line = "";

	foreach my $fastq_file (@$p_fastq_files){
		my $fastqp;
		if($fastq_file =~ /\.gz$/){
			open ($fastqp, "gunzip -c $fastq_file | ") or die "cannot open $fastq_file\n";
		}
		else{
			open ($fastqp, "$fastq_file") or die "cannot open $fastq_file\n";
		}
		my $out_file = "$fastq_file.$outfile_suffix";
		open (my $outfilep, " | gzip -f - > $out_file") or die "cannot open $out_file\n";
		
		print STDERR "Reducing into out_file: $out_file\n";
		
		my $fastq_lcnt = 1;
		my $fastq_reccnt = 0;	# in fastq, 1 origin
		
		while (<$fastqp>) {
			chomp;
			if (/^$/) { next; }		# skip
			if ($fastq_lcnt > 4) { $fastq_lcnt = 1; }
			if ($fastq_lcnt == 1) {
				if (/^\@/) {
					$fastq_reccnt++;
					$write_line = "$_\n";
				}
				else {
					die "invalid fastq $fastq_file";
				}	
			}
			elsif ($fastq_lcnt == 3) {
				$write_line .= "$_\n";
			}
			elsif ($fastq_lcnt == 4) {
				$write_line .= "$_\n";
				if ($ReduceNum) {
					if ($Reduce_Mode eq "writeout") {
						# if defined, I write
						# what means "defined"
						if (defined($Random_WriteOut_ReadNo{$fastq_reccnt}) and
							$Random_WriteOut_ReadNo{$fastq_reccnt} == $fastq_reccnt) {
							print $outfilep $write_line;
						} else {
							;
						}
					} else {
						# if not defined, I write
						if (!defined($Random_WriteOut_ReadNo{$fastq_reccnt})) {
							# print STDERR "!defined $fastq_reccnt\n";
							print $outfilep $write_line;
						} else {
							# print STDERR "\tdefined $fastq_reccnt\n";
							# if ($Random_WriteOut_ReadNo{$fastq_reccnt} != $fastq_reccnt) {
							# 	print STDERR "!= fastq_reccnt $fastq_reccnt\n";
							# 	print $outfilep $write_line;
							# }
						}
					}	
				}
				else{
					print $outfilep $write_line;	
				}
			}
			else {
				$write_line .= "$_\n";
			}
			$fastq_lcnt++;
			
		}
		
		close($fastqp);
		close($outfilep);
		print STDERR "Done.\n";
	}
}

# ==================================================
# usage
# ==================================================
sub usage() {
	my($msg) = @_;          
	if ($msg) { print STDERR "$msg\n"; }

	print STDERR "\n";
	print STDERR "Usage: $Command <fastq1.gz> <fastq2.gz> <num>\n";
	print STDERR "    <fastq1.gz> <fastq2.gz> : paired-end\n";
	print STDERR "    <num>                   : reduce number\n";
	print STDERR "\n";
    
    exit(0);
}
