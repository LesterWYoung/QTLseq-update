# ##################################################			
# global environment			
# ##################################################			
Set_TOPPATH_SCRIPTS(){			
	TOPPATH_SCRIPTS="./ibrc_scripts"		
	echo "${TOPPATH_SCRIPTS}"		
}			
			
Set_TOPPATH_COVAL(){			
	TOPPATH_COVAL="./execute_Coval/Coval-1.4"		
	echo "${TOPPATH_COVAL}"		
}			
			
# ##################################################			
# for 1.qualify_read			
# ##################################################			
# ==================================================			
# Bat_fastq_quality_filter.sh			
# ==================================================			
Set_READ_QVAL(){			
	READ_QVAL=30		
	echo "${READ_QVAL}"		
}			
			
Set_READ_PVAL(){			
	READ_PVAL=90		
	echo "${READ_PVAL}"		
}			
			
Set_READ_QOPT(){			
	READ_QOPT=""		
	echo "${READ_QOPT}"		
			
	# READ_QOPT="-Q 33"	in case of CASAVA 1.8 later	
	# READ_QOPT=""		in case of the previous Illumina-specific offset value of 64
}			
#--                              #2014/04/08 kikuchi			
Set_READ_QVAL_MY_CULTIVAR(){			
	READ_QVAL=30		
	echo "${READ_QVAL}"		
}			
			
Set_READ_PVAL_MY_CULTIVAR(){			
	READ_PVAL=90		
	echo "${READ_PVAL}"		
}			
			
Set_READ_QOPT_MY_CULTIVAR(){			
	READ_QOPT=""		
	echo "${READ_QOPT}"		
			
	# READ_QOPT="-Q 33"		in case of CASAVA 1.8 later
	# READ_QOPT=""		in case of the previous Illumina-specific offset value of 64
}			
			
# ##################################################			
# for 2.make_consensus 			
# ##################################################			
# ==================================================			
# Bat_make_symbolic_link_of_qualified_fastq.sh.sh			
# ==================================================			
Set_SRC_READ_PATH_MY_CULTIVAR(){			
	SRC_READ_PATH_MY_CULTIVAR="/Users/LWY/Unix/QTLseq-update-local/QTLseq-update/original_pipeline/1.qualify_read/Adelie/q30p90/sep_pair"		
	echo "${SRC_READ_PATH_MY_CULTIVAR}"		
}			
			
Set_MY_CULTIVAR_NAME(){			
	MY_CULTIVAR_NAME="Adelie"		
	echo "${MY_CULTIVAR_NAME}"		
}			
			
# ==================================================			
# Bat_bwa2bam.sh			
# ==================================================			
Set_PUBLIC_REF_FASTA(){			
	REF_FASTA="/Users/LWY/Unix/QTLseq-update-local/QTLseq-update/original_pipeline/downloaded_fasta/INSERT_ref_seq_fasta_filename_here"		
	echo "${REF_FASTA}"		
}			
			
# ==================================================			
# Bat_run_coval-refine-bam.pl.sh			
# ==================================================			
Set_MIS_MATCH_FOR_MAKE_CONSENSUS(){			
	MIS_MATCH=10		
	echo "${MIS_MATCH}"		
}			
			
Set_DONTREALIGN_FOR_MAKE_CONSENSUS(){			
	DONTREALIGN="yes"		
	echo "${DONTREALIGN}"		
			
	# DONTREALIGN=""		
	# DONTREALIGN="--disalign"		
}			
			
# ==================================================			
# Bat_run_coval-call-pileup.pl.sh			
# ==================================================			
Set_COVAL_CALL_MINNUM_FOR_MAKE_CONSENSUS()			
{			
	COVALCALL_MINNUM=3		
	echo "${COVALCALL_MINNUM}"		
}			
			
Set_COVAL_CALL_MAXR_FOR_MAKE_CONSENSUS()			
{			
	COVALCALL_MAXR=10000		
	echo "${COVALCALL_MAXR}"		
}			
			
Set_COVAL_CALL_MINFREQ_FOR_MAKE_CONSENSUS()			
{			
	COVALCALL_MINFREQ=0.5		
	echo "${COVALCALL_MINFREQ}"		
}			
			
Set_COVAL_CALL_MINTNUM_FOR_MAKE_CONSENSUS()			
{			
	COVALCALL_MINTNUM=3		
	echo "${COVALCALL_MINTNUM}"		
}			
			
Set_COVAL_CALL_MINQUALBASE_FOR_MAKE_CONSENSUS()			
{			
	COVALCALL_MINQUALBASE=3		
	echo "${COVALCALL_MINQUALBASE}"		
}			
			
Set_COVAL_CALL_MINQUALAVE_FOR_MAKE_CONSENSUS()			
{			
	COVALCALL_MINQUALAVE=20		
	echo "${COVALCALL_MINQUALAVE}"		
}			
			
Set_COVAL_CALL_CALLTYPE_FOR_MAKE_CONSENSUS()			
{			
	COVALCALL_CALLTYPE="illumina"		
	echo "${COVALCALL_CALLTYPE}"		
			
	# COVALCALL_CALLTYPE="sanger"		
	# COVALCALL_CALLTYPE="illumina"		
}			
			
# ##################################################			
# for 3.alignment			
# ##################################################			
Set_BULK_NAME(){			
	BULK_NAME="mybulk"		
	echo "${BULK_NAME}"		
}			
			
Set_BULK_NAME_IDA(){			
	BULK_NAME_ID="A"		
	echo "${BULK_NAME_ID}"		
}			
			
Set_BULK_NAME_IDB(){			
	BULK_NAME_ID="B"		
	echo "${BULK_NAME_ID}"		
}			
			
Set_REF_FASTA()			
{			
	REF_FASTA="/Users/LWY/Unix/QTLseq-update-local/QTLseq-update/original_pipeline/2.make_consensus/90.align_to_this_fasta/00.reference/Adelie_q30p90_MSR_Cov_10_S-snp_RYKMSWBDHV2ACGT.fa"		
	echo "${REF_FASTA}"		
}			
			
# --------------------------------------------------			
# for Bat_make_symbolic_link_of_qualified_fastq.sh			
# --------------------------------------------------			
Set_SRC_READ_PATH(){			
	SRC_READ_PATH="/Users/LWY/Unix/QTLseq-update-local/QTLseq-update/original_pipeline/1.qualify_read/equalized"		
			
	echo "${SRC_READ_PATH}"		
}			
			
# --------------------------------------------------			
# for Bat_bwa2bam.sh			
# --------------------------------------------------			
Set_BWA_CPU()			
{			
	BWA_CPU=5		
	echo "${BWA_CPU}"		
}			
			
# --------------------------------------------------			
# Bat_run_coval-refine-bam.pl.sh			
# --------------------------------------------------			
Set_DONTREALIGN()			
{			
	DONTREALIGN=""		
	echo "${DONTREALIGN}"		
			
	# DONTREALIGN=""		
	# DONTREALIGN="--disalign"		
}			
			
Set_MIS_MATCH_DEFAULT()			
{			
	MIS_MATCH_DEFAULT=10		
	echo "${MIS_MATCH_DEFAULT}"		
}			
			
# --------------------------------------------------			
# Bat_run_coval-call-pileup.pl.sh			
# --------------------------------------------------			
Set_COVAL_CALL_MINNUM()			
{			
	COVALCALL_MINNUM=3		
	echo "${COVALCALL_MINNUM}"		
}			
			
Set_COVAL_CALL_MAXR()			
{			
	COVALCALL_MAXR=75		
	echo "${COVALCALL_MAXR}"		
}			
			
Set_COVAL_CALL_MINFREQ()			
{			
	COVALCALL_MINFREQ=0.0		
	echo "${COVALCALL_MINFREQ}"		
}			
			
Set_COVAL_CALL_MINTNUM()			
{			
	COVALCALL_MINTNUM=3		
	echo "${COVALCALL_MINTNUM}"		
}			
			
Set_COVAL_CALL_MINQUALBASE()			
{			
	COVALCALL_MINQUALBASE=3		
	echo "${COVALCALL_MINQUALBASE}"		
}			
			
Set_COVAL_CALL_MINQUALAVE()			
{			
	COVALCALL_MINQUALAVE=20		
	echo "${COVALCALL_MINQUALAVE}"		
}			
			
Set_COVAL_CALL_CALLTYPE()			
{			
	COVALCALL_CALLTYPE="illumina"		
	echo "${COVALCALL_CALLTYPE}"		
			
	# COVALCALL_CALLTYPE="sanger"		
	# COVALCALL_CALLTYPE="illumina"		
}			
			
# --------------------------------------------------			
# Bat_exclude_common_snps.pl.sh			
# --------------------------------------------------			
Set_MISs()			
{			
	MISs="2 3 4"		
	echo "${MISs}"		
}			
			
Set_PILEUPDB_PATH()			
{			
	PILEUPDB_PATH="/Users/LWY/Unix/QTLseq-update-local/QTLseq-update/original_pipeline/2.make_consensus/90.align_to_this_fasta/30.coval_call"		
	echo "${PILEUPDB_PATH}"		
}			
			
Set_PILEUPDB_NAME()			
{			
	PILEUPDB_NAME="Adelie_q30p90"		
	echo "${PILEUPDB_NAME}"		
}			
			
Set_PILEUPDB_MIS_FIXED()			
{			
	PILEUPDB_MIS_FIXED=10 # if 0, each mismatch value		
	echo "${PILEUPDB_MIS_FIXED}"		
}			
			
			
			
			
			
			
			
			
			
			
# --------------------------------------------------			
# Bat_awk_custom.sh			
# --------------------------------------------------			
			
Set_MINDEPTH()			
{			
	MINDEPTH=3		
	echo "${MINDEPTH}"		
}			
			
			
			
			
Set_MIN_CONSENSUS_QUALITY()			
{			
	MIN_CONSENSUS_QUALITY=20		
	echo "${MIN_CONSENSUS_QUALITY}"		
}			
			
Set_MIN_SNP_INDEX()			
{			
	MIN_SNP_INDEX=0.0		
	echo "${MIN_SNP_INDEX}"		
}			
			
# ##################################################			
# for 4.search_for_pair			
##################################################			
			
			
# ##################################################			
# for 5.compare			
# ##################################################			
# --------------------------------------------------			
# Bat_cbind_confidence_interval.sh			
# --------------------------------------------------			
			
Set_CONFINTRVL_MISs()			
{			
	MISs="2 3 4"		
	echo "${MISs}"		
}			
			
Set_CONFINTRVL_CALC_MODE()			
{			
	CALC_MODE=1	    # 1:DoCalculateNow  0:UsePreviousOne	
	echo "${CALC_MODE}"		
}			
			
Set_CONFINTRVL_INDIVIDUALS()			
{			
	INDIVIDUALS=20		
	echo "${INDIVIDUALS}"		
}			
			
Set_CONFINTRVL_NUM_OF_TRIALS()			
{			
	NUM_OF_TRIALS=10000		
	echo "${NUM_OF_TRIALS}"		
}			
			
Set_CONFINTRVL_CUTOFF()			
{			
	CUTOFF=0.3		
	echo "${CUTOFF}"		
}			
			
Set_CONFINTRVL_PREVIOUS_PATH()			
{			
	PREVIOUS_PATH=./5.compare/10.cbind_confidence_interval		
	echo "${PREVIOUS_PATH}"		
}			
			
Set_CONFINTRVL_POPULATION()			
{			
	POPULATION=RIL		
	echo "${POPULATION}"		
}			
			
# --------------------------------------------------			
# Bat_awk_custom2.sh			
# --------------------------------------------------			
			
Set_AWK_CUSTOM2_DEPTHs()			
{			
	DEPTHs="5 7"		
	echo "${DEPTHs}"		
}			
			
Set_AWK_CUSTOM2_MIN_SNPINDEX()			
{			
	MIN_SNPINDEX=0.3		
	echo "${MIN_SNPINDEX}"		
}			
			
Set_AWK_CUSTOM2_MAX_SNPINDEX()			
{			
	MAX_SNPINDEX=1		
	echo "${MAX_SNPINDEX}"		
}			
			
			
