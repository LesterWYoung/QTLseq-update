#!/bin/bash			
			
			
# ==================================================			
# load config.txt 			
# ==================================================			
if [ -z $1 ] ; then			
	. ./config.txt		
else			
	. $1		
fi			
			
# ==================================================			
# noname directory -> your samples			
# ==================================================			
			
my_cultivar_name="1.qualify_read/${Key1_My_cultivar_sample_name}"			
bulk_name_ida="1.qualify_read/${Key1_Bulked_sample_name}_${Key1_Bulked_sample_Type_A}"			
bulk_name_idb="1.qualify_read/${Key1_Bulked_sample_name}_${Key1_Bulked_sample_Type_B}"			
qdir_my_cultivar="q${Key1_Phred_quality_score_for_my_cultivar}p${Key1_Percentage_of_above_score_for_my_cultivar}"			
qdir="q${Key1_Phred_quality_score_for_bulked}p${Key1_Percentage_of_above_score_for_bulked}"			
			
MKDIR="mkdir -p"			
						
if [ ! -e "${bulk_name_ida}" ] ; then			
	CMD="mv 1.qualify_read/mybulk_A ${bulk_name_ida}"		
	echo ${CMD}		
	eval ${CMD}		
fi			
			
if [ ! -e "${bulk_name_idb}" ] ; then			
	CMD="mv 1.qualify_read/mybulk_B ${bulk_name_idb}"		
	echo ${CMD}		
	eval ${CMD}		
fi			

#Creation of secondary reference moved to "Determine if secondary reference is being created" 			

#Don't need the following direcories now
#echo "${MKDIR} ${bulk_name_ida}/${qdir}/sep_pair"			
#eval "${MKDIR} ${bulk_name_ida}/${qdir}/sep_pair"			
			
#echo "${MKDIR} ${bulk_name_idb}/${qdir}/sep_pair"			
#eval "${MKDIR} ${bulk_name_idb}/${qdir}/sep_pair"			
			
#echo "${MKDIR} ${my_cultivar_name}/${qdir_my_cultivar}/sep_pair"			
#eval "${MKDIR} ${my_cultivar_name}/${qdir_my_cultivar}/sep_pair"			
			
# ==================================================			
# convert some variables 			
# ==================================================			
# --------------------------------------------------			
# relative path name --> absolute path name			
# --------------------------------------------------			
relativepath_to_absolutepath(){			
	if [ $1 = "." ] ; then		
		abspath=`pwd`	
	else		
		abspath=$(cd $(dirname "$1") && pwd)/$(basename "$1")	
	fi		
	echo ${abspath}		
}			
			
# --------------------------------------------------			
#abspath=`relativepath_to_absolutepath ${Key0_TopPath_Scripts}`			
#Key0_TopPath_Scripts="${abspath}"			
			
# --------------------------------------------------			
#abspath=`relativepath_to_absolutepath ${Key0_TopPath_Coval}`			
#Key0_TopPath_Coval="${abspath}"			
			
# --------------------------------------------------			
abspath=`relativepath_to_absolutepath ${Key0_TopPath_Work}`			
Key0_TopPath_Work="${abspath}"			
			
# --------------------------------------------------			
abspath=`relativepath_to_absolutepath ${Key2_Path_public_reference_FASTA}`			
Key2_Path_public_reference_FASTA="${abspath}"

# --------------------------------------------------			
# Determine if secondary reference is being created
# --------------------------------------------------			
if [[ ${Key0_Make_secondary_reference} = "yes/no" ]] ; then
	echo "Change Key0_Make_secondary_reference in ./config.txt to yes or no"
	exit 1
else			
	if [[ ${Key0_Make_secondary_reference} = "yes" ]]; then
		if [[ ${my_cultivar_name} = "INSERT-secondary-ref-name-OR-leave-empty" ]]; then
			echo "Change Key1_My_cultivar_sample_name to secondary reference name"
			exit 2
		else
			if [[ ${my_cultivar_name} = "" ]]; then
				echo "Insert value for Key1_My_cultivar_sample_name in config.txt"
				exit 3
			else
			CMD="mv 1.qualify_read/${my_cultivar_name}"		
			echo ${CMD}		
			eval ${CMD}
			fi
		fi
	fi
fi

#Don't need the following as not using now			
#if [ ${Key2_Coval_realign_for_my_cultivar} = "yes" ] ; then			
#	Key2_Coval_realign_for_my_cultivar=""		
#else			
#	Key2_Coval_realign_for_my_cultivar="--disalign"		
#fi			
			
# --------------------------------------------------			
if [ ${Key3_Coval_realign_for_bulked} = "yes" ] ; then			
	Key3_Coval_realign_for_bulked=""		
else			
	Key3_Coval_realign_for_bulked="--disalign"		
fi			
			
# --------------------------------------------------			
if [ ${Key1_Score_type_for_bulked} = "sanger" ] ; then			
	QOPT_for_bulked="-Q 33"		
else			
	QOPT_for_bulked=""		
fi			
			
# --------------------------------------------------			
if [ ${Key1_Score_type_for_my_cultivar} = "sanger" ] ; then			
	QOPT_for_my_cultivar="-Q 33"		
else			
	QOPT_for_my_cultivar=""		
fi			
			
# --------------------------------------------------			
if [ ${Key3_Mode_reference_FASTA} -eq 0 ] ; then			
	Key3_Path_reference_FASTA="${Key0_TopPath_Work}"		
	Key3_Path_reference_FASTA="${Key3_Path_reference_FASTA}/2.make_consensus"		
	Key3_Path_reference_FASTA="${Key3_Path_reference_FASTA}/90.align_to_this_fasta"		
	Key3_Path_reference_FASTA="${Key3_Path_reference_FASTA}/00.reference"		
	Key3_Path_reference_FASTA="${Key3_Path_reference_FASTA}/${Key1_My_cultivar_sample_name}_${qdir_my_cultivar}_MSR_Cov_${Key2_Maximum_number_mismatch_default_for_my_cultivar}_S-snp_RYKMSWBDHV2ACGT.fa"		
else			
    abspath=`relativepath_to_absolutepath ${Key3_Path_reference_FASTA}`			
    Key3_Path_reference_FASTA="${abspath}"			
fi			
			
# --------------------------------------------------			
if [ ${Key3_Mode_PileupDB} -eq 0 ] ; then			
	Key3_Path_PileupDB="${Key0_TopPath_Work}"		
	Key3_Path_PileupDB="${Key3_Path_PileupDB}/2.make_consensus"		
	Key3_Path_PileupDB="${Key3_Path_PileupDB}/90.align_to_this_fasta"		
	Key3_Path_PileupDB="${Key3_Path_PileupDB}/30.coval_call"		
	Key3_PileupDB_head_name="${Key1_My_cultivar_sample_name}_${qdir_my_cultivar}"		
	Key3_PileupDB_mismatch_fixed=${Key2_Maximum_number_mismatch_default_for_my_cultivar}		
else			
    abspath=`relativepath_to_absolutepath ${Key3_Path_PileupDB}`			
    Key3_Path_PileupDB="${abspath}"			
fi			
			
# --------------------------------------------------			
# comma delimiter --> space delimiter			
# --------------------------------------------------			
# Key3_Maximum_mismatches=2,3,4			
arry=`echo ${Key3_Maximum_mismatches} | tr -s ',' ' '`			
Key3_Maximum_mismatches=${arry}			
			
# --------------------------------------------------			
# Key5_Depths=5,7			
arry=`echo ${Key5_Depths} | tr -s ',' ' '`			
Key5_Depths=${arry}			
			
# --------------------------------------------------			
# Key5_Miss=5,7			
arry=`echo ${Key5_Miss} | tr -s ',' ' '`			
Key5_Miss=${arry}			
			
OUTNAME="0.common/common.fnc"			
			
			
# ##################################################			
# here document			
# ##################################################			
			
cat << EOT > ${OUTNAME}			
# ##################################################			
# global environment			
# ##################################################			
Set_TOPPATH_SCRIPTS(){			
	TOPPATH_SCRIPTS="${Key0_TopPath_Scripts}"		
	echo "\${TOPPATH_SCRIPTS}"		
}			
			
Set_

Set_TOPPATH_COVAL(){			
	TOPPATH_COVAL="${Key0_TopPath_Coval}"		
	echo "\${TOPPATH_COVAL}"		
}			
			
# ##################################################			
# for 1.qualify_read			
# ##################################################			
# ==================================================			
# Bat_fastq_quality_filter.sh			
# ==================================================			
Set_READ_QVAL(){			
	READ_QVAL=${Key1_Phred_quality_score_for_bulked}		
	echo "\${READ_QVAL}"		
}			
			
Set_READ_PVAL(){			
	READ_PVAL=${Key1_Percentage_of_above_score_for_bulked}		
	echo "\${READ_PVAL}"		
}			
			
Set_READ_QOPT(){			
	READ_QOPT="${QOPT_for_bulked}"		
	echo "\${READ_QOPT}"		
			
	# READ_QOPT="-Q 33"	in case of CASAVA 1.8 later	
	# READ_QOPT=""		in case of the previous Illumina-specific offset value of 64
}			
#--                              #2014/04/08 kikuchi			
Set_READ_QVAL_MY_CULTIVAR(){			
	READ_QVAL=${Key1_Phred_quality_score_for_my_cultivar}		
	echo "\${READ_QVAL}"		
}			
			
Set_READ_PVAL_MY_CULTIVAR(){			
	READ_PVAL=${Key1_Percentage_of_above_score_for_my_cultivar}		
	echo "\${READ_PVAL}"		
}			
			
Set_READ_QOPT_MY_CULTIVAR(){			
	READ_QOPT="${QOPT_for_my_cultivar}"		
	echo "\${READ_QOPT}"		
			
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
	SRC_READ_PATH_MY_CULTIVAR="${Key0_TopPath_Work}/1.qualify_read/${Key1_My_cultivar_sample_name}/q${Key1_Phred_quality_score_for_my_cultivar}p${Key1_Percentage_of_above_score_for_my_cultivar}/sep_pair"		
	echo "\${SRC_READ_PATH_MY_CULTIVAR}"		
}			
			
Set_MY_CULTIVAR_NAME(){			
	MY_CULTIVAR_NAME="${Key1_My_cultivar_sample_name}"		
	echo "\${MY_CULTIVAR_NAME}"		
}			
			
# ==================================================			
# Bat_bwa2bam.sh			
# ==================================================			
Set_PUBLIC_REF_FASTA(){			
	REF_FASTA="${Key2_Path_public_reference_FASTA}"		
	echo "\${REF_FASTA}"		
}			
			
# ==================================================			
# Bat_run_coval-refine-bam.pl.sh			
# ==================================================			
Set_MIS_MATCH_FOR_MAKE_CONSENSUS(){			
	MIS_MATCH=${Key2_Maximum_number_mismatch_default_for_my_cultivar}		
	echo "\${MIS_MATCH}"		
}			
			
Set_DONTREALIGN_FOR_MAKE_CONSENSUS(){			
	DONTREALIGN="${Key2_Coval_realign_for_my_cultivar}"		
	echo "\${DONTREALIGN}"		
			
	# DONTREALIGN=""		
	# DONTREALIGN="--disalign"		
}			
			
# ==================================================			
# Bat_run_coval-call-pileup.pl.sh			
# ==================================================			
Set_COVAL_CALL_MINNUM_FOR_MAKE_CONSENSUS()			
{			
	COVALCALL_MINNUM=${Key2_Minimum_number_of_depth_for_my_cultivar}		
	echo "\${COVALCALL_MINNUM}"		
}			
			
Set_COVAL_CALL_MAXR_FOR_MAKE_CONSENSUS()			
{			
	COVALCALL_MAXR=${Key2_Maximum_number_of_depth_for_my_cultivar}		
	echo "\${COVALCALL_MAXR}"		
}			
			
Set_COVAL_CALL_MINFREQ_FOR_MAKE_CONSENSUS()			
{			
	COVALCALL_MINFREQ=${Key2_Minimum_frequency_for_my_cultivar}		
	echo "\${COVALCALL_MINFREQ}"		
}			
			
Set_COVAL_CALL_MINTNUM_FOR_MAKE_CONSENSUS()			
{			
	COVALCALL_MINTNUM=${Key2_Minimum_number_of_reads_heterozygous_for_my_cultivar}		
	echo "\${COVALCALL_MINTNUM}"		
}			
			
Set_COVAL_CALL_MINQUALBASE_FOR_MAKE_CONSENSUS()			
{			
	COVALCALL_MINQUALBASE=${Key2_Minimum_basecall_quality_for_my_cultivar}		
	echo "\${COVALCALL_MINQUALBASE}"		
}			
			
Set_COVAL_CALL_MINQUALAVE_FOR_MAKE_CONSENSUS()			
{			
	COVALCALL_MINQUALAVE=${Key2_Minimum_average_basecall_quality_for_my_cultivar}		
	echo "\${COVALCALL_MINQUALAVE}"		
}			
			
Set_COVAL_CALL_CALLTYPE_FOR_MAKE_CONSENSUS()			
{			
	COVALCALL_CALLTYPE="${Key1_Score_type_for_my_cultivar}"		
	echo "\${COVALCALL_CALLTYPE}"		
			
	# COVALCALL_CALLTYPE="sanger"		
	# COVALCALL_CALLTYPE="illumina"		
}			
			
# ##################################################			
# for 3.alignment			
# ##################################################			
Set_BULK_NAME(){			
	BULK_NAME="${Key1_Bulked_sample_name}"		
	echo "\${BULK_NAME}"		
}			
			
Set_BULK_NAME_IDA(){			
	BULK_NAME_ID="${Key1_Bulked_sample_Type_A}"		
	echo "\${BULK_NAME_ID}"		
}			
			
Set_BULK_NAME_IDB(){			
	BULK_NAME_ID="${Key1_Bulked_sample_Type_B}"		
	echo "\${BULK_NAME_ID}"		
}			
			
Set_REF_FASTA()			
{			
	REF_FASTA="${Key3_Path_reference_FASTA}"		
	echo "\${REF_FASTA}"		
}			
			
# --------------------------------------------------			
# for Bat_make_symbolic_link_of_qualified_fastq.sh			
# --------------------------------------------------			
Set_SRC_READ_PATH(){			
	SRC_READ_PATH="${Key0_TopPath_Work}/1.qualify_read/equalized"		
			
	echo "\${SRC_READ_PATH}"		
}			
			
# --------------------------------------------------			
# for Bat_bwa2bam.sh			
# --------------------------------------------------			
Set_BWA_CPU()			
{			
	BWA_CPU=${Key2_BWA_CPU}		
	echo "\${BWA_CPU}"		
}			
			
# --------------------------------------------------			
# Bat_run_coval-refine-bam.pl.sh			
# --------------------------------------------------			
Set_DONTREALIGN()			
{			
	DONTREALIGN="${Key3_Coval_realign_for_bulked}"		
	echo "\${DONTREALIGN}"		
			
	# DONTREALIGN=""		
	# DONTREALIGN="--disalign"		
}			
			
Set_MIS_MATCH_DEFAULT()			
{			
	MIS_MATCH_DEFAULT=${Key3_Maximum_number_mismatch_default_for_bulked}		
	echo "\${MIS_MATCH_DEFAULT}"		
}			
			
# --------------------------------------------------			
# Bat_run_coval-call-pileup.pl.sh			
# --------------------------------------------------			
Set_COVAL_CALL_MINNUM()			
{			
	COVALCALL_MINNUM=${Key3_Minimum_number_of_depth_for_bulked}		
	echo "\${COVALCALL_MINNUM}"		
}			
			
Set_COVAL_CALL_MAXR()			
{			
	COVALCALL_MAXR=${Key3_Maximum_number_of_depth_for_bulked}		
	echo "\${COVALCALL_MAXR}"		
}			
			
Set_COVAL_CALL_MINFREQ()			
{			
	COVALCALL_MINFREQ=${Key3_Minimum_frequency_for_bulked}		
	echo "\${COVALCALL_MINFREQ}"		
}			
			
Set_COVAL_CALL_MINTNUM()			
{			
	COVALCALL_MINTNUM=${Key3_Minimum_number_of_reads_heterozygous_for_bulked}		
	echo "\${COVALCALL_MINTNUM}"		
}			
			
Set_COVAL_CALL_MINQUALBASE()			
{			
	COVALCALL_MINQUALBASE=${Key3_Minimum_basecall_quality_for_bulked}		
	echo "\${COVALCALL_MINQUALBASE}"		
}			
			
Set_COVAL_CALL_MINQUALAVE()			
{			
	COVALCALL_MINQUALAVE=${Key3_Minimum_average_basecall_quality_for_bulked}		
	echo "\${COVALCALL_MINQUALAVE}"		
}			
			
Set_COVAL_CALL_CALLTYPE()			
{			
	COVALCALL_CALLTYPE="${Key1_Score_type_for_bulked}"		
	echo "\${COVALCALL_CALLTYPE}"		
			
	# COVALCALL_CALLTYPE="sanger"		
	# COVALCALL_CALLTYPE="illumina"		
}			
			
# --------------------------------------------------			
# Bat_exclude_common_snps.pl.sh			
# --------------------------------------------------			
Set_MISs()			
{			
	MISs="${Key3_Maximum_mismatches}"		
	echo "\${MISs}"		
}			
			
Set_PILEUPDB_PATH()			
{			
	PILEUPDB_PATH="${Key3_Path_PileupDB}"		
	echo "\${PILEUPDB_PATH}"		
}			
			
Set_PILEUPDB_NAME()			
{			
	PILEUPDB_NAME="${Key3_PileupDB_head_name}"		
	echo "\${PILEUPDB_NAME}"		
}			
			
Set_PILEUPDB_MIS_FIXED()			
{			
	PILEUPDB_MIS_FIXED=${Key3_PileupDB_mismatch_fixed} # if 0, each mismatch value		
	echo "\${PILEUPDB_MIS_FIXED}"		
}			
			
			
			
			
			
			
			
			
			
			
# --------------------------------------------------			
# Bat_awk_custom.sh			
# --------------------------------------------------			
			
Set_MINDEPTH()			
{			
	MINDEPTH=${Key3_Mindepth}		
	echo "\${MINDEPTH}"		
}			
			
			
			
			
Set_MIN_CONSENSUS_QUALITY()			
{			
	MIN_CONSENSUS_QUALITY=${Key3_Min_consensus_quality}		
	echo "\${MIN_CONSENSUS_QUALITY}"		
}			
			
Set_MIN_SNP_INDEX()			
{			
	MIN_SNP_INDEX=${Key3_Min_snp_index}		
	echo "\${MIN_SNP_INDEX}"		
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
	MISs="${Key5_Miss}"		
	echo "\${MISs}"		
}			
			
Set_CONFINTRVL_CALC_MODE()			
{			
	CALC_MODE=${Key5_Calc_mode}	    # 1:DoCalculateNow  0:UsePreviousOne	
	echo "\${CALC_MODE}"		
}			
			
Set_CONFINTRVL_INDIVIDUALS()			
{			
	INDIVIDUALS=${Key5_Individuals}		
	echo "\${INDIVIDUALS}"		
}			
			
Set_CONFINTRVL_NUM_OF_TRIALS()			
{			
	NUM_OF_TRIALS=${Key5_Num_of_trials}		
	echo "\${NUM_OF_TRIALS}"		
}			
			
Set_CONFINTRVL_CUTOFF()			
{			
	CUTOFF=${Key5_Cutoff}		
	echo "\${CUTOFF}"		
}			
			
Set_CONFINTRVL_PREVIOUS_PATH()			
{			
	PREVIOUS_PATH=${Key5_Previous_path}		
	echo "\${PREVIOUS_PATH}"		
}			
			
Set_CONFINTRVL_POPULATION()			
{			
	POPULATION=${Key5_Population}		
	echo "\${POPULATION}"		
}			
			
# --------------------------------------------------			
# Bat_awk_custom2.sh			
# --------------------------------------------------			
			
Set_AWK_CUSTOM2_DEPTHs()			
{			
	DEPTHs="${Key5_Depths}"		
	echo "\${DEPTHs}"		
}			
			
Set_AWK_CUSTOM2_MIN_SNPINDEX()			
{			
	MIN_SNPINDEX=${Key5_Min_snpindex}		
	echo "\${MIN_SNPINDEX}"		
}			
			
Set_AWK_CUSTOM2_MAX_SNPINDEX()			
{			
	MAX_SNPINDEX=${Key5_Max_snpindex}		
	echo "\${MAX_SNPINDEX}"		
}			
			
			
EOT
