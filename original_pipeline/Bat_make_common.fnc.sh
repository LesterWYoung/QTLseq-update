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
# change generic directory names to user configured names			
# ==================================================			
			
my_cultivar_name="1.qualify_read/${Key1_My_cultivar_sample_name}"			
bulk_name_ida="1.qualify_read/${Key1_Bulked_sample_name}_${Key1_Bulked_sample_Type_A}"			
bulk_name_idb="1.qualify_read/${Key1_Bulked_sample_name}_${Key1_Bulked_sample_Type_B}"			
			
MKDIR="mkdir -p"			
						
if [ ! -e "${bulk_name_ida}" ] ; then			
	CMD="MKDIR ${bulk_name_ida}"		
	echo ${CMD}		
	eval ${CMD}		
fi			
			
if [ ! -e "${bulk_name_idb}" ] ; then			
	CMD="MKDIR ${bulk_name_idb}"		
	echo ${CMD}		
	eval ${CMD}		
fi			

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
		abspath=$(\cd $(dirname "$1") && pwd)/$(basename "$1")	
	fi		
	echo ${abspath}		
}			
			
# --------------------------------------------------			
abspath=`relativepath_to_absolutepath ${Key0_TopPath_Scripts}`			
Key0_TopPath_Scripts="${abspath}"			
			
# --------------------------------------------------			
abspath=`relativepath_to_absolutepath ${Key0_TopPath_Work}`			
Key0_TopPath_Work="${abspath}"			
			
# --------------------------------------------------
# Path to trimmomatic: created March 2019 LY
# --------------------------------------------------
abspath=`relativepath_to_absolutepath ${Key1_Path_to_trimmomatic}`
Key1_Path_to_trimmomatic="${abspath}"

# --------------------------------------------------
# Path to secondary reference directories: created March 2019 LY
# --------------------------------------------------
abspath=`relativepath_to_absolutepath 1.qualify_read/secondary_readfiles`
trimmed_parental_readfile_dir="${abspath}"
abspath=`relativepath_to_absolutepath 2.make_consensus/secondary_reference`
secondary_ref_dir="${abspath}"

# --------------------------------------------------			
# Determine if secondary reference is being created: LY
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
			else		# create diectory to put secondary reference
			CMD="MKDIR ${trimmed_parental_readfile_dir}"
			echo ${CMD}
			eval ${CMD}
			CMD="MKDIR ${secondary_ref_dir}"
			echo ${CMD}
			eval ${CMD}
			
			fi
		fi
	fi
fi

# --------------------------------------------------
# Reference sequence location
# --------------------------------------------------
abspath=`relativepath_to_absolutepath ${Key2_Path_public_reference_FASTA}`			
Key2_Path_public_reference_FASTA="${abspath}"


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
			
Set_TOPPATH_COVAL(){			
	TOPPATH_COVAL="${Key0_TopPath_Coval}"		
	echo "\${TOPPATH_COVAL}"		
}			

Set_MAKE_SECONDARY_REF(){
	MAKE_SECONDARY_REF="${Key0_Make_secondary_reference}"
	echo "\${MAKE_SECONDARY_REF}"
}
			
# ##################################################			
# for 1.qualify_read
# now using trimmomatic instead 13 mar 2019			
# ##################################################			
Set_TRIMMOMATIC(){			
	TRIMMOMATIC="${Key1_Path_to_trimmomatic}"		
	echo "\${TRIMMOMATIC}"		
}
# ==================================================			
# Bat_fastq_quality_filter.sh			
# ==================================================			
Set_READ_QVAL(){			
	READ_QVAL=${Key1_Phred_quality_score_for_bulked}		
	echo "\${READ_QVAL}"		
}			

# Function to setting $MIN_LEN_BULKED_READS. Trimmomatic will remove reads shorter that this number (and it's mate): LY
Set_MIN_LEN_BULKED_READS(){
	MIN_LEN_BULKED_READS=${Key1_Min_length_bulked_reads}
	echo "\${MIN_LEN_BULKED_READS}"
}			

Set_READ_QVAL_MY_CULTIVAR(){			
	READ_QVAL=${Key1_Phred_quality_score_for_my_cultivar}		
	echo "\${READ_QVAL}"		
}

# Function to setting $MIN_LUN_MY_CULTIVAR_READS. Trimmomatic will remove reads shorter that this number (and it's mate)
Set_MIN_LEN_MY_CULTIVAR_READS(){
	MIN_LEN_MY_CULTIVAR_READS=${Key1_Min_length_my_cultivar_reads}
	echo "\${MIN_LEN_MY_CULTIVAR_READS}"
}			
			
# ##################################################			
# for 2.make_consensus 			
# ##################################################			
# ==================================================			
# If secondary ref being made, define these functions			
# ==================================================			
Set_PATH_TO_TRIMMED_PARENTAL_READS(){			
	PATH_TO_TRIMMED_PARENTAL_READS="${trimmed_parental_readfile_dir}"
	echo "\${PATH_TO_TRIMMED_PARENTAL_READS}"		
}			

Set_PATH_TO_SECONDARY_REF(){			
	PATH_TO_SECONDARY_REF="${secondary_ref_dir}"		
	echo "\${PATH_TO_SECONDARY_REF}"		
}			

Set_MY_CULTIVAR_NAME(){			
	MY_CULTIVAR_NAME="${Key1_My_cultivar_sample_name}"		
	echo "\${MY_CULTIVAR_NAME}"		
}			
			
# ==================================================			
# Bat_bowtie2_secondaryref.sh			
# ==================================================			
Set_PUBLIC_REF_FASTA(){			
	REF_FASTA="${Key2_Path_public_reference_FASTA}"		
	echo "\${REF_FASTA}"		
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
			
Set_REF_FASTA(){			
	REF_FASTA="${Key3_Path_reference_FASTA}"		
	echo "\${REF_FASTA}"		
}			
			
# --------------------------------------------------			
# for Bowtie2 alignment		
# --------------------------------------------------			
Set_BOWTIE2_CPU(){			
	BOWTIE2_CPU=${Key2_Bowtie2_CPU}		
	echo "\${BOWTIE2_CPU}"	
}			

Set_BOWTIE2_OPTIONS(){
	BOWTIE2_OPTIONS='${Key2_BT2op1} ${Key2_BT2op2} ${Key2_BT2op3} ${Key2_BT2op4} ${Key2_BT2op5}'
	echo "\${BOWTIE2_OPTIONS}"
}

# --------------------------------------------------
# for bcftools processing
# --------------------------------------------------
Set_BCFT_MPILEUP(){
	BCFT_MPILEUP='${Key2_BCFT_mpileup}'
	echo "\${BCFT_MPILEUP}"
}

Set_BCFT_CALL(){
        BCFT_CALL='${Key2_BCFT_call}' 
        echo "\${BCFT_CALL}"
}

Set_BCFT_NORM(){
        BCFT_NORM='${Key2_BCFT_norm}' 
        echo "\${BCFT_NORM}"
}

Set_BCFT_FILTER(){
        BCFT_FILTER='${Key2_BCFT_filter}' 
        echo "\${BCFT_FILTER}"
}

Set_BCFT_CONSENSUS(){
        BCFT_CONSENSUS='${Key2_BCFT_consensus}' 
        echo "\${BCFT_CONSENSUS}"
}
			
EOT
