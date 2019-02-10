#! /bin/sh			
			
# ==================================================			
# Those keys enable users to customize this script			
# ==================================================			
BULK_NAME=""			
IDs=""			
			
MISs=""			
			
PILEUPDB_PATH=""			
PILEUPDB_NAME=""			
PILEUPDB_MIS_FIXED=""	
QVAL=""			
PVAL=""			
			
# --------------------------------------------------			
# load common functions			
. ../0.common/common.fnc			
# --------------------------------------------------			
			
if [ -z ${BULK_NAME} ]; then			
	BULK_NAME=`Set_BULK_NAME`		
fi			
			
if [ -z ${IDs} ]; then			
	IDA=`Set_BULK_NAME_IDA`		
	IDB=`Set_BULK_NAME_IDB`		
	IDs=(${IDA} ${IDB})		
fi			

if [ -z ${QVAL} ]; then			
	QVAL=`Set_READ_QVAL`		
fi			
			
if [ -z ${PVAL} ]; then			
	PVAL=`Set_READ_PVAL`		
fi			
	
			
if [ -z ${MISs} ]; then			
	MISs=`Set_MISs`		
fi			
			
if [ -z ${PILEUPDB_PATH} ]; then			
	PILEUPDB_PATH=`Set_PILEUPDB_PATH`		
fi			
			
if [ -z ${PILEUPDB_NAME} ]; then			
	PILEUPDB_NAME=`Set_PILEUPDB_NAME`		
fi			
			
if [ -z ${PILEUPDB_MIS_FIXED} ]; then			
	PILEUPDB_MIS_FIXED=`Set_PILEUPDB_MIS_FIXED`		
fi			
			
# ==================================================			
# environment			
# ==================================================			
TOPPATH_SCRIPTS=""			
TOPPATH_COVAL=""			
			
if [ -z ${TOPPATH_SCRIPTS} ]; then			
	TOPPATH_SCRIPTS=`Set_TOPPATH_SCRIPTS`		
fi			
			
if [ -z ${TOPPATH_COVAL} ]; then			
	TOPPATH_COVAL=`Set_TOPPATH_COVAL`		
fi			
			
# ==================================================			
			
if [ -z $1 ];then			
	echo "missing args!"		
	echo "[usage]"		
	echo "    $0 <INT1>"		
	echo "        <INT1> : 0 or 1 (0:${IDs[0]} / 1:${IDs[1]})"		
	exit 0		
else			
	if [ $1 -eq 0 -o $1 -eq 1 ];then		
		myID=$1	
	else		
		echo "invalid args!"	
		echo "[usage]"	
		echo "    $0 <INT1>"	
		echo "        <INT1> : 0 or 1 (0:${IDs[0]} / 1:${IDs[1]})"	
		exit 0	
	fi		
fi			
			
NAME="${BULK_NAME}_${IDs[${myID}]}"			
QNAME="${NAME}_q${QVAL}p${PVAL}"			
# --------------------------------------------------			
			
SRCPATH="30.coval_call"			
			
OUTPATH="40.exclude_common_snps"			
			
for mis in ${MISs}			
do			
	SNPPILEUP="${SRCPATH}/${QNAME}_MSR_Cov_${mis}_S-snp.pileup"		
	PREFIX0=`basename ${SNPPILEUP} .pileup`		
			
	# --------------------------------------------------		
	pileupdb_dir=${PILEUPDB_PATH}		
	mismatch=${PILEUPDB_MIS_FIXED}		
	if [ ${PILEUPDB_MIS_FIXED} -eq 0 ];then		
		mismatch=${mis}	
	fi		
	# DBPILEUP1="${pileupdb_dir}/Dummy/Null/Null_MSR_Cov_${mismatch}_S_HT.pileup"		
	# DBPILEUP2="${pileupdb_dir}/Dummy/Null/Null_MSR_Cov_${mismatch}_S_HT.pileup"		
	# DBPILEUP3="${pileupdb_dir}/Dummy/Null/Null_MSR_Cov_${mismatch}_S_HT.pileup"		
	# ...		
	# DBPILEUP12="${pileupdb_dir}/Dummy/Null/Null_MSR_Cov_${mismatch}_S_HT.pileup"		
	DBPILEUP1="${pileupdb_dir}/${PILEUPDB_NAME}_MSR_Cov_${mismatch}_S-snp.pileup"		
			
	CMD1="${TOPPATH_SCRIPTS}/3./extract_common_snp_pos.pl"		
	# --------------------------------------------------		
	echo "Now runnning extract_common_snp_pos.pl................"		
	# --------------------------------------------------		
	echo "${CMD1} ${SNPPILEUP} ${DBPILEUP1}"		
	eval "${CMD1} ${SNPPILEUP} ${DBPILEUP1}"		
	# --------------------------------------------------		
	echo "extract_common_snp_pos.pl were finished!!"		
			
	# --------------------------------------------------		
	# created files automatically		
	# 	${PREFIX0}-common-pos.txt	
	# 	${DBPILEUP1}-common-pos.txt	
	# 	${DBPILEUP2}-common-pos.txt	
	# 	...	
	# --------------------------------------------------		
			
	RENAME="mv"		
	COMMONS="${OUTPATH}/${PREFIX0}-common-pos.pileup"		
	echo "${RENAME} ${PREFIX0}-common-pos.txt ${COMMONS}" 		
	eval "${RENAME} ${PREFIX0}-common-pos.txt ${COMMONS}" 		
			
	DELETE="rm -f"		
	echo "${DELETE} *-common-pos.txt"		
	eval "${DELETE} *-common-pos.txt"		
			
			
	CMD2="${TOPPATH_SCRIPTS}/3./exclude_common_snps.pl"		
	# --------------------------------------------------		
	echo "Now runnning exclude_common_snps.pl................"		
	# --------------------------------------------------		
	echo "${CMD2} ${SNPPILEUP} ${COMMONS}"		
	eval "${CMD2} ${SNPPILEUP} ${COMMONS}"		
	# --------------------------------------------------		
	echo "exclude_common_snps.pl were finished!!"		
			
	echo "${RENAME} ${PREFIX0}-rmc2snp.pileup ${OUTPATH}/" 		
	eval "${RENAME} ${PREFIX0}-rmc2snp.pileup ${OUTPATH}/" 		
done

