#! /bin/sh			
			
# ==================================================			
# Those keys enable users to customize this script			
# ==================================================			
BULK_NAME=""			
IDs=""			
MISs=""			
			
QVAL=""			
PVAL=""			
			
# --------------------------------------------------			
MIN_DEPTH=""			
MIN_CONSENSUS_QUALITY=""			
MIN_SNP_INDEX=""			
			
# ==================================================			
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
			
if [ -z ${MIN_DEPTH} ]; then			
	MIN_DEPTH=`Set_MINDEPTH`		
fi			
			
if [ -z ${MIN_CONSENSUS_QUALITY} ]; then			
	MIN_CONSENSUS_QUALITY=`Set_MIN_CONSENSUS_QUALITY`		
fi			
			
if [ -z ${MIN_SNP_INDEX} ]; then			
	MIN_SNP_INDEX=`Set_MIN_SNP_INDEX`		
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
			
NAME="${BULK_NAME}"			
			
			
NAME="${NAME}_${IDs[${myID}]}"			
QNAME="${NAME}_q${QVAL}p${PVAL}"			
			
DEPTHs=${MIN_DEPTH}			
			
# --------------------------------------------------			
			
SRCPATH="40.exclude_common_snps"			
OUTPATH="50.awk_custom"			
			
INNAMEHEAD="${SRCPATH}/${QNAME}_MSR_Cov"			
INNAMETAIL="S-snp-rmc2snp.pileup"			
			
			
CMD0="cat"			
CMD1="awk '\$8>=${dep}'"			
CMD2="awk '\$5 >= ${MIN_CONSENSUS_QUALITY}'"			
CMD3="awk '\$9 >= ${MIN_SNP_INDEX}'"			
			
for mis in ${MISs}			
do			
			
	INPILEUP="${INNAMEHEAD}_${mis}_${INNAMETAIL}"		
			
	DESTPATH="mut_index_${mis}"		
	echo "mkdir -p ${OUTPATH}/${DESTPATH}"		
	eval "mkdir -p ${OUTPATH}/${DESTPATH}"		
			
	for dep in ${DEPTHs}		
	do		
		CMD1="awk '\$8>=${dep}'"	
		OUTPILEUP="${OUTPATH}/${DESTPATH}/${QNAME}_cov${mis}_co${dep}.txt"	
		echo "${CMD0} ${INPILEUP} | ${CMD1} | ${CMD2} | ${CMD3} > ${OUTPILEUP}"	
		eval "${CMD0} ${INPILEUP} | ${CMD1} | ${CMD2} | ${CMD3} > ${OUTPILEUP}"	
			
	done		
done			
	
