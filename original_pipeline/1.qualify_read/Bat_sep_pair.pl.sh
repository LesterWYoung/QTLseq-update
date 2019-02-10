#! /bin/sh			
			
Command=`basename $0`			
Logf="log.$Command.txt"			
			
# ==================================================			
# Those keys enable users to customize this script			
# ==================================================			
MY_CULTIVAR_NAME=""			
BULK_NAME=""			
IDs=""			
			
QVAL=""			
PVAL=""			
QOPT=""			
			
# --------------------------------------------------			
# load common functions			
. ../0.common/common.fnc			
# --------------------------------------------------			
if [ -z ${MY_CULTIVAR_NAME} ]; then			
	MY_CULTIVAR_NAME=`Set_MY_CULTIVAR_NAME`		
fi			
			
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
	if [ $1 -eq 9 ]; then		
		QVAL=`Set_READ_QVAL_MY_CULTIVAR`	
	fi		
fi			
			
if [ -z ${PVAL} ]; then			
	PVAL=`Set_READ_PVAL`		
	if [ $1 -eq 9 ]; then		
		PVAL=`Set_READ_PVAL_MY_CULTIVAR`	
	fi		
fi			
			
if [ -z ${QOPT} ]; then			
	QOPT=`Set_READ_QOPT`		
	if [ $1 -eq 9 ]; then		
		QOPT=`Set_READ_QOPT_MY_CULTIVAR`	
	fi		
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
	echo "        <INT1> : 0 or 1 (0:${IDs[0]} / 1:${IDs[1]} / 9:${MY_CULTIVAR_NAME})"		
	exit 0		
else			
	if [ $1 -eq 0 -o $1 -eq 1 -o $1 -eq 9 ];then		
		myID=$1	
	else		
		echo "invalid args!"	
		echo "[usage]"	
		echo "    $0 <INT1>"	
		echo "        <INT1> : 0 or 1 (0:${IDs[0]} / 1:${IDs[1]} / 9:${MY_CULTIVAR_NAME})"	
		exit 0	
	fi		
fi			
			
# ==================================================			
if [ $1 -eq 0 -o $1 -eq 1 ];then			
	NAME="${BULK_NAME}_${IDs[${myID}]}"		
else			
	NAME="${MY_CULTIVAR_NAME}"		
fi			
# ==================================================			
			
OUTDESCRPT="RnFq${QVAL}p${PVAL}"			
			
CMD1="${TOPPATH_SCRIPTS}/1./sep_pair3.pl"			
CMD2="gunzip -c"			
CMD3="fastx_quality_stats ${QOPT} -o"			
			
INPATH="${NAME}/q${QVAL}p${PVAL}"			
INFASTQs=`ls ${INPATH}/${NAME}_[0-9]*.${OUTDESCRPT}.1.gz`			
			
OUTPATH="${NAME}/q${QVAL}p${PVAL}/sep_pair"			
OUTFASTQs=""			
			
for myfastq in ${INFASTQs}			
do			
	BASEFASTQ=`basename ${myfastq} .1.gz`		
	infastqgz1="${INPATH}/${BASEFASTQ}.1.gz"		
	infastqgz2="${INPATH}/${BASEFASTQ}.2.gz"		
			
	CMD="${CMD1} ${infastqgz1} ${infastqgz2} ${OUTPATH}"		
	echo ${CMD}		
	eval ${CMD}		
			
	OUTFASTQs="${OUTFASTQs} ${BASEFASTQ}S.1.gz ${BASEFASTQ}S.2.gz ${BASEFASTQ}S.s.gz"		
done			
			
for myfastq in ${OUTFASTQs}			
do			
	CMD="${CMD2} ${OUTPATH}/${myfastq} | ${CMD3} ${OUTPATH}/${myfastq}_stats.txt"		
	echo ${CMD}		
	eval ${CMD}		
done

