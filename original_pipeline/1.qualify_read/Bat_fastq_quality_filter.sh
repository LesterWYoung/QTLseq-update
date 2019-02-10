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
if [ -z $1 ];then			
	echo "missing args!"		
	echo "[usage]"		
	echo "    $0 <INT1>"		
	echo "        <INT1> : 0 or 1 or 9(0:${IDs[0]} / 1:${IDs[1]} / 9:${MY_CULTIVAR_NAME})"		
	exit 0		
else			
	if [ $1 -eq 0 -o $1 -eq 1 -o $1 -eq 9 ];then		
		myID=$1	
	else		
		echo "invalid args!"	
		echo "[usage]"	
		echo "    $0 <INT1>"	
		echo "        <INT1> : 0 or 1 or 9(0:${IDs[0]} / 1:${IDs[1]} / 9:${MY_CULTIVAR_NAME})"	
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
INPATH=${NAME}			
INFASTQs=`ls ${INPATH}/${NAME}_[0-9]*.Rn.[12].gz`			
echo "${INFASTQs}"			
			
CMD1="gunzip -c"			
CMD2="fastq_quality_filter ${QOPT} -v -q ${QVAL} -p ${PVAL} -z -o"			
CMD3="fastx_quality_stats ${QOPT} -o"			
OUTDESCRPT="Fq${QVAL}p${PVAL}"			
			
OUTPATH="${NAME}/q${QVAL}p${PVAL}"			
			
OUTFASTQs=""			
for myfastq in ${INFASTQs}			
do			
	BASEFASTQ=`basename ${myfastq}`		
	OUTFASTQ=""		
#-                              #${NAME}_[0-9]*.Rn.[12].gz    			
	BASEFASTQ=${BASEFASTQ%.*}      # -> ${NAME}_[0-9]*.Rn.[12]		
	OUTFASTQ=${BASEFASTQ%.*}       # -> ${NAME}_[0-9]*.Rn		
	my1or2=${BASEFASTQ##*.}        # -> [12]		
			
    OUTFASTQ="${OUTFASTQ}${OUTDESCRPT}.${my1or2}.gz"			
			
	CMD="${CMD1} ${myfastq} | ${CMD2} ${OUTPATH}/${OUTFASTQ}"		
	echo ${CMD}		
	eval ${CMD}		
			
	OUTFASTQs="${OUTFASTQs} ${OUTFASTQ}"		
done			
			
for myfastq in ${OUTFASTQs}			
do			
	CMD="${CMD1} ${OUTPATH}/${myfastq} | ${CMD3} ${OUTPATH}/${myfastq}_stats.txt"		
	echo ${CMD}		
	eval ${CMD}		
done

