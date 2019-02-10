#! /bin/sh			
			
# ==================================================			
# Those keys enable users to customize this script			
# ==================================================			
BULK_NAME=""			
IDs=""			
SRC_PATH=""			
			
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
			
if [ -z ${SRC_PATH} ]; then			
	SRC_PATH=`Set_SRC_READ_PATH`		
fi			

if [ -z ${QVAL} ]; then			
	QVAL=`Set_READ_QVAL`		
fi			
			
if [ -z ${PVAL} ]; then			
	PVAL=`Set_READ_PVAL`		
fi			
			
# --------------------------------------------------			
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
			
# --------------------------------------------------			
			
OUTPATH="10.bwa2bam"			
SRC_TAIL="*_[0-9]*.RnFq${QVAL}p${PVAL}S.[12].gz"			
FQ_HEAD="p_q${QVAL}p${PVAL}"			
FQ_TAIL="sequence.txt.gz"			
			
echo "mkdir -p ${OUTPATH}/${NAME}"			
eval "mkdir -p ${OUTPATH}/${NAME}"			
			
			
READs=`ls ${SRC_PATH}/${NAME}/${SRC_TAIL}`			
for rd in ${READs}			
do			
#                                   #${SRC_PATH}/*_[0-9]*.RnFq30p90S.[12].gz			
	rd_basename=`basename ${rd}`    # -> *_[0-9]*.RnFq30p90S.[12].gz		
	mylane=${rd_basename##*_}        # -> [0-9]*.RnFq30p90S.[12].gz		
	mylane=${mylane%%.*}            # -> [0-9]*		
	my1or2=${rd_basename%.*}        # -> *_[0-9]*.RnFq30p90S.[12]		
	my1or2=${my1or2##*.}            # -> [12]		
			
	echo "ln -s ${rd} ${OUTPATH}/${NAME}/${FQ_HEAD}_${mylane}_${my1or2}_${FQ_TAIL}"		
	eval "ln -s ${rd} ${OUTPATH}/${NAME}/${FQ_HEAD}_${mylane}_${my1or2}_${FQ_TAIL}"		
done
