#! /bin/sh			
			
# ==================================================			
# Those keys enable users to customize this script			
# ==================================================			
MY_CULTIVAR_NAME=""			
SRC_PATH=""			
			
QVAL=""			
PVAL=""			
# ==================================================			
			
# --------------------------------------------------			
# load common functions			
. ../0.common/common.fnc			
# --------------------------------------------------			
if [ -z ${MY_CULTIVAR_NAME} ]; then			
	MY_CULTIVAR_NAME=`Set_MY_CULTIVAR_NAME`		
fi			
			
			
if [ -z ${SRC_PATH} ]; then			
	SRC_PATH=`Set_SRC_READ_PATH_MY_CULTIVAR`		
fi			
			
if [ -z ${QVAL} ]; then			
	QVAL=`Set_READ_QVAL_MY_CULTIVAR`		
fi			
			
if [ -z ${PVAL} ]; then			
	PVAL=`Set_READ_PVAL_MY_CULTIVAR`		
fi			
			
			
			
# --------------------------------------------------			
			
NAME=${MY_CULTIVAR_NAME}			
OUTPATH="10.bwa2bam"			
			
# ==================================================			
SRC_NAME="*_[0-9]*.RnFq${QVAL}p${PVAL}S.[12].gz"			
SEQ_HEAD="p_q${QVAL}p${PVAL}"			
			
SEQ_TAIL="sequence.txt.gz"			
			
echo "mkdir -p ${OUTPATH}/${NAME}"			
eval "mkdir -p ${OUTPATH}/${NAME}"			
			
READs=`ls ${SRC_PATH}/${SRC_NAME}`			
			
for rd in ${READs}			
do			
			
#                                   #${SRC_PATH}/*_[0-9]*.RnFq30p90S.[12].gz			
	rd_basename=`basename ${rd}`    # -> *_[0-9]*.RnFq30p90S.[12].gz		
	mylane=${rd_basename##*_}        # -> [0-9]*.RnFq30p90S.[12].gz		
	mylane=${mylane%%.*}            # -> [0-9]*		
	my1or2=${rd_basename%.*}        # -> *_[0-9]*.RnFq30p90S.[12]		
	my1or2=${my1or2##*.}            # -> [12]		
			
	echo "ln -s ${rd} ${OUTPATH}/${NAME}/${SEQ_HEAD}_${mylane}_${my1or2}_${SEQ_TAIL}"		
	eval "ln -s ${rd} ${OUTPATH}/${NAME}/${SEQ_HEAD}_${mylane}_${my1or2}_${SEQ_TAIL}"		
			
done

