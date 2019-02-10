#! /bin/sh			
			
# ==================================================			
# Those keys enable users to customize this script			
# ==================================================			
MY_CULTIVAR_NAME=""			
			
# --------------------------------------------------			
# load common functions			
. ../../0.common/common.fnc			
# --------------------------------------------------			
if [ -z ${MY_CULTIVAR_NAME} ]; then			
	MY_CULTIVAR_NAME=`Set_MY_CULTIVAR_NAME`		
fi			
			
NAME=${MY_CULTIVAR_NAME}

SRC_PATH="../10.bwa2bam/${NAME}"			
SEQ_TAIL="sequence.txt.gz"			
			
OUTPATH="10.bwa2bam"			
			
echo "mkdir -p ${OUTPATH}/${NAME}"			
eval "mkdir -p ${OUTPATH}/${NAME}"			
			
READs=`ls ${SRC_PATH}/*${SEQ_TAIL}`			
			
for rd in ${READs}			
do			
	echo "ln -s ../../${rd} ${OUTPATH}/${NAME}/"		
	eval "ln -s ../../${rd} ${OUTPATH}/${NAME}/"		
done
