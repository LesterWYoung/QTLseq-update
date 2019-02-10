#! /bin/sh			
			
Command=`basename $0`			
Logf="log.$Command.txt"			
			
# ==================================================			
# Those keys enable users to customize this script			
# ==================================================			
BULK_NAME=""			
IDs=""			
MISs=""			
MIN_DEPTH=""
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
			
if [ -z ${MIN_DEPTH} ]; then			
	MIN_DEPTH=`Set_MINDEPTH`		
fi			
# ==================================================			
# environment			
# ==================================================			
TOPPATH_SCRIPTS=""			
			
if [ -z ${TOPPATH_SCRIPTS} ]; then			
	TOPPATH_SCRIPTS=`Set_TOPPATH_SCRIPTS`		
fi			
			
NAME="${BULK_NAME}"			
			
for mis in ${MISs}			
do			
	SRCPATH="../3.alignment/50.awk_custom/mut_index_${mis}"		
	OUTPATH="10.paired_or_unpaired"		
	OUTPATH="${OUTPATH}/mut_index_${mis}"		
			
	echo "mkdir -p ${OUTPATH}"		
	eval "mkdir -p ${OUTPATH}"		
			
    CMD="${TOPPATH_SCRIPTS}/4./paired_or_unpaired.pl"			
	NAME_A="${NAME}_${IDs[0]}_q${QVAL}p${PVAL}"		
	NAME_B="${NAME}_${IDs[1]}_q${QVAL}p${PVAL}"		
			
	A_SNP="${NAME_A}_cov${mis}_co${MIN_DEPTH}.txt"		
	B_SNP="${NAME_B}_cov${mis}_co${MIN_DEPTH}.txt"		
	#--		
	A_SNP_UNPAIRED="${OUTPATH}/unpaired_${A_SNP}"		
	B_SNP_UNPAIRED="${OUTPATH}/unpaired_${B_SNP}"		
			
	A_SNP_PAIRED="${OUTPATH}/paired_${A_SNP}"		
	B_SNP_PAIRED="${OUTPATH}/paired_${B_SNP}"		
			
	echo "${CMD} ${SRCPATH}/${A_SNP} ${SRCPATH}/${B_SNP} ${A_SNP_PAIRED} ${B_SNP_PAIRED} ${A_SNP_UNPAIRED} ${B_SNP_UNPAIRED}"		
	eval "${CMD} ${SRCPATH}/${A_SNP} ${SRCPATH}/${B_SNP} ${A_SNP_PAIRED} ${B_SNP_PAIRED} ${A_SNP_UNPAIRED} ${B_SNP_UNPAIRED}"		
			
done
