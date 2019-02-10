#! /bin/sh			
			
# ==================================================			
# Those keys enable users to customize this script			
# ==================================================			
MY_CULTIVAR_NAME=""			
MIS_MATCH_DEFAULT=""		
QVAL=""			
PVAL=""			
# --------------------------------------------------			
# load common functions			
. ../0.common/common.fnc			
# --------------------------------------------------			
if [ -z ${MY_CULTIVAR_NAME} ]; then			
	MY_CULTIVAR_NAME=`Set_MY_CULTIVAR_NAME`		
fi			
			
			
if [ -z ${QVAL} ]; then			
	QVAL=`Set_READ_QVAL_MY_CULTIVAR`		
fi			
			
if [ -z ${PVAL} ]; then			
	PVAL=`Set_READ_PVAL_MY_CULTIVAR`		
fi			
			
if [ -z ${MIS_MATCH_DEFAULT} ]; then			
	MIS_MATCH_DEFAULT=`Set_MIS_MATCH_FOR_MAKE_CONSENSUS`		
fi			
# --------------------------------------------------			
if [ -z $1 ];then			
	MIS_MATCH=${MIS_MATCH_DEFAULT}		
else			
	MIS_MATCH=$1		
fi			
# ==================================================			
NAME=${MY_CULTIVAR_NAME}			
			
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
QNAME="${NAME}_q${QVAL}p${PVAL}"			
			
INPATH="30.coval_call"			
INPILEUP="${INPATH}/${QNAME}_MSR_Cov_${MIS_MATCH}_S-snp.pileup"			
			
OUTPATH="40.RYKMSWBDHV_to_ACGT"			
OUTPILEUP=`basename ${INPILEUP} .pileup`			
OUTPILEUP="${OUTPATH}/${OUTPILEUP}_RYKMSWBDHV2ACGT.pileup"			
			
CMD="${TOPPATH_SCRIPTS}/2./RYKMSWBDHV_to_ACGT.pl"			
			
echo "${CMD} ${INPILEUP} ${OUTPILEUP}"			
eval "${CMD} ${INPILEUP} ${OUTPILEUP}"
