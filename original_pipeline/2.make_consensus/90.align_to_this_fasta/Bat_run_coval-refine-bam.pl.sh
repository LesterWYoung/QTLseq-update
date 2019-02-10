#! /bin/sh		
		
# ==================================================		
# Those keys enable users to customize this script		
# ==================================================		
REF=""		
MY_CULTIVAR_NAME=""		
MIS_MATCH_DEFAULT=""		
DONTREALIGN=""		
		
QVAL=""		
PVAL=""		
# --------------------------------------------------		
# load common functions		
. ../../0.common/common.fnc		
# --------------------------------------------------		
if [ -z ${REF} ]; then		
	REF=`Set_REF_FASTA`	
fi		
		
if [ -z ${MY_CULTIVAR_NAME} ]; then		
	MY_CULTIVAR_NAME=`Set_MY_CULTIVAR_NAME`	
fi		
		
if [ -z ${QVAL} ]; then		
	QVAL=`Set_READ_QVAL_MY_CULTIVAR`	
fi		
		
if [ -z ${PVAL} ]; then		
	PVAL=`Set_READ_PVAL_MY_CULTIVAR`	
fi		
		
if [ -z ${DONTREALIGN} ]; then		
	DONTREALIGN=`Set_DONTREALIGN_FOR_MAKE_CONSENSUS`	
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
		
# --------------------------------------------------		
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
Command=`basename $0`		
Logf="log.${Command}.${QNAME}.${MIS_MATCH}.txt"		
		
# --------------------------------------------------		
INPATH="10.bwa2bam/${NAME}"		
INBAM="${INPATH}/${QNAME}_MSR.bam"		
		
OUTPATH="20.coval_refine"		
OUTPREF=`basename ${INBAM} .bam`		
OUTPREF="${OUTPATH}/${OUTPREF}_Cov_${MIS_MATCH}_S"		
		
OPTION=""		
OPTION="${OPTION} -ref ${REF}"		
OPTION="${OPTION} -pref ${OUTPREF}"		
OPTION="${OPTION} -n ${MIS_MATCH}" 		
OPTION="${OPTION} ${DONTREALIGN}"		
		
CMD="${TOPPATH_COVAL}/coval refine"		
CMD="${CMD} ${OPTION} ${INBAM}"		
		
# start=`stop_watch.pl start`		
		
echo "${CMD}"		
${CMD} 2>&1 | tee -a ${OUTPATH}/$Logf		
		
# stop_watch.pl $start
