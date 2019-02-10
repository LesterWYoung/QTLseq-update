#! /bin/sh			
			
Command=`basename $0`			
Logf="log.${Command}.txt"			
			
# ==================================================			
# Those keys enable users to customize this script			
# ==================================================			
REF=""			
MY_CULTIVAR_NAME=""			
MIS_MATCH_DEFAULT=""
QVAL=""			
PVAL=""			
# --------------------------------------------------			
# load common functions			
. ../0.common/common.fnc			
# --------------------------------------------------			
if [ -z ${REF} ]; then			
	REF=`Set_PUBLIC_REF_FASTA`		
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
			
INPATH="40.RYKMSWBDHV_to_ACGT"			
INPILEUP="${INPATH}/${QNAME}_MSR_Cov_${MIS_MATCH}_S-snp_RYKMSWBDHV2ACGT.pileup"			
			
OUTPATH="50.make_consensus"			
OUTFA=`basename ${INPILEUP} .pileup`			
OUTFA="${OUTPATH}/${OUTFA}.fa"			
			
CMD="${TOPPATH_SCRIPTS}/2./make_consensus.pl"			
CMD="${CMD} -ref ${REF} ${INPILEUP} > ${OUTFA}"			
			
# start=`stop_watch.pl start`			
			
echo ${CMD}			
eval ${CMD}			
			
# stop_watch.pl $start
