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
. ../../0.common/common.fnc				
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
QNAME="${NAME}_q${QVAL}p${PVAL}"				
				
SRC_PATH="../50.make_consensus"				
SRC_NAME="${QNAME}_MSR_Cov_${MIS_MATCH}_S-snp_RYKMSWBDHV2ACGT.fa"				
OUTPATH="00.reference"				
# ==================================================				
				
echo "ln -s ../${SRC_PATH}/${SRC_NAME} ${OUTPATH}/"				
eval "ln -s ../${SRC_PATH}/${SRC_NAME} ${OUTPATH}/"

