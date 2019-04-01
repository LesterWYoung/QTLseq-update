#! /bin/sh		
		
# ==================================================		
# Those keys enable users to customize this script		
# ==================================================		
BULK_NAME=""		
IDs=""		
MISs=""		
MIN_DEPTH=""		
		
MODE=""		
INDIVIDUALS=""		
NUM_OF_TRIALS=""		
CUTOFF=""		
PREVIOUS_PATH=""		
POPULATION=""		
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
	MISs=`Set_CONFINTRVL_MISs`	
fi		
		
if [ -z ${MIN_DEPTH} ]; then		
	MIN_DEPTH=`Set_MINDEPTH`	
fi		
		
if [ -z ${MODE} ]; then		
	MODE=`Set_CONFINTRVL_CALC_MODE`	
fi		
		
if [ -z ${INDIVIDUALS} ]; then		
	INDIVIDUALS=`Set_CONFINTRVL_INDIVIDUALS`	
fi		
		
if [ -z ${NUM_OF_TRIALS} ]; then		
	NUM_OF_TRIALS=`Set_CONFINTRVL_NUM_OF_TRIALS`	
fi		
		
if [ -z ${CUTOFF} ]; then		
	CUTOFF=`Set_CONFINTRVL_CUTOFF`	
fi		
		
if [ -z ${PREVIOUS_PATH} ]; then		
	PREVIOUS_PATH=`Set_CONFINTRVL_PREVIOUS_PATH`	
fi		
		
if [ -z ${POPULATION} ]; then		
	POPULATION=`Set_CONFINTRVL_POPULATION`	
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
NAME="${BULK_NAME}"		
NAME_A="${NAME}_${IDs[0]}_q${QVAL}p${PVAL}"		
NAME_B="${NAME}_${IDs[1]}_q${QVAL}p${PVAL}"		
NAME_AB="${NAME}_${IDs[0]}${IDs[1]}_q${QVAL}p${PVAL}"		
NAME="merge_${NAME_AB}_paired"		
		
DEPTHs=${MIN_DEPTH}		
		
# --------------------------------------------------		
# internal Constants		
# --------------------------------------------------		
UsePreviousOne=0		
DoCalculateNow=1		
# --------------------------------------------------		
		
SRCPATH="../4.search_for_pair/40.merge_paired"		
OUTPATH="10.cbind_confidence_interval"		
		
CMD1="Rscript"		
CMD1="${CMD1} ${TOPPATH_SCRIPTS}/5./qtl_seq_sumilation_v4.R"		
if [ ${MODE} -eq ${DoCalculateNow} ]; then		
	echo "${CMD1} ${INDIVIDUALS} ${NUM_OF_TRIALS} ${CUTOFF} ${POPULATION}"	
	eval "${CMD1} ${INDIVIDUALS} ${NUM_OF_TRIALS} ${CUTOFF} ${POPULATION}"	
		
	echo "mv ${INDIVIDUALS}individuals.txt ${OUTPATH}/"	
	eval "mv ${INDIVIDUALS}individuals.txt ${OUTPATH}/"	
	PREVIOUS_PATH="${OUTPATH}"	
fi		
		
CMD2="${TOPPATH_SCRIPTS}/5./cbind_simulated_interval_to_pileup.pl"		
CMD2="${CMD2} ${PREVIOUS_PATH}/${INDIVIDUALS}individuals.txt"		
		
for mis in ${MISs}		
do		
	DESTPATH="${OUTPATH}/mut_index_${mis}"	
	echo "mkdir -p ${DESTPATH}"	
	eval "mkdir -p ${DESTPATH}"	
		
	for dep in ${DEPTHs}	
	do	
		INPILEUP="${SRCPATH}/mut_index_${mis}/${NAME}_cov${mis}_co${dep}.txt"
		OUTPILEUP="${DESTPATH}/${NAME}_pvalue_cov${mis}_co${dep}.txt"
		
		echo "${CMD2} ${INPILEUP} ${OUTPILEUP}"
		eval "${CMD2} ${INPILEUP} ${OUTPILEUP}"
		
	done	
done

