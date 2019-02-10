#! /bin/sh				
				
# ==================================================				
# Those keys enable users to customize this script				
# ==================================================				
BULK_NAME=""				
IDs=""				
MISs=""				
DEPTHs=""				
				
QVAL=""				
PVAL=""				
				
MIN_DEPTH=""				
MIN_SNPINDEX=""				
MAX_SNPINDEX=""				
				
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
				
if [ -z ${MIN_DEPTH} ]; then				
	MIN_DEPTH=`Set_MINDEPTH`			
fi				
				
if [ -z ${MISs} ]; then				
	MISs=`Set_CONFINTRVL_MISs`			
fi				
				
if [ -z ${DEPTHs} ]; then				
	DEPTHs=`Set_AWK_CUSTOM2_DEPTHs`			
fi				
				
if [ -z ${MIN_SNPINDEX} ]; then				
	MIN_SNPINDEX=`Set_AWK_CUSTOM2_MIN_SNPINDEX`			
fi				
				
if [ -z ${MAX_SNPINDEX} ]; then				
	MAX_SNPINDEX=`Set_AWK_CUSTOM2_MAX_SNPINDEX`			
fi				
				
# ==================================================				
				
				
# --------------------------------------------------				
NAME="${BULK_NAME}"				
NAME_A="${NAME}_${IDs[0]}_q${QVAL}p${PVAL}"				
NAME_B="${NAME}_${IDs[1]}_q${QVAL}p${PVAL}"				
NAME_AB="${NAME}_${IDs[0]}${IDs[1]}_q${QVAL}p${PVAL}"				
NAME="merge_${NAME_AB}_paired_pvalue"				
				
# --------------------------------------------------				
				
MKDIR="mkdir -p"				
SRCPATH="10.cbind_confidence_interval"				
OUTPATH="20.awk_custom2"				
				
CMD0="cat"				
CMD1="awk '(\$9 >= ${MIN_SNPINDEX}) || (\$20 >= ${MIN_SNPINDEX})'"				
CMD2="awk '(\$9 <= ${MAX_SNPINDEX}) || (\$20 <= ${MAX_SNPINDEX})'"				
				
for mis in ${MISs}				
do				
	echo "${MKDIR} ${OUTPATH}/mut_index_${mis}/"			
	eval "${MKDIR} ${OUTPATH}/mut_index_${mis}/"			
				
	for dep in ${DEPTHs}			
	do			
		inpileup="${SRCPATH}/mut_index_${mis}/${NAME}_cov${mis}_co${MIN_DEPTH}.txt"		
		outpileup=`basename ${inpileup} _co${MIN_DEPTH}.txt`		
		outpileup="filtered_${outpileup}_co${dep}.txt"		
				
		AWKDEPTH="awk '\$8>=${dep} && \$19>=${dep}'"		
				
		echo "${CMD0} ${inpileup} | ${AWKDEPTH} | ${CMD1} | ${CMD2} > ${OUTPATH}/mut_index_${mis}/${outpileup}"		
		eval "${CMD0} ${inpileup} | ${AWKDEPTH} | ${CMD1} | ${CMD2} > ${OUTPATH}/mut_index_${mis}/${outpileup}"		
	done			
done
