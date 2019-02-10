#! /bin/sh			
			
Command=`basename $0`			
Logf="log.$Command.txt"			
			
# ==================================================			
# Those keys enable users to customize this script			
# ==================================================			
BULK_NAME=""			
IDs=""			
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
			
REF_FASTA=""			
			
if [ -z ${REF_FASTA} ]; then			
	REF_FASTA=`Set_REF_FASTA`		
fi			
			
			
# ==================================================			
NAME="${BULK_NAME}"			
NAME_A="${NAME}_${IDs[0]}_q${QVAL}p${PVAL}"			
NAME_B="${NAME}_${IDs[1]}_q${QVAL}p${PVAL}"			
NAME_AB="${NAME}_${IDs[0]}${IDs[1]}_q${QVAL}p${PVAL}"			
			
MIS=4			
if [ -z $1 ];then			
	echo "missing args!"		
	echo "[usage]"		
	echo "    $0 <INT1>"		
	echo "        <INT1> : MISMATCH"		
	exit 0		
else			
	MIS=$1		
fi			
# ==================================================			
			
SRCPATH="10.paired_or_unpaired/mut_index_${MIS}"			
OUTPATH="40.merge_paired"			
OUTPATH="${OUTPATH}/mut_index_${MIS}"			
			
echo "mkdir -p ${OUTPATH}"			
eval "mkdir -p ${OUTPATH}"			
			
PAIRED_A="${SRCPATH}/paired_${NAME_A}_cov${MIS}_co${MIN_DEPTH}.txt"			
PAIRED_B="${SRCPATH}/paired_${NAME_B}_cov${MIS}_co${MIN_DEPTH}.txt"			
			
PAIRED_AB="${OUTPATH}/paste_paired_A_and_paired_B_cov${MIS}_co${MIN_DEPTH}.txt"			
echo "paste ${PAIRED_A} ${PAIRED_B} > ${PAIRED_AB}"			
eval "paste ${PAIRED_A} ${PAIRED_B} > ${PAIRED_AB}"			
			
# --------------------------------------------------			
			
			
SRCPATH="30.unpaired_to_paired/mut_index_${MIS}"			
OUTPATH="40.merge_paired"			
OUTPATH="${OUTPATH}/mut_index_${MIS}"			
			
CMD="${TOPPATH_SCRIPTS}/4./merge_paired.pl"			
			
REF_FAI="${REF_FASTA}.fai"			
			
UNPAIRED_A="${SRCPATH}/paste_unpaired_${NAME_A}_and_${NAME_B}_cov${MIS}_co${MIN_DEPTH}.txt"			
UNPAIRED_B="${SRCPATH}/paste_${NAME_A}_and_unpaired_${NAME_B}_cov${MIS}_co${MIN_DEPTH}.txt"			
			
AWK_PRINTF="awk '{printf(\"%s\\t%d\\t%s\\t%s\\t%d\\t%d\\t%d\\t%d\\t%.2f\\t%s\\t%s\\t%s\\t%d\\t%s\\t%s\\t%d\\t%d\\t%d\\t%d\\t%.2f\\t%s\\t%s\\n\",\$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$11,\$9,\$10,\$12,\$13,\$14,\$15,\$16,\$17,\$18,\$19,\$22,\$20,\$21)}'"			
			
MERGED_PAIRED_AB_awk="${OUTPATH}/merge_${NAME_AB}_paired_cov${MIS}_co${MIN_DEPTH}.awk.txt"			
MERGED_PAIRED_AB="${OUTPATH}/merge_${NAME_AB}_paired_cov${MIS}_co${MIN_DEPTH}.txt"			
			
echo "${CMD} ${REF_FAI} ${PAIRED_AB} ${UNPAIRED_A} ${UNPAIRED_B} ${MERGED_PAIRED_AB_awk}"			
eval "${CMD} ${REF_FAI} ${PAIRED_AB} ${UNPAIRED_A} ${UNPAIRED_B} ${MERGED_PAIRED_AB_awk}"			
			
echo "cat ${MERGED_PAIRED_AB_awk} | ${AWK_PRINTF} > ${MERGED_PAIRED_AB}"			
eval "cat ${MERGED_PAIRED_AB_awk} | ${AWK_PRINTF} > ${MERGED_PAIRED_AB}"			
			
			
rm ${PAIRED_AB}			
rm ${MERGED_PAIRED_AB_awk}

