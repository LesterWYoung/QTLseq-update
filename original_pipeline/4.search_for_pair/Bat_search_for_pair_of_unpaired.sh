#! /bin/sh				
				
Command=`basename $0`				
Logf="log.$Command.txt"				
				
# ==================================================				
# Those keys enable users to customize this script				
# ==================================================				
BULK_NAME=""				
IDs=""				
REF=""		
MINDEPTH=""				
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
				
if [ -z ${REF} ]; then				
	REF=`Set_REF_FASTA`			
fi				
				
if [ -z ${MINDEPTH} ]; then				
	MINDEPTH=`Set_MINDEPTH`			
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
	
				
MIS=4				
if [ -z $1 ];then				
	echo "missing args!"			
	echo "[usage]"			
	echo "    $0 <INT1> [<INT2>]"			
	echo "        <INT1> : 0 or 1 (0:${IDs[0]} / 1:${IDs[1]})"			
	echo "        <INT2> : MISMATCH (if null, set ${MIS})"			
	exit 0			
else				
				
	if [ $1 -eq 0 -o $1 -eq 1 ];then			
		myID=$1		
	else			
		echo "invalid args!"		
		echo "[usage]"		
		echo "    $0 <INT1> [<INT2>]"		
		echo "        <INT1> : 0 or 1 (0:${IDs[0]} / 1:${IDs[1]})"		
		echo "        <INT2> : MISMATCH (if null, set ${MIS})"		
		exit 0		
	fi			
				
fi				
				
if [ -z $2 ];then				
	MIS=${MIS}			
else				
	MIS=$2			
fi				
# ==================================================				
				
SNPPATH="10.paired_or_unpaired/mut_index_${MIS}"				
BAMPATH="../3.alignment/20.coval_refine"				
OUTPATH="20.search_for_pair_of_unpaired"				
OUTPATH="${OUTPATH}/mut_index_${MIS}"				
				
echo "mkdir -p ${OUTPATH}"				
eval "mkdir -p ${OUTPATH}"				
				
# --------------------------------------------------				
# separate common/unigue				
# --------------------------------------------------				
EXE="${TOPPATH_SCRIPTS}/4./comp_snpA_to_bamB_v2.0.pl"				
				
if [ ${myID} -eq 0 ]; then				
	PILEUP1="${SNPPATH}/unpaired_${NAME_B}_cov${MIS}_co${MINDEPTH}.txt"			
	BAM2="${BAMPATH}/${NAME_A}_MSR_Cov_${MIS}_S.bam"			
				
	echo "${EXE} ${PILEUP1} ${BAM2} ${REF} ${OUTPATH}"			
	eval "${EXE} ${PILEUP1} ${BAM2} ${REF} ${OUTPATH}"			
else				
	PILEUP2="${SNPPATH}/unpaired_${NAME_A}_cov${MIS}_co${MINDEPTH}.txt"			
	BAM1="${BAMPATH}/${NAME_B}_MSR_Cov_${MIS}_S.bam"			
				
	echo "${EXE} ${PILEUP2} ${BAM1} ${REF} ${OUTPATH}"			
	eval "${EXE} ${PILEUP2} ${BAM1} ${REF} ${OUTPATH}"			
fi

