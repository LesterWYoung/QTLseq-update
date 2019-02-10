#! /bin/sh			
			
# ==================================================			
# Those keys enable users to customize this script			
# ==================================================			
REF=""			
BULK_NAME=""			
IDs=""			
MIS_MATCH_DEFAULT=""	
QVAL=""			
PVAL=""			
# --------------------------------------------------			
# for coval call			
MINNUM=""			
MAXR=""			
MINFREQ=""			
MINTNUM=""			
MINQUALBASE=""			
MINQUALAVE=""			
CALLTYPE=""			
			
# --------------------------------------------------			
# load common functions			
. ../0.common/common.fnc			
# --------------------------------------------------			
			
if [ -z ${REF} ]; then			
	REF=`Set_REF_FASTA`		
fi			
			
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
	
			
if [ -z ${MIS_MATCH_DEFAULT} ]; then			
	MIS_MATCH_DEFAULT=`Set_MIS_MATCH_DEFAULT`		
fi			
			
if [ -z ${MINNUM} ]; then			
	MINNUM=`Set_COVAL_CALL_MINNUM`		
fi			
			
if [ -z ${MAXR} ]; then			
	MAXR=`Set_COVAL_CALL_MAXR`		
fi			
			
if [ -z ${MINFREQ} ]; then			
	MINFREQ=`Set_COVAL_CALL_MINFREQ`		
fi			
			
if [ -z ${MINTNUM} ]; then			
	MINTNUM=`Set_COVAL_CALL_MINTNUM`		
fi			
			
if [ -z ${MINQUALBASE} ]; then			
	MINQUALBASE=`Set_COVAL_CALL_MINQUALBASE`		
fi			
			
if [ -z ${MINQUALAVE} ]; then			
	MINQUALAVE=`Set_COVAL_CALL_MINQUALAVE`		
fi			
			
if [ -z ${CALLTYPE} ]; then			
	CALLTYPE=`Set_COVAL_CALL_CALLTYPE`		
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
			
			
# ==================================================			
			
if [ -z $1 ];then			
	echo "missing args!"		
	echo "[usage]"		
	echo "    $0 <INT1> [<INT2>]"		
	echo "        <INT1> : 0 or 1 (0:${IDs[0]} / 1:${IDs[1]})"		
	echo "        <INT2> : MIS_MATCH (if null, set ${MIS_MATCH})"		
	exit 0		
else			
	if [ $1 -eq 0 -o $1 -eq 1 ];then		
		myID=$1	
	else		
		echo "invalid args!"	
		echo "[usage]"	
		echo "    $0 <INT1> [<INT2>]"	
		echo "        <INT1> : 0 or 1 (0:${IDs[0]} / 1:${IDs[1]})"	
		echo "        <INT2> : MIS_MATCH (if null, set ${MIS_MATCH})"	
		exit 0	
	fi		
fi			
			
NAME="${BULK_NAME}_${IDs[${myID}]}"			
			
if [ -z $2 ];then			
	MIS_MATCH=${MIS_MATCH_DEFAULT}		
else			
	MIS_MATCH=$2		
fi			
			
# --------------------------------------------------			
QNAME="${NAME}_q${QVAL}p${PVAL}"			
			
INPATH="20.coval_refine"			
INBAM="${INPATH}/${QNAME}_MSR_Cov_${MIS_MATCH}_S.bam"			
			
OUTPATH="30.coval_call"			
OUTSNPINDEL=`basename ${INBAM} .bam`			
OUTSNPINDEL="${OUTPATH}/${OUTSNPINDEL}.snpindel"			
			
if [ -f ${OUTSNPINDEL} ]; then			
	echo "${OUTSNPINDEL} exists!"		
else			
	CMD0="samtools pileup -vcf"		
	echo "${CMD0} ${REF} ${INBAM} > ${OUTSNPINDEL}"		
	${CMD0} ${REF} ${INBAM} > ${OUTSNPINDEL}		
fi			
# --------------------------------------------------			
			
OUTPREF=`basename ${INBAM} .bam`			
OUTPREF="${OUTPREF}"			
			
OPTION=""			
OPTION="${OPTION} --pref ${OUTPATH}/${OUTPREF}"			
OPTION="${OPTION} --num ${MINNUM}" 			
OPTION="${OPTION} --maxr ${MAXR}"			
OPTION="${OPTION} --freq ${MINFREQ}" 			
OPTION="${OPTION} --tnum ${MINTNUM}" 			
OPTION="${OPTION} --qual_base ${MINQUALBASE}"			
OPTION="${OPTION} --qual_ave ${MINQUALAVE}"			
OPTION="${OPTION} --calltype ${CALLTYPE}"			
			
CMD1="${TOPPATH_COVAL}/coval call"			
TEELOG="${OUTPATH}/log.${OUTPREF}.txt"			
			
# start=`stop_watch.pl start`			
			
echo "${CMD1} ${OPTION} ${OUTSNPINDEL}"			
${CMD1} ${OPTION} ${OUTSNPINDEL} 2>&1 | tee -a ${TEELOG}			
			
# stop_watch.pl $start			