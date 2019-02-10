#! /bin/sh			
			
# ==================================================			
# Those keys enable users to customize this script			
# ==================================================			
REF=""			
MY_CULTIVAR_NAME=""			
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
			
if [ -z ${MINNUM} ]; then			
	MINNUM=`Set_COVAL_CALL_MINNUM_FOR_MAKE_CONSENSUS`		
fi			
if [ -z ${MAXR} ]; then			
	MAXR=`Set_COVAL_CALL_MAXR_FOR_MAKE_CONSENSUS`		
fi			
if [ -z ${MINFREQ} ]; then			
	MINFREQ=`Set_COVAL_CALL_MINFREQ_FOR_MAKE_CONSENSUS`		
fi			
if [ -z ${MINTNUM} ]; then			
	MINTNUM=`Set_COVAL_CALL_MINTNUM_FOR_MAKE_CONSENSUS`		
fi			
if [ -z ${MINQUALBASE} ]; then			
	MINQUALBASE=`Set_COVAL_CALL_MINQUALBASE_FOR_MAKE_CONSENSUS`		
fi			
if [ -z ${MINQUALAVE} ]; then			
	MINQUALAVE=`Set_COVAL_CALL_MINQUALAVE_FOR_MAKE_CONSENSUS`		
fi			
if [ -z ${CALLTYPE} ]; then			
	CALLTYPE=`Set_COVAL_CALL_CALLTYPE_FOR_MAKE_CONSENSUS`		
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
			
#--         #2014/04/22 kikuchi			
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