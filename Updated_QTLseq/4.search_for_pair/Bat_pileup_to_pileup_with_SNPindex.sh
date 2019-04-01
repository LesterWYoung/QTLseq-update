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
		
MISs=(${MIS})		
		
for mis in ${MISs}		
do		
	SRCPATH="20.search_for_pair_of_unpaired/mut_index_${MIS}"	
	OUTPATH="mut_index_${mis}"	
	OUTPATH="30.unpaired_to_paired/${OUTPATH}"	
	PASTE="${OUTPATH}/work_${myID}_${MIS}.txt"	
		
	echo "mkdir -p ${OUTPATH}"	
	eval "mkdir -p ${OUTPATH}"	
		
	EXE="${TOPPATH_SCRIPTS}/4./pileup_to_pileup_with_SNPindex.pl"	
		
	if [ ${myID} -eq 0 ]; then	
		my_name0=${NAME_B}
		my_name1=${NAME_A}
	else	
		my_name0=${NAME_A}
		my_name1=${NAME_B}
	fi	
		
                                                                            ##pileup		
	inpileup1="${SRCPATH}/common_pos_${my_name1}_MSR_Cov_${mis}_S.bam.pileup"	
	outpileup1=`basename ${inpileup1}` 	
	outpileup1="${OUTPATH}/${outpileup1}"	
	echo "${EXE} ${inpileup1} ${outpileup1}"	
	eval "${EXE} ${inpileup1} ${outpileup1}"	
                                                                            ##unpaired		
	inpileup0="${SRCPATH}/common_pos_unpaired_${my_name0}_cov${MIS}_co${MIN_DEPTH}.txt"      	
	outpileup0=`basename ${inpileup0}`	
	outpileup0="${OUTPATH}/${outpileup0}"	
	echo "ln -s ../../${inpileup0} ${outpileup0}"	
	eval "ln -s ../../${inpileup0} ${outpileup0}"	
		
		
	if [ ${myID} -eq 0 ]; then	
		echo "paste ${outpileup1} ${outpileup0} > ${PASTE}"
		eval "paste ${outpileup1} ${outpileup0} > ${PASTE}"
		pastepileup="paste_${my_name1}_and_unpaired_${my_name0}_cov${MIS}_co${MIN_DEPTH}.txt"
	else	
		echo "paste ${outpileup0} ${outpileup1} > ${PASTE}"
		eval "paste ${outpileup0} ${outpileup1} > ${PASTE}"
		pastepileup="paste_unpaired_${my_name0}_and_${my_name1}_cov${MIS}_co${MIN_DEPTH}.txt"
	fi	
		
	pastepileup="${OUTPATH}/${pastepileup}"	
		
	DEPTH3="awk '\$8>=${MIN_DEPTH} && \$19>=${MIN_DEPTH}'"	
		
	CMD="cat ${PASTE} | ${DEPTH3} > ${pastepileup}"	
	echo ${CMD}	
	eval ${CMD}	
		
	echo "rm ${PASTE}"	
	eval "rm ${PASTE}"	
		
done

