#! /bin/sh			
			
Command=`basename $0`			
Logf="log.$Command.txt"			
			
# ==================================================			
# Those keys enable users to customize this script			
# ==================================================			
REF=""			
CPU=""			
MY_CULTIVAR_NAME=""			
			
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
	
if [ -z ${CPU} ]; then			
	CPU=`Set_BWA_CPU`		
fi			
			
NAME=${MY_CULTIVAR_NAME}			
			
#==================================================			
Enabled=1			
Disabled=0			
			
remove_sais=${Enabled}			
remove_sams=${Enabled}			
remove_bam0s=${Enabled}			
add_bai=${Enabled}			
add_unmap=${Enabled}			
			
#==================================================			
QNAME="${NAME}_q${QVAL}p${PVAL}"			
			
OUTPATH="10.bwa2bam"			
#--------------------------------------------------			
			
#--------------------------------------------------			
BWA_INDEX="bwa index"			
CHK_EXIST="${REF}.bwt"			
if [ -e ${CHK_EXIST} ];then			
	echo "${CHK_EXIST} exists! continue."		
else			
	PREFIX=${REF}		
			
	echo "${BWA_INDEX} -p ${PREFIX} -a is ${REF}"		
	eval "${BWA_INDEX} -p ${PREFIX} -a is ${REF}"		
			
fi			
			
#--------------------------------------------------			
			
#--------------------------------------------------			
BWA_ALN="bwa aln"			
OPTIONS=""			
OPTIONS="${OPTIONS} -n 0.04"			
OPTIONS="${OPTIONS} -o 1"			
OPTIONS="${OPTIONS} -k 2"			
OPTIONS="${OPTIONS} -t ${CPU}"			
			
READS=`ls ${OUTPATH}/${NAME}/*_[12]_sequence.txt.gz`			
			
_i_=0			
myREADs=()			
mySAIs=()			
for myread in ${READS}			
do			
			
	mysai=`basename ${myread}`		
	mysai="${OUTPATH}/${NAME}/${mysai}.sai"		
			
	if [ -e ${mysai} ];then		
		echo "${mysai} exists! continue."	
	else		
		echo "${BWA_ALN} ${OPTIONS} ${REF} ${myread} > ${mysai}"	
		eval "${BWA_ALN} ${OPTIONS} ${REF} ${myread} > ${mysai}"	
	fi		
			
	myREADs[${_i_}]="${myread}"		
	mySAIs[${_i_}]="${mysai}"		
			
	_i_=`expr ${_i_} + 1`		
done			
			
num_of_read=`expr ${_i_} - 1`			
num_of_pair=`expr ${_i_} / 2 - 1`			
			
#--------------------------------------------------			
			
#--------------------------------------------------			
BWA_SAMPE="bwa sampe"			
OPTIONS=""			
OPTIONS="${OPTIONS} -a 1000"			
OPTIONS="${OPTIONS} -n 3"			
OPTIONS="${OPTIONS} -N 10"			
			
mySAM0s=()			
_p_=0			
while [ ${_p_} -le ${num_of_pair} ]			
do			
	mySAM0s[${_p_}]="${OUTPATH}/${NAME}/${QNAME}_${_p_}.sam"		
	_j0_=`expr ${_p_} \* 2 + 0`		
	_j1_=`expr ${_p_} \* 2 + 1`		
	mypairread="${myREADs[${_j0_}]} ${myREADs[${_j1_}]}"		
	mypairsai="${mySAIs[${_j0_}]} ${mySAIs[${_j1_}]}"		
			
	echo "${BWA_SAMPE} ${OPTIONS} ${REF} ${mypairsai} ${mypairread} > ${mySAM0s[${_p_}]}"		
	eval "${BWA_SAMPE} ${OPTIONS} ${REF} ${mypairsai} ${mypairread} > ${mySAM0s[${_p_}]}"		
			
	_p_=`expr ${_p_} + 1`		
done			
			
#--------------------------------------------------			
			
#--------------------------------------------------			
SAMTOOLS_VIEW="samtools view"			
OPTIONS=""			
OPTIONS="-Sb"			
			
myBAM0s=()			
_p_=0			
while [ ${_p_} -le ${num_of_pair} ]			
do			
	mySAM=${mySAM0s[${_p_}]}		
	myBAM=`basename ${mySAM} .sam`		
	myBAM="${myBAM}.bam"		
	myBAM="${OUTPATH}/${NAME}/${myBAM}"		
			
	myBAM0s[${_p_}]=${myBAM}		
			
	echo "${SAMTOOLS_VIEW} ${OPTIONS} ${mySAM} > ${myBAM}"		
	eval "${SAMTOOLS_VIEW} ${OPTIONS} ${mySAM} > ${myBAM}"		
			
	_p_=`expr ${_p_} + 1`		
done			
			
#--------------------------------------------------			
			
#--------------------------------------------------			
myBAMall="${OUTPATH}/${NAME}/${QNAME}_all.bam"			
			
if [ ${num_of_pair} -eq 0 ]; then			
	_p_=0		
	echo "mv ${myBAM0s[${_p_}]} ${myBAMall}"		
	eval "mv ${myBAM0s[${_p_}]} ${myBAMall}"		
else			
	SAMTOOLS_MERGE="samtools merge"		
	OPTIONS=""		
			
	myBAMall="${OUTPATH}/${NAME}/${QNAME}_all.bam"		
	myBAMs=""		
	_p_=0		
	while [ ${_p_} -le ${num_of_pair} ]		
	do		
		myBAMs="${myBAMs} ${myBAM0s[${_p_}]}"	
		_p_=`expr ${_p_} + 1`	
	done		
			
	echo "${SAMTOOLS_MERGE} ${OPTIONS} ${myBAMall} ${myBAMs}"		
	eval "${SAMTOOLS_MERGE} ${OPTIONS} ${myBAMall} ${myBAMs}"		
fi			
			
#--------------------------------------------------			
			
#--------------------------------------------------			
SAMTOOLS_VIEW="samtools view"			
OPTIONS=""			
OPTIONS="-F 0x0004"			
OPTIONS="${OPTIONS} -bu"			
			
myBAMmap="${OUTPATH}/${NAME}/${QNAME}_mapped.bam"			
			
echo "${SAMTOOLS_VIEW} ${OPTIONS} ${myBAMall} > ${myBAMmap}"			
eval "${SAMTOOLS_VIEW} ${OPTIONS} ${myBAMall} > ${myBAMmap}"			
			
#--------------------------------------------------			
			
#--------------------------------------------------			
SAMTOOLS_SORT="samtools sort"			
OPTIONS=""			
			
myBAMums="${OUTPATH}/${NAME}/${QNAME}_UMS.bam"			
myBAMumsPrfx=`basename ${myBAMums} .bam`			
myBAMumsPrfx="${OUTPATH}/${NAME}/${myBAMumsPrfx}"			
			
echo "${SAMTOOLS_SORT} ${OPTIONS} ${myBAMmap} ${myBAMumsPrfx}"			
eval "${SAMTOOLS_SORT} ${OPTIONS} ${myBAMmap} ${myBAMumsPrfx}"			
			
#--------------------------------------------------			
			
#--------------------------------------------------			
SAMTOOLS_RMDUP="samtools rmdup"			
OPTIONS=""			
			
myBAMmsr="${OUTPATH}/${NAME}/${QNAME}_MSR.bam"			
			
echo "${SAMTOOLS_RMDUP} ${OPTIONS} ${myBAMums} ${myBAMmsr}"			
eval "${SAMTOOLS_RMDUP} ${OPTIONS} ${myBAMums} ${myBAMmsr}"			
			
#--------------------------------------------------			
			
#--------------------------------------------------			
echo "rm -f ${myBAMmap}"			
eval "rm -f ${myBAMmap}"			
echo "rm -f ${myBAMums}"			
eval "rm -f ${myBAMums}"			
			
#--------------------------------------------------			
			
#--------------------------------------------------			
if [ ${add_bai} -eq ${Enabled} ]; then			
	SAMTOOLS_INDEX="samtools index"		
	OPTIONS=""		
	echo "${SAMTOOLS_INDEX} ${OPTIONS} ${myBAMmsr}"		
	eval "${SAMTOOLS_INDEX} ${OPTIONS} ${myBAMmsr}"		
fi			
			
#--------------------------------------------------			
if [ ${add_unmap} -eq ${Enabled} ]; then			
	SAMTOOLS_VIEW="samtools view"		
	OPTIONS=""		
	OPTIONS="-f 0x0004"		
	OPTIONS="${OPTIONS} -b"		
			
	myBAMunmap="${OUTPATH}/${NAME}/${QNAME}_unmp.bam"		
			
	echo "${SAMTOOLS_VIEW} ${OPTIONS} ${myBAMall} > ${myBAMunmap}"		
	eval "${SAMTOOLS_VIEW} ${OPTIONS} ${myBAMall} > ${myBAMunmap}"		
			
	echo "rm -f ${myBAMall}"		
	eval "rm -f ${myBAMall}"		
fi			
			
#--------------------------------------------------			
if [ ${remove_sais} -eq ${Enabled} ]; then			
	_p_=0		
	while [ ${_p_} -le ${num_of_pair} ]		
	do		
		_j0_=`expr ${_p_} \* 2 + 0`	
		_j1_=`expr ${_p_} \* 2 + 1`	
		mypairsai="${mySAIs[${_j0_}]} ${mySAIs[${_j1_}]}"	
			
		echo "rm -f ${mypairsai}"	
		eval "rm -f ${mypairsai}"	
			
		_p_=`expr ${_p_} + 1`	
	done		
fi			
			
#--------------------------------------------------			
if [ ${remove_sams} -eq ${Enabled} ]; then			
	_p_=0		
	while [ ${_p_} -le ${num_of_pair} ]		
	do		
		echo "rm -f ${mySAM0s[${_p_}]}"	
		eval "rm -f ${mySAM0s[${_p_}]}"	
		_p_=`expr ${_p_} + 1`	
	done		
fi			
			
#--------------------------------------------------			
if [ ${remove_bam0s} -eq ${Enabled} ]; then			
	_p_=0		
	while [ ${_p_} -le ${num_of_pair} ]		
	do		
		echo "rm -f ${myBAM0s[${_p_}]}"	
		eval "rm -f ${myBAM0s[${_p_}]}"	
		_p_=`expr ${_p_} + 1`	
	done		
fi
