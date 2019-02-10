#!/bin/sh				
				
# ==================================================				
# Those keys enable users to customize this script				
# ==================================================				
BULK_NAME=""				
IDs=""				
				
QVAL=""				
PVAL=""				
QOPT=""				
				
AUTO_EQUALIZE=1				
# if you set AUTO_EQUALIZE=0, then 				
LANE_A=(5 6)				
LANE_B=(8 7)				
WHICHLARGE=(0 0)			# set 0 if you want reduce A / set 1 if you want reduce B	
SMALLERSIZE=(20964328 21025805)				
				
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
				
if [ -z ${QOPT} ]; then				
	QOPT=`Set_READ_QOPT`			
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
				
# ==================================================				
# if AUTO_EQUALIZE=1				
# ==================================================				
if [ ${AUTO_EQUALIZE} -eq 1 ]; then				
				
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~			
	# update ${LANE_A[@]}, ${LANE_B[@]} and get each size			
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~			
				
	NAME_TAIL="_[0-9]*.RnFq${QVAL}p${PVAL}S.1.gz"			
#-	NAME_TAIL="_[0-9]*.RnFq30p90S.1.gz"			
				
	LANE_A=()			
	LANE_B=()			
	SIZE_A=()			
	SIZE_B=()			
				
	_i_=0			
	while [ ${_i_} -le 1 ]			
	do			
		READs=`ls ${NAME}_${IDs[${_i_}]}/q${QVAL}p${PVAL}/sep_pair/${NAME}_${IDs[${_i_}]}${NAME_TAIL}`		
#-		READs=`ls ${NAME}_${IDs[${_i_}]}/q30p90/sep_pair/${NAME}_${IDs[${_i_}]}${NAME_TAIL}`		
				
		echo ${READs}		
		for rd in ${READs}		
		do		
				
				
				
#                                           #${NAME}_${IDs[${_i_}]}/q${QVAL}p${PVAL}/sep_pair/${NAME}_${IDs[${_i_}]}_[0-9]*.RnFq${QVAL}p${PVAL}S.1.gz				
			rd_basename=`basename ${rd}`    # -> ${NAME}_${IDs[${_i_}]}_[0-9]*.RnFq${QVAL}p${PVAL}S.1.gz	
			mylane=${rd_basename##*_}        # -> [0-9]*.RnFq${QVAL}p${PVAL}S.1.gz	
			mylane=${mylane%%.*}            # -> [0-9]*	
			echo ${mylane}	
			sizewc=`gunzip -c ${rd} | wc -l`	
			sizewc=`expr ${sizewc} / 4`	
			echo ${sizewc}	
			if [ ${_i_} -eq 0 ]; then	
				LANE_A=(${LANE_A[@]} ${mylane})
				SIZE_A=(${SIZE_A[@]} ${sizewc})
			else	
				LANE_B=(${LANE_B[@]} ${mylane})
				SIZE_B=(${SIZE_B[@]} ${sizewc})
			fi	
		done		
				
		_i_=`expr ${_i_} + 1`		
	done			
				
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~			
	# update ${WHICHLARGE[@]}, ${SMALLERSIZE[@]}			
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~			
	WHICHLARGE=()			
	SMALLERSIZE=()			
	_j_=0			
	while [ ${_j_} -lt ${#LANE_A[@]} ]			
	do			
		echo ${LANE_A[${_j_}]}		
		echo ${LANE_B[${_j_}]}		
				
		echo ${SIZE_A[${_j_}]}		
		echo ${SIZE_B[${_j_}]}		
				
		if [ ${SIZE_A[${_j_}]} -lt ${SIZE_B[${_j_}]} ]; then		
			WHICHLARGE=(${WHICHLARGE[@]} 1)	
			SMALLERSIZE=(${SMALLERSIZE[@]} ${SIZE_A[${_j_}]})	
		else		
			WHICHLARGE=(${WHICHLARGE[@]} 0)	
			SMALLERSIZE=(${SMALLERSIZE[@]} ${SIZE_B[${_j_}]})	
		fi		
				
		_j_=`expr ${_j_} + 1`		
	done			
fi				
				
echo "==== LANE_A ===="				
_j_=0				
while [ ${_j_} -lt ${#LANE_A[@]} ]				
				
do				
	echo ${LANE_A[${_j_}]}			
	_j_=`expr ${_j_} + 1`			
done				
				
echo "==== LANE_B ===="				
_j_=0				
while [ ${_j_} -lt ${#LANE_B[@]} ]				
				
do				
	echo ${LANE_B[${_j_}]}			
	_j_=`expr ${_j_} + 1`			
done				
				
echo "==== WHICHLARGE ===="				
_j_=0				
while [ ${_j_} -lt ${#WHICHLARGE[@]} ]				
				
do				
	echo ${WHICHLARGE[${_j_}]}			
	_j_=`expr ${_j_} + 1`			
done				
				
echo "==== SMALLERSIZE ===="				
_j_=0				
while [ ${_j_} -lt ${#SMALLERSIZE[@]} ]				
				
do				
	echo ${SMALLERSIZE[${_j_}]}			
	_j_=`expr ${_j_} + 1`			
done				
				
# ==================================================				
# Equalizing process				
# ==================================================				
				
OUTPATH="equalized"				
				
for idid in ${IDs[@]}				
				
do				
	echo "mkdir -p ${OUTPATH}/${NAME}_${idid}"			
	eval "mkdir -p ${OUTPATH}/${NAME}_${idid}"			
done				
				
_i_=0				
for mywhich in ${WHICHLARGE[@]}				
				
do				
	if [ ${mywhich} -eq 0 ]; then			
		MY_ID=${IDs[0]}		
		my_id=${IDs[1]}		
				
		MY_LANE=${LANE_A[${_i_}]}		
		my_lane=${LANE_B[${_i_}]}		
	else			
		MY_ID=${IDs[1]}		
		my_id=${IDs[0]}		
				
		MY_LANE=${LANE_B[${_i_}]}		
		my_lane=${LANE_A[${_i_}]}		
	fi			
				
	# --------------------------------------------------			
	# upper char : reduce read size 			
	# lower char : not reduce read size 			
	# --------------------------------------------------			
	MY_PATH="${NAME}_${MY_ID}/q${QVAL}p${PVAL}/sep_pair"			
	MY_NAME="${MY_PATH}/${NAME}_${MY_ID}"			
	my_path="${NAME}_${my_id}/q${QVAL}p${PVAL}/sep_pair"			
	my_name="${my_path}/${NAME}_${my_id}"			
				
	INPAIREDFQ="${MY_NAME}_${MY_LANE}.RnFq${QVAL}p${PVAL}S.1.gz ${MY_NAME}_${MY_LANE}.RnFq${QVAL}p${PVAL}S.2.gz"			
	inpairedfq="${my_name}_${my_lane}.RnFq${QVAL}p${PVAL}S.1.gz ${my_name}_${my_lane}.RnFq${QVAL}p${PVAL}S.2.gz"			
				
	NEWSIZE=${SMALLERSIZE[${_i_}]}			
	# --------------------------------------------------			
	CMD="${TOPPATH_SCRIPTS}/1./reduce_read.pl"			
	CMD="${CMD} ${INPAIREDFQ} ${NEWSIZE}"			
				
	echo ${CMD}			
	eval ${CMD}			
	# --------------------------------------------------			
				
	for infile in ${INPAIREDFQ}			
	do			
		LINKFILE="${infile}.reduced-${NEWSIZE}.gz"		
		LINKFILE="../../${LINKFILE}"		
		DESTPATH=`basename ${LINKFILE} .reduced-${NEWSIZE}.gz`		
		DESTPATH="${NAME}_${MY_ID}/equalized_${DESTPATH}"		
				
		echo "ln -s ${LINKFILE} ${OUTPATH}/${DESTPATH}"		
		eval "ln -s ${LINKFILE} ${OUTPATH}/${DESTPATH}"		
	done			
				
	for infile in ${inpairedfq}			
	do			
		linkfile="${infile}"		
		linkfile="../../${linkfile}"		
		destpath=`basename ${linkfile}`		
		destpath="${NAME}_${my_id}/equalized_${destpath}"		
				
		echo "ln -s ${linkfile} ${OUTPATH}/${destpath}"		
		eval "ln -s ${linkfile} ${OUTPATH}/${destpath}"		
	done			
				
	_i_=`expr ${_i_} + 1`			
done

