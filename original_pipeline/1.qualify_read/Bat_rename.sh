#! /bin/sh			
			
# ==================================================			
# Those keys enable users to customize this script			
# ==================================================			
MY_CULTIVAR_NAME=""			
BULK_NAME=""			
IDs=""			
			
# --------------------------------------------------			
# load common functions			
 . ../0.common/common.fnc			
# --------------------------------------------------			
if [ -z ${MY_CULTIVAR_NAME} ]; then			
	MY_CULTIVAR_NAME=`Set_MY_CULTIVAR_NAME`		
fi			
			
if [ -z ${BULK_NAME} ]; then			
	BULK_NAME=`Set_BULK_NAME`		
fi			
			
if [ -z ${IDs} ]; then			
	IDA=`Set_BULK_NAME_IDA`		
	IDB=`Set_BULK_NAME_IDB`		
	IDs=(${IDA} ${IDB})		
fi			
			
# ==================================================			
if [ -z $1 ];then			
	echo "missing args!"		
	echo "[usage]"		
	echo "    $0 <INT1>"		
	echo "        <INT1> : 0 or 1 or 9(0:${IDs[0]} / 1:${IDs[1]} / 9:${MY_CULTIVAR_NAME})"		
	exit 0		
else			
	if [ $1 -eq 0 -o $1 -eq 1 -o $1 -eq 9 ];then		
		myID=$1	
	else		
		echo "invalid args!"	
		echo "[usage]"	
		echo " $0 <INT1>"	
		echo " <INT1> : 0 or 1 or 9(0:${IDs[0]} / 1:${IDs[1]} / 9:${MY_CULTIVAR_NAME})"	
		exit 0	
	fi		
fi			
			
# ==================================================			
if [ $1 -eq 0 -o $1 -eq 1 ];then			
	NAME="${BULK_NAME}_${IDs[${myID}]}"		
else			
	NAME="${MY_CULTIVAR_NAME}"		
fi			
			
# --------------------------------------------------			
			
INPATH=${NAME}			
OUTPATH=${NAME}			
			
			
# --------------------------------------------------			
# in case of demultiplexing fastq which has a name like as p_1_GACTTG_1_sequence.txt.gz			
# --------------------------------------------------			
NAME_TAIL="_[0-9]*_[12]_sequence.txt.gz"			
READs=`ls ${INPATH}/*${NAME_TAIL}`			
			
for rd in ${READs}			
do			
#                                   #${INPATH}/*_[0-9]*_[12]_sequence.txt.gz			
    rd_basename=`basename ${rd}`    # -> *_[0-9]*_[12]_sequence.txt.gz			
    mylane=${rd_basename#*_}        # -> [0-9]*_[12]_sequence.txt.gz			
    mylane=${mylane%%_*}            # -> [0-9]*			
	my1or2=${rd_basename%_*}        # -> *_[0-9]*_[12]		
    my1or2=${my1or2##*_}            # -> [12]			
			
	echo "ln -s ${rd_basename} ${OUTPATH}/${NAME}_${mylane}.Rn.${my1or2}.gz"		
	eval "ln -s ${rd_basename} ${OUTPATH}/${NAME}_${mylane}.Rn.${my1or2}.gz"		
done

