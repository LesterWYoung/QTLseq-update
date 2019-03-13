#! /bin/sh

# --------------------------------------------------
# load common functions
. ../0.common/common.fnc
# --------------------------------------------------

MY_CULTIVAR_NAME=`Set_MY_CULTIVAR_NAME`

IDA=`Set_BULK_NAME_IDA`
IDB=`Set_BULK_NAME_IDB`
IDs=(${IDA} ${IDB})

if [ 

if [ -z $1 ];then				
	echo "missing args!"			
	echo "[usage]"			
	echo "    $0 <INT>"			
	echo "        <INT> : 0 or 1 or 9 (0:${IDs[0]} / 1:${IDs[1]} / 9:${MY_CULTIVAR_NAME})"			
	exit 0			
fi				
				
myid=$1				
				
				
if  [ -e fin${myid}.txt ];then				
    echo "----------------------------------------"				
    echo "rm fin${myid}.txt"				
    echo "----------------------------------------"				
    rm -f fin${myid}.txt				
fi				
				
				
echo "----------------------------------------"				
echo "Run Bat_rename.sh ${myid} ..."				
echo "----------------------------------------"				
./Bat_rename.sh ${myid}				
				
echo "----------------------------------------"				
echo "Run Bat_fastq_quality_filter.sh ${myid} ..."				
echo "----------------------------------------"				
./Bat_fastq_quality_filter.sh ${myid}				
				
echo "----------------------------------------"				
echo "Run Bat_sep_pair.pl.sh ${myid} ..."				
echo "----------------------------------------"				
./Bat_sep_pair.pl.sh ${myid}				
				
# --------------------------------------------------				
if [ ${myid} -eq 0 -o ${myid} -eq 1 ]; then				
	touch fin${myid}.txt			
fi				
				
# --------------------------------------------------				
equalize_ready=0				
if [ ${myid} -eq 0 -o ${myid} -eq 1 ]; then				
	if [ ${myid} -eq 0 ]; then			
		if [ -e fin1.txt ]; then		
			equalize_ready=1	
		fi		
	else			
		if [ -e fin0.txt ]; then		
			equalize_ready=1	
		fi		
	fi			
				
fi				
				
# --------------------------------------------------				
if [ ${equalize_ready} -eq 1 ]; then				
	echo "----------------------------------------"			
	echo "Run Bat_equalize_read.sh ..."			
	echo "----------------------------------------"			
	./Bat_equalize_read.sh			
fi
