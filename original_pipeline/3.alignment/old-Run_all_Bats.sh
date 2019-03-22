#! /bin/sh			
			
# --------------------------------------------------			
# load common functions			
. ../0.common/common.fnc			
# --------------------------------------------------			
IDA=`Set_BULK_NAME_IDA`			
IDB=`Set_BULK_NAME_IDB`			
IDs=(${IDA} ${IDB})			
			
MISs=`Set_MISs`			
			
if [ -z $1 ];then			
	echo "missing args!"		
	echo "[usage]"		
	echo "    $0 <INT>"		
	echo "        <INT> : 0 or 1 (0:${IDs[0]} / 1:${IDs[1]})"		
	exit 0		
fi			
			
myid=$1			
			
			
# ==================================================			
./Bat_make_symbolic_link_of_qualified_fastq.sh ${myid}			
			
# ==================================================			
./Bat_bwa2bam.sh ${myid}			
			
# ==================================================			
for mis in ${MISs}			
do			
	./Bat_run_coval-refine-bam.pl.sh ${myid} ${mis}		
	./Bat_run_coval-call-pileup.pl.sh ${myid} ${mis}		
done			
			
# ==================================================			
./Bat_exclude_common_snps.pl.sh ${myid}			
			
# ==================================================			
./Bat_awk_custom.sh ${myid}			

