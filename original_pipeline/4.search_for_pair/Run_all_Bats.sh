#! /bin/sh				
				
# --------------------------------------------------				
# load common functions				
. ../0.common/common.fnc				
# --------------------------------------------------				
MISs=`Set_MISs`				
				
echo "----------------------------------------"				
echo "Run Bat_paired_or_unpaired.sh ..."				
echo "----------------------------------------"				
./Bat_paired_or_unpaired.sh				
				
# ==================================================				
for mis in ${MISs}				
do				
	for myid in 0 1			
	do			
		echo "----------------------------------------"		
		echo "Run Bat_search_for_pair_of_unpaired.sh ${myid} ${mis} ..."		
		echo "----------------------------------------"		
		./Bat_search_for_pair_of_unpaired.sh ${myid} ${mis}		
				
		echo "----------------------------------------"		
		echo "Run Bat_pileup_to_pileup_with_SNPindex.sh ${myid} ${mis} ..."		
		echo "----------------------------------------"		
		./Bat_pileup_to_pileup_with_SNPindex.sh ${myid} ${mis}		
	done			
				
	echo "----------------------------------------"			
	echo "Run Bat_merge_paired.sh ${mis} ..."			
	echo "----------------------------------------"			
	./Bat_merge_paired.sh ${mis}			
done
