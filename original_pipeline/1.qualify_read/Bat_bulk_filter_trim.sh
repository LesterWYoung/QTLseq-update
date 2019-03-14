#! /bin/sh

#This script uses trimmomatic to produce filtered and trimmed files
#This version is for trimming the bulk readfiles from mybulk_A and mybulk_B
#The reads are in ../mybulk_[AB]/name_[0-9]*_R[12].fastq.gz

source ../0.common/common.fnc	#load common variables including BULK_NAME_IDA, BULK_NAME_IDB
				# and $Key1_Phred_quality_score_for_bulked (default=30)

#----------------------------------------------------
# Load variable names by calling Set_xxxxxx function.
# values for these variables can be user-altered in config.txt and set by running Bat_common_function.sh
#----------------------------------------------------
echo "\n----------\nTrimming and filtering bulk readfiles\n----------"
printf "Path to trimmomatic: "; Set_TRIMMOMATIC	
					# TRIMMOMATIC=${Key1_Path_to_trimmomatic}. Probably .ibrc_scripts/1./trimmomatic
printf "Min Phred value for bulked reads: "; Set_READ_QVAL	
					# READ_QVAL=${Key1_Phred_quality_score_for_bulked} in config.txt
printf "Min length of reads after trimming: "; Set_MIN_LEN_BULKED_READS	
					# MIN_LEN_BULKED_READS=90 set in config.txt and used by trimmomatic
printf "Diretory prefix: "; Set_BULK_NAME	
					# BULK_NAME="${Key1_Bulked_sample_name}" This is the directory the readfiles are located

#----------------------------------------------------
# Trim and filter mybulk_A and mybulk_B readfiles
#----------------------------------------------------
for bulk_AB in A B; do

	CMD="Set_BULK_NAME_ID${bulk_AB}"
	printf "Procesing readfiles from bulk: "
	eval ${CMD}					#gives $BULK_NAME_ID value of BULK_NAME_IDA or $BULK_NAME_IDB 

	CMD="find ../${BULK_NAME}_${bulk_AB}/${BULK_NAME_ID}_[0-9]*_R[12].fastq.gz > readfilelist_AB"
	eval ${CMD}					#assumes unchanged ../mybulk_A and ../mybulk_B directory names

	numfiles=`wc -l < readfilelist_AB`
	if [ $numfiles -eq 0 ]; then      	# if no readfiles are present then give error message and exit
		echo "Readfiles for Bulk ${BULK_NAME_ID} are missing."
		echo "Put them in ../${BULK_NAME}_${BULK_NAME_ID} using name_[0-9*]_R[12].fastq.gz format"
		exit 1
	fi

	number_of_R1=`grep -c "R1" readfilelist_AB`		#get number of forward readfiles
	number_of_R2=`grep -c "R2" readfilelist_AB`		#get number of reverse readfiles
	if [ $number_of_R1 -ne $number_of_R2 ]; then		#exit if number of R1 and R2 readfiles is not the same
		echo "There are an uneven number of forward and reverse readfiles in ../${BULK_NAME}_${BULK_NAME_ID}"
		echo "Both forward (R1) and reverse (R2) readfiles for each individual should be in the same directory"
		cat readfilelist_AB
		exit 3
	fi

	echo "Running trimmomatic on the following ../${BULK_NAME}_${BULK_NAME_ID} readfiles:"; cat readfilelist_AB

	while read -r readfile; do
		readfilename=`basename $readfile`		#gets file name, removing path
								# echo "filename being used = $readfilename"
		shortreadfilename=`cut -f 1,2 -d "_" <<< $readfilename`
								#gets name_[0-9]* from readfilename
								# echo "short name = $shortreadfilename"
	        read_direction=`cut -f 3 -d "_" <<< $readfilename | cut -f 1 -d "."`
        							#searches for "-R[12]." to get read direction
								#used cut -f 3 as bulk reads have additional number in name
								# echo "read dir = $read_direction"

		if [[ $read_direction == "R1" ]]; then		#if foward readfile, then assign both forward and reverse variable names
        							#echo "found R1"
			readfileR1=$readfile
			outfileR1paired="${BULK_NAME}_${BULK_NAME_ID}/${shortreadfilename}-paired_R1.fastq.gz"
			outfileR1unpaired="${BULK_NAME}_${BULK_NAME_ID}/${shortreadfilename}-unpaired_R1.fastq.gz"
								# echo "first readfile = $readfileR1  $outfileR1paired   $outfileR1unpaired"
        							#echo "found R2"
			readfileR2=`echo $readfile | sed 's/R1/R2/'`	#changes R1 in $readfile to R2
			outfileR2paired="${BULK_NAME}_${BULK_NAME_ID}/${shortreadfilename}-paired_R2.fastq.gz"
			outfileR2unpaired="${BULK_NAME}_${BULK_NAME_ID}/${shortreadfilename}-unpaired_R2.fastq.gz"

			if [ ! -f ${readfileR2} ]; then		# checks to see if reverse file exists. If not, then exit
				echo "Warning! Expected to find ${readfileR2} to match ${readfileR1}, but it doesn't exist!"
				echo "Both R1 and R2 readfiles for an individual should have the same name and number"
				exit 4
			fi
								# echo "second readfile = $readfileR2  $outfileR2paired   $outfileR2unpaired"
		else						# assumes forward readfile was read in last iteration, so just continue while loop
			continue
		fi

		CMD="java -jar ${TRIMMOMATIC} PE"
		CMD="$CMD $readfileR1 $readfileR2"			#fwd and rev readfiles to be trimmed
		CMD="$CMD $outfileR1paired $outfileR1unpaired"		#fwd output reads (paired with rev reads and unpaired)
		CMD="$CMD $outfileR2paired $outfileR2unpaired"		#rev output reads (paired with fwd reads and unpaired)
		CMD="$CMD LEADING:4 TRAILING:4"			#trims bases from each end if qual<15. Change if desired
		CMD="$CMD SLIDINGWINDOW:4:${READ_QVAL}"			# cut if average quality score for sliding window < set value
		CMD="$CMD MINLEN:${MIN_LEN_BULKED_READS}"		#filter out reads < set value after trimming
		CMD="$CMD AVGQUAL:${READ_QVAL}"				#remove read if average quality is less than set value
		
#		echo ${CMD}
#		eval ${CMD}
	
	done < readfilelist_AB					# end of while loop. Goes through all readfiles in bulk
done								# end of for loop. Does bulk A dnthe bulk B

echo "Finished trimming and filtering bulk readfiles"
