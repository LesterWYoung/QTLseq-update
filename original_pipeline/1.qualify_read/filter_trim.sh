#!/bin/sh

#This script uses trimmomatic to produce filtered and trimmed files
#This version is for trimming the secondary reference (assuming it's being generated by this pipeline)
#The reads are in ../parental_readfiles/name_R[12].fastq.gz

#. ../0.common/common.fnc	#load common variables, including Key1_My_cultivar_sample_name
				# and $Key1_Phred_quality_score_for_my_cultivar (default=30)

# if parental readfiles are missing and only README present
if [ `find ../parental_readfiles/* | wc -l` -eq 1 ]
	then echo "Parental readfiles are missing. Put them in ../parental_readfiles using name_R[12].fastq.gz format"
	exit 1
fi

find ../parental_readfiles/*.fastq.gz > readfilelist
if [ wc -l readfilelist -gt 2 ]; then	#if more than two readfiles are present, then exit
	echo "There are more than two readfiles in ../parental_readfiles/"
	exit 2
fi

if [ wc -l readfilelist -eq 1 ]; then	#if there is only one readfile present, also exit
	echo "Two readfiles are required (one forward and one reverse) are reqiuired in ../parental_readfiles/"
	exit 3
fi

while read -r readfile
	do
	readfilename=`cut -f 3 -d "/" <<< $readfile`	#gets file name, removing path
	echo "filename being used = $readfilename"
							#can maybe change this to $Key1_My_cultivar_sample_name

        read_direction=`cut -f 2 -d "_" <<< $readfilename | cut -f 1 -d "."`
        						#searches for "-R[12]." to get read direction
	echo "read dir = $read_direction"

	if [[ $read_direction == "R1" ]]; then		#if foward readfile, then assign forward variable names
        	echo "found R1"
		readfileR1=$readfile
		outfileR1paired="secondary_readfiles/${readfilename}-paired_R1.fastq.gz"
		outfileR1unpaired="secondary_readfiles/${readfilename}-unpaired_R1.fastq.gz"
        else

        if [[ $read_direction == "R2" ]]
        then
        	echo "found R2"
		readfileR2=$readfile
		outfileR2paired="secondary_readfiles/${readfilename}-paired_R2.fastq.gz"
		outfileR2unpaired="secondary_readfiles${readfilename}-unpaired_R2.fastq.gz"
echo "$readfileR1  $outfileR1paired   $outfileR1unpaired"
echo "$readfileR2  $outfileR2paired   $outfileR2unpaired"



done < readfilelist
exit 4
CMD="java -jar trimmomatic-0.35.jar PE"
CMD="$CMD $readfile1 $readfile2 $outfile1paired"
CMD="$CMD $outfile1unpaired $outfile2paired $outfile2unpaired"
CMD="$CMD LEADING:15 TRAILING:15"			#trims bases from each end if qual<15. Change if desired
CMD="$CMD SLIDINGWINDOW:4:${Key1_Phred_quality_score_for_my_cultivar}"		#cutoff if quality score for sliding window
CMD="$CMD MINLEN:110"								#filter out reads<110 bp after trimming.

#echo ${CMD}
#eval ${CMD}

