#!/bin/bash

#load, from 0.common/common.fnc, the location abs of the parental reads, the mybulk_A and mybulk_B, and the reference sequence directories and names
#also loads values for options in bowtie2 alignment and bcftools 
#if generating a secondary reference, does things different from aligning pools

 
#Get list of readfiles to be aligned. These are in name_R[12].fastq.gz format for parents and name_[0-9]*_R[12].fastq.gz format for the individuals in the pools

#for this testscript start assuming secondary reference is being aligned here. Put parental readfiles in ../readfiles
 
find ../readfiles/* > readfilelist			#use find instead of list as directory path is included in output

while read -r readfile
do
	readfilename=`cut -f 3 -d "/" <<< $readfile`	#gets file name, removing path
							#echo "filename being used = $readfilename"
							#can maybe change this to $Key1_My_cultivar_sample_name
       	if [[ $readfilename == "readme.txt" ]]
       	then
       		echo "readme.txt found"
        	continue				#skips ahead while loop if $readfilename is readme.txt
        fi

        read_direction=`cut -f 2 -d "_" <<< $readfilename | cut -f 1 -d "."`	
        						#searches for "-R[12]." to get read direction
							#echo "read dir = $read_direction"
        if [[ $read_direction == "R2" ]]
        then 
        	echo "found R2"
        	continue				#skips ahead while loop if $readfilename = 2 as previous iteration aligned readfiles
        fi

        forward_readfile=${readfile}
        reverse_readfile=`echo ${readfile} | sed 's/_R1./_R2./'`
	outfilename=`cut -f 1 -d "_" <<< $readfilename`

	echo "fwd = $forward_readfile"
	echo "rev = $reverse_readfile"
	echo "outfilename = $outfilename"
	
done < readfilelist

#now align reads against reference sequence in ../reference_sequence
#first, determine if bowtie2 library is present in ../reference_sequence. Build library if not present

ref_seq_name="Beth_13_4-6.fa"
if ! [[ `find ../reference_sequence/ref_seq.1.bt2` ]]		#if ref_seq.1.bt2 is not present then build library
then
	echo "Building bowtie2 library using $ref_seq_name" 
	bowtie2-build "../reference_sequence/${ref_seq_name}" ref_seq
fi

#build bowtie2 alignment command
CMD="bowtie2 --no-unal  --no-discordant"		#only paired concordant reads aligning to ref are kept. add --no-mixed if you want alignment to be more stringent
CMD="$CMD -p 4"						#runs bowtie2 multithreaded *** replace 4 with ${NUM_THREADS} in pipeline ***
CMD="$CMD --rg-id $outfilename"				#add RG tag to results
CMD="$CMD -x ../reference_sequence/ref_seq"		#library name for reference sequence 
CMD="$CMD -1 $forward_readfile -2 $reverse_readfile"	#readfile locations and names
CMD="$CMD -S $outfilename.sam"

echo $CMD
eval $CMD

#bowtie2 -p 4 --no-unal --no-mixed --no-discordant --rg-id RajaLG5 -x LG5-1.150..1.375Mbp \
#	-q -1 RxBF2_Raja_1.fq.gz -2 RxBF2_Raja_2.fq.gz | \
#	samtools view -@ 4 -Shb - | samtools sort -@4 -o $outname -

#outname="BisonLG5-1.15Mbp.sort.bam"

#bowtie2 -p 4 --no-unal --no-mixed --no-discordant --rg-id BisonLG5 -x LG5-1.150..1.375Mbp \
#	-q -1 RxBF2_Bison_1.fq.gz -2 RxBF2_Bison_2.fq.gz | \
#	samtools view -@ 4 -Shb - | samtools sort -@4 -o $outname -

#note: use trimmomatic filtered reads (performed on raw fastq reads from PBI)
# RxBF2_Bison_1_fq.gz is the 2P output from trimmomatic
