#!/bin/bash

# --------------------------------------------------
# load common functions
. ../0.common/common.fnc
# --------------------------------------------------

echo "----------\nGenerating secondary reference using trimmed readfiles in 1.qualify_read/secondary_readfiles"
	printf "Path to trimmed parental readfiles: "; Set_PATH_TO_TRIMMED_PARENTAL_READS
		short_trimmed_reads_dir="../`echo ${PATH_TO_TRIMMED_PARENTAL_READS} | rev | cut -f 1,2 -d "/" | rev`"
		echo "Shortened to ${short_trimmed_reads_dir}"
	printf "Directory to save secondary reference files:"; Set_PATH_TO_SECONDARY_REF
		short_sec_ref_dir=`basename ${PATH_TO_SECONDARY_REF}`; echo "Shortened to ${short_sec_ref_dir}"
	printf "Name of parental cultivar used for secondary reference: "; Set_MY_CULTIVAR_NAME
 
	CMD="find ${short_trimmed_reads_dir}/${MY_CULTIVAR_NAME}-paired_R[12].fastq.gz > readfilelist"
	eval ${CMD}	
	echo "Parental readfiles being used to generate secondary reference"; cat readfilelist

	numfiles=`wc -l < readfilelist`

	if [ $numfiles -ne 2 ]; then      		# Return an error if there aren't two readfiles present
		echo "The trimmed readfiles required to generate a secondary reference are are missing or not in the correct format."
		echo "Generate them using 1.qualify_read/Bat_secondary_filter_trim.sh or put them in ${short_trimmed_reads_dir}"
		exit 20
	fi

	while read -r readfile; do
        	readfilename=`basename $readfile`               #gets file name, removing path
                                                        	 echo "filename being used = $readfilename"

        	read_direction=`cut -f 2 -d "_" <<< $readfilename | cut -f 1 -d "."`
                                                        	#searches for "-R[12]." to get read direction
                                                        	 echo "read dir = $read_direction"

		if [[ ${read_direction} == "R1" ]]; then
		        forward_readfile=${readfile}
			else 
        		reverse_readfile=${readfile}
`		fi
		outfilename="${short_sec_ref_dir}/${MY_CULTIVAR_NAME}_secondaryref"

			echo "fwd = $forward_readfile"
			echo "rev = $reverse_readfile"
			echo "outfilename = $outfilename"
	done < readfilelist

exit

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
