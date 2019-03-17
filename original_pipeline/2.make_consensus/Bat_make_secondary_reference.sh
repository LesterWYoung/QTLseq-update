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
		short_sec_ref_dir="./`basename ${PATH_TO_SECONDARY_REF}`"; echo "Shortened to ${short_sec_ref_dir}"
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
                readfilename=`basename $readfile`	#gets file name, removing path
                                                        # echo "filename being used = $readfilename"

                read_direction=`cut -f 2 -d "_" <<< $readfilename | cut -f 1 -d "."`
                                                        #searches for "-R[12]." to get read direction
							#echo "Read direction: $read_direction"
        	if [[ ${read_direction} == "R1" ]]; then
			forward_readfile=${readfile}
							#echo ${read_direction} ${forward_readfile}
		else
		if [[ ${read_direction} == "R2" ]]; then
			reverse_readfile=${readfile}
							#echo ${read_direction} ${reverse_readfile}	
		else
			echo "Cannot parse directions on readfiles in ${short_trimmed_reads}"
			exit
		fi
		fi
	done < readfilelist
	outfilename="${short_sec_ref_dir}/${MY_CULTIVAR_NAME}_secondaryref.bam"

#now align reads against reference sequence in Set_PUBLIC_REF_FASTA

#----------
# Determine if bowtie2 library is present. Build library if not present
#----------
printf "Reference sequence being used: "; Set_PUBLIC_REF_FASTA
ref_seq_dir=`dirname ${REF_FASTA}`
ref_seq_filename=`basename ${REF_FASTA}`
ref_seq_basename=`cut -f 1 -d "." <<< $ref_seq_filename`
bowtie_library_name="${ref_seq_dir}/${ref_seq_basename}.1.bt2"

if [ ! -f "${bowtie_library_name}" ]; then           #if ${ref_seq_basename}.1.bt2 is not present then build library
	echo "Building bowtie2 library using '${ref_seq_basename}' as library name"
	CMD="bowtie2-build ${REF_FASTA} ${ref_seq_dir}/${ref_seq_basename}"
	echo ${CMD}
	eval ${CMD}
	
else
	echo "Bowtie2 library for ${ref_seq_filename} exists"
fi

# ----------
# Perform bowtie2 alignment
# ----------
echo "\n----------\nPerforming bowtie2 alignment to generate secondary reference\n----------"
printf "Bowtie2 options used: "; Set_BOWTIE2_OPTIONS
printf "Number of threads to be used: "; Set_BOWTIE2_CPU

CMD="bowtie2 ${BOWTIE2_OPTIONS}"			# change options in ../config.txt.
CMD="$CMD -p ${BOWTIE2_CPU}"				# runs bowtie2 multithreaded
CMD="$CMD --rg-id ${MY_CULTIVAR_NAME}_secondaryref"	# add RG tag to results
CMD="$CMD -x ${ref_seq_dir}/${ref_seq_basename}"	# library name for reference sequence 
CMD="$CMD -1 $forward_readfile -2 $reverse_readfile"	# readfile locations and names
							# CMD="$CMD -S $outfilename.sam"
CMD="$CMD | samtools view -@ 4 -Shb - | samtools sort -@ 4 -o $outfilename"

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
	
