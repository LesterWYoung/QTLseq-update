#!/bin/bash

# --------------------------------------------------
# load common functions
. ../0.common/common.fnc
# --------------------------------------------------

echo "\n--------------------"
echo "Generating secondary reference using trimmed readfiles in 1.qualify_read/secondary_readfiles"
echo "--------------------"
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
#--------------------
# Determine if bowtie2 library is present. Build library if not present
#--------------------
echo "\n----------\nLooking for bowtie2 library for reference sequence"
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
	echo "Bowtie2 library for ${ref_seq_filename} is present in ${ref_seq_dir}"
fi

# --------------------
# Perform bowtie2 alignment
# ---------------------
echo "\n--------------------\nPerforming bowtie2 alignment to generate secondary reference"
echo "--------------------"
printf "Using these bowtie2 options: "; Set_BOWTIE2_OPTIONS
printf "Number of threads to use: "; Set_BOWTIE2_CPU

CMD="bowtie2 ${BOWTIE2_OPTIONS}"			# change options in ../config.txt.
CMD="$CMD -p ${BOWTIE2_CPU}"				# runs bowtie2 multithreaded
CMD="$CMD --rg-id ${MY_CULTIVAR_NAME}_secondaryref"	# add RG tag to results
CMD="$CMD -x ${ref_seq_dir}/${ref_seq_basename}"	# library name for reference sequence 
CMD="$CMD -1 $forward_readfile -2 $reverse_readfile"	# readfile locations and names
							# CMD="$CMD -S $outfilename.sam"
CMD="$CMD | samtools view -@ 4 -Shb - | samtools sort -@ 4 -o ${outfilename}"

echo ${CMD}
#eval ${CMD}

# --------------------
# Call SNPs and generate consensus sequence
# --------------------
# This part of the script takes the .bam file generated above calls SNPs 
# and generates a consensus sequence. It puts the consensus sequence into
# ./secondary_reference/ with the name ${MY_CULTIVAR_NAME}_secondaryref.fasta

echo "\n--------------------\nFiltering bowtie2 alignment with bcftools"
echo "--------------------" 
printf "bam file used: ${outfilename}\n"
printf "bcftools mpileup options used: "; Set_BCFT_MPILEUP
printf "bcftools call options used: "; Set_BCFT_CALL
printf "bcftools normalize options used: "; Set_BCFT_NORM
printf "bcftools filter options used: "; Set_BCFT_FILTER

secondary_vcffile="${short_sec_ref_dir}/${MY_CULTIVAR_NAME}_secondaryref.vcf.gz"

CMD="${BCFT_MPILEUP} $outfilename | ${BCFT_NORM} | ${BCFT_CALL} | ${BCFT_FILTER}"
echo ${CMD}
#eval ${CMD}

CMD="bcftools index ${secondary_vcffile}"
#eval ${CMD}

printf "bcftools consensus options used: "; Set_BCFT_CONSENSUS
secondary_ref="${short_sec_ref_dir}/${MY_CULTIVAR_NAME}_secondaryref.fa"
CMD="${BCFT_CONSENSUS}"
echo ${CMD}
#eval ${CMD}
