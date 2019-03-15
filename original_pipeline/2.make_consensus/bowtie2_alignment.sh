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
	outfilename="${short_sec_ref_dir}/${MY_CULTIVAR_NAME}_secondaryref"

#now align reads against reference sequence in Set_PUBLIC_REF_FASTA
#first, determine if bowtie2 library is present in Set_PUBLIC_REF_FASTA. Build library if not present

printf "Location of reference sequence: "; Set_PUBLIC_REF_FASTA
exit
ref_seq_name="Beth_13_4-6.fa"
if ! [[ `find ../reference_sequence/ref_seq.1.bt2` ]]           #if ref_seq.1.bt2 is not present then build library
then
        echo "Building bowtie2 library using $ref_seq_name"
        bowtie2-build "../reference_sequence/${ref_seq_name}" ref_seq
fi
