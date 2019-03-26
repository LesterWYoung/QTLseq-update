#! /bin/sh			
			
# --------------------------------------------------			
# load common functions			
. ../0.common/common.fnc			
# --------------------------------------------------			
echo "\n--------------------\nSetting reference sequence"
echo "---------------------\n"

printf "Was a secondary reference made? "; Set_MAKE_SECONDARY_REF
if [[ ${MAKE_SECONDARY_REF} == "yes" ]]; then
	printf "A secondary reference was made using readfiles from "; Set_MY_CULTIVAR_NAME
	ref_to_use="../2.make_consensus/secondary_reference/${MY_CULTIVAR_NAME}_secondaryref.fa"
	echo "Using ${ref_to_use} as reference sequence"
else
	printf "Using reference sequence "; Set_PUBLIC_REF_FASTA
	ref_to_use="${REF_FASTA}"
fi

#--------------------
# Determine if bowtie2 library is present. Build library if not present
#--------------------
echo "\n----------\nLooking for bowtie2 library for reference sequence ${ref_to_use}"
ref_seq_dir=`dirname ${ref_to_use}`
ref_seq_filename=`basename ${ref_to_use}`
ref_seq_basename=`cut -f 1 -d "." <<< $ref_seq_filename`
bowtie_library_name="${ref_seq_dir}/${ref_seq_basename}.1.bt2"

if [ ! -f "${bowtie_library_name}" ]; then           #if ${ref_seq_basename}.1.bt2 is not present then build library
	echo "Building bowtie2 library using '${ref_seq_basename}' as library name"
       	CMD="bowtie2-build ${ref_to_use} ${ref_seq_dir}/${ref_seq_basename}"
       		echo ${CMD}
       		eval ${CMD}

	else
       		echo "Bowtie2 library for ${ref_seq_filename} is present in ${ref_seq_dir}"
fi

printf "Directory to get trimmed readfiles: "; Set_BULK_NAME

for bulk_AB in A B; do
        CMD="Set_BULK_NAME_ID${bulk_AB}"
        printf "Aligning trimmed readfiles from bulk: "
        eval ${CMD}                                     #gives $BULK_NAME_ID value of BULK_NAME_IDA or $BULK_NAME_IDB 
	
	output_directory="${BULK_NAME}_${bulk_AB}_alignments"
	CMD="\mkdir ${output_directory}"
	eval ${CMD}
        
	CMD="find ../1.qualify_read/${BULK_NAME}_${bulk_AB}/${BULK_NAME_ID}_[0-9]*-paired_R[12].fastq.gz > readfilelist_AB"
        eval ${CMD}                                     # gets list of mybulk_A or mybulk_B readfiles filtered by trimmomatic

        number_of_R1=`grep -c "R1" readfilelist_AB`             #get number of forward readfiles
        number_of_R2=`grep -c "R2" readfilelist_AB`             #get number of reverse readfiles
        if [ $number_of_R1 -ne $number_of_R2 ]; then            #exit if number of R1 and R2 readfiles is not the same
                echo "There are an uneven number of forward and reverse readfiles in ../${BULK_NAME}_${BULK_NAME_ID}"
                echo "Both forward (R1) and reverse (R2) readfiles for each individual should be in the same directory"
                cat readfilelist_AB
                exit 3
        fi

        echo "Using bowtie2 to align the following readfiles against ${ref_to_use}:"; cat readfilelist_AB

        while read -r readfile; do
                readfilename=`basename $readfile`               #gets file name, removing path
                                                                #echo "filename being used = $readfilename"
                shortreadfilename=`cut -f 1,2 -d "_" <<< $readfilename | cut -f 1 -d "-"`
                                                                #gets name_[0-9]* from readfilename
                                                                #echo "short name = $shortreadfilename"
                read_direction=`cut -f 3 -d "_" <<< $readfilename | cut -f 1 -d "."`
                                                                #searches for "-R[12]." to get read direction
                                                                #used cut -f 3 as bulk reads have additional number in name
                                                                #echo "read dir = $read_direction"

                if [[ ${read_direction} == "R1" ]]; then
                        forward_readfile=${readfile}
                        reverse_readfile=`echo $readfile | sed 's/R1/R2/'`
			                                	#echo ${read_direction} ${forward_readfile}
                else
	                if [[ ${read_direction} == "R2" ]]; then
				continue			#skip this readfile as it was just done
                	fi
                fi

		outfilename="${output_directory}/${shortreadfilename}.bam"

		echo "Aligning ${forward_readfile} and ${reverse_readfile}"
		printf "Using these bowtie2 options: "; Set_BOWTIE2_OPTIONS
		printf "Number of threads to use: "; Set_BOWTIE2_CPU

		CMD="bowtie2 ${BOWTIE2_OPTIONS}"                        # change options in ../config.txt.
		CMD="$CMD -p ${BOWTIE2_CPU}"                            # runs bowtie2 multithreaded
		CMD="$CMD --rg-id ${MY_BULK_ID}_${shortreadfilename}"     # add RG tag to results
		CMD="$CMD -x ${ref_seq_dir}/${ref_seq_basename}"        # library name for reference sequence 
		CMD="$CMD -1 $forward_readfile -2 $reverse_readfile"    # readfile locations and names
          		                                              # CMD="$CMD -S $outfilename.sam"
		CMD="$CMD | samtools view -@ 4 -Shb - | samtools sort -@ 4 -o ${outfilename}"

		echo ${CMD}
		#eval ${CMD}
		
		echo ${outfilename} >> bulk_alignment_filenames

	continue				#skips next file (presumably the R2 file) in the while loop

	done < readfilelist_AB						# go onto next readfile		
done									# go onto mybulk_B


# --------------------
# Call SNPs and determine genotypes
# --------------------
echo "\n--------------------\nSNP and INDEL calling using bcftools" 
printf "bam files used:\n"; cat bulk_alignment_filenames
printf "bcftools mpileup options used: "; Set_BCFT_MPILEUP_BULK
printf "bcftools call options used: "; Set_BCFT_CALL_BULK
printf "bcftools normalize options used: "; Set_BCFT_NORM_BULK
printf "bcftools filter options used: "; Set_BCFT_FILTER_BULK

bulk_vcffile="${output_directory}/${shortreadfilename}.vcf.gz"

CMD="${BCFT_MPILEUP} | ${BCFT_NORM} | ${BCFT_CALL} | ${BCFT_FILTER}"
echo ${CMD}
#eval ${CMD}

CMD="bcftools index -f ${bulk_vcffile}"
echo ${CMD}
#eval ${CMD}						#echo ${outfilename}
			




exit
			
# ==================================================			
./Bat_exclude_common_snps.pl.sh ${myid}			
			
# ==================================================			
./Bat_awk_custom.sh ${myid}			

