# ##################################################			
# global environment			
# ##################################################			
Set_TOPPATH_SCRIPTS(){			
	TOPPATH_SCRIPTS="/Users/LWY/Unix/QTLseq-update-local/QTLseq-update/original_pipeline/ibrc_scripts"		
	echo "${TOPPATH_SCRIPTS}"		
}			
			
Set_TOPPATH_COVAL(){			
	TOPPATH_COVAL=""		
	echo "${TOPPATH_COVAL}"		
}			

Set_MAKE_SECONDARY_REF(){
	MAKE_SECONDARY_REF="yes"
	echo "${MAKE_SECONDARY_REF}"
}
			
# ##################################################			
# for 1.qualify_read
# now using trimmomatic instead 13 mar 2019			
# ##################################################			
Set_TRIMMOMATIC(){			
	TRIMMOMATIC="/Users/LWY/Unix/QTLseq-update-local/QTLseq-update/original_pipeline/ibrc_scripts/1./trimmomatic-0.38.jar"		
	echo "${TRIMMOMATIC}"		
}
# ==================================================			
# Bat_fastq_quality_filter.sh			
# ==================================================			
Set_READ_QVAL(){			
	READ_QVAL=30		
	echo "${READ_QVAL}"		
}			

# Function to setting . Trimmomatic will remove reads shorter that this number (and it's mate): LY
Set_MIN_LEN_BULKED_READS(){
	MIN_LEN_BULKED_READS=90
	echo "${MIN_LEN_BULKED_READS}"
}			

Set_READ_QVAL_MY_CULTIVAR(){			
	READ_QVAL=30		
	echo "${READ_QVAL}"		
}

# Function to setting . Trimmomatic will remove reads shorter that this number (and it's mate)
Set_MIN_LEN_MY_CULTIVAR_READS(){
	MIN_LEN_MY_CULTIVAR_READS=90
	echo "${MIN_LEN_MY_CULTIVAR_READS}"
}			
			
# ##################################################			
# for 2.make_consensus 			
# ##################################################			
# ==================================================			
# If secondary ref being made, define these functions			
# ==================================================			
Set_PATH_TO_TRIMMED_PARENTAL_READS(){			
	PATH_TO_TRIMMED_PARENTAL_READS="/Users/LWY/Unix/QTLseq-update-local/QTLseq-update/original_pipeline/1.qualify_read/secondary_readfiles"
	echo "${PATH_TO_TRIMMED_PARENTAL_READS}"		
}			

Set_PATH_TO_SECONDARY_REF(){			
	PATH_TO_SECONDARY_REF="/Users/LWY/Unix/QTLseq-update-local/QTLseq-update/original_pipeline/2.make_consensus/secondary_reference"		
	echo "${PATH_TO_SECONDARY_REF}"		
}			

Set_MY_CULTIVAR_NAME(){			
	MY_CULTIVAR_NAME="Adelie"		
	echo "${MY_CULTIVAR_NAME}"		
}			
			
# ==================================================			
# Bat_bowtie2_secondaryref.sh			
# ==================================================			
Set_PUBLIC_REF_FASTA(){			
	REF_FASTA="/Users/LWY/Unix/QTLseq-update-local/QTLseq-update/original_pipeline/downloaded_fasta/LG5frag.fa"		
	echo "${REF_FASTA}"		
}			
			


# ##################################################			
# for 3.alignment			
# ##################################################			
Set_BULK_NAME(){			
	BULK_NAME="mybulk"		
	echo "${BULK_NAME}"		
}			
			
Set_BULK_NAME_IDA(){			
	BULK_NAME_ID="A"		
	echo "${BULK_NAME_ID}"		
}			
			
Set_BULK_NAME_IDB(){			
	BULK_NAME_ID="B"		
	echo "${BULK_NAME_ID}"		
}			
			
Set_REF_FASTA(){			
	REF_FASTA=""		
	echo "${REF_FASTA}"		
}			
			
# --------------------------------------------------			
# for Bowtie2 alignment		
# --------------------------------------------------			
Set_BOWTIE2_CPU(){			
	BOWTIE2_CPU=4		
	echo "${BOWTIE2_CPU}"	
}			

Set_BOWTIE2_OPTIONS(){
	BOWTIE2_OPTIONS="--no-discordant"
	BOWTIE2_OPTIONS=" --no-unal"
	BOWTIE2_OPTIONS=" --no-mixed"
	BOWTIE2_OPTIONS=" --sensitive"
	BOWTIE2_OPTIONS=" --sensitive"
	echo "${BOWTIE2_OPTIONS}"
}
			
