#! /bin/sh

# --------------------------------------------------
# load common functions
. ../0.common/common.fnc
# --------------------------------------------------

echo "\n----------\nGenerating a secondary reference, if this is being made\n----------"

printf "Is secondary reference being made? "; Set_MAKE_SECONDARY_REF
if [[ ${MAKE_SECONDARY_REF} == "yes" ]]; then
	echo "Generating secondary reference using trimmed readfiles in 1.qualify_read/secondary_readfiles"
	printf "Name of parental cultivar used for secondary reference: "
	Set_MY_CULTIVAR_NAME		# MY_CUTIVAR_NAME=${Key1_My_cultivar_sample_name} in config.txt 
	CMD="find ../1.qualify_read/secondary_readfiles/${MY_CULTIVAR_NAME}-paired_R[12].fastq.gz > readfilelist"
	eval ${CMD}	
	echo "Parental readfiles being used to generate secondary reference"; cat readfilelist

	numfiles=`wc -l < readfilelist`
	if [ $numfiles -ne 2 ]; then      	# Return an error if there aren't two readfiles present
		echo "The trimmed readfiles required to generate a secondary reference are are missing or not in the correct format."
		echo "Generate them using 1.qualify_read/Bat_secondary_filter_trim.sh or put them in 1.qualify_read/secondary_readfiles"
		exit 20
	fi

	

fi



