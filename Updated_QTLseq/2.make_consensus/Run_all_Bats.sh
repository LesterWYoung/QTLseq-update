#! /bin/sh

# --------------------------------------------------
# load common functions
. ../0.common/common.fnc
# --------------------------------------------------

echo "\n----------\nGenerating a secondary reference, if this is being made\n----------"

printf "Is secondary reference being made? "; Set_MAKE_SECONDARY_REF
if [[ ${MAKE_SECONDARY_REF} == "yes" ]]; then
	. ./Bat_make_secondary_reference.sh
#change Set_REF_FASTA so that ref = ${MY_CULTIVAR_NAME}_secondaryref.fa
else
	echo "Proceed to bulk alignment against "; Set_REF_FASTA
fi



