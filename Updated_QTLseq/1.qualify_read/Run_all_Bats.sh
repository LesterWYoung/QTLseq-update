#! /bin/sh

# --------------------------------------------------
# load common functions
. ../0.common/common.fnc
# --------------------------------------------------

echo "\n----------\nTrimming and filtering reads for bulks and secondary reference (if this is being made)\n----------"

printf "Is secondary reference being made? "; Set_MAKE_SECONDARY_REF
if [[ ${MAKE_SECONDARY_REF} == "yes" ]]; then
	sh Bat_secondary_filter_trim.sh
fi

sh Bat_bulk_filter_trim.sh

	

