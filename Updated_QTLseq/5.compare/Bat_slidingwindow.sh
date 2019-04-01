#! /bin/sh					
					
# ==================================================					
# Those keys enable users to customize this script					
# ==================================================					
REF=""					
BULK_NAME=""					
IDs=""					
MISs=""					
DEPTHs=""				
QVAL=""					
PVAL=""					
					
#Average=(0 1 2 3)					
Average=(2 3)					
AverageOption=("2000000 200000" "4000000 400000" "2000000 50000" "4000000 50000")					
AverageDescrpt=("2M200K" "4M400K" "2M50K" "4M50K")					
					
MinCount=10					
					
# --------------------------------------------------					
# load common functions					
. ../0.common/common.fnc					
# --------------------------------------------------					
if [ -z ${REF} ]; then					
	REF=`Set_REF_FASTA`				
fi					
					
if [ -z ${BULK_NAME} ]; then					
	BULK_NAME=`Set_BULK_NAME`				
fi					
					
if [ -z ${IDs} ]; then					
	IDA=`Set_BULK_NAME_IDA`				
	IDB=`Set_BULK_NAME_IDB`				
	IDs=(${IDA} ${IDB})				
fi					

if [ -z ${QVAL} ]; then					
	QVAL=`Set_READ_QVAL`				
fi					
					
if [ -z ${PVAL} ]; then					
	PVAL=`Set_READ_PVAL`				
fi					

if [ -z ${MISs} ]; then					
	MISs=`Set_CONFINTRVL_MISs`				
fi					
					
if [ -z ${DEPTHs} ]; then					
	DEPTHs=`Set_AWK_CUSTOM2_DEPTHs`				
fi					
					
# ==================================================					
# environment					
# ==================================================					
TOPPATH_SCRIPTS=""					
TOPPATH_COVAL=""					
					
if [ -z ${TOPPATH_SCRIPTS} ]; then					
	TOPPATH_SCRIPTS=`Set_TOPPATH_SCRIPTS`				
fi					
					
if [ -z ${TOPPATH_COVAL} ]; then					
	TOPPATH_COVAL=`Set_TOPPATH_COVAL`				
fi					
					
# ==================================================					

NAME="${BULK_NAME}"					
NAME_A="${NAME}_${IDs[0]}_q${QVAL}p${PVAL}"					
NAME_B="${NAME}_${IDs[1]}_q${QVAL}p${PVAL}"					
NAME_AB="${NAME}_${IDs[0]}${IDs[1]}_q${QVAL}p${PVAL}"					
	
FASTA_FAI="${REF}.fai"					
INPILEUP_AB="filtered_merge_${NAME_AB}_paired_pvalue"					
					
WithUnmaskedGraph=1					
## 1 : enabled					
## 0 : disabled					
# ==================================================					
					
# --------------------					
# 					
# --------------------					
MKDIR="mkdir -p"					
SRCPATH="20.awk_custom2"					
DSTPATH="90.slidingwindow"					
					
SCRIPT1="${TOPPATH_SCRIPTS}/5./pileup_to_slidingwindow.R"					
CMD="Rscript ${SCRIPT1}"					
					
SCRIPT2="${TOPPATH_SCRIPTS}/5./make_png.R"					
CONFIGs=("config_delta_SNPindex.ini" \					
		 "config_SNPindex_A.ini" \			
		 "config_SNPindex_B.ini" \			
		 "config_depth.ini" \			
		 "config_density.ini" \			
		 "config_out_of_conf_intrvl.ini")			
CMD2="Rscript ${SCRIPT2}"					
					
for mis in ${MISs}					
do					
	echo "${MKDIR} ${DSTPATH}/mut_index_${mis}/"				
	eval "${MKDIR} ${DSTPATH}/mut_index_${mis}/"				
					
	for dep in ${DEPTHs}				
	do				
					
		PILEUP="${SRCPATH}/mut_index_${mis}/${INPILEUP_AB}_cov${mis}_co${dep}.txt"			
		for ave in ${Average[@]}			
		do			
			OUTPATH="${DSTPATH}/mut_index_${mis}"		
			OUTNAME="${INPILEUP_AB}_sldwnd${AverageDescrpt[$ave]}_cov${mis}_co${dep}.txt"		
			OUTNAME2="mask${MinCount}_${OUTNAME}"		
			echo "${CMD} ${PILEUP} ${FASTA_FAI} ${AverageOption[$ave]} ${MinCount} ${OUTPATH} ${OUTNAME}"		
			eval "${CMD} ${PILEUP} ${FASTA_FAI} ${AverageOption[$ave]} ${MinCount} ${OUTPATH} ${OUTNAME}"		
					
			OUTNAME="${OUTPATH}/${OUTNAME}"		
			OUTNAME2="${OUTPATH}/${OUTNAME2}"		
					
			OUTPNGPATH="${DSTPATH}/pngs/${AverageDescrpt[$ave]}/mut_index_${mis}"		
			echo "mkdir -p ${OUTPNGPATH}"		
			eval "mkdir -p ${OUTPNGPATH}"		
			OUTPNGPATH2="${OUTPNGPATH}/mask${MinCount}"		
			echo "mkdir -p ${OUTPNGPATH2}"		
			eval "mkdir -p ${OUTPNGPATH2}"		
					
			for mycfg in ${CONFIGs[@]}		
			do		
				if [ ${mycfg} = "config_delta_SNPindex.ini" ]; then	
					myheader=${NAME_AB}
				fi	
				if [ ${mycfg} = "config_SNPindex_A.ini" ]; then	
					myheader=${NAME_A}
				fi	
				if [ ${mycfg} = "config_SNPindex_B.ini" ]; then	
					myheader=${NAME_B}
				fi	
				if [ ${mycfg} = "config_depth.ini" ]; then	
					myheader="depth_${NAME_AB}"
				fi	
				if [ ${mycfg} = "config_density.ini" ]; then	
					myheader="density_${NAME_AB}"
				fi	
				if [ ${mycfg} = "config_out_of_conf_intrvl.ini" ]; then	
					myheader="ratio_${NAME_AB}"
				fi	
					
				OUTPNG="${OUTPNGPATH}/${myheader}_sldwnd${AverageDescrpt[$ave]}_cov${mis}_co${dep}.png"	
				OUTPNG2="${OUTPNGPATH2}/mask${MinCount}_${myheader}_sldwnd${AverageDescrpt[$ave]}_cov${mis}_co${dep}.png"	
					
				mycfg="${TOPPATH_SCRIPTS}/5./${mycfg}"	
					
				if [ ${WithUnmaskedGraph} -eq 1 ]; then	
					echo "${CMD2} ${mycfg} ${PILEUP} ${OUTNAME} ${OUTPNG}"
					eval "${CMD2} ${mycfg} ${PILEUP} ${OUTNAME} ${OUTPNG}"
				fi	
				echo "${CMD2} ${mycfg} ${PILEUP} ${OUTNAME2} ${OUTPNG2}"	
				eval "${CMD2} ${mycfg} ${PILEUP} ${OUTNAME2} ${OUTPNG2}"	
			done		
		done			
	done				
					
done
