#!/usr/bin/env bash
set -x
BINDIR=$(dirname "$0")
INPUT=$1
OUTPUT=$2
TOTAL=$3
check_errs() 
{ 
  # Function. 
  if [ "${1}" != "0" ]; then 
    exit ${2} 
  fi 
} 

split() {
	INPUT1=$1
	OUTPUT1=$2
	INDEX1=$3
	TOTAL1=$4
	TABLENAME1=$5
	SET_NAME=$6
	awk 'NR==1 || NR%'"${TOTAL1}"' == '"${INDEX1}" ${INPUT1}_${TABLENAME1}.tsv > ${OUTPUT1}_test_${SET_NAME}_${TABLENAME1}.tsv
	awk 'NR==1 || NR%'"${TOTAL1}"' != '"${INDEX1}" ${INPUT1}_${TABLENAME1}.tsv > ${OUTPUT1}_train_${SET_NAME}_${TABLENAME1}.tsv
}

split_all() {
	INPUT2=$1
	OUTPUT2=$2
	TOTAL2=$3
	TABLENAME2=$4
	for ((ID = 0; ID < $TOTAL2; ID++))
	do
		split $INPUT2 $OUTPUT2 $ID $TOTAL2 $TABLENAME2 $ID
	done
}

split_allfiles() {
	INPUT3=$1
	OUTPUT3=$2
	TOTAL3=$3
	split ${INPUT3} ${OUTPUT3} 4 $TOTAL3 FeatureMatrix Split
        split ${INPUT3} ${OUTPUT3} 4 $TOTAL3 Outcome Split
        split ${INPUT3} ${OUTPUT3} 4 $TOTAL3 Phenotype Split

	split_all ${OUTPUT3}_train_Split ${OUTPUT3} $TOTAL3 FeatureMatrix
        split_all ${OUTPUT3}_train_Split ${OUTPUT3} $TOTAL3 Outcome
        split_all ${OUTPUT3}_train_Split ${OUTPUT3} $TOTAL3 Phenotype
}

split_allfiles ${INPUT}/sc1_Phase1_GE ${OUTPUT}/sc1_Phase1_GE $TOTAL
split_allfiles ${INPUT}/sc2_Phase1_CN ${OUTPUT}/sc2_Phase1_CN $TOTAL
split_allfiles ${INPUT}/sc3_Phase1_CN_GE ${OUTPUT}/sc3_Phase1_CN_GE $TOTAL
