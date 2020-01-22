#!/usr/bin/env bash
set -x
BINDIR=$(dirname "$0")
BINARY=$1
IN=$2
IN_TEST=$3
OUT=$4
TOTAL=$5
TRAIN=$6

check_errs() 
{ 
  # Function. 
  if [ "${1}" != "0" ]; then 
    exit ${2} 
  fi 
} 


${BINDIR}/Sentieon_BioM_SPC.sh $BINARY $IN $IN_TEST ${OUT}_1  $TOTAL $TRAIN 0.5 2   0 50 0   15 3 2   4  10.0 1 0.1 1 10 10 30
${BINDIR}/Sentieon_BioM_SPC.sh $BINARY $IN $IN_TEST ${OUT}_2  $TOTAL $TRAIN 2.0 4   0 50 20  15 3 2   3  10.0 1 0.1 1 10 5  30
${BINDIR}/Sentieon_BioM_SPC.sh $BINARY $IN $IN_TEST ${OUT}_3  $TOTAL $TRAIN 1.5 2.5 0 50 0   15 3 2   3  10.0 1 0.1 1 10 5  30
${BINDIR}/Sentieon_BioM_SPC.sh $BINARY $IN $IN_TEST ${OUT}_4  $TOTAL $TRAIN 1.0 2   0 50 0   15 3 2   4  10.0 1 0.1 1 10 10 30
${BINDIR}/Sentieon_BioM_SPC.sh $BINARY $IN $IN_TEST ${OUT}_5  $TOTAL $TRAIN 0.5 3   0 50 0   15 4 2   4  10.0 1 0.1 1 10 5  30
${BINDIR}/Sentieon_BioM_SPC.sh $BINARY $IN $IN_TEST ${OUT}_6  $TOTAL $TRAIN 1.0 3   1 50 0   15 4 2   4  10.0 1 0.1 1 10 5  30
${BINDIR}/Sentieon_BioM_SPC.sh $BINARY $IN $IN_TEST ${OUT}_7  $TOTAL $TRAIN 0.5 2   0 50 10  15 4 2   4  10.0 1 0.1 1 10 0  30
${BINDIR}/Sentieon_BioM_SPC.sh $BINARY $IN $IN_TEST ${OUT}_8  $TOTAL $TRAIN 1.0 2   0 50 0   15 3 2   3  11.0 1 0.1 1 10 5  30
${BINDIR}/Sentieon_BioM_SPC.sh $BINARY $IN $IN_TEST ${OUT}_9  $TOTAL $TRAIN 0.5 3   0 50 0   15 4 2   4   9.0 1 0.1 1 10 5  30
${BINDIR}/Sentieon_BioM_SPC.sh $BINARY $IN $IN_TEST ${OUT}_10 $TOTAL $TRAIN 1.0 2   0 50 0   15 3 2   3  10.0 1 0.1 1 10 5  30
${BINDIR}/Sentieon_BioM_SPC.sh $BINARY $IN $IN_TEST ${OUT}_11 $TOTAL $TRAIN 1.0 2   0 50 0   15 3 2   4  10.0 1 0.1 1 10 5  30
${BINDIR}/Sentieon_BioM_SPC.sh $BINARY $IN $IN_TEST ${OUT}_12 $TOTAL $TRAIN 1.0 2   0 50 0   15 4 2   10 10.0 1 0.1 1 10 5  30
${BINDIR}/Sentieon_BioM_SPC.sh $BINARY $IN $IN_TEST ${OUT}_13 $TOTAL $TRAIN 1.0 2   0 50 0   15 4 2   4  10.0 1 0.1 1 10 5  30
${BINDIR}/Sentieon_BioM_SPC.sh $BINARY $IN $IN_TEST ${OUT}_14 $TOTAL $TRAIN 1.0 2   0 50 0   10 4 2   4  10.0 1 0.1 1 10 5  30
${BINDIR}/Sentieon_BioM_SPC.sh $BINARY $IN $IN_TEST ${OUT}_15 $TOTAL $TRAIN 0.5 3   0 50 0   10 4 2   4  10.0 1 0.1 1 10 5  30
${BINDIR}/Sentieon_BioM_SPC.sh $BINARY $IN $IN_TEST ${OUT}_16 $TOTAL $TRAIN 1.0 3   0 50 20  15 3 1.5 4  10.0 1 0.1 1 10 5  30
${BINDIR}/Sentieon_BioM_SPC.sh $BINARY $IN $IN_TEST ${OUT}_17 $TOTAL $TRAIN 0.5 3   0 50 0   10 4 2   4  10.0 1 0.1 1 10 10 30
${BINDIR}/Sentieon_BioM_SPC.sh $BINARY $IN $IN_TEST ${OUT}_18 $TOTAL $TRAIN 2.5 2   0 50 0   10 4 2   4  10.0 1 0.1 1 10 10 30
${BINDIR}/Sentieon_BioM_SPC.sh $BINARY $IN $IN_TEST ${OUT}_19 $TOTAL $TRAIN 2.5 2   0 50 10  10 4 2   4  10.0 1 0.1 1 10 10 30


