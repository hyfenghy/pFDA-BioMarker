#!/usr/bin/env bash
set -x
BINDIR=$(dirname "$0")
BINARY=$1
IN=$2
IN_TEST=$3
OUT=$4
TOTAL=$5
TRAIN=$6
TREEMARGIN=$7
DISTWT=$8
ALLFLAG=$9
ITER=${10}
MARGINTHRES=${11}
REALFEAT_CNT=${12}
MINNEIGH=${13}
MINWT=${14}
MAXNEIGH=${15}
MAXWT=${16}
NEIGHSTEP=${17}
WTSTEP=${18}
DISTSCALE=${19}
RANDOMCNT=${20}
TOPFIXED=${21}
FEAT_CNT=${22}

check_errs() 
{ 
  # Function. 
  if [ "${1}" != "0" ]; then 
    exit ${2} 
  fi 
} 

TREE_MODEL=${OUT}_Tree.model
KNN_MODEL=${OUT}_KNN.model
$BINDIR/Sentieon_BioM_Train.sh $BINARY $TRAIN $TREE_MODEL $KNN_MODEL $TREEMARGIN $DISTWT $ALLFLAG $ITER ${OUT}.cnt

OUT_CV=${OUT}_CV
LOG=${OUT}_CV_PERF
rm $LOG

MAXCNT=`cat ${OUT}.cnt`
for CNT in $(seq 0.2 0.1 $MAXCNT )
do
    rm ${OUT_CV}*perf
    TRAINTEST_ARG="$TREEMARGIN	$DISTWT	$ALLFLAG	$ITER	$MARGINTHRES	$CNT	$REALFEAT_CNT	$MINNEIGH	$MINWT	$MAXNEIGH	$MAXWT	$NEIGHSTEP	$WTSTEP	$DISTSCALE	$RANDOMCNT	$TOPFIXED	$FEAT_CNT"
    for ((ID = 0; ID < $TOTAL; ID++))
    do
	${BINDIR}/Sentieon_BioM_TrainTestEval.sh "$BINARY" ${IN}_${ID} ${IN_TEST}_${ID} $OUT_CV ${OUT_CV}_${ID}_perf $TRAINTEST_ARG
    done
	
    echo -n -e "CV" '\t' "$TRAINTEST_ARG" '\t' >> $LOG
    $BINARY 7 $OUT_CV $TOTAL >> $LOG

    OUT_SPLIT=${OUT}_SPLIT
    echo -n -e "SPLIT" "$TRAINTEST_ARG" '\t' >> $LOG
    ${BINDIR}/Sentieon_BioM_TrainTestEval.sh "$BINARY" ${IN}_Split ${IN_TEST}_Split $OUT_SPLIT $LOG $TRAINTEST_ARG

    $BINDIR/Sentieon_BioM_Test.sh $BINARY ${IN}_Split ${IN}_Split ${OUT_SPLIT}_Train ${OUT_SPLIT}_Tree.model ${OUT_SPLIT}_KNN.model $MARGINTHRES $CNT $ALLFLAG $REALFEAT_CNT $MINNEIGH $MINWT $MAXNEIGH $MAXWT $NEIGHSTEP $WTSTEP $DISTSCALE $RANDOMCNT $TOPFIXED $FEAT_CNT
    echo -n -e "SPLIT_TRAIN" "$TRAINTEST_ARG" '\t' >> $LOG
    $BINARY 5 ${IN}_Split_Outcome.tsv ${OUT_SPLIT}_Train_Final.prob $CNT >> $LOG
done

