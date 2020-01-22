#!/usr/bin/env bash
set -x
BINDIR=$(dirname "$0")
BINARY=$1
INPUT=$2
TREE_MODEL=$3
KNN_MODEL=$4
TREEMARGIN=$5
DISTWT=$6
ALLFLAG=$7
ITER=$8
RATIOF=$9

$BINARY 0 ${INPUT}_Outcome.tsv ${INPUT}_Phenotype.tsv ${INPUT}_FeatureMatrix.tsv 0 $TREEMARGIN > $TREE_MODEL

$BINARY 2 ${INPUT}_Outcome.tsv ${INPUT}_Phenotype.tsv ${INPUT}_FeatureMatrix.tsv $DISTWT $ALLFLAG $ITER $RATIOF |sort -k2,2gr - > $KNN_MODEL

