#!/bin/bash


GRAPH_PATH=$1
for NUM_ESTIMATORS in 100 200 500 1000 2000 5000
do
  echo $NUM_ESTIMATORS
  sbt "run ${GRAPH_PATH} ${NUM_ESTIMATORS} --numEPart=1" >> 3motif.log 2>&1
done

for NUM_ESTIMATORS in 10000 20000 50000 100000 200000 500000
do
  echo $NUM_ESTIMATORS
  sbt "run ${GRAPH_PATH} ${NUM_ESTIMATORS} --numEPart=1" >> 3motif.log 2>&1
done

NUM_ESTIMATORS=1000000
while [ $NUM_ESTIMATORS -le 10000000 ]
do
  echo $NUM_ESTIMATORS
  sbt "run ${GRAPH_PATH} ${NUM_ESTIMATORS} --numEPart=1" >> 3motif.log 2>&1
  NUM_ESTIMATORS=$((NUM_ESTIMATORS+1000000))
done

cat output3Motif.txt | awk -F ' ' '{print $2}'

