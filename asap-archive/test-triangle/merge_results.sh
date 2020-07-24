#!/bin/bash

GRAPH_NAME=$1
YYYYMMDD=$2

MERGED_FILE=output/${YYYYMMDD}_triangle_${GRAPH_NAME}_merged.txt
cat output/${YYYYMMDD}_triangle_${GRAPH_NAME}_1.txt | awk -F ' ' '{print $1}' > $MERGED_FILE
TMP_FILE=output/tmp_results.txt
INTMD_FILE=output/intmd_merge.txt

for RUN_ID in {1..10}
do
  cat output/${YYYYMMDD}_triangle_${GRAPH_NAME}_${RUN_ID}.txt | awk -F ' ' '{print $2}' > $TMP_FILE
  paste -d'|' $MERGED_FILE $TMP_FILE > $INTMD_FILE
  #paste -d'|' $MERGED_FILE $TMP_FILE
  echo $RUN_ID
  mv $INTMD_FILE $MERGED_FILE
  cat $MERGED_FILE
done;

rm $TMP_FILE
