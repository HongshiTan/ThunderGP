#!/bin/bash


GRAPH_PATH=$1
GRAPH_NAME=$(cut -d'/' -f3 <<< ${GRAPH_PATH} | cut -d'.' -f1)
YYMMDD=$(date '+%Y%m%d')
LOGFILE=3motif.log
OUTPUTFILE=output3Motif.txt

echo "Processing graph: " $GRAPH_NAME

for RUN_ID in {1..10}
do
  if [ -f $LOGFILE ]; then
    rm $LOGFILE
  fi

  if [ -f $OUTPUTFILE ]; then
    rm $OUTPUTFILE
  fi
  echo "Iteration " $RUN_ID
  ./run.sh $GRAPH_PATH
  mv $LOGFILE log/${YYMMDD}_3motif_${GRAPH_NAME}_${RUN_ID}.log
  mv $OUTPUTFILE output/${YYMMDD}_3motif_${GRAPH_NAME}_${RUN_ID}.txt
done

./merge_results.sh ${GRAPH_NAME}
