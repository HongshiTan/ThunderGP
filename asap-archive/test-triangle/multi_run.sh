#!/bin/bash


GRAPH_PATH=$1
GRAPH_NAME=$(cut -d'/' -f3 <<< ${GRAPH_PATH} | cut -d'.' -f1)
YYMMDD=$(date '+%Y%m%d')
LOGFILE=triangle.log
OUTPUTFILE=output_triangle.txt

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
  sbt "run ${GRAPH_PATH} 0 10 --numEPart=2" >> $LOGFILE 2>&1
  mv $LOGFILE log/${YYMMDD}_triangle_${GRAPH_NAME}_${RUN_ID}.log
  mv $OUTPUTFILE output/${YYMMDD}_triangle_${GRAPH_NAME}_${RUN_ID}.txt
done

./merge_results.sh ${GRAPH_NAME}
