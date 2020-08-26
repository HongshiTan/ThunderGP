#!/bin/bash

tmp_string=`date +%Y%m%d%T`
date_str=${tmp_string//:}
log_path=test_log_${date_str}
mkdir -p ${log_path}

test() {
var=$3
for (( id=1; id<=${var}; id++ ))
do
	echo " $1  $2  ${id}"
	./approximation_triangle_scheme_3 0 $1  > ./${log_path}/$2_${id}_t3.log
	./approximation_motifs_scheme_3 0 $1    > ./${log_path}/$2_${id}_m3.log
done;
}

test /graph_data/citeseer.txt citeseer 1

test /graph_data/youtube.ungraph.txt  youtube 2

test /graph_data/LiveJournal.txt  livejournal 1
