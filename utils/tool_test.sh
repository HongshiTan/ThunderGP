#!/bin/sh

date_str=`date +%Y%m%d%T`
log_path=test_log_${date_str}
mkdir -p ${log_path}

test() {

./approximation_triangle_scheme_1 0 $1  > ./${log_path}/$2_$3_t1.log
./approximation_triangle_scheme_2 0 $1  > ./${log_path}/$2_$3_t2.log
./approximation_motifs_scheme_1   0 $1  > ./${log_path}/$2_$3_m1.log
./approximation_motifs_scheme_2   0 $1  > ./${log_path}/$2_$3_m2.log

}

#test /graph_data/citeseer.txt citeseer 1
#test /graph_data/citeseer.txt citeseer 2
#test /graph_data/citeseer.txt citeseer 3
#test /graph_data/citeseer.txt citeseer 4
#test /graph_data/citeseer.txt citeseer 5

 
#test /graph_data/youtube.ungraph.txt  youtube 1
#test /graph_data/youtube.ungraph.txt  youtube 2
#test /graph_data/youtube.ungraph.txt  youtube 3
#test /graph_data/youtube.ungraph.txt  youtube 4
#test /graph_data/youtube.ungraph.txt  youtube 5

test /graph_data/LiveJournal.txt  livejournal 1
test /graph_data/LiveJournal.txt  livejournal 2
test /graph_data/LiveJournal.txt  livejournal 3
test /graph_data/LiveJournal.txt  livejournal 4
test /graph_data/LiveJournal.txt  livejournal 5
