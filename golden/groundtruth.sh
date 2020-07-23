#!/bin/bash

test() {

if [ ! -d "/graph_data/peregrine/$1" ];
	then  
		mkdir -p /graph_data/peregrine/$1 
		./bin/convert_to_binary.sh /graph_data/$2 1134217728 /graph_data/peregrine/$1/
fi

bin/count  /graph_data/peregrine/$1/ 3-motifs 8 > $1.log


}

#test citseer  citeseer.txt 

test LiveJournal  LiveJournal.txt

test youtube com-youtube.ungraph.txt


