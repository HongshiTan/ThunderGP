#!/bin/sh
echo "r19"
cat $1/r19.log  | grep $2
echo "r21"
cat $1/r21.log 	| grep $2
echo "r24"
cat $1/r24.log 	| grep $2
echo "mg"
cat $1/mg.log 	| grep $2
echo "wt"
cat $1/wt.log 	| grep $2
echo "gg"
cat $1/gg.log 	| grep $2
echo "pk"
cat $1/pk.log 	| grep $2
echo "lj"
cat $1/lj.log 	| grep $2
echo "tw"
cat $1/tw.log 	| grep $2



#[[:blank:]]