#! /bin/bash

enatool=../tools/enaBrowserTools-master/python3/enaDataGet
metadata=../data/PRJEB14695.txt


acc=$(grep "whole" $metadata | cut -f4)



for a in $acc; do
	echo $a
	$enatool -m -f fastq -d ../data $a
done
