#! /bin/bash

samples=$(find /scratch2/data/weilguny/data/samples -name "*.fastq.gz")


for s in $samples; do
	name=$(basename $s)
	linking="ln -s $s $name"
	echo $linking
	eval $linking
done
