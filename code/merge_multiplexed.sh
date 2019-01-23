#! /bin/bash

dirs=$(ls)

for d in $dirs; do
	echo $d
	forward=$(ls $d/*_1.fastq.gz)
	reverse=$(ls $d/*_2.fastq.gz)

	forward_merger="zcat $forward | gzip -c >$d/sample_1.fastq.gz"
	echo $forward_merger

	reverse_merger="zcat $reverse | gzip -c >$d/sample_2.fastq.gz"
	echo $reverse_merger
	
	#eval $forward_merger
	eval $reverse_merger


done
