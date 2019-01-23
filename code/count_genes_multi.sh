bams=$(ls ../04_map2/*/Aligned.out.bam)
ann=../../data/smansoni.gtf
dir=multi_map

mkdir $dir
for b in $bams; do
	o="$dir/$(echo $b | cut -f3 -d'/').counts"
	#echo $b
	#echo $o
	/proj/rpz/scripts/tmpfile_env --in $b --out $o --lock=L "featureCounts -p -M --fraction -a $ann -o %o1 %i1"
done