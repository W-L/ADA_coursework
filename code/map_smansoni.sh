
forward=$(ls ../00_links/*_1.fastq.gz)

for f in $forward; do
	r=$(ls ../00_links/$(basename $f _1.fastq.gz)_2.fastq.gz)
	/proj/rpz/scripts/tmpfile_env --in ../01_index/smansoni_index_fixed --in $f --in $r --dir $(basename $f _1.fastq.gz) --lock=L 'STAR --genomeDir %i1  --readFilesIn %i2 %i3 --runThreadN 8 --readFilesCommand zcat --outFilterMultimapNmax 20 --outFileNamePrefix %d1/ --outFilterType BySJout --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMtype BAM Unsorted --outTmpDir $TMPDIR/tmp --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3'
done

