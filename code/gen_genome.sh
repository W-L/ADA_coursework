#! /bin/bash

#STAR --runMode genomeGenerate --genomeDir smansoni_index --sjdbGTFfile ../../data/smansoni.gtf --genomeFastaFiles ../../data/schistosoma_mansoni.PRJEA36577.WBPS10.genomic_masked.fa --runThreadN 8 

/proj/rpz/scripts/tmpfile_env --in ../../data/smansoni.gtf --in ../../data/schistosoma_mansoni.PRJEA36577.WBPS10.genomic_masked.fa --dir smansoni_index 'STAR --runMode genomeGenerate --genomeDir %d1 --sjdbGTFfile %i1 --genomeFastaFiles %i2 --runThreadN 8'

STAR --runMode genomeGenerate --genomeDir smansoni_index_fixed --sjdbGTFfile ../../data/smansoni.gtf --genomeFastaFiles ../../data/schistosoma_mansoni.PRJEA36577.WBPS10.genomic_masked.fa --runThreadN 4

/proj/rpz/scripts/tmpfile_env --in ../../data/smansoni.gtf --in ../../data/schistosoma_mansoni.PRJEA36577.WBPS10.genomic_masked.fa --dir smansoni_index 
STAR --runMode genomeGenerate --genomeDir smansoni_index_nogtf --genomeFastaFiles ../../data/schistosoma_mansoni.PRJEA36577.WBPS10.genomic_masked.fa --runThreadN 4 --genomeChrBinNbits 15



# working
STAR --runMode genomeGenerate --genomeDir smansoni_index_fixed --sjdbGTFfile ../../data/smansoni.gtf --genomeFastaFiles ../../data/schistosoma_mansoni.PRJEA36577.WBPS10.genomic_masked.fa --runThreadN 4 --genomeSAindexNbases 6