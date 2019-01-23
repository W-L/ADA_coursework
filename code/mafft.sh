# align seqs

prots=$(ls extension/*)
cat $prots >multi_fasta

mafft-einsi --ep 0 --genafpair --maxiterate 1000 multi_fasta >multi_align_out

~/software/iqtree-1.6.6-MacOSX/bin/iqtree -s multi_align_out -m PMB+F+R3