blastn -task blastn-short -query <reads> \
-subject ./primer.fa -out <output> \
-num_threads 1 -outfmt 7 -evalue 1.0e-5
