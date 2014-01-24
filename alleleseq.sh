cat *cap.contigs > all_contigs.fa
blastall -p blastn -m 7 -i all_contigs.fa -d st15genes.fa > all_contigs.fa.xml
python ~/poppin/alleleseq.py ../st15genes.fa all_contigs.fa.xml > alleleseq.fa
