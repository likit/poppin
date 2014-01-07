Poppin
======

Poppin uses BLAT to align genes of interest to contigs from assembly.
Then it select all contigs that match genes and save them in separate files.
Sequences in each file are aligned to crate a consensus sequence(s).

Consensus sequences from assemblies from different strains can
compared or used to build a phylogenetic tree.

Input
-----
*contigs
  Poppin accepts contigs/sequences in FASTA format.
  The input file should contain multiple contigs.
*genes
  Sequences of genes of interest in FASTA format.
  
  
Prerequisite software
---------------------

Poppin uses BLAT and CAP3. Please install them first.
