Poppin
======

Poppin uses BLAT to align a gene of interest to contigs from assembly.

Then it extract parts of a contig that match the gene and reassemble them to
create a consensus sequence using CAP3.

Consensus sequences from assemblies from different strains can
compared or used to build a phylogenetic tree.

Input
-----

  Poppin accepts contigs/sequences and gene sequences in FASTA format.
  The input files should contain multiple sequences.
  
Prerequisite software
---------------------

Poppin uses BLAT and CAP3. Please install them first.
