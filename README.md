## GENOM project: Evolution of structural signatures in families of soluble protein domains

This repo was made for the GENOM project in the context of the M2 Bioinformatique et Modélisation.
This project was done by Liam LG, Alexis T, and Yann Z.

The purpose of the project was the creation of a BLOSUM type matrix built from Pfam alignments of soluble domains obtained from SCOPe.
In order to create such a matrix, HCA (Hydrophobic Cluster Analysis) was used, and raw amino acid sequences were binarized according to a strong hydrophobic amino acid alphabet.

The test notebook contains the code for the creation of the full RP15 substitution matrix made on 3200+ Pfam alignments.
The tutorial notebook contains similar code but only on part of the data, with 2 Pfam alignments (that matched soluble domains with known 3D structures).
