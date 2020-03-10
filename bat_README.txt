# created:
# 	Joe Bemister
# 	28/08/2015
#
# updated:
#	Joseph Bemister-Buffington
#	10/03/2020
#

Bfactor Alignment Tool (BAT) requires the input of 3 files:

1. (-a) an alignment file prepared from the parsable data retrieved from a DALI pairwise alignment.
2. (-t) a target pdb file that corresponds to the first position (target position) of the alignment file.
3. (-q) a query pdb file that corresponds to the second position (query position) of the alignment file.

The script takes the residue number, chain ID, residue ID, and bfactor value of each alpha carbon present in the target and query pdb files.
It then outputs these values in the residue positions identified in the alignment files. The output file is a csv created in the directory
of the alignment file and is entitled "<name of alignment file>_bat.csv".
