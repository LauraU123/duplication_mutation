# Duplication mutation scripts

To reconstruct the mutations in the G gene, a series of steps are needed.
This includes reconstructing all branches in the tree, followed by finding all

1. Reconstructing from root

The first step of the workflow includes reconstructing all of the branch sequences
from the input tree file.
It accomplishes this by recursively adding mutations on each branch to the root sequence.
Terminal sequences are copied from an input fasta file. 

Inputs:
	- root sequence
	- tree file (with mutation annotations for each branch)
	- terminal sequences (fasta)
	
Output:
	- fasta file of all reconstructed branches and terminal sequences
	


