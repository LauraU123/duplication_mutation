# Duplication mutation scripts

To reconstruct the mutations in the G gene, a series of steps are needed.
This includes reconstructing all branches in the tree, followed by finding all

1.**Reconstructing from root**

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
	

2. **Pairwise alignment to reference**

All the reconstructed sequences are aligned to each other.

Input:
	- reference fasta (a or b)
	- all reconstructed branches and terminal sequences (Step 1. output)
Output:
	- pairwise aligned sequences
	

3. **Cut out duplication**

The G duplication is cut out of the alignment. The indices must be manually provided.

Input:
	- pairwise aligned branches and sequences
Params:
	- duplication locations (start, end)
Output:
	- fasta file containing only the relevant part of the alignment


4. **Alignment of the duplication**

The cut out G duplicated sequences are Multiple Sequence Aligned to each other

Input:
	- duplicated G region
Output:
	- MSA duplicated G region
