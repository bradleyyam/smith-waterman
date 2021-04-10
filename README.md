# Smith-Waterman Algorithm

This implementation uses the optimized wavefront method and runs in O(n) time. References can be found here: 
See reference paper here: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btaa777/5904262, some nice illustrations here: http://cs.rhodes.edu/welshc/COMP465_F16/Lecture10.pdf and the basic idea here: https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm 

To install with pip, call pip3 install git+git://github.com/bradleyyam/smith-waterman.git

Otherwise, clone directly into your machine and follow the usage norms: 

```
### Usage: python src/main.py -i <input file> -s <score file>
### Example: python src/main.py -i input.txt -s blosum62.txt
```

# Folder Structure

The scripts are in /src, including another version of the same algorithm that doesn't allow double extensions in two directions.
The sample data are in /data, including score matrices and input sequences.
The output results are in /results. These can be compared with sample-output1.txt and sample-output2.txt in the /data folder.
The other scripts for setup are in /swalign for pip install compatibility.
Other reference files are in /ref.

# API

There is one method: runSW(inputFile, scoreFile, openGap, extGap):
    It takes the string of the input file location,
    the string of the score file location,
    a negative penalty for opening a gap
    and a negative penalty for extending a gap.

# Data Structures

```
## SW Class Matrix

The SWMatrix holds all the information and methods needed to find the best local alignments. For efficiency, it uses seven matrices (four alignment and three traceback) and runs in O(n) time. 
The SWMatrix is broken down into the following parts:
### score: contains the score matrix which determines the affinity between any two symbols
### o: is the penalty for opening a gap
### e: is the penalty for extending a gap
### seq1: is the first sequence prepended with a tab.
### seq2: is the second sequence prepended with a tab.
### F: This matrix contains all the final scores for every position, including matches and extensions. This is the matrix we print to the display at the very end. At each step, it saves the max of the other three alignment matrices at the same coordinates.
### M: This matrix contains all the alignment scores for matches at every position. This means that every number in this matrix represents a match between two symbols.
### Ix: This matrix contains all the alignment scores for extensions along any column (or along any x value). Every number in this matrix represents either extending an open gap or creating a new one along a column.
### Iy: This matrix contains all the alignment scores for extensions along any row (or along any y value). Every number in this matrix represents either extending an open gap or creating a new one along a row.
### TM: This is the first traceback matrix for matches. For every given match in the M matrix, it contains a number from 0-3 indicating where the previous value came from, 0 represents the M matrix, 1 represents the Ix matrix, and 2 represents the Iy matrix, 3 represents that the max was the arbitrary 0 value and hence a halt to the traceback.
### TIx: This is the second traceback matrix for extensions along any column, the numbers and their meanings are the same as the M matrix.
### TIy: This is the third traceback matrix for extensions along any row, the numbers and their meanings are the same as the previous two matrices.
```

```
## Methods

### __init__ intializes the matrix with the score matrix, the two sequences, and the gap and extension penalties.
### fillMatrix uses a dynamic programming algorithm to calculate the scores in all seven matrices. This is done as an optimization of the basic Smith-Waterman Algorithm using a wavefront.
### In the basic Smith-Waterman Algorithm, each cell in the final matrix is computed by taking the max of the the cell in the [i-1, j-1] position + affinity score, any cell in [i, k] for k < j or any cell in [l, j] for l < i + opening gap + extension gap * (i-l-1) or (j-k-1), or 0 (to get the local alignment).
### The optimized Smith-Waterman Algorithm effectively does the same thing, but using thre matrices we already know the optimal value of any match or extension up till the [i, j] cell, therefore, the algorithm doesn't need to look across all l or j.
### The optimized algorithm only needs to look for matches in [i-1, j], [i-1, j-1], [i, j-1], and gap openings or extensions in [i-1, j] and [j-1, i]. The maximum values have been precomputed.
### The tradeoff of such an algorithm is the extra space needed to hold three traceback matrices, as each alignment matrix needs to know which other matrix it came from. This is commonly known as the wavefront method. See reference paper here: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btaa777/5904262, some nice illustrations here: http://cs.rhodes.edu/welshc/COMP465_F16/Lecture10.pdf
### The optimized time is O(n) from a O(n^2) runtime.
### Additionally, the algorithm is also configurable to subtle changes such as: do we allow consecutive extensions in orthogonal directions, or must a match happen before this happens? Is this scored as a double extension, or as another gap? 

### getMax returns the maximum value in the F matrix, i.e. the best alignment score.

### getMaxCoord returns the coordinates of the maximum value in the matrix. This is used in the traceback function. 

### completeFront is a helper function in the traceback function to add the rest of the symbols on the front of the completed alignment

### completeBack is a helper function in the traceback function to add the rest of the symbols on the back of the completed alignment

### traceback enables us to reconstruct the optimal alignment by tracing our steps back through the three matrices. 0-3 represents states, but also matrices, each state corresponds to a matrix, hence the state tells us which matrix to reference for any given [i, j] coordinate.
### Since the states also tell us if it is a match or an extension and in which direction, we can use state information to directly rebuild the sequence alignment.
```

