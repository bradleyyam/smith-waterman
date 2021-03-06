#!/usr/bin/python

__author__ = "Bradley Yam"
__email__ = "bradley.yam@yale.edu"
__copyright__ = "Copyright 2021"
__license__ = "GPL"
__version__ = "1.0.0"

### Usage: python src/main.py -i <input file> -s <score file>
### Example: python src/main.py -i input.txt -s blosum62.txt
### Requirements: pandas
### Run command in root folder of the project.
### Note: Smith-Waterman Algorithm (Wavefront Optimization)

import argparse
import pandas as pd

# SW Class Matrix

# The SWMatrix holds all the information and methods needed to find the best local alignments. For efficiency, it uses seven matrices (four alignment and three traceback) and runs in O(n) time. 
# The SWMatrix is broken down into the following parts:
## score: contains the score matrix which determines the affinity between any two symbols
## o: is the penalty for opening a gap
## e: is the penalty for extending a gap
## seq1: is the first sequence prepended with a tab.
## seq2: is the second sequence prepended with a tab.
## F: This matrix contains all the final scores for every position, including matches and extensions. This is the matrix we print to the display at the very end. At each step, it saves the max of the other three alignment matrices at the same coordinates.
## M: This matrix contains all the alignment scores for matches at every position. This means that every number in this matrix represents a match between two symbols.
## Ix: This matrix contains all the alignment scores for extensions along any column (or along any x value). Every number in this matrix represents either extending an open gap or creating a new one along a column.
## Iy: This matrix contains all the alignment scores for extensions along any row (or along any y value). Every number in this matrix represents either extending an open gap or creating a new one along a row.
## TM: This is the first traceback matrix for matches. For every given match in the M matrix, it contains a number from 0-3 indicating where the previous value came from, 0 represents the M matrix, 1 represents the Ix matrix, and 2 represents the Iy matrix, 3 represents that the max was the arbitrary 0 value and hence a halt to the traceback.
## TIx: This is the second traceback matrix for extensions along any column, the numbers and their meanings are the same as the M matrix.
## TIy: This is the third traceback matrix for extensions along any row, the numbers and their meanings are the same as the previous two matrices.

# Methods

## __init__ intializes the matrix with the score matrix, the two sequences, and the gap and extension penalties.
## fillMatrix uses a dynamic programming algorithm to calculate the scores in all seven matrices. This is done as an optimization of the basic Smith-Waterman Algorithm using a wavefront.
### In the basic Smith-Waterman Algorithm, each cell in the final matrix is computed by taking the max of the the cell in the [i-1, j-1] position + affinity score, any cell in [i, k] for k < j or any cell in [l, j] for l < i + opening gap + extension gap * (i-l-1) or (j-k-1), or 0 (to get the local alignment).
### The optimized Smith-Waterman Algorithm effectively does the same thing, but using thre matrices we already know the optimal value of any match or extension up till the [i, j] cell, therefore, the algorithm doesn't need to look across all l or j.
### The optimized algorithm only needs to look for matches in [i-1, j], [i-1, j-1], [i, j-1], and gap openings or extensions in [i-1, j] and [j-1, i]. The maximum values have been precomputed.
### The tradeoff of such an algorithm is the extra space needed to hold three traceback matrices, as each alignment matrix needs to know which other matrix it came from. This is commonly known as the wavefront method. See reference paper here: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btaa777/5904262, some nice illustrations here: http://cs.rhodes.edu/welshc/COMP465_F16/Lecture10.pdf
### The optimized time is O(n) from a O(n^2) runtime.
### Additionally, the algorithm is also configurable to subtle changes such as: do we allow consecutive extensions in orthogonal directions, or must a match happen before this happens? Is this scored as a double extension, or as another gap? 

## getMax returns the maximum value in the F matrix, i.e. the best alignment score.

## getMaxCoord returns the coordinates of the maximum value in the matrix. This is used in the traceback function. 

## completeFront is a helper function in the traceback function to add the rest of the symbols on the front of the completed alignment

## completeBack is a helper function in the traceback function to add the rest of the symbols on the back of the completed alignment

## traceback enables us to reconstruct the optimal alignment by tracing our steps back through the three matrices. 0-3 represents states, but also matrices, each state corresponds to a matrix, hence the state tells us which matrix to reference for any given [i, j] coordinate.
### Since the states also tell us if it is a match or an extension and in which direction, we can use state information to directly rebuild the sequence alignment.

class SWMatrix:
    def __init__(self, score, seq1, seq2, openGap=-2, extGap=-1):
        # Store Metadata
        self.score = score
        self.o = openGap
        self.e = extGap
        self.seq1 = "\t" + seq1
        self.seq2 = "\t" + seq2

        # Initialize the final score matrix
        self.F = pd.DataFrame(0, index = [char for char in self.seq2], columns = [char for char in self.seq1])

        # Initialize the Match Matrix such that M(0,0) = 0
        self.M = pd.DataFrame(0, index = [char for char in self.seq2], columns = [char for char in self.seq1])

        #Initialize the Gap X Matric such that Ix(i, 0) = open penalty + gap penalty * i (in this instance we floor everything to 0)
        self.Ix = pd.DataFrame(0, index = [char for char in self.seq2], columns = [char for char in self.seq1])

        #Initialize the Gap Y Matric such that Ix(i, 0) = open penalty + gap penalty * j (in this instance we floor everything to 0)
        self.Iy = pd.DataFrame(0, index = [char for char in self.seq2], columns = [char for char in self.seq1])
        
        #Initialize traceback matrices which will hold a number that refers to which path it took from the previous matrix
        self.TM = pd.DataFrame(0, index = [char for char in self.seq2], columns = [char for char in self.seq1])
        self.TIx = pd.DataFrame(0, index = [char for char in self.seq2], columns = [char for char in self.seq1])
        self.TIy = pd.DataFrame(0, index = [char for char in self.seq2], columns = [char for char in self.seq1])

    def fillMatrix(self):
        for i in range(1, len(self.seq2)):
            for j in range(1, len(self.seq1)):
                #compute the score for the match matrix
                score = self.score.at[self.seq1[j], self.seq2[i]]
                a = self.M.iat[i-1, j-1] + score #match xi with yj
                b = self.Ix.iat[i-1, j-1] + score #insertion in x
                c = self.Iy.iat[i-1, j-1] + score #insertion in y
                ret1 = max([a, b, c, 0])
                self.M.iat[i,j] = ret1 
                #a number from 0-3 will tell us where the sequence came from
                self.TM.iat[i,j] = [a,b,c,0].index(ret1) 

                #compute the score for the Ix matrix
                a = self.M.iat[i-1,j] + self.o #open gap in x
                b = self.Ix.iat[i-1,j] + self.e #extend gap in x
                c = self.Iy.iat[i-1,j] + self.o #open gap after existing gap in y
                ret2 = max([a, b, c, 0]) #magic to make state match with matrix
                self.Ix.iat[i,j] = ret2
                #a number from 0-2 will tell us where the sequence came from
                self.TIx.iat[i,j] = [a,b,c,0].index(ret2) 

                #compute the score for the Iy matrix
                a = self.M.iat[i,j-1] + self.o #open gap in y
                b = self.Ix.iat[i,j-1] + self.o #open gap after existing gap in x
                c = self.Iy.iat[i,j-1] + self.e#extend gap in y
                ret3 = max([a, b, c, 0])
                self.Iy.iat[i,j] = ret3
                #a number from 0-2 will tell us where the sequence came from
                self.TIy.iat[i,j] = [a,b,c,0].index(ret3) 

                #get best score
                self.F.iat[i,j] = max([ret1, ret2, ret3])

    def getMax(self):
        return self.F.max().max()

    def getMaxCoord(self):
        i = self.F.max(axis=1).argmax()
        j = self.F.iloc[i].argmax()
        return i, j

    def completeFront(self, i, j):
        self.matchStr1.append('(')
        self.matchStr2.append('(')
        self.matchLine.append(' ')
        while i > 0 or j > 0:
            if j > 0:
                self.matchStr1.append(self.seq1[j])
            else:
                self.matchStr1.append(' ')
            self.matchLine.append(' ')
            if i > 0:
                self.matchStr2.append(self.seq2[i])
            else:
                self.matchStr2.append(' ')
            i -= 1
            j -= 1
    
    def completeBack(self, i, j):
        self.matchStr1.append(')')
        self.matchStr2.append(')')
        self.matchLine.append(' ')
        while i < len(self.seq2) or j < len(self.seq1):
            if j < len(self.seq1):
                self.matchStr1.append(self.seq1[j])
            else:
                self.matchStr1.append(' ')
            self.matchLine.append(' ')
            if i < len(self.seq2):
                self.matchStr2.append(self.seq2[i])
            else:
                self.matchStr2.append(' ')
            i += 1
            j += 1

    def traceback(self, x, y):
        # 0 = M matrix, 1 = Ix matrix, 2 = Iy matrix
        state = 0
        self.matchStr1 = []
        self.matchStr2 = []
        self.matchLine = []
        i = x
        j = y
        i, j = self.getMaxCoord()
        currVal = self.M.iat[i,j]
        while currVal != 0 and state != 3:
            if self.seq1[j] == self.seq2[i]:
                self.matchLine.append('|')
            else:
                self.matchLine.append(' ')
            if state == 0:
                nstate = self.TM.iat[i,j]
                self.matchStr1.append(self.seq1[j])
                self.matchStr2.append(self.seq2[i])
                i -= 1
                j -= 1
                state = nstate
                currVal = self.M.iat[i,j]
            elif state == 1:
                nstate = self.TIx.iat[i,j]
                self.matchStr1.append('-')
                self.matchStr2.append(self.seq2[i])
                i -= 1
                state = nstate
                currVal = self.Ix.iat[i,j]
            elif state == 2:
                nstate = self.TIy.iat[i,j]
                self.matchStr1.append(self.seq1[j])
                self.matchStr2.append('-')
                j -= 1
                state = nstate
                currVal = self.Iy.iat[i,j]
        self.completeFront(i, j)
        self.matchStr1.reverse()
        self.matchLine.reverse()
        self.matchStr2.reverse()
        i, j = self.getMaxCoord()
        self.completeBack(i+1, j+1)

### This is one way to read in arguments in Python.
parser = argparse.ArgumentParser(description='Smith-Waterman Affine Gap Algorithm')
parser.add_argument('-i', '--input', help='input file', required=True)
parser.add_argument('-s', '--score', help='score file', required=True)
parser.add_argument('-o', '--opengap', help='open gap', required=False,
default=-2)
parser.add_argument('-e', '--extgap', help='extension gap', required=False,
default=-1)
args = parser.parse_args()

### Implement your Smith-Waterman Algorithm
def runSW(inputFile, scoreFile, openGap, extGap):
    s = pd.read_csv(scoreFile, delimiter="\s+")

    with open(inputFile) as f:
        seq = f.read().splitlines()     

    m = SWMatrix(s, seq[0], seq[1], openGap, extGap)

 ### calculation
    m.fillMatrix()
    m.traceback(0,0)

 ### write output
    seqHeader = "-----------\n|Sequences|\n-----------"
    scoreHeader = "--------------\n|Score Matrix|\n--------------"
    alignHeader = "----------------------\n|Best Local Alignment|\n----------------------"
##seq
    print(seqHeader)
    print("sequence1")
    print(seq[0])
    print("sequence2")
    print(seq[1])
##score
    print(scoreHeader)
    for j in range(len(m.seq1)):
        print(m.seq1[j], end = '\t')
    print("\n", end = '')
    for i in range(len(m.seq2)):
        if i > 0: print(m.seq2[i], end = '\t')
        else: print('\t', end = '')
        for j in range(len(m.seq1)):
            print(m.F.iloc[i,j], end = "\t")
        print("\n", end = '')
#alignment
    print(alignHeader)
    print("Alignment Score:", end='')
    print(m.getMax())
    print("Alignment Results:")
    print(''.join(m.matchStr1))
    print(''.join(m.matchLine))
    print(''.join(m.matchStr2))


### Run your Smith-Waterman Algorithm
runSW(args.input, args.score, int(args.opengap), int(args.extgap))