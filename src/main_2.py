#!/usr/bin/python

__author__ = "Bradley Yam"
__email__ = "bradley.yam@yale.edu"
__copyright__ = "Copyright 2021"
__license__ = "GPL"
__version__ = "1.0.0"

### Usage: python main.py -i <input file> -s <score file>
### Example: python main.py -i input.txt -s blosum62.txt
### Note: Smith-Waterman Algorithm

import argparse
import pandas as pd

# SW Class Matrix

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
                # if i == 46 and j == 37: print(i, j, a, b, c, score)
                self.M.iat[i,j] = ret1 
                #a number from 0-3 will tell us where the sequence came from
                self.TM.iat[i,j] = [a,b,c,0].index(ret1) 

                #compute the score for the Ix matrix
                a = self.M.iat[i-1,j] + self.o #open gap in x
                b = self.Ix.iat[i-1,j] + self.e #extend gap in x
                ret2 = max([a, b, -1, 0]) #magic to make state match with matrix
                # if i == 45 and j == 37: print(i, j, a, b)
                self.Ix.iat[i,j] = ret2
                #a number from 0-2 will tell us where the sequence came from
                self.TIx.iat[i,j] = [a,b,-1,0].index(ret2) 

                #compute the score for the Iy matrix
                a = self.M.iat[i,j-1] + self.o #open gap in y
                b = self.Iy.iat[i,j-1] + self.e #extend gap in y
                ret3 = max([a, -1, b, 0])
                self.Iy.iat[i,j] = ret3
                # if i == 46 and j == 37: print(i, j, a, b)
                #a number from 0-2 will tell us where the sequence came from
                self.TIy.iat[i,j] = [a,-1,b,0].index(ret3) 

                #get best score
                self.F.iat[i,j] = max([ret1, ret2, ret3])

    #returns the highest score in the matrix
    def getMax(self):
        return self.F.max().max()

    #returns the coordinates of the highest score in the matrix
    def getMaxCoord(self):
        i = self.F.max(axis=1).argmax()
        j = self.F.iloc[i].argmax()
        return i, j

    #fills out the front of the alignemnt
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
    
    #fills out the back of the alignemnt
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

    #beginning from the max coord, traces back based on the trace matrices to find the best local alignment
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
        #the algorithm will terminate when the score hits 0, i.e. state == 3 or currVal = 0.
        #states represent the matching state, 0 = match, 1 = skip on the x-axis, 2 = skip on the y-axis, 3 = halt, traced from 0 minimum.
        while currVal != 0 and state != 3:
            if self.seq1[j] == self.seq2[i]:
                self.matchLine.append('|')
            else:
                self.matchLine.append(' ')
            # print(state, i, j, self.seq1[j], self.seq2[i])
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

    # print(m.M.iloc[44:47, 35:38])
    # print(m.Ix.iloc[44:47, 35:38])
    # print(m.Iy.iloc[44:47, 35:38])

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