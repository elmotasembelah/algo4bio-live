from math import ceil
import numpy 
from numpy import zeros
from .pre_algorithms import Boyermoore, k_mer_indexing , Seqalignment

class Algorithms:
    
    @staticmethod
    def localAlignment(seqs, costs):

        x = seqs[0]
        y = seqs[1]

        resultdict = {}
        ''' Calculate local alignment values of sequences x and y using
            dynamic programming.  Return maximal local alignment value. '''

        V = numpy.zeros((len(x)+1, len(y)+1), dtype=int)
        for i in range(1, len(x)+1):
            for j in range(1, len(y)+1):
                V[i, j] = max(V[i-1, j-1] + Seqalignment.localalignmentcost(x[i-1], y[j-1], costs), # diagonal
                            V[i-1, j  ] + Seqalignment.localalignmentcost(x[i-1], '-', costs),    # vertical
                            V[i  , j-1] + Seqalignment.localalignmentcost('-',    y[j-1], costs), # horizontal
                            0)                               # empty

        argmax = numpy.where(V == V.max())
        print (numpy.max(argmax))
        traceback = Seqalignment.localalignmenttraceback(V, x, y, costs)

        vector = numpy.vectorize(numpy.float)
        best = numpy.max(argmax)


        matrix = str(V) 
        resultMatrix = []
        start = 1
        while len(matrix) > 0:

 
            resultMatrix.append(matrix[start: matrix.find(']')+1])
            matrix = matrix[matrix.find(']')+2:]


        resultdict['matrix'] = resultMatrix
        resultdict['score'] = "Best score=%d, in cell %s" % (best , numpy.unravel_index(numpy.argmax(V), V.shape))
        resultdict['align3'] = traceback

        return resultdict


    @staticmethod
    def globalAlignment(seqs, costs):

        x = seqs[0]
        y = seqs[1]

        resultdict = {}
        """ Calculate global alignment value of sequences x and y using
            dynamic programming.  Return global alignment value. """
        D = zeros((len(x)+1, len(y)+1), dtype=int)
        for j in range(1, len(y)+1):
            D[0, j] = D[0, j-1] + Seqalignment.globalalignmentcost('-', y[j-1], costs)
        for i in range(1, len(x)+1):
            D[i, 0] = D[i-1, 0] + Seqalignment.globalalignmentcost(x[i-1], '-', costs)
        for i in range(1, len(x)+1):
            for j in range(1, len(y)+1):
                D[i, j] = min(D[i-1, j-1] + Seqalignment.globalalignmentcost(x[i-1], y[j-1], costs), # diagonal
                            D[i-1, j  ] + Seqalignment.globalalignmentcost(x[i-1], '-', costs),    # vertical
                            D[i  , j-1] + Seqalignment.globalalignmentcost('-',    y[j-1], costs)) # horizontal

        matrix = str(D) 
        resultMatrix = []
        start = 1
        while len(matrix) > 0:

 
            resultMatrix.append(matrix[start: matrix.find(']')+1])
            matrix = matrix[matrix.find(']')+2:]


        resultdict['matrix'] = D
        resultdict['value'] = D[len(x), len(y)]

        return resultdict
 

    @staticmethod
    def hashmapIndex(seqs, index):
        p = seqs[0]
        t = seqs[1]
        resultdict = {}
        count = 0
        k = index.k
        offsets = []
        for i in index.query(p):
            count = count + 1 
            if p[k:] == t[i+k:i+len(p)]:  # verify that rest of P matches
                offsets.append(i)

        resultdict['positions'] = offsets
        resultdict['count'] = count

        return resultdict 


    @staticmethod
    def binarySearchKmerIndex(seqs, index):
        
        p = seqs[0]
        t = seqs[1]
        resultdict = {}
        count = 0

        k = index.k
        offsets = []
        for i in index.query(p):
            count = count + 1
            if p[k:] == t[i+k:i+len(p)]:  # verify that rest of P matches
                offsets.append(i)

        resultdict['positions'] = offsets
        resultdict['count'] = count
        return resultdict


    @staticmethod
    def editDistance(seqs):
        seq1 = seqs[0]
        seq2 = seqs[1]
        resultdict = {}
        differences = []
        """ Calculate edit distance between sequences x and y using
            matrix dynamic programming.  Return distance. """
        D = zeros((len(seq1)+1, len(seq2)+1), dtype=int)
        D[0, 1:] = range(1, len(seq2)+1)
        D[1:, 0] = range(1, len(seq1)+1)
        for i in range(1, len(seq1)+1):
            for j in range(1, len(seq2)+1):
                delt = 1 if seq1[i-1] != seq2[j-1] else 0
                D[i, j] = min(D[i-1, j-1]+delt, D[i-1, j]+1, D[i, j-1]+1)


        resultdict['count'] = D[len(seq1), len(seq2)]
        resultdict['positions'] = D

        return resultdict


    @staticmethod
    def hammingDistance(seqs):
        seq1 = seqs[0]
        seq2 = seqs[1]
        resultdict = {}
        differences = [] 

        if len(seq1) != len(seq2):
            resultdict['alig1'] = "The sequences that has been entered are not of the same length, Please review the entered sequences"
            return resultdict

        nmm = 0
        for i in range(0, len(seq1)):
            if seq1[i] != seq2[i]:
                differences.append(i)
                nmm += 1

        resultdict['count'] = nmm
        resultdict['positions'] = differences
        return resultdict


    @staticmethod
    def boyerMooreApproximatematching(self, seqs, maxdist):

        p = seqs[0]
        t = seqs[1]
        resultdict = {}
        comparisionsCount = 0

        segment_length = int( ceil(  len(p) / (maxdist+1)) ) 
        all_matches = set()
        
        for i in range (maxdist+1):
            start = i*segment_length
            end = min((i+1)*segment_length, len(p))
            p_bm = Boyermoore.BoyerMoore(p[start:end], 'ACGT')
            tempseqs = [p[start:end], t]
            matches, cc = self.boyer_MooreExactMatching(tempseqs, p_bm)
            comparisionsCount = comparisionsCount + cc

            for m in matches:
                if m < start or m-start+len(p) > len(t):
                    continue

                mismatches = 0 

                for j in range (0, start):
                    comparisionsCount = comparisionsCount + 1
                    if not p[j] == t[m-start+j]:
                        mismatches += 1
                        if mismatches > maxdist:
                            break

                for j in range (end, len(p)):
                    comparisionsCount = comparisionsCount + 1
                    if not p[j] == t[m-start+j]:
                        mismatches += 1
                        if mismatches > maxdist:
                            break 

                if mismatches <= maxdist:    
                    all_matches.add(m - start)
        

        resultdict['positions'] = list(all_matches)    
        resultdict['count'] = comparisionsCount
        return resultdict


    @staticmethod
    def naiveApproximateMatching(seqs, maxdist = 2):

        p = seqs[0]
        t = seqs[1]
        resultdict = {}
        comparisionsCount = 0
    
        occurrences = []
        for i in range(0, len(t) - len(p) + 1): # for all alignments
            nmm = 0
            
            for j in range(0, len(p)):          # for all characters
                comparisionsCount = comparisionsCount + 1
                if t[i+j] != p[j]:               # does it match?
                    nmm += 1                     # mismatch
                    if nmm > maxdist:
                        break                    # exceeded maximum distance
            if nmm <= maxdist:
                occurrences.append((i, nmm))
        resultdict['positions'] = occurrences
        resultdict['count'] = comparisionsCount
        return resultdict


    @staticmethod
    def boyer_MooreExactMatching(seqs, p_bm = 0):
        p = seqs[0]
        t = seqs[1]
        comparisionsCount = 0

        check = False
        if p_bm == 0 :
            p_bm = Boyermoore.BoyerMoore(p, alphabet='ACGT')
            check = True
        """ Do Boyer-Moore matching """
        resultdict = {}
        i = 0
        occurrences = []
        while i < len(t) - len(p) + 1:
            shift = 1
            mismatched = False
            for j in range(len(p)-1, -1, -1):
                comparisionsCount = comparisionsCount + 1
                if p[j] != t[i+j]:
                    skip_bc = p_bm.bad_character_rule(j, t[i+j])
                    skip_gs = p_bm.good_suffix_rule(j)
                    shift = max(shift, skip_bc, skip_gs)
                    mismatched = True
                    break
            if not mismatched:
                occurrences.append(i)
                skip_gs = p_bm.match_skip()
                shift = max(shift, skip_gs)
            i += shift

        if check :
            resultdict['positions'] = occurrences    
            resultdict['count'] = comparisionsCount
            return resultdict
        return occurrences, comparisionsCount 


    @staticmethod
    def naiveExactMatching(seqs):
        p = seqs[0]
        g = seqs[1]
        resultdict = {}
        comparisionsCount = 0
        occurrences = []
        for i in range(len(g) - len(p) + 1):  # loop over alignments
            match = True
            for j in range(len(p)):  # loop over characters
                comparisionsCount = comparisionsCount + 1
                if g[i+j] != p[j]:  # compare characters
                    match = False  # mismatch; reject alignment
                    break
            if match:
                occurrences.append(i)  # all chars matched; record
                
        resultdict['positions'] = occurrences
        resultdict['count'] = comparisionsCount
        return resultdict