import numpy
from numpy import zeros

# Pre local alignment
def localalignmentcost(xc, yc, costs):
    
    if xc == yc: return costs[0] # match
    if xc == '-' or yc == '-': return costs[2] # gap
    return costs[1]

def localalignmenttraceback(V, x, y, costs):
    """ Trace back from given cell in local-alignment matrix V """
    # get i, j for maximal cell
    i, j = numpy.unravel_index(numpy.argmax(V), V.shape)
    xscript, alx, aly, alm = [], [], [], []
    while (i > 0 or j > 0) and V[i, j] != 0:
        diag, vert, horz = 0, 0, 0
        if i > 0 and j > 0:
            diag = V[i-1, j-1] + localalignmentcost(x[i-1], y[j-1], costs)
        if i > 0:
            vert = V[i-1, j] + localalignmentcost(x[i-1], '-', costs)
        if j > 0:
            horz = V[i, j-1] + localalignmentcost('-', y[j-1], costs)
        if diag >= vert and diag >= horz:
            match = x[i-1] == y[j-1]
            xscript.append('M' if match else 'R')
            alm.append('|' if match else ' ')
            alx.append(x[i-1]); aly.append(y[j-1])
            i -= 1; j -= 1
        elif vert >= horz:
            xscript.append('D')
            alx.append(x[i-1]); aly.append('-'); alm.append(' ')
            i -= 1
        else:
            xscript.append('I')
            aly.append(y[j-1]); alx.append('-'); alm.append(' ')
            j -= 1
    xscript = (''.join(xscript))[::-1]
    alignment = '\n'.join(map(lambda x: ''.join(x), [alx[::-1], alm[::-1], aly[::-1]]))

    
    return  alignment

# Pre global alignment
def globalalignmentcost (xc, yc, costs):
    """ Cost function assigning 0 to match, 2 to transition, 4 to
        transversion, and 8 to a gap """
    if xc == yc: return costs[0] # match
    if xc == '-' or yc == '-': return costs[3] # gap
    minc, maxc = min(xc, yc), max(xc, yc)
    if minc == 'A' and maxc == 'G': return costs[1] # transition
    elif minc == 'C' and maxc == 'T': return costs[1] # transition
    return costs[2] # transversion

