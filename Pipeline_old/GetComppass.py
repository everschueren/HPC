#!/usr/bin python

import numpy, copy, sys
from operator import itemgetter


# --- read in the input
def ReadInput(file):
    '''
    Reads in matrix file and returns the dictionary as well as pointer lists:
       OUT = {Bait1:{Exp1:[Prey1,Prey2,...],
                     Exp2:[Prey1,Prey2,...],
                     ...
                     }
              Bait2:{Exp1:[Prey1,Prey2,...],
                     Exp2:[Prey1,Prey2,...],
                     ...
                     }
              ...
              }
	Returns: 
	   - M3D matrix
	   - list of preys
	   - list of baits
	   - list of experiments
	   - list of experiment clusters (named decoys)
    '''
    
    data = open(file,'r')
    D = data.readlines()
    data.close()

    M3D = {}
    experiments = D[0].strip().split('\t')[4:] #
    if len(experiments) - len(set(experiments)) != 0:
        print 'MatrixFormatingError: not all of the experiments has unique ID'
        raise
    baits = D[1].strip().split('\t')[4:] #
    preys = [i.strip().split('\t')[0] for i in D[3:]] 
    ProteinLengths = numpy.array([int(i.strip().split('\t')[2]) for i in D[3:]]) #
    
    decoys = D[2].strip().split('\t')[4:] #
    Decoys = {}
    for n,d in enumerate(baits):
        if d not in Decoys:
            composed = decoys[n].split('|')
            if len(composed) > 1: Decoys[d] = [i for i in composed if i != d]
            if len(composed) == 1 and composed[0]==d: Decoys[d] = []
            if len(composed) == 1 and composed[0]!=d: Decoys[d] = composed
            #if d in composed: Decoys[d] = []
            #else: Decoys[d] = composed
    
    # create matrix
    matrix = numpy.zeros((len(preys),len(baits)))
    for n,d in enumerate(D[3:]):
        d = d.strip().split('\t')
        
        PeptideCounts = numpy.array([float(i) for i in d[4:]]) #
        matrix[n,:] = PeptideCounts
   
    # fill in the M3D with the SIN scores
    for e in xrange(len(experiments)):
        if baits[e] not in M3D:
            M3D[baits[e]] = {experiments[e]:list(matrix[:,e] / ProteinLengths / numpy.sum(matrix[:,e]))}
        else:
            M3D[baits[e]][experiments[e]]  = list(matrix[:,e] / ProteinLengths / numpy.sum(matrix[:,e]))
    
    return (M3D,preys,list(set(baits)),experiments,Decoys)
       

# --- Calculate CompPass scores
def GetCompPASS(M3D,ACCNs,B,ordered = 1, normalize_p = 0):
    '''
    INPUT:
       M3D  = {
               BaitIndex1:
                   {
                    ExpID1:[SIN1,SIN2...],
                    EXPID2:[SIN1,SIN2...],
                    ...
                    }
                BaitIndex2:
                ...
                }
       ACCNs = [PreyID1,PreyID2,...] prey index list for SIN scores in M3D
       B = [BaitID1,BaitID2,...] bait index list for BaitIndex in M3D
    
    OUTPUT:
       ALL = [
              (BaitID1,PreyID1,Dr),
              ...
              ]
       if ordered=1, output list will be sorted
    '''

    # --- calculate 2D matrix with average SIN scores
    M2D = numpy.zeros( (len(M3D),len(ACCNs)) )

    for prey in xrange(len(ACCNs)):
        for bait in M3D:
            
            totalSIN = 0.
            for exp in M3D[bait]:
                totalSIN += M3D[bait][exp][prey]

            M2D[B.index(bait),prey] = totalSIN / len(M3D[bait])
            
    
    # --- prepare for final list
    ALL = [] 

    for bait in M3D:
        for prey in xrange(len(ACCNs)):
            
            p = 0 # number of replicates where specific
                  # prey was detected
            p_total = 0
            for exp in M3D[bait]:
                if M3D[bait][exp][prey] > 0.0: p += 1
                p_total += 1.0
                
            if normalize_p==1: p /= p_total # p normalization

            f = 0 # number of times specific prey was
                  # detected with different baits
            for bm in M2D:
                if bm[prey] > 0.0: f += 1
            if f==0: continue # happens, when SIN was too low
            #print prey, len(M3D), f, p
	    #Dr = ( ((len(M3D) / f)**p)*M2D[B.index(bait),prey] )**0.5
	    Dr = ( ( len(M3D)**(0.5*p) ) / ( f**(0.5*p) ) ) * (M2D[B.index(bait),prey] **0.5)

            ALL.append( (bait,ACCNs[prey],Dr) )

    if ordered==1:
        ALL = sorted(ALL,key=itemgetter(0,2))
    
    return ALL

def CompPASS(file, ord, norm):
	M3D, Preys, Baits, Exps, BaitClusters = ReadInput(file)
	return GetCompPASS(M3D, Preys, Baits, ordered = ord, normalize_p = norm)

if 1:
	args = sys.argv
	if len(args) < 2 or len(args) > 4:
		print """
Correct usage: python GetComppass.py <InputFile> [<sort (0/1)> <normalize p (0/1)>]
			  """
	else:
		ord,norm = 1,0
		if len(args)>=3: ord = int(args[2])
		if len(args)==4: norm = int(args[3])
		
		Scores = CompPASS(args[1],ord,norm)
		for s in Scores: print '\t'.join([str(i) for i in s])


