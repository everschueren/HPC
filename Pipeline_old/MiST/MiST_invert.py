########################################################################################
#                  MiST - Mass spectrometry interaction STatistics                     #
#                                Peter Cimermancic                                     #
#                                   Krogan Lab                                         #
#                                    May 2010                                          #
########################################################################################

import numpy, copy, sys, mdp
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
       

def ThreeMetrics(M3D,ACCNs,Decoys,filter=0):
    '''
    Returns the three metrics
    '''

    # --- Normalize experiments
    M3D_normal = {}
    for bait in M3D:
        if bait == None: continue
        for exp in M3D[bait]:
            S = numpy.array(M3D[bait][exp])
            Z = S / numpy.sum(S)
        
            if bait not in M3D_normal:
                M3D_normal[bait] = {exp:Z}
            else:
                if exp not in M3D_normal[bait]:
                    M3D_normal[bait][exp] = Z
                else:
                    print 'ExperimentError: repeated experiment %i when normalizing data' % exp
                    raise

    M3D = None # delete hash

    # --- Reproducibility and Averages
    L = len(ACCNs)

    AvgReplicates = {}
    EntReplicates = {} # keep entropies for replicates
    Baits = []
    
    for bait in M3D_normal:
        Baits.append(bait)
        averages = numpy.zeros(L)
        entropies = numpy.zeros(L)
        TEMP = []
        for exp in M3D_normal[bait]:
            TEMP.append(M3D_normal[bait][exp])
        TEMP = numpy.array(TEMP)
        shp = numpy.shape(TEMP)
    
        for x in xrange(shp[1]):
            sumx = numpy.sum(TEMP[:,x])
            if(sumx>0):
                protx = TEMP[:,x] / sumx    # normalize per protein in replicates
            else:
                protx = TEMP[:,x]
            for y in protx:
                if y > 0.0 and y != 1.:   # if TEMP[y,x] == 0: H == 0 by definition
                    entropies[x] += y * numpy.log2(y)
                    averages[x] += y
                elif y == 1.:             # to avoid confusion to all-zeros case
                    entropies[x] += (y-1e-10) * numpy.log2(y-1e-10)
            averages[x]  = sumx
       
        if shp[0] != 1:
            entropies /= (numpy.log2(1./shp[0])) # devide entropies by max entropy
        elif shp[0] == 1:
            #entropies *= (-1)
            entropies = numpy.ones(shp[1])
        averages /= shp[0]

        if bait not in EntReplicates:
            EntReplicates[bait] = entropies
            AvgReplicates[bait] = averages
        else:
            print 'BaitError: bait %s appears multiple times' % bait
            raise

    Reproducibility = []
    Abundancy = []
    for bait in Baits:
        Reproducibility.append(numpy.array(EntReplicates[bait]))
        Abundancy.append(numpy.array(AvgReplicates[bait]))
    Reproducibility = numpy.array(Reproducibility)
    Abundancy = numpy.array(Abundancy)

    
    # --- Specificity (from AvgReplicates)
    '''
        Calculated as a fraction of abundancies
    '''
    SHP = numpy.shape(Abundancy)
    summ = numpy.zeros(SHP[1]) # sum over baits
    Specificity = numpy.zeros(SHP) # probabilities per protein over all baits
    
    for i in xrange(SHP[1]):         # over preys

        for j in xrange(SHP[0]):     # over baits
            
            M = Decoys[Baits[j]]                 # keep the list of decoys
            MP = None
            if len(M) > 0:

                MP = [Baits.index(z) for z in M]  # keep the pointers for decoys
                MP = list(set(range(len(Baits))) ^ set(MP))
            else: MP = range(len(Baits))
            
            summ[i] = numpy.sum(Abundancy[MP,i])
                
            if summ[i] > 0.:
                Specificity[j,i] = Abundancy[j,i] / summ[i]
            else: Specificity[j,i] = 0.

            if filter == 1:
                if Specificity[j,i] == 1.0 and Reproducibility[j,i] < 10e-7: Specificity[j,i] = 0.0
    return (1-Reproducibility,1-Abundancy,1-Specificity,Baits) #changed by erik v
    # return (Reproducibility,Abundancy,Specificity,Baits)
        
# --- Write output for the three metrics
def OutputMetrics(R,A,S,B,P,out=1,FileName=None):
    '''
    Write the output file for the three metrics
    '''

    if out==1:
        output = open(FileName+'_metrics.out','w')
        output.write('Bait\tPrey\tReproducibility\tAbundance\tSpecificity\n')
    X,Y = numpy.shape(R)
    Matrix = numpy.zeros((X*Y,3))
    Pairs = []
    #Matrix1 = []
    c = 0
    for y in xrange(Y):
        prey = P[y]
        for x in xrange(X):
            bait = B[x]
            Pairs.append((bait,prey))
            #if R[x,y]+A[x,y]+S[x,y] > 0: Matrix1.append(numpy.array([R[x,y],A[x,y],S[x,y]]))
            Matrix[c,0] = R[x,y]
            Matrix[c,1] = A[x,y]            
            Matrix[c,2] = S[x,y]
            c += 1
            if out==1:
                line = '%s\t%s\t%.4f\t%.4f\t%.4f\n' % (bait,prey,R[x,y],A[x,y],S[x,y])
                output.write(line)
    #Matrix1 = numpy.array(Matrix1)
    #print len(Matrix1)
    #print len(Matrix)
    if out==1:
        output.close()
    
    return (Matrix,Pairs)


# --- No training of the three metrics
def NoTraining(R,A,S,B,P):
    '''
    Calculate MiST scores from the three metrics withouth training
    '''

    X,Y = numpy.shape(R)
    c = 0
    Scores = []
    for y in xrange(Y):
        prey = P[y]
        for x in xrange(X):
            bait = B[x]

            r = R[x,y]
            a = A[x,y]            
            s = S[x,y]
            
            mist = 0.68551*s + 0.30853*r + 0.00596*a
            Scores.append(mist)
            
    return Scores

# --- Perform PCA
def PCA(matrix,pairs,filter=0):
    '''
    Calculate PCA
    '''

    pca = mdp.nodes.PCANode(output_dim=3)
    pca.train(matrix)
    
    P = pca(matrix)
    PcaScore = P[:,0]

    Eigens = pca.d
    Eigens /= numpy.sum(Eigens)

    EIV = numpy.transpose(pca.v[:,0]) / numpy.sum(pca.v[:,0])

    Variance = [Eigens[0],Eigens[0]+Eigens[1],numpy.sum(Eigens)]
   
    scores = []
    for i in matrix:
        scores.append(numpy.sum(Eigens*i))
       
    scores = PcaScore #numpy.array(scores)
    Scores = (scores - numpy.min(scores)) / (numpy.max(scores) - numpy.min(scores))
    
    if filter == 2:
        for x in xrange(len(Scores)):
            if matrix[x,0] < 10e-7 and matrix[x,2] == 1: Scores[x] = 0.1

    return (Scores,Variance,EIV)
    
# --- Output final PCA score
def OutputPCA(Pairs,scores,filename):
    
    output = open(filename+'_mist.out','w')
    output.write('Bait\tPrey\tMiST score\n')
    
    Z = zip([i[0] for i in Pairs],[i[1] for i in Pairs],scores)
    sortedPCA = sorted(Z,key=itemgetter(0,2))

    for line in sortedPCA:
        l = '%s\t%s\t%.5f\n' % line
        output.write(l)

    output.close()


# --- MAIN ---
def __main__():

#     print '''
# \t*****************************************************************
# \t*  Wellcome to MiST - Mass spectormetry Interaction STatistics  *
# \t*       written by Peter Cimermancic (Nevan Krogan Lab)         *
# \t*                           May 2010                            *
# \t*****************************************************************
# 
#           '''

    if len(sys.argv) < 5:
        explain = '''
Correct Usage: python MiST.py <input>  <output>  <filter (0/1)> <training (0/1)>
Vignette for input file:
\tInput file should contain tab separated information in columns
\tas in this example:\n
\t\t#	#\t#\t#	A1	A2	B1	B2	C1	C2
\t\tPreys\t#	Length\t#	A	A	B	B	C	C
\t\tBC\t#\t#\t#\tA\tA\tA\tA\tC\tC
\t\tprot1\t#	188\t#	12	4	7	24	16	21
\t\tprot2\t#	157\t#	0	0	0	1	0	2
\t\tprot3\t#	723\t#	9	21	18	57	24	0
\t\tprot4\t#	186\t#	6	10	7	14	15	21
\t\tprot5\t#	988\t#	0	0	0	0	0	12\n
\tTwo files will be generated - *_metrics.out and *_mist.out, both
\tcontaing bait,prey pair and either the three metrics (reproducibility,
\tabundance,specificity) or MiST score, respectively.

\tThe third line (BC - Bait Composition) is useful when one would like omit
\tcertain baits in a calculation of the specificity of the given bait. For 
\texample, pull-downs can be done with full-length protein as well as with
\tits domains. In this case, one can expect to find similar preys in both cases.
\tThe field in this line is a bait that should be omitted to calculate specificity,
\tOR bait itself (this is mandatory). If there are multiple baits to be omitted,
\tseparate them with '|' and use no white spaces in-between.

\tFilter argument filters out the preys detected only ones if 1.
    
                  '''
        print explain
        
    else:
    
		training = int(sys.argv[-1])
	    
		A,B,C,E,D = ReadInput(sys.argv[-4])
		print '\tNumber of Preys: %i' % len(B)
		print '\tNumber of Experiments: %i' % len(E)
		print '\tNumber of Baits: %i' % len(C)
		R,A,S,Baits = ThreeMetrics(A,B,D,int(sys.argv[-2]))
		Matrix,Pairs = OutputMetrics(R,A,S,Baits,B,out=1,FileName=sys.argv[-3])
        
		if training == 1: 
			score,variance,eigens = PCA(Matrix,Pairs,int(sys.argv[-2]))
			print '''\tPercentage of variance described (cumulatively):\n\t\tPC1: %.5f\n\t\tPC2: %.5f\n\t\tPC3: %.5f''' % tuple(variance)

		 	print '''\tEigenvector - weights:\n\t\tReproducibility: %.5f\n\t\tAbundance: %.5f\n\t\tSpecificity: %.5f''' % tuple(eigens)

		elif training == 0:
			score = NoTraining(R,A,S,Baits,B)	
    
		OutputPCA(Pairs,score,sys.argv[-3])
		# print 'Thank you for using MiST!\n'


__main__()
