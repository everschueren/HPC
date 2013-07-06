import sys

if len(sys.argv) != 2: 
	print 'Correct Usage: python FileCheck.py <input file>'
	sys.exit(1)

file = sys.argv[-1]
data = open(file,'r')
Lines = [i.strip().split('\t') for i in data.readlines()]
data.close()


# --- check lenghts
Lengths = [len(i) for i in Lines]
if len(set(Lengths)) == 1: print "\tLine lengths:      OK"
else:
	L = len(Lines[0])
	for x,l in enumerate(Lines):
		if len(l) != L: print "\tERROR: Line %i disagrees with #Exp (%i)!" % (x+1,L)

# --- check if values
Mistakes = []
for x,l in enumerate(Lines):
	if x<3: continue
	for y,i in enumerate(l):
		if y==0 or y==1 or y==3: continue
		try: V = float(i)
		except ValueError: Mistakes.append((x+1,y+1))

if len(Mistakes)==0: print "\tValues in lines:   OK"
else:
	for m in Mistakes: print "\tERROR: Element in line %i, row %i is not a value!" % m		

# --- check experiment uniqueness
if len(Lines[0][3:]) == len(set(Lines[0][3:])): print "\tExperiments:       OK"
else:
	for i in set(Lines[0][3:]):
		if Lines[0][3:].count(i)>1: print "\tERROR: experiment %s repeated!" % i
		
# --- check prey uniqueness
Preys = [i[0] for i in Lines[3:]]
if len(Preys) == len(set(Preys)): print "\tPreys:             OK"
else:
	for i in set(Preys):
		if Preys.count(i)>1: print "\tERROR: prey %s repeated!" % i
		
