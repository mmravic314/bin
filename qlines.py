## Input a MASTER output.seq file. convert sequence to single residues, count number of unique sequences under an RMSD input

# input arg 1: path to .seq file
# input arg 2: rmsd cut off (assumes this value is first column of .seq file, so could be CA or fullbackbone RMSD)

import sys, os
from PDButil import UnNatAA

store = []
rmsdCut = float( sys.argv[2] )
with open( sys.argv[1] ) as file:
	for i in file: 
		if len(i )< 1: continue
		seq = ''.join( [ UnNatAA[res] for res in i.split()[2:] ] )
		rmsd= float( i.split()[0] )
		
		if seq not in store and rmsd <= rmsdCut:
			store.append( seq )
		print seq, '\n'

print os.path.basename(sys.argv[1]), len( store )