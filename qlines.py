## Input a MASTER output.seq file. convert sequence to single residues, count number of unique sequences under an RMSD input
#  also report length of overall sequence searched and number of floating fragments found

# input arg 1: path to .seq file
# input arg 2: rmsd cut off (assumes this value is first column of .seq file, so could be CA or fullbackbone RMSD)
# input arg 3: path to match file

# example command line
# for i in *.seq; do python ~/bin/qlines.py $i ${i::-3}match ; done


import sys, os, math
from PDButil import UnNatAA

def score( f, l, m):
	return (2**f)*l*m/48


# find number of separate fragments found
mFile = open( sys.argv[2], 'rU' )
try:
	m = mFile.readlines()[0]
	fragNum = len( [ tuple( x.strip(',').strip('()').split(',') ) for x in m.split('[')[-1].split(']')[0].split() ] )
except:
	fragNum = 0

rmsdDict = { 4:3.0, 3:2.0, 2:1.25, 0:0 }

if fragNum == 0:
	print os.path.basename(sys.argv[1]), 0 
	sys.exit()


store = []
rmsdCut = rmsdDict[ fragNum ]
with open( sys.argv[1] ) as file:
	for i in file: 
		if len(i )< 1: continue
		seq = ''.join( [ UnNatAA[res] for res in i.split()[2:] ] )
		rmsd= float( i.split()[0] )
		
		if seq not in store and rmsd <= rmsdCut:
			store.append( seq )
		#print seq, '\n'

		if len(store) == 100:
			break



print os.path.basename(sys.argv[1]), len( store ), len( store[0] ), fragNum, score( fragNum, len( store[0] ), len( store ) )

