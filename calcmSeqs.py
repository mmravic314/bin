import sys,os
from collections import defaultdict

mDict = defaultdict(list)
for i in os.listdir( '.' ):
	 if i[-6:] == '8.mseq':
	 	
	 	intF 	= i[:2]
	 	pdb 	= i[3:19]
	 	numLines= len( open(i, 'rU').readlines() )

	 	#print intF,pdb,numLines
	 	mDict[pdb].append( ( intF, numLines ) )

for k,v in mDict.items():
	print k + '\t' + '\t'.join( [ str( x[1] ) for x in sorted(v) ] )