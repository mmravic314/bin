import sys,os
from PDButil import natAA

# Turn each match seq file into a multiple sequence alignment
# Hash it by alignment, print as multiple sequence alignment Fasta file together (and separate)

# 
fragList = ['BC', 'CD']
for f in os.listdir( sys.argv[1] ):
	if f[:2] == 'AB':
		pathList = [ f, fragList[0] + f[2:], fragList[1] + f[2:]  ]
		pathOutL = [ k[:-4] + '.mseq' for k in pathList ]
		for p in pathList:
			inF = open( p, 'rU' )

		print
	print 
	sys.exit()  


