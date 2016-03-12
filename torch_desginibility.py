# make term-like fragments of model, search database for structural matches
# 

from prody import *
import sys, os, numpy as np

def selStr( lst ):			# return string from list
	return ' '.join( [ str(x) for x in lst] )


inPDB 	= parsePDB( sys.argv[1] )
intF 	= [10,12,14,16,18,21,23]


intFset	= inPDB.select( '(not element H) and chain A B C D E F G H I J K resnum ' + selStr( intF ) )

## Dump model if any backbone clashes, 



# Extract each helical interface, search, save
helices = [ inPDB.select( 'chain %s' % (c) ) for c in ['X', 'Y', 'Z'] ]

# make helix-sheet fragments, search
hSels	= [ np.arange( 7, 15 ), np.arange( 12, 20 ), np.arange( 17, 25 ) ]

hFrags 	= [ helices[1].select( 'resnum ' + selStr(h) ).copy() for h in hSels  ] 


# Make search fragments
for fr in hFrags:
	print fr
	# find three CB's closest to a interface residue 
	seeds = []

	for a in fr.select( 'name CB' ):

		for k in sorted( iterNeighbors( a, 7, intFset), key=lambda x: calcDistance( x[0], x[1] ) ) : 
			
			print a.getResnum(), k[1].getResnum(), calcDistance( a, k[1] )
			d = calcDistance( a, k[1] )
			if d < 7:
				seeds.append(  (a, d)  )
			break
	
	print 	
	if len(seeds) >= 3:
		sds = [ a[0] for a in sorted( seeds[:3], key=lambda x: x[1] )  ]
		
		for j in sds:	print j.getResnum()

		break



