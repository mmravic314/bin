## Marco Mravic UCSF Biophysics DeGrado Lab August 2015 ##
# From an input directory of match files from MASTER, make pariwise RMSD matrix
# Will likely be sparse, so initialize a cap at 1. All other values will be scalred from 0 - RMSD/3.5 (max distance is 1)
# ouput a pickle file hash indexed by pairs of PDB name with user given pathname (bulky but user friendly)

# python 

import sys,os, cPickle as pic, numpy as np


pdbList = []
rmsdD 	= {}

for f in os.listdir( sys.argv[1] ):
	
	if os.path.splitext(f)[-1] != '.m':
		continue
	pdb = os.path.splitext(f)[0]
	pdbList.append( pdb )
	print 'Entering', f
	with open( os.path.abspath( os.path.join( sys.argv[1], f)  ) ) as file:
		for i in file:
			j 			= i.split()
			bbRMSD 		= round( float( j[0] ) / 4.0, 4 ) 
			matePdb 	= os.path.splitext( os.path.basename( j[1] ) )[0]
			pairID 		= tuple( sorted( [ matePdb, pdb ] ) )

			# This should take the closer of two matches if two parts of one fragment matches to the query fragment
			try:
				rmsdD[ pairID ]
			except KeyError: 
				rmsdD[ pairID ] = bbRMSD


# Fill out distance matrix and hash the pdb pair by indices (indices as tuple)

lookupDict 	= {}											
dMatrix 	= np.ones ( ( len( pdbList), len(pdbList) ) )

for i in enumerate( sorted( pdbList ) )  :
	for j in enumerate( sorted( pdbList ) ) :
			pairID 		= tuple( sorted( [ i[-1], j[-1] ] ) )

			try:
				dMatrix[ i[0], j[0] ] = rmsdD[ pairID ]
				dMatrix[ j[0], i[0] ] = rmsdD[ pairID ]
				lookupDict[ str( (i[0], j[0]) ) ] = pairID
				lookupDict[ str( (j[0], i[0]) ) ] = pairID

			except KeyError:
				rmsdD[ pairID ] = 1.0
				

pic.dump( dMatrix, 		open( sys.argv[2], 'wb' ) )
pic.dump( lookupDict, 	open( sys.argv[2] + '__lookUpHash' , 'wb' ) )

