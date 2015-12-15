# Marco Mravic 2015 DeGrado Lab
# input text file with number of subangstrom matches to each of 3 interfaces in a given pdb
# Check the
# output a histogram of each 
# output an list of pathes for the "filtered" bundle designs, to determine helix-heilx geometry (distance & angle).

import os, sys, numpy as np
from collections import defaultdict,Counter
from prody import *

limit, maxN	= 100, 180
designDict 	= defaultdict(list)

selDict = { 
			'AB': ( 'chain A resnum 14 to 25', 'chain B resnum 14 to 25' ) ,
			'BC': ( 'chain B resnum 14 to 25', 'chain C resnum 14 to 25' ) ,
			'CD': ( 'chain C resnum 14 to 25', 'chain D resnum 14 to 25' )
			}



def filterCoords():

	return 

header = '''###SUMMARY OF DESIGNABILITY SCORE AND INTERFACE DISTANCES
#path AB_subA_matches BC_subA_matches CD_subA_matches AB_CAcontact BC_CAcontact CD_CAcontact Average_AD_DIST'''
print header
with open( sys.argv[ 1 ] ) as file:
	for i in file:

		# Filter any non-matches or any interfaces with less than 100 sub-angstrom matches
		if len( i.split() )  < 4  or min( [ int(x) for x in i.split()[1:]] ) < limit:
			continue

		path = os.path.abspath( '/'.join( i.split()[0].split('/')[5:7] ) )

		if not os.path.exists( path ):
			continue
		### look at c-alpha distance matrix & record the average Ca-Ca distance between a-d-e residues for three closest distances
		pdb 	= parsePDB( path, subset = 'ca' )
		cadist 	= []
		for k,v in sorted( selDict.items() ):

			h1, h2 	= pdb.select( v[0] ), pdb.select( v[1] )
			dMat 	= buildDistMatrix( h1, h2 )
			caD 	= [ min(d) for d in dMat ]	
			
			cadist.append( round( np.mean( sorted( caD )[:3] ), 1) )
	
		# Filter any interfaces with larger a-d-e c-alpha distances of intended cluster geometry e.g. > "6,8,6" Angstroms
		if cadist[0] > 6 or cadist[1] > 7.5 or cadist[2] > 6 or min( cadist ) < 4: 
			continue

		# find average closest distance from residues  across CA
		ha, hd = pdb.select( 'chain A resnum 5 to 30' ), pdb.select( 'chain D resnum 5 to 30' )
		
		dMatSC 		= buildDistMatrix( ha, hd )
		caDsc 		= [ min(d) for d in dMatSC ]
		superCore 	= round( np.mean( caDsc ), 2 )


		
		print i.rstrip(), ' '.join( [ str(x) for x in cadist] ), superCore

