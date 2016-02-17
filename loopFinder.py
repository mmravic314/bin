#For match file from master, look up helices and find if they are connected... find if loop is connecting segment

import sys, os, subprocess as sp, numpy as np
from prody import *
from operator import itemgetter
from itertools import groupby

inPDB = parsePDB( sys.argv[1], subset = 'bb' )

# python ~/bin/loopFinder.py 56_helices.pdb 56_helices.m ~/termanal/support.default/151218_masterDB_parsedPDB/ ~/tertBuilding/CMP_bobo/56_helices/



# look into file
lineNum = 0
with open( sys.argv[2] ) as file:
	for m in file:
		
		rmsd = float(  m[:7].strip() )
		if rmsd > 0.8: 
			break

		fragMatch = [ tuple( x.strip(',').strip('()').split(',') ) for x in m.split('[')[-1].split(']')[0].split() ]
		
		gapL = int ( fragMatch[1][0] ) - int( fragMatch[0][1] ) 

		if gapL > 14 or gapL < 0:
			lineNum +=1
			continue

		#Convert match residue index ranges to array of all ( for prody selection, e.g. [38,39,40,41,59,60,61,62] )
		all_match = []
		for k in fragMatch:
			rng = np.arange( int( k[0] ), 1 + int( k[1] )  )
			all_match.extend( rng )

		resis = ' '.join( [ str( x ) for x in sorted( all_match[:] )] )


		selStr 		= 'resnum ' + resis
		selRng 		= 'resnum %d to %d' % ( all_match[0], all_match[-1] )
		look_upPath = os.path.join( sys.argv[3], os.path.basename( m.split()[1] ).split('.')[0] + '.pdb' )	
		mPdb 		= parsePDB( look_upPath, subset = 'bb' )

		r = 0
		for res in mPdb.iterResidues():
			res.setResnum( r )
			r += 1
		print gapL#, rmsd, 'match%d' % ( lineNum )# fragMatch, selRng, selStr
		#print


		wholeset 	= mPdb.select( selRng ).copy()
		ends		= mPdb.select( selStr ).copy()

		wholeset.setTitle( str(rmsd) )
		ends.setTitle( str(rmsd) )
		loopPath = os.path.join( sys.argv[4], 'loopMatch' + str(lineNum) + '.pdb' )
		endsPath = os.path.join( sys.argv[4], 'endsMatch' + str(lineNum) + '.pdb' )

		
		# find transformation matrix from original PDB (aligned to target) to native PDB
		betterMat	 = superpose( ends, inPDB )[1]
		wholeset 	 = applyTransformation( betterMat, wholeset.copy())


		writePDB( loopPath, wholeset )
		writePDB( endsPath, ends )

		## Align loop file to input PDB coords


		lineNum +=1