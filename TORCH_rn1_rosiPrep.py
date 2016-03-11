##  given input list of good models
#   unzip & move PDBs their own folder
#   Write resfile for system

# python ~/bin/TORCH_rn1_rosiPrep.py ../rnd1_models2Design.txt ~/peptideAmyloid/RND1_designTrj030515/ ~/peptideAmyloid/OsakaModels/

import sys, os, shutil, gzip
from prody import *
from collections import defaultdict

if not os.path.exists( sys.argv[2] ):
	os.mkdir( sys.argv[2] )

IntF=[ 1,0,1,0,1,0,1,0,1,0,0,1, 1 ]

print
with open( sys.argv[1] ) as file:
	for i in file:

		if i[0] == '#': continue

		oldPath = os.path.join( sys.argv[3], 'model_%s.pdb.gz' % (i.split()[0] ) )
 		
 		newDir 	= os.path.join( sys.argv[2], 'model_%s/' % (i.split()[0] ) )
 		
 		newPath = os.path.join( newDir, 'model_%s.pdb' % (i.split()[0] ) )

 		print oldPath


 		if not os.path.exists( oldPath ): 
 			print 'not found'
 			continue
 		

 		if not os.path.exists( newDir ):
 			os.mkdir( newDir )

		with gzip.open(oldPath, 'r') as f_in, open(newPath, 'w') as f_out:
			shutil.copyfileobj(f_in, f_out)

		# format new file for rosetta
		inPdb 	= parsePDB( newPath, chain='BCDEFGHIJXYZ' )
		ind 	= 1
		chResi 	= defaultdict(list)

		for r in inPdb.iterResidues():
			r.setResnum( ind )
			
			chResi[ r.getChid() ].append( ind )

			ind += 1 

		writePDB( newPath, inPdb )
		print 'wrote', newPath, '\n'

		# write resfile just once, since all are the same
		txt = 'start\n'
		for k,v in chResi.items():
			if k in ['X', 'Y', 'Z']:
				for r in v:
					txt += '%d %s ALLAAxc\n' % ( r, k )
			else:
				for r, c in zip( v, IntF ):
					if c: txt += '%d %s NATAA\n' % ( r, k )
					else: txt += '%d %s NATRO\n' % ( r, k )
		
		resF = open(  'resfile' , 'w' )
		resF.write( txt )
		resF.close()

		sys.exit()

