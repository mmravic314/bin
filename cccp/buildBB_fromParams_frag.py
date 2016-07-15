import sys, os, subprocess as sp, numpy as np
from collections import defaultdict
from prody import *

# python ~/bin/cccp/buildBB_fromParams_frag.py generateCrickBB tmp.par ~/bin/40_allAlaIdealHelix.pdb ~/bin/cccp/tests/trial.pdb



path2Funct 	= sys.argv[1]
paramsF 	= sys.argv[2]	
frag 		= parsePDB( sys.argv[3] ).select( 'resnum 1 2 3 4 bb' ).copy()
fragCA		= frag.select( 'ca' ).copy()

params 		= []


prefix 		= sys.argv[4]
dirPath 	= os.path.dirname( prefix )
base		= os.path.basename( prefix ).split('.pdb')[0]
if not os.path.exists( dirPath ):
	os.mkdir(dirPath)


pSet = 0
with open( paramsF ) as fin:
	for p in fin:
		if p[0] == '#': continue

		# parse parameters, format for octave function call
		par 		= p.rstrip().split(',')
		chN, chL 	= tuple( [ int( n ) for n in par[:2] ] )
		call 		= path2Funct + '(%s)' % p.rstrip()
		print p.rstrip()

		# call octave with given parameters, record output. Quit if function error found
		cmd 	= ['octave', '--eval', call]
		out 	= sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
		stdout, err 	= out.communicate()
		if len(err) > 0:
			print 'SOME ERROR FOUND... quitting\n', err
			sys.exit()


		# parse coordinates, save into chains
		coords 	= defaultdict(list)
		atm		= 0 
		ch 		= 65 		# Ascii value for chain to start as A
		for k in stdout.split('\n'):
			if len( k.rstrip() ) < 8: continue

			x, y, z = tuple( [ round( float( n ), 3) for n in k.rstrip().split() ] )
			#print atm, x,y,z, chr(ch)
			coords[ chr(ch) ].append( np.array( [x,y,z] ) )
			
			# make new chain if atom named is multiple of chain length
			atm += 1
			if atm%chL == 0:
				ch +=1

		# With coordinates, write full backbone atoms for chain
		helices = []
		for c in sorted( coords.keys() ):
			stp = 1

			# from N-terminus, align  3 residue fragment by CA superposition, 
			for i in coords[c][:-3]: 

				target 	= coords[c][stp - 1:stp + 3]
				mat		= calcTransformation( fragCA.copy(), np.array( target ) )
				obj 	= applyTransformation( mat, frag.copy() ).select( 'resnum 1' ).copy()
				obj.setResnums( [stp for x in obj.iterAtoms()] )
				if stp == 1:
					helix = obj
				else:
					helix += obj
				stp += 1

			# Add last two residues to C-terminus
			target 	= coords[c][-4:]
			mat		= calcTransformation( fragCA.copy(), np.array( target ) )
			obj 	= applyTransformation( mat, frag.copy() ).select( 'resnum 2 3 4' ).copy()
			for r in obj.iterResidues():
				r.setResnum( stp )
				stp += 1

			helix += obj

			helix.setChids( [ c for r in helix.iterAtoms() ] )
			if c == 'A':
				helices = helix
			else:
				helices += helix

		pathout = os.path.join( dirPath, '%s_%d.pdb' % ( base, pSet )  )
		helices.setTitle( p.rstrip() )
		writePDB( pathout, helices )

		## Add stuff here to spit out NAMD files to do minimization

		##

		pSet += 1




			