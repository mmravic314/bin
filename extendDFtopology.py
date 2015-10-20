## From helical pairs with one helix aligned to DF1, and the other helix is splayed relative to original 4-helix bundle
## Perform MASTER search to find close backbone matches in PDB to helical interface. Use to fuse backbone of matches, extending helices

# python ~/bin/extendDFtopology

from prody import *
import numpy as np, sys, subprocess as sp, os

# input text file containing name of PDB and residue ranges (referred to as mobile)
# input directory containing PDB's on the list (must all be in same dir)
# input path to master directory (including MASTER binary and create PDB directory)
# input path to template pdb file of DF1 dimer (referred to as target)

searchDir = os.path.join( os.path.dirname( sys.argv[4] ), 'extSearch/' )
if not os.path.exists( searchDir ):
	os.mkdir( searchDir )

# Load target PDB
target = parsePDB( sys.argv[4],  subset = 'bb' )

# Store filepath as index that looks up tuple of the 'prody' style selection strings for mobile and target backbone atoms to make in query PDB
hashMobile = []
with open( sys.argv[1] ) as file:
	for i in file:
		if i[0] == '#': continue

		basname 	= i.split()[1] 

		selStrMob 	= 'chain %s resnum %s' %  (  i.split()[2] , ' '.join( [ str(x) for x in np.arange(  int( i.split()[3].split('-')[0] ), 1 + int( i.split()[3].split('-')[-1] )  ) ]  ) )  
		selStrTar 	= 'chain %s resnum %s' %  (  i.split()[4] , ' '.join( [ str(x) for x in np.arange(  int( i.split()[5].split('-')[0] ), 1 + int( i.split()[5].split('-')[-1] )  ) ] )  )

		inPath		= os.path.join( sys.argv[2], basname )
		oPath		= os.path.join( searchDir, basname[:-4] + '-QUERY_%s.pdb' % ( i.split()[3] ) )
		mobStrct 	= parsePDB( inPath, subset = 'bb' ).select( selStrMob )
		
		# Rename mobile chain if it's chain ID is 'B'
		if mobStrct.getChids()[0] == target.select( selStrTar   ).getChids()[0]:
			mobStrct.setChids( [ 'A' for x in mobStrct.getChids() ] )

		writePDB( 'tmpT.pdb' , target.select( selStrTar  ) )
		writePDB( 'tmpM.pdb' , mobStrct )

		mobileSeg	= parsePDB( 'tmpM.pdb' )
		targetSeg 	= parsePDB( 'tmpT.pdb' )
		oStruct		= mobileSeg + targetSeg
	
		writePDB( oPath, oStruct )

		# write PDS query file
		cmd = [ os.path.join( sys.argv[3], 'createPDS' ), '--type', 'query', '--pdb', oPath ]
		print
		sp.call( cmd )
		print 

		# C2 symmetric operations on object (via DF1 constant chain superposition)
		oStructSym	= oStruct.copy()
		selStrTarSym= 'chain %s resnum %s' %  (  "A" , ' '.join( [ str(x) for x in np.arange(  int( i.split()[5].split('-')[0] ), 1 + int( i.split()[5].split('-')[-1] )  ) ] )  )
		symTargSeg 	= target.select( selStrTarSym )

		
		

		print 
		#break
