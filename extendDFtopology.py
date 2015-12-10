## From helical pairs with one helix aligned to DF1, and the other helix is splayed relative to original 4-helix bundle
## (DOESNT DO THIS) Perform MASTER search to find close backbone matches in PDB to helical interface. Use to fuse backbone of matches, extending helices

# python ~/bin/extendDFtopology.py ~/splayBundle/c6-splayfit-queries.txt ~/splayBundle/pairs_Cluster-006 ~/bin/  ~/splayBundle/df1L13G-DIMERCenteredX.pdb  ~/bin/protein-comparison-tool_20130326/runCE.sh

from prody import *
import numpy as np, sys, subprocess as sp, os

# input text file containing name of PDB and residue ranges (referred to as mobile)
# input directory containing PDB's on the list (must all be in same dir)
# input path to master directory (including MASTER binary and create PDB directory) OBSOLETE
# input path to template pdb file of DF1 dimer (referred to as target)
# input path to octave/mathlab function for fitting coiled-coil backbones

searchDir = os.path.join( os.path.dirname( sys.argv[4] ), 'tm%s_extSearch/' % ( sys.argv[2][-1:].strip() ) )
if not os.path.exists( searchDir ):
	os.mkdir( searchDir )

# Load target PDB
target 			= parsePDB( sys.argv[4],  subset = 'bb' )


# Store filepath as index that looks up tuple of the 'prody' style selection strings for mobile and target backbone atoms to make in query PDB
hashMobile = []



with open( sys.argv[1] ) as file:
	for i in file:

		if i[0] == '#': continue

		basname 	= i.split()[1] 
		print basname

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
	
		#writePDB( oPath, oStruct )

		# write PDS query file
		#cmd = [ os.path.join( sys.argv[3], 'createPDS' ), '--type', 'query', '--pdb', oPath ]
		#print cmd
		#sp.call( cmd[0] )
		print 

		### C2 symmetric operations on object (via DF1 constant chain superposition)
		# replicate atom group, Change chain IDs (from AB to CD) set outfile path
		oStructSym	= oStruct.copy()
		oStructSym.setChids
		oPathSym 	= os.path.join( searchDir, basname[:-4] + '-SYM_%s.pdb' % ( i.split()[3] ) )

		# load atom group containing the C2 symmetry mate of DF1 helix (resi 27-48) segment 
		selStrTarSym	= 'chain %s resnum %s' %  (  "A" , ' '.join( [ str(x) for x in np.arange(  int( i.split()[5].split('-')[0] ), 1 + int( i.split()[5].split('-')[-1] )  ) ] )  )
		symTargSeg 		= target.select( selStrTarSym )
		#symPath 		= os.path.join( os.path.dirname( sys.argv[4] ), 'tmpTSym.pdb' )
		#writePDB( symPath , symTargSeg )

		# Calculate transformation matrix for C2 symmetry mate, minimizing RMSD
		# then apply matrix to splay(mobile) + target
		aligned		= superpose(  targetSeg, target.select( selStrTarSym ) )
		# print calcRMSD( target.select( selStrTarSym ).getCoords() , target=aligned[0].getCoords() )
		oStructSym 	= applyTransformation( aligned[1], oStructSym )

		for x in oStructSym:
			if x.getChid() == 'A':
				x.setChid( 'D' )
			else:
				x.setChid( 'C' )


		
		# Manipulate/join prody atom groups to have C2 mate of splay(mobile) + target object; 
		# Exit each segment to have unqique chain ID for coiled-coil fitting; save file
		writePDB( 'tmpOsym.pdb' , oStructSym ) 
		oStructSym	= parsePDB( 'tmpOsym.pdb' )
		c2_coil 	= oStruct + oStructSym


		writePDB( oPathSym , c2_coil )
		print 'wrote',  oPathSym
		

		## perform fits of coiled coil backbone with octave (DOESN't WORK!)
		#outTxtPath 	= oPathSym[:-4]
		#args 		= ','.join( [ '\'%s\'' % (oPathSym), '4', '\'GENERAL\'', '1', '\'%s\'' % ('tmp'), '[]', '[]', '[]', '[]' ] )
		#fitCMD		= [ 'octave', '--silent', '--eval', '\"fcrick(%s)\"' % ( args ) ]
		#print fitCMD 
		#print
		#sp.call( fitCMD )

		

		print 

try:
	os.remove( 'tmpOsym.pdb' )
	os.remove( 'tmpM.pdb' )
	os.remove( 'tmpT.pdb' )
except OSError:
	pass
#os.remove( oPath )

## perform fits
#fcrick("df1L13G-DIMERcccp.pdb",4,'GENERAL',1,0.5,[],[],[],[])
###
##
# from #~/bin/cccp/:
# octave --silent --eval "fcrick('/Users/mmravic/splayBundle/tm1_extSearch/1m0l-008_007-0040_0062_B_0008_0030_B_aligned-SYM_34-46.pdb',4,'GENERAL',1,'tmp',[],[],[],[])"


