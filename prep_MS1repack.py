# with an input ca-only coiled coil, generate input files for rosetta 
# Example:
# python ~/bin/prep_MS1repack.py ~/tertBuilding/MS1/ ~/rosetta/ ~/tertBuilding/MS1/Xms1_Spread.pdb  ~/tertBuilding/MS1/Xms1_alignMEM.pdb helix_Relax-MS1mini.xml 



import sys, os, subprocess as sp, time, numpy as np
from prody import *
from PDButil import NatAA_3
from collections import defaultdict



### I/O

workingDir	= sys.argv[1]


rosiBase 	= sys.argv[2]
rosiDB 		= os.path.join( rosiBase, 'database/')
rosiScrps	= os.path.join( rosiBase, 'source/bin/rosetta_scripts.linuxgccrelease' )
rosiSpanGen = os.path.join( rosiBase, 'source/bin/spanfile_from_pdb.linuxgccrelease')

spreadPDB 	= parsePDB( sys.argv[3] )
alignPDB 	= parsePDB( sys.argv[4] ).select(' resnum 1 to 28')

designXML 	= sys.argv[5]

FNULL 		= open(os.devnull, 'w')		# Place to dump stdout

###

## Make span file for the spread out case
os.chdir( workingDir )
spanPath = sys.argv[3][:-4] + '.span'
if not os.path.exists( spanPath ):

		cmdSpan = [ rosiSpanGen, 
		'-database', rosiDB, 
		'-in:file:s', sys.argv[3]
		]	
		sp.call( cmdSpan, stdout=FNULL, stderr=sp.STDOUT  )
##


### Prepare the resfile for a position starting at residue 9, 12, and 13
## also prepare the two spread PDB files in these formats
#     AQLLIA VLLLIAT NLILLIA VARLRYL VG
#  K AQWLLIA VLLLIAT NLILLIA VARLRYL VG	# 9, 31 (C-term die)
# KK AQWLLIA VLLLIAT NLILLIA VARLRYL V 	# 12, 31 (N-term die)
#RKK AQWLLIA VLLLIAT NLILLIA VARLRYL  	# 13, 31 (N-term die) ignore this

hepDir = {}
hepDir[9] 	= ('KLLIAVLLLIATNLILLIAVAFLRYLVG', [1,2,3,32])
hepDir[8] 	= ('KLLIAVLLLIATNLILLIAVAFLRYLVG', [1,2,31,32])
#hepDir[12] 	= ('KAQWLLIAVLLLIATNLILLIAVARLRYLVG', 32)
hepDir[6] 	= ('KLLIAVLLLIATNLILLIAVAFLRYLVG', [29,30,31,32])
hepDir[5] 	= ('KLIAVLLLIATNLILLIAVAFLRYLVGG', [29,30,31,32])


for k, v in hepDir.items():
	
	# resile writing for the two unique sequences
	if k not in [5,6]: continue

	res_file 	= os.path.join( workingDir, 'resfile%d' % k ) 
	txt 		= 'start\n'

	rnum = 1
	for res in v[0]:
		txt +=  '%d A PIKAA %s\n%d B PIKAA %s\n%d C PIKAA %s\n' % ( rnum, res, rnum, res, rnum, res )
		rnum += 1
	oF = open( res_file, 'w' )
	oF.write( txt )
	oF.close()

	# spread pdb rewrite, including CA constraints for repack
	pdbFile = os.path.join( workingDir, os.path.basename( sys.argv[3] )[:-4] + '_%d.pdb' % k )
	
	subset 	= spreadPDB.copy()
	# cst_txt += 'CoordinateConstraint %s %d N 1 %f %f %f HARMONIC 0.0 0.67\n' % ( name, resnum, x,y,z  )


	txt = ''
	overallResi = 1
	for chain in subset.iterChains():
			r 		= 1
			for resi, name in zip( chain.iterResidues(), v[0] ):

				resi.setResnum( r )
				resi.setResname( NatAA_3[name] )
				# Skip first residue constraint
				if  overallResi == 1: 
					overallResi += 1
					r += 1
					continue

				# not need for restraints
				#x,y,z, 	= tuple( [ '%0.3f' % round( x, 3 ) for x in resi.select( 'ca' ).getCoords()[0]  ] ) 
				#txt 	+= 'CoordinateConstraint CA %d CA 1 %s %s %s HARMONIC 0.0 0.4\n' % ( overallResi, x,y,z )

				r += 1
				overallResi += 1


	#cst_path =  os.path.join( workingDir, 'bbCA_%d.cst' % k )
	#oF = open( cst_path, 'w' )
	#oF.write( txt )
	#oF.close()

	#writePDB( pdbFile, subset )


	## Do design here, if it's the first tim doing it
	spreadDir = os.path.join( workingDir, 'spread_outputs/' )
	if not os.path.exists( spreadDir ):
		os.mkdir( spreadDir ) 

		# Skip rosetta designs if you dont need to
	continue
	n = '10'
	cmdSp = [ rosiScrps, 
'-parser:protocol', designXML, 					# Path to Rosetta script (see above)
'-in:file:s', sys.argv[3],						# Input PDB structure
'-nstruct', n, 								# Generate 1000 models
'-mp:setup:spanfiles', spanPath,				# Input spanfile
'-mp:scoring:hbond', 						# Turn on membrane hydrogen bonding
'-relax:jump_move', 'true', 				# Allow jumps to move during relax
'-out:prefix',spreadDir,
'-out:suffix','_%d' % k,
'-packing:resfile', res_file,
'-out:overwrite',
'-packing:pack_missing_sidechains', '0' ]

	start = time.time()
	#sp.call( cmdSp ) #, stdout=FNULL, stderr=sp.STDOUT  )
	print 'Elapsed:', time.time() - start
	print cmdSp

#sys.exit()

###
sc = defaultdict(list)
for f in os.listdir( spreadDir ):

	label = f.split('_')[2]	
	with open( os.path.join( spreadDir, f ) ) as fin:
			for i in fin:

				if i[:5] == 'score':
					sc[label].append( float( i.split()[1] ) ) 


for k,v in sc.items():
	print k, min( v )



### Prep all regular files for design
for f in [ x for x in os.listdir( workingDir ) if '_trial' in x and x[-3:] == 'pdb' ]:
	

	# create new working directory for this designs and it's outputs
	innerDir 	= os.path.join( workingDir, f[:-4] )
	print 'entering...', innerDir
	if not os.path.exists( innerDir ):
		os.mkdir( innerDir )

	
	pdbPath 	= os.path.join( workingDir, f )
	pdbObj 		= parsePDB( pdbPath )
	# make versions of backbone for heptad with 'a' starting on residue 9 or 12
	# also write backbone constraint file
	txt = ''
	for k, v in hepDir.items():
	
		pdbFile = os.path.join( innerDir, f[:-4] + '_%d.pdb' % k )
		subset 	= pdbObj.select( 'not resnum %s' % ' '.join( [ str(x) for x in v[1]] ) ).copy()
		superpose( subset, alignPDB )
		txt = ''
		for chain in subset.iterChains():
			r 		= 1
			for resi, name in zip( chain.iterResidues(), v[0] ):
				resi.setResnum( r )
				resi.setResname( NatAA_3[name] )
				#print r, name, chain, v[1]
				# Skip first residue constraint

				r += 1
				overallResi += 1
		#cst_Filename =  os.path.join( innerDir, 'bbCA_%d.cst' % k )
		#oF = open( cst_Filename, 'w' )
		#oF.write( txt )
		#oF.close()
		'''
		superpose( subset, alignPDB )
		'''


		
		writePDB( pdbFile, subset )
	#sys.exit()



