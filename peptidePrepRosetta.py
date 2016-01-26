# look through a directory tree and find PDB files, rewrite copy fit for Rosetta
# Write a res file for each and a coordinate constriant file for Calpha to fix during minimization

# input path to directory tree
# python ~/bin/peptidePrepRosetta.py ~/peptideAmyloid/searchFrags/2helix/ ~/peptideAmyloid/rosettaFixBB/ > inputKey.txt
interface 	= [10,12,14,16,18,21,23]
renumNot	= []    # for non-interface residues in the renumbered scheme

import os,sys
from prody import *
from PDButil import UnNatAA

index 		= 1
storeNames 	= []

for dirN in sorted( os.listdir( sys.argv[1] ) ):
	if dirN[0] == ".": continue
	dirPath = os.path.join( sys.argv[1], dirN )
	for f in sorted( os.listdir( dirPath ) ):
		if f[-4:] != '.pdb': continue

		oldPath = os.path.join( dirPath, f )
		treePath = os.path.join( sys.argv[2], 'input%d' % (index) )
		newPath = os.path.join( treePath, 'input%d.pdb' % (index) )

		storeNames = ( os.path.basename( newPath ), oldPath ) 

		# Enter file and rewrite for Rosetta format
		# Remember few things 
		remark 	= open( oldPath, 'rU' ).readlines()[0]
		resi	= []
		ch 		= 'X'
		for s in remark.split()[3:]:
			if s == 'ChY':
				ch = 'Y'
				continue
			resi.append( ( ch, s ) )

		inPDB 			= parsePDB( oldPath )
		# find all the amyloid chains within 25 angstroms from any helix atom ... to reduce simulation box
		helixSet 		= inPDB.select( 'chain X Y' ).copy()
		closeChs		= inPDB.select( 'exwithin 7 of helixSet', helixSet = helixSet ).getChids()
		closeChs		= [ c for c in list( set( closeChs )) if c not in ['X', 'Y'] ] 
		inPDB			= inPDB.select( 'chain X Y %s' % ( ' '.join( closeChs ) ) ).copy()

		# renumber residues
		ca 				= []
		count 			= 1
		startAmyloid 	= []
		chains 			= []
		res_str			= ''
		for r in inPDB.iterResidues():

			res = ( r.getChid(), str( r.getResnum() ) )							# find the residues indicated in the REMARK line to be Ca fixed
			if res in resi:
				ca.append( ( count, list( r.select('ca').getCoords()[0] ) ) )

				# OSAKA Interfacial residues can repack
			if r.getResnum() in interface and r.getChid() not in ['X', 'Y']:
				res_str += '%d %s PIKAA %s\n' % ( count, r.getChid(), UnNatAA[r.getResname()] )

				# OSAKA non-interfacial residues just stick to natural rotamer
			if r.getChid() not in ['X', 'Y'] and r.getResnum() not in interface:
				res_str += '%d %s NATRO\n' % ( count, r.getChid() )

			r.setResnum( count )  # set residue number 

			count  += 1

		# make residue list for amyloids
		resSet = [ (x.getResnum(), x.getCoords()[0] ) for x in inPDB.select( 'calpha chain %s' % ( ' '.join( closeChs ) ) ).copy().iterResidues() ] 

		# write new tree for this input PDB
		if not os.path.exists( treePath ):
			os.mkdir( treePath )
		writePDB( newPath, inPDB )

		# find e-g pairs on surface-face interface. call 1 E and other K (NAH maybe do this later)
		
		# Write constraint file
		cstStr = ''
		# Constraint alpha helices pretty loosly
		for r in ca:
			if ca[0] == r: continue # skip the first residue, since it is the reference coordinate for constraints
			cstStr += 'CoordinateConstraint CA %d CA %d %f %f %f HARMONIC 0.0 0.3\n' % ( r[0], ca[0][0], r[1][0], r[1][1], r[1][2] )
		# tightly constraint beta sheet Ca's
		for r in resSet:
			if resSet[0] == r: 
				cstStr += 'CoordinateConstraint CA %d CA %d %f %f %f HARMONIC 0.0 0.15\n' % ( r[0], resSet[1][0], r[1][0], r[1][1], r[1][2]  )
				continue
			cstStr += 'CoordinateConstraint CA %d CA %d %f %f %f HARMONIC 0.0 0.15\n' % ( r[0], resSet[0][0], r[1][0], r[1][1], r[1][2]  )

		cstPath = open( os.path.join( treePath, 'input%d.cst' % (index) ), 'w' )
		cstPath.write( cstStr )
		cstPath.close()

		# Write res file
		resStr = 'start\n* X NOTAA WCP\n* Y NOTAA WCP\n' + res_str
		resPath = open( os.path.join( treePath, 'input%d.resfile' % (index) ), 'w' )
		resPath.write( resStr )
		resPath.close()


		print storeNames[0], storeNames[1]
		index 	+= 1
