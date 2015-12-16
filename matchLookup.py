# given a match file output from master, 
# Download the pdb, if necessary and select residues within 5A of sidechains
# look for segments of helical residues nearby to matches in original PDB 
# (selected residue within match based on query) 

import sys, os, subprocess as sp, numpy as np
from prody import *
from operator import itemgetter
from itertools import groupby

## Example


# look through match file, hash residue numbers by pbd name + chain (allow matches to more than 1 pdb)
# *note: only search surface residues of input pdb (hard coded?)

### On any chain of input structures, find interaction partners in normal pdb to THESE residues
interface = [10,12,14,16,18,21,23]
# convert an input pdb file to index of eligible residues in fragments (e.g. 1:[0,2,4,7,8]; 2:[0,2,4,7,8]  )
# Which can be converted to selection string of a matched pdb
def queryInterfaceConvert( pdbPath ):
	fragmentDict = {}

	pdb = parsePDB( pdbPath )

	frag 			= 0
	for c in pdb.iterChains():
		resiInterface 	= []
		index 			= 0

		for r in c.iterResidues():
			
			if r.getResnum() in interface:
				resiInterface.append( index )
			index += 1

		fragmentDict[ frag ] = resiInterface
		frag += 1



	return fragmentDict

def continutityCheck( array ):
	ind = 1
	while ind < len( array ):
			if array[ind-1] + 1 != array[ind]:
				return False
			ind +=1
	return True


fragmentDict = queryInterfaceConvert( sys.argv[1] )

#print fragmentDict.items()

queryName = os.path.basename( sys.argv[1] ).split('.')[0]
print "Entering matches for", queryName, '\n'
mNum = 0
with open( sys.argv[2] ) as file:
	for i in file:
		mNum += 1

		# for each match index, convert to fragment 
		fragMatch = [ tuple( x.strip(',').strip('()').split(',') ) for x in i.split('[')[-1].split(']')[0].split() ]
		
		#Convert match residue index ranges to array of all ( for prody selection, e.g. [38,39,40,41,59,60,61,62] )
		all_match = []
		for k in fragMatch:
			rng = np.arange( int( k[0] ), 1 + int( k[1] )  )
			all_match.extend( rng )

		all_match = ' '.join( [ str( x ) for x in all_match[:]] )

		chID 	= os.path.basename( i.split()[1] ).split('_')[1][0]
		# Prepare to store residueList
		resiInterface 	= []
		disclude		= []
		discludeMini 	= []	# don't exlcude residues +/- 2, in case at termini
		

		# Select matches' residue numbers corresponding to interface residues on query fragment
		fIndex = 0
		for f in fragMatch:
			resIndex 	= 0
			# begin to disclude +/- 2 residues of each fragment (check if these overlap with other fragments in query)  
			disclude.extend( [ int(f[0]) - 1 , int(f[0]) - 2 , 1 + int( f[-1] ), 2 + int( f[-1] ) ] ) 
			for f in np.arange( int(f[0]), 1+ int( f[-1] ) ):

				if resIndex in fragmentDict[fIndex]:
					resiInterface.append( f )
				
				disclude.append( f )
				discludeMini.append( f )

				resIndex += 1

			fIndex += 1

		
		print os.path.basename( i.split()[1] ), i.split()[0], 

		# Download whole pdb if necesasary to self-named directory in working dir
		pdbPath = os.path.join( queryName, os.path.basename( i.split()[1] ).split('.')[0] + '.pdb' )
		pdbM 	= os.path.basename( i.split()[1] ).split('_')[0].upper()
		if not os.path.exists( pdbPath ):

			p 		=  os.path.basename( i.split()[1] ).split('_')[0].upper()
			# download whole pdb, take only specified chain, renumber residues?
			sp.call( [ 'wget', 'www.rcsb.org/pdb/files/%s.pdb' % (p) , '-P', queryName ])
			pact 	=  os.path.join( queryName, '%s.pdb' % (p) )
			

			writePDB( pdbPath, parsePDB( pact, chain = chID ).select( 'protein' ) )

		# look for residues within chain that are within 4 angstroms of side chain atoms
		# disclude non-interfacial and +/- 2 from ends of fragment
		selStr 		= 'resnum %s not bb' % ' '.join( [ str( x ) for x in resiInterface ] ) 

		subset 		= parsePDB( pdbPath )

		#Renumber residues in this single chain protein 
		r = 1
		for res in subset.iterResidues():
			res.setResnum( r )
			r += 1

		pdbInterface= subset.select( selStr )

		#print resiInterface, selStr


		# Search pdb for neighbors by each atom in interface set
		pdbCont 	= Contacts( subset )
		contacts 	= []
		for a in pdbInterface.iterAtoms():
			for atom in [ atm for atm in pdbCont.select( 4.2, a.getCoords() ).iterAtoms() if atm.getResnum() not in disclude ]: 
				if atom.getResnum() not in contacts:
					
					contacts.append( atom.getResnum() )
					# try to add +/- 2 for each contacted residue (check helicity later)
					more = [ atom.getResnum() -2, atom.getResnum() -1, atom.getResnum() +1, atom.getResnum() +2 ]
					for k in more:
						if subset[chID, k]:
							if k not in contacts:
								contacts.append( k )

		# print match seaction pdb... for debugging reasons
		#resiInterface.extend( contacts )
		#selStrOut 	= 'resnum %s' % ' '.join( [ str( x ) for x in resiInterface ] ) 
		#writePDB( pdbPath[:-4]  + '_m.pdb', subset.select( selStrOut ) )
		#sys.exit()

		# for residues within 4.2 angstroms of target interface, find continuous helical fragments at least 4 residues
		helix = []
		for c in contacts:
			resi 	= subset[chID, c]
			try:
				phi, psi = calcPhi(resi), calcPsi(resi)
			except ValueError: ## Skip initial/final residues with missing dihedral
				continue
			if -130 < phi < -20 and -90 < psi < 30:	
				helix.append( c )

		# evaluate continuity of helix, first by checking if there were any found 
		if len( helix ) < 4:
			print 'zero or <4 helical residues; SKIP\n'
			continue

		helix = sorted( helix )
		
		# random pertrbations to helix array.... to check that this break-up is working

		helix = [323,350,351,352,353,354] + helix
		# break into continuous groups of residues, get all the different helical segments >4 long
		helixGroups = []
		for k, g in groupby(enumerate(helix), lambda (i,x):i-x):
			catch =  map(itemgetter(1), g)
			if len( catch ) >= 4:
				helixGroups.append( catch ) 

		# Quit this match if no valid groups of helical residues
		if len(helixGroups) < 1:
			print "no valid helical segments found\n"
			continue

		# try to entend each helix further if possible.
		## Then save different copies of match and extended helical regions
		exthelixGroups= []
		for h in helixGroups:
			tmp = h[:]
			# Work from the N-terminal residue of the helix ad
			nStr = h[0] - 1
			while nStr > subset[0].getResnum():
				# check helical segments for each group
				resi = subset[ chID, nStr ] 
				if resi == None: break
				try:
					
					phi, psi = calcPhi(resi), calcPsi(resi)
					if -90 < phi < -35 and -70 < psi < 0:		# Strict helical 
						tmp.append( nStr )

					else:		# Stop extending helix if residue not helical
						break

				except ValueError, TypeError: ## Skip if residue not there or phi/psi can't be calculated
					break

				nStr -= 1

			# back
			cStr = h[0] + 1
			while nStr < subset[-1].getResnum():
				# check helical segments for each group
				resi 	= subset[ chID, cStr ]
				if resi == None: break
				try:

					phi, psi = calcPhi(resi), calcPsi(resi)
					if -90 < phi < -35 and -70 < psi < 0:		# Strict helical 
						tmp.append( cStr )

					else:		# Stop extending helix if residue not helical
						break

				except ValueError, TypeError: ## Skip if residue not there or phi/psi can't be calculated
					break

				cStr += 1

			# store the extended residue set
			exthelixGroups.append( sorted( tmp ) )

		print "FOUND!"
		#Save files of the helices (both macthed and extended) in a directory in the working directory
		outDir = queryName + '_helixMatches/'
		if not os.path.exists( outDir ):
			os.mkdir( outDir )


		index = 0
		for h,ext in zip( helixGroups, exthelixGroups ):

			hPath 		= os.path.join( outDir,  'match%d_%d_%s.pdb' 		% ( mNum, index, pdbM) )
			extPath 	= os.path.join( outDir,  'match%d_%d_%s_ext.pdb' 	% ( mNum, index, pdbM) )
			hStr 		= ' '.join( [ str( n ) for n in h ] )
			extStr 		= ' '.join( [ str( n ) for n in ext ] )


			print hStr

			hsel 		= subset.select( 'resnum %s %s' % ( hStr, all_match  ) )
			writePDB( hPath, hsel )

			if hStr != extStr:
				extsel 		= subset.select( 'resnum %s %s' % ( extStr, all_match  ) )
				writePDB( extPath, hsel )
				
			index += 1	

		print
		




