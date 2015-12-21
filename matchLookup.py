# given a match file output from master, 
# Download the pdb, if necessary and select residues within 5A of sidechains
# look for segments of helical residues nearby to matches in original PDB 
# (selected residue within match based on query) 

import sys, os, subprocess as sp, numpy as np
from prody import *
from operator import itemgetter
from itertools import groupby

## Example
# > python ~/bin/matchLookup.py n9_2strand.pdb n9_2strand.pds.m  n9_2strand.pds.struct
## path2originalPDB(assumed as path to .pds version as well), path2Matchfile path2StructuresDir path2OriginalPDbs(outputs of createPDS) path2masterDir
## NEWEST > python ~/bin/matchLookup.py n9_2strand.pdb n9_2strand.pds.m  n9_2strand.pds.struct ~/peptideAmyloid/masterSubDir/ ~/termanal/



# look through match file, hash residue numbers by pbd name + chain (allow matches to more than 1 pdb)
# *note: only search surface residues of input pdb (hard coded?)

### On any chain of input structures, find interaction partners in normal pdb to THESE residues
interface = [10,12,14,16,18,21,23]
# convert an input pdb file to index of eligible residues in fragments (e.g. 1:[0,2,4,7,8]; 2:[0,2,4,7,8]  )
# Which can be converted to selection string of a matched pdb
def queryInterfaceConvert( pdb ):
	
	fragmentDict = {}

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

pdbIn 		 = parsePDB( sys.argv[1] )
fragmentDict = queryInterfaceConvert( pdbIn )
pdsPath 	 = sys.argv[1][:-4] + '.pds'

#check = parsePDB( sys.argv[5] )

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

		all_match = ' '.join( [ str( x ) for x in sorted( all_match[:] )] )
		#print i
		#print all_match

		#chID 	= os.path.basename( i.split()[1] ).split('_')[1][0]
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

		# Look for 
		pdbPath = os.path.join( sys.argv[4], os.path.basename( i.split()[1] ).split('.')[0] + '.pdb' )
		pdbM 	= os.path.basename( i.split()[1] ).split('_')[0].upper()
		print pdbPath

		#if not os.path.exists( pdbPath ):

		#	p 		=  os.path.basename( i.split()[1] ).split('_')[0].upper()
			
			# extract createPDS parsed chain
		#	pact 	=  os.path.join( sys.argv[4], '%s_.pdb' % (p) )

			#NAH# download whole pdb, take only specified chain, renumber residues?
			#sp.call( [ 'wget', 'www.rcsb.org/pdb/files/%s.pdb' % (p) , '-P', queryName ])
			#pact 	=  os.path.join( queryName, '%s.pdb' % (p) )
			#writePDB( pdbPath, parsePDB( pact).select( 'protein' ) )

		# look for residues within chain that are within 4 angstroms of side chain atoms
		# disclude non-interfacial and +/- 2 from ends of fragment
		selStr 		= 'resnum %s not bb' % ' '.join( [ str( x ) for x in resiInterface ] ) 

		subset 		= parsePDB( pdbPath )

		#Renumber residues 
		r = 0
		for res in subset.iterResidues():
			res.setResnum( r )
			r += 1

		pdbInterface= subset.select( selStr )

		# Search pdb for neighbors by each atom in interface set
		pdbCont 	= Contacts( subset )
		contacts 	= []
		for a in pdbInterface.iterAtoms():
			for atom in [ atm for atm in pdbCont.select( 4.2, a.getCoords() ).iterAtoms() if atm.getResnum() not in disclude ]: 
				if atom.getResnum() not in contacts:
					
					contacts.append( atom.getResnum() )
					chID = subset.select( 'resnum %d' % ( atom.getResnum() )   ).getChids()[0]
					# try to add +/- 2 for each contacted residue (check helicity later)
					more = [ atom.getResnum() -2, atom.getResnum() -1, atom.getResnum() +1, atom.getResnum() +2 ]
					for k in more:

						if subset[chID, k] and k >= 0:
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
			chID 	= subset.select( 'resnum %d' % (c) ).getChids()[0]
			resi 	= subset[chID, c]
			
			try:
				phi, psi = calcPhi(resi), calcPsi(resi)
			except ValueError, TypeError: ## Skip initial/final residues with missing dihedral
				continue
			if -130 < phi < -20 and -90 < psi < 30:	
				helix.append( c )

		# evaluate continuity of helix, first by checking if there were any found 
		if len( helix ) < 4:
			print 'zero or <4 helical residues; SKIP\n'
			continue

		helix = sorted( helix )
		
		# random pertrbations to helix array.... to check that this break-up is working

		#helix = [323,350,351,352,353,354] + helix
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
				chID 	= subset.select( 'resnum %d' % (nStr) ).getChids()[0]
				resi 	= subset[chID, nStr] 
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
				chID 	= subset.select( 'resnum %d' % (cStr) ).getChids()[0]
				resi 	= subset[chID, cStr] 
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

		if not os.path.exists( 'tmp' ):
			os.mkdir( 'tmp' )

		############### obsolete section to align expanded match (+helix) with the original template
		# Before saving helical segments, find transformation matrix for best RMSD for aligning match and original pdb segment 
		# Apply to the entire section 
		#match 		= subset.select( 'bb resnum %s' % ( all_match ) ).copy()				# matched regions in original PDB
		#indexMatch 	= str( mNum )
		#if len( indexMatch ) == 1:
		#	indexMatch = '00' + indexMatch
		#elif len( indexMatch ) == 2:
		#	indexMatch = '0' + indexMatch

		#mPath 		= os.path.join( sys.argv[3], 'match%s.pdb' % indexMatch  ) 		# aligned match
		#alignedM	= parsePDB( mPath , subset = 'bb')
		
		#print calcRMSD( match, pdbIn.select( 'bb' ) )
		#transMat 	= calcTransformation( match, alignedM )
		#superpose( match, alignedM.copy() )	
		#print calcRMSD( match, pdbIn.select( 'bb' ) )

		#quitFlag = False
		# to use master to find a match from the objects ( found match as target and full unit  )
		#if calcRMSD( match, pdbIn.select( 'bb' ) )	> 1.2:
		#	quitFlag = True
		#print transMat
		##################################

########## Actually save match files, aligned to original query match site
		index = 0
		for h,ext in zip( helixGroups, exthelixGroups ):
			#transMat 	= calcTransformation( match, pdbIn.select( 'ca' ))	

			hPath 		= os.path.join( outDir,  'match%d_%d_%s.pdb' 		% ( mNum, index, pdbM[:-4]) )
			extPath 	= os.path.join( outDir,  'match%d_%d_%s_ext.pdb' 	% ( mNum, index, pdbM[:-4]) )
			hStr 		= ' '.join( [ str( n ) for n in h ] )
			extStr 		= ' '.join( [ str( n ) for n in ext ] )

			print hPath
			

			## re-find match,except with helix this time. Save as tmp files. move files into output dir and rename 
			if hStr != extStr:
				extsel 		= subset.copy().select( 'resnum %s %s' % ( extStr, all_match  ) )
				writePDB( extPath, extsel )

				# create target PDS of this match in tmp directory
				opdsPath 	= 'tmp/matchEXT.pds'
				cpdsCMD		= [ os.path.join( sys.argv[5], 'createPDS' ), '--type', 'target', '--pdb', extPath, '--pds', opdsPath ]
				sp.call( cpdsCMD )

				# run master search with this pds and the input pds(query), save top match in 
				mstrCMD 	= [ os.path.join( sys.argv[5], 'master' ), '--query', pdsPath, '--target', opdsPath, '--topN',
								'1', '--rmsdCut', '3.2', '--outType', 'full', '--structOut', 'tmp'
							  ] 
				outMatchPath = 'tmp/full1.pdb'
				sp.call( mstrCMD )

				sp.call( [ 'mv', outMatchPath, extPath ] )


			else:
				hsel 		= subset.copy().select( 'resnum %s %s' % ( hStr, all_match  ) )
				writePDB( hPath, hsel )

				# create target PDS of this match in tmp directory
				opdsPath 	= 'tmp/match.pds'
				cpdsCMD		= [ os.path.join( sys.argv[5], 'createPDS' ), '--type', 'target', '--pdb', extPath, '--pds', opdsPath ]
				sp.call( cpdsCMD )

				# run master search with this pds and the input pds(query), save top match in 
				mstrCMD 	= [ os.path.join( sys.argv[5], 'master' ), '--query', pdsPath, '--target', opdsPath, '--topN',
								'1', '--rmsdCut', '3.2', '--outType', 'full', '--structOut', 'tmp'
							  ] 
				outMatchPath = 'tmp/full1.pdb'
				sp.call( mstrCMD )

				sp.call( [ 'mv', outMatchPath, hPath ] )
				
			

			index += 1	
			print

		print
sp.call( ['rm', '-r', 'tmp'] )




######### 
'''
to fix problem of poor rmsd by Prody, try to align match+helix(target) to original(query)

>  ~/termanal/createPDS --type target --pdb searchFrags/n9_2strand_helixMatches/match49_0_3EOJ_ext.pdb

> ~/termanal/master --query searchFrags/n9_2strand.pds --target searchFrags/n9_2strand_helixMatches/match49_0_3EOJ_ext.pds --outType full --rmsdCut 10 --structOut tmp_rematch/ --topN 1

WORKS!!
'''