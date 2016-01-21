# input 1: directory of .dssp and .pdb files
# input 2: pdb to be docked on ( pdb helix is aligned already)
# input 3: output directory
# input 4: pre-approved files
# input 5: filename to store cPickle with hashed resnums of docked helix indexed by filename

#find residue/chain info from helix. Print to file with full amyloid to input dir
# Go through files, send to

'python ~/bin/parseDssp.py n9_2strand_helixMatches/ FullInterface.pdb 1helix/ helix_hits_visual/  docked1helix.pkl'

import os,sys, cPickle as pic
from prody import *

# Make list of files from the visually approved lit of hits, and only work on these files
approved = []
for k in os.listdir( sys.argv[4] ):
	approved.append(k.split('.')[0])

# set up pickled hash. Load if it exists and update. 
cpath 	= sys.argv[5]
if os.path.exists( cpath ):
	picDict = pic.load( open( cpath , 'rb' ) )
else:
	picDict = {}

# import 
surface = parsePDB( sys.argv[2] )
for f in os.listdir( sys.argv[1] ):

	if f[-4:] == 'dssp':

		dic 		= {}
		fullPath 	= os.path.join( sys.argv[1], f )
		basename 	= f.split('.')[0]
		helix 		= []

		if basename not in approved:
			continue

		# parse dssp file
		inF = open( fullPath, 'rU' )

		flg 		= 0		#start reading after first '#'
		seg 		= []
		current 	= []
		strandFlg 	= 0
		for i in inF:
			if flg:
				# save list of residues (segment) if, then clear current
				if i[13] == '!':

					seg.append( current )
					current = []
					continue

				current.append( (i[7:10].strip(), i[11].strip(), i[16].strip()) )

			if i[2] == '#':  flg += 1

		seg.append( current )
		# save residues, this catches end of file
		strandFlg = 0

		for k in seg:
			strandFlg = 0
			for d in k:
				if d[2] in ['B', 'E']: # ignore strands for strands
					strandFlg =1
					break
			if not strandFlg:
				helix = k	

		if len(  [x[0] for x in helix] ) < 8: continue			# Skip too short helices
		inF.close()

		selstr = 'chain %s resnum ' % helix[0][1] + ' '.join( [x[0] for x in helix] )   
		#print basename
		#print selstr

		helixPDB = parsePDB( fullPath[:-5] ).select( selstr )  ## File path string slicing here assumes xyz.pdb.dssp
		# Change chain ID of target helix to 'X'
		helixPDB.setChids( ['X' for x in helixPDB.iterAtoms() ] )
		#print helixPDB ################ stopped here

		writePDB( 'tmp.pdb', helixPDB)
		helixPDB = parsePDB( os.path.join(  os.getcwd(), 'tmp.pdb' ), subset = 'bb' )

		# Write helix BB coords into file with full target interface
		docked 	= surface + helixPDB
		outPath = os.path.join( sys.argv[3], basename[5:] + '%stemp.pdb' % os.path.basename( sys.argv[1])[:2] )
		writePDB( outPath, docked )

		# hash list of resnums by filename, store in cPickle file
		if outPath not in picDict.keys():
			picDict[outPath] = [x for x in helix] 




pic.dump( picDict , open( cpath , 'wb') )
os.remove( 'tmp.pdb' )
 
