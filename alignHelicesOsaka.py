
'''
Marco Mravic DeGrado Lab Jan 2016
# input 1: input directory
# input 2: output directory
# input 3: cluster of helical pairs
# input 4: path to TMalign binary (compiled from fortran)
# input 5: path to cPickle hashing of helix resnums for each helix-osaka match 
# input 6: path to full interface of OSAKA PDB

#for i in *helix_matches/; do cd $i; for f in `ls`; do ~/bin/dsspcmbi $f $f.dssp; done; cd ..; done
#python ~/bin/parseDssp.py n9_2strand_helixMatches/ FullInterface.pdb 1helix/ helix_hits_visual/  docked1helix.pkl
#python ~/bin/parseDssp.py c9_2strand_helixMatches/ FullInterface.pdb 1helix/ helix_hits_visual/  docked1helix.pkl


python ~/bin/alignHelicesOsaka.py ~/peptideAmyloid/searchFrags/1helix/ ~/peptideAmyloid/searchFrags/2helix ~/splayBundle/Cluster-004/ ~/bin/TMalign docked1helix.pkl FullInterface.pdb 
'''

import sys, os, subprocess as sp, re, cPickle as pic, shutil, numpy as np 
from prody import *
from operator import itemgetter
from itertools import groupby
from PDButil import ext_aHelixBOTH

########## Global Defs ################
# renaming elements in TMalgin outpur files
eleD = {  'CA':'C', 'N':'N', 'O':'O', 'C':'C' }
#aligned chain (flg=0), the floating chain (flg=1), target helix (flg=2)
chD  = {  0:'Z', 1:'Y', 2:'X' }
# resnums for interfacial residues of OSAKA target surface
interface = [10,12,14,16,18,21,23]
# Full interface pdb of osaka fragment
osakaPDB = parsePDB( sys.argv[6] )

picDict = pic.load( open( sys.argv[5], 'rb' ) )

# path to ideal 10 residue alpha helix based on GCN4
idealH = os.path.join( os.path.dirname( sys.argv[0] ), 'idealHelix10ala.pdb' )
print idealH

#######################################



###########   Main   ##################

# work through all helices docked on osaka frag
for f in sorted( os.listdir( sys.argv[1] ) ):
	if f[-4:] != '.pdb': continue
	if f[-7:] == 'EXT.pdb' or f == 'tmp.pdb': continue

	targetPDB = os.path.join( sys.argv[1], f )	

	# Begin 2helix/ tree
	if not os.path.exists( sys.argv[2] ):
		os.mkdir( sys.argv[2] )
		

	# Save only helical part of target PDB file as tmp for alignment of cluster helial pairs
	print f
	helix 	= picDict[  '1helix/' + f  ]
	selstr 	= 'chain X resnum %s' %  (' '.join( [x[0] for x in helix] ) ) 
	pathTmp = os.path.join(  sys.argv[1] , 'tmp.pdb' )
	writePDB( pathTmp,  parsePDB( targetPDB).select( selstr ) )
	targetPDB = parsePDB( pathTmp, subset = 'bb' )

	# Add 4 residues to each end based off ideal helix. Hand editted so adding 4 doesn't go into space away from interface
	# hardcoded the path to indela 10 residue helix from GCN4 coords Grigoryan et al 2011. see PDButil module

	extPDB 		= ext_aHelixBOTH( targetPDB, parsePDB( idealH, subset = 'bb' ) )
	newPath 	= os.path.join(  sys.argv[1] , f[:-4] + 'EXT.pdb' ) 
	writePDB( newPath,  extPDB )
	targetPDB 	= parsePDB( newPath )

	# Evaluate if this docked helix (after exending) has a >12 residue segment
	# if any backbone atom less than 4.3 A from any Amyloid C-Beta or backbone less more than 14 to the closest C-beta: kill that residue
	# Of those residues that are valid, find continuous segs on chain. abort further work if valid segment is <12 residues long 

	matrix 	= buildDistMatrix(  targetPDB, osakaPDB.select( 'name CB resnum %s' % ( ' '.join( [str(x) for x in interface ] ) )  ) )
	minD 	= [ min( x ) for x in matrix ]
	index 	= 0
	H1bad 	= []
	# look for "too close" and "not itneracting"
	for a in minD:
		if a > 13.5 or a < 4.2:
				resi = targetPDB.getResnums()[0] + int( index ) / 4
				if resi not in H1bad:
					H1bad.append( resi )
		index += 1
	H1good 	= [ x for x in targetPDB.select( 'calpha' ).getResnums() if x not in H1bad ]
	# find longest continuous segment
	segs = []
	for k, g in groupby(enumerate( H1good ), lambda (i,x):i-x):
			catch =  map(itemgetter(1), g)
			if len( catch ) > len( segs ): segs = catch
	if len(segs) < 15:
		print 'skipped', f, segs
		os.remove( newPath )
		continue
	
	print "FOUND eligible:", f, segs
	print "... NOW: aligning cluster of pairs to helix, writing in 2helix/ tree"
	writePDB( newPath, targetPDB.select( 'resnum %s' % ' '.join( [ str(x) for x in segs]  )  ) )
	# make sub dir for each each 1helix docking
	outPath = os.path.join( sys.argv[2], f.split('.')[0][:-8] )
	if not os.path.exists( outPath ):
		os.mkdir(outPath)		


	#inner loop through input cluster pairs, find best alignment to single helix
	for pdb in sorted( os.listdir(sys.argv[3]) ):
		
		mobilePDB 	= os.path.join( sys.argv[3], pdb )
		outF 		= os.path.join( outPath, pdb[:-4] + 'X.pdb' )

		# Full command to output a PDB file and shortened file to just read TMalign output text
		#cmd  	= ['bash', sys.argv[3], '-file1', targPath, '-file2', f.rstrip(), '-outFile', outF,  '-outputPDB', '-printCE']
		miniCmd = [sys.argv[4], mobilePDB, newPath]
		#print miniCmd, '\n'

		out = sp.Popen( miniCmd  , stdout = sp.PIPE ).communicate()[0]			#  full command
		#print out
		# parse output and
		print os.path.basename( mobilePDB ), 'to', os.path.basename( newPath )
		#### Aligned length=   15, RMSD=   0.58, Seq_ID=n_identical/n_aligned= 0.000  ###
		match = re.search( r'Aligned\slength=\s+(\d+),\sRMSD=\s+(\d\.\d\d),' , out)
		if not match:
			AlignLength, RMSD = 0.0, 0.0
			print "Null alignment... skip"
		#	continue
		else:
			AlignLength, RMSD = int( match.group(1) ) , float( match.group(2) )
			print AlignLength, RMSD
			if RMSD < 1.00 and AlignLength >= 12:
				outCmd = [sys.argv[4], mobilePDB, newPath, '-o', outF ]
				# output the superposition
				sp.Popen( outCmd  , stdout = sp.PIPE ).communicate()[0]	

		#	break
				print
	if len( os.listdir( outPath ) ) <= 1:
		shutil.rmtree( outPath )
		continue

	## Clear the extraneous files
	for pdb in os.listdir( outPath ):
		
		potPath = os.path.join( outPath, pdb ) 
		if pdb[-3:] != 'lig':
			os.remove( potPath )
			continue			
	
	# Load original template/target helix after it was extended
	extPDB = parsePDB( newPath )

	# For all the aligned files created, dump all the non-necessary files, reformat file and keep only closely interacting helical pair
	for pdb in sorted( os.listdir( outPath ) ):
		
		potPath = os.path.join( outPath, pdb ) 
		if pdb[-3:] != 'lig':
			continue

		# Reformatting
		inPDB = parsePDB( potPath, subset = 'bb' )

		# rewrite occupany = 1, beta-factor = 10.00, element = autodetect
		# track which continuous segment this is, the aligned chain (flg=0, ch Z ), the floating chain (flg=1, ch Y), target helix (flg=2, ch X)
		flg 		= 0
		prvRes 		= 0
		prvChain 	= 'A'
		for a in inPDB.iterAtoms():
			a.setBeta(10.00)
			a.setOccupancy(1.00)
			a.setElement( eleD[ a.getName() ] )
			if a.getResnum() < prvRes or prvChain != a.getChid() :
				flg += 1
			prvRes 		= a.getResnum()
			prvChain 	= a.getChid()
			a.setChid( chD[ flg ] )

		print "\n", 'checking aligned', pdb, 

		# If helical pair too far away...just skip and delete alignmetn file
		if not inPDB.select( 'calpha chain Y and within 10 of chain X' ): 
			os.remove( potPath )
			print
			continue

		
		# tricky heiristic to get only interacting portion of aligned floating helix 
		subsetStr 	=  list( inPDB.select( 'calpha chain Y and within 10 of chain X' ).getResnums() )

		for r in subsetStr[:]:
			for k in [ r-2, r-1, r+1, r+2 ]:
				if k > 0 and k not in subsetStr:
					subsetStr.append(k)
		print len( subsetStr ),

		if len( subsetStr ) < 15: 
			os.remove( potPath )
			print
			continue 

		## pick out only "good" segment of the aligned floating chain.. to validate
		subset		=   inPDB.select( 'chain Y resnum %s' % ' '.join( [ str(x) for x in subsetStr] ) ).copy() 

		# Now look at floating aligned helix with respect to their alignment on the OSAKA structure
		# if any backbone atom less than 4.3 A from any Amyloid C-Beta or backbone less more than 16 to the closest C-beta: kill that residue
		# Of those residues that are valid, find continuous segs on each chain and only save if it's >=12 long
		matrix 	= buildDistMatrix(  subset, osakaPDB.select( 'name CB' ) )
		minD 	= [ min( x ) for x in matrix ]
		index 	= 0
		H1bad 	= []
		for a in minD:
			if a > 14 or a < 4.2:
				resi = subset.getResnums()[0] + int( index ) / 4
				if resi not in H1bad:
					H1bad.append( resi )
			index += 1
		H1good 	= [ x for x in subset.select('calpha').getResnums() if x not in H1bad ]
		# find longest continuous segment, for the set of valid residues
		if len( H1good ) < 15:
			print "Segments found not too small, or nothing found!!!" 
			os.remove( potPath )
			continue
			
		segs = []
		for k, g in groupby(enumerate( H1good ), lambda (i,x):i-x):
			catch =  map(itemgetter(1), g)
			if len( catch ) > len( segs ): segs = catch

		print segs
		if len( segs ) < 15:				## Skip and delete alignment file if non-valid segment found
			os.remove( potPath )
			continue
		print "FOUND", pdb
		# If good floating helix found,
		## Join the "good" segment of the aligned floating chain and the original target helix
		
		finalSeg	= subset.select( 'resnum %s' % ( ' '.join( [ str(x) for x in segs ] ) ) ).copy() 
		outPDB 		= extPDB + finalSeg + osakaPDB
		newOutPath 	= os.path.join( outPath, pdb.split('.')[0] + '-screened.pdb' )  
		os.remove( potPath )

		# Make sure to save the residue numbers of the calpha chains to keep fixed during minimization
		# Chain X: those from original search match; Chain Y: 4 closest CA's to Chain X
		Xconst		= ' '.join( [ str( x ) for x in list( extPDB.select('calpha').getResnums() ) if x in [int( j[0]) for j in helix] ] )
		print list( extPDB.select('calpha').getResnums() ), [j[0] for j in helix]

		xyMatrix 	= buildDistMatrix( finalSeg.select('calpha')  ,extPDB.select('calpha') )
		store		= []
		for row, index in zip( xyMatrix, np.arange(  len( xyMatrix )  ) ):
			store.append( ( min(row), index ) )
		Yconst = ' '.join( sorted( [ str( finalSeg.select('calpha').copy()[ i[1] ].getResnum() ) for i in sorted( store )[:4] ] ) )

		title = 'ChX %s ChY %s' % ( Xconst, Yconst )
		print "WROTE:   ", title, newOutPath
		outPDB.setTitle( title )
		writePDB( newOutPath, outPDB )
		print
	

os.remove( pathTmp )

'''
-file1  /Users/mmravic/peptideAmyloid/searchFrags/tmp.pdb -file2 /Users/mmravic/splayBundle/Cluster-004/1AH7-007_005-0206_0242_A_0172_0189_A.pdb  -printCE

'''