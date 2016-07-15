
## following pairIntructions.txt, stop when homology reduction is needed. used TM-only database, homology reduced
## HERE: search PDB structures, relabel by chain  and find potential trimers with a loose distance map

import os, sys, cPickle as pic, numpy as np, time
from prody import *
from collections import defaultdict, Counter
from PDButil import UnNatAA
import networkx as nx 


## example command line
# python ~/bin/extractTrimerTM.py ~/jan2016TMpdb/non_redund.txt ~/jan2016TMpdb/TM_redundant/ ~/jan2016TMpdb/3_helix/ ~/jan2016TMpdb/extendedDefinitions.txt segLabeled_nonRedun-TM/

# in 1: list of PDBs eligible
# in 2:	path to fdirectory with TM-only version of PDBs
# in 3: target directory to make a pre-lim Trimer subdir & parent (renumbered/chained)
# in 4:


## I/O
tmPDB_dir 	= sys.argv[2]
outDir 		= sys.argv[3]
def_file 	= open( sys.argv[4], 'rU' )
targDir 	= sys.argv[5]

### Read in list of eligible PDBs:
inFile 	= open( sys.argv[1], 'rU' )
pdbList = []

for i in inFile:
	if len(i) > 1:
		pdb = i.rstrip()[:4]
	pdbList.append( pdb )

# locate Tm segments, collable as ProDy atom groups
storeDict 		= defaultdict(list)	
prvPDB 			= '' 
for i in def_file:
	pdb = i[:4]
	if pdb not in pdbList: continue
	# skip short segments
	if int( i[18:23].strip() ) - int( i[13:18].strip() ) < 11: continue

	storeDict[ i[:4] ].append( 'chain %s resnum %s to %s'% ( i[5], i[13:18].strip(), i[18:23].strip() ) )


##
def compare_Trimers( t1, t2 ):

			seqMat = np.zeros( (3,3) )
			c,r = 0, 0

			t1_seqs = [   ''.join( [ UnNatAA[ q ] for q in n.getResnames() ])   for n in t1.iterSegments() ]
			t2_seqs = [   ''.join( [ UnNatAA[ q ] for q in n.getResnames() ])   for n in t2.iterSegments() ]
			
			## Sequence comparions should be super negative,
			for k in t1_seqs:
				c=0
				for j in t2_seqs:

					result 			= pairwise2.align.globalxs(j, k, -1000, -1000 )[0][2]
					if result > 0:
						perc 		= result/min( len(k), len(j) )
					else:
						perc = 0
					seqMat[c][r] 	= perc
					

					c +=1
				r += 1

			mini = [ max(comp) for comp in seqMat ]
			if np.mean( mini ) > 0.95:
				return 1
			else:
				return 0
##

# load each pdb with at least 3 helices, rename/number TM segement, and log if close to any two helices 

def relabel_TMsegs( storeDict, targDir):
	prvPDB = ''
	trimersByPDB = {}
	for pdb, segs in storeDict.items():


		pdbPath = os.path.join( tmPDB_dir, pdb + '.pdb' )
		print pdbPath
		newPDBpth = os.path.join( targDir, pdb + '.pdb' )
		inPDB = parsePDB( pdbPath )
		prvCh = ''
		ascii = 97
		neighborDict = defaultdict(list)
	# relabel proteins into segments
		for h in segs:
		
			ch = h[6]
			if ch == prvCh:
				ascii += 1
				newCh = ch + chr(ascii)
			else:
				ascii = 97
				newCh = ch + chr(ascii)
	#	print h, newCh
			prvCh 	= ch
			inPDB.select( h ).setSegnames( newCh )
		writePDB( newPDBpth,  inPDB )
	return

## Relabel segments of TM helices by chain, then by segment
if len( os.listdir( targDir ) ) < 200: 
	relabel_TMsegs( storeDict, targDir )
	tmPDB_dir = targDir
else:
	tmPDB_dir = targDir


def pairwise_contacts( storeDict ):
	prvPDB = ''
	trimersByPDB = {}
	for pdb in storeDict.keys():

		pdbPath 		= os.path.join( tmPDB_dir, pdb + '.pdb' )
		inPDBFull 		= parsePDB( pdbPath )
		inPDB 			= inPDBFull.select( 'ca' ).copy()
		pairList 		= []

		print pdbPath
		segments = []
		segNames = []
		for k in inPDB.iterSegments():
			if len( k.getSegnames()[0] ) <2: 
				continue
			segments.append(k)
			segNames.append( k.getSegnames()[0] )

		linkMat = np.zeros( (len( segments ), len( segments ) ) )

		# find all pairwise chain interactions, log in matrix
		seg = 0
		for h in segments:

			mates = inPDB.select( 'exwithin 9 of hel', hel = h )
			if mates:
				mates = set( mates.getSegnames()  )
			else:
				continue
			
			for p in mates:

				if len(p) != 2: continue		## avoiding undocumented helices

				name = tuple( sorted( [p, h.getSegnames()[0]] ) )
				dMat 	= buildDistMatrix( h, inPDB.select( 'segment  %s' % p ) )
				cnt 	= 0

				for d in dMat:
					if min(d) < 8: 
						cnt += 1
						if cnt == 5: 
							linkMat[seg][segNames.index( p )]+=1
							break
			seg += 1

		## Analyze matrix for 3 body contacts
		# Grab all k=3 cliques
		net = nx.Graph(linkMat)
		trimers = [ x for x in nx.enumerate_all_cliques( net ) if len( x ) == 3 ]

		if len( trimers ) == 0:
			continue

		## Extract and print these prelim trimers

		saved = []

		# remove redundants via sequence,
		# input c-alpha prody objects


		for tri in trimers:

			segs 	= [ segNames[ k ] for k in tri ] 
			call 	= '%s_%s.pdb' % ( pdb, '-'.join( segs ) )
			newPath = os.path.join( outDir, call )
			atoms 	= ''
			for hel in segs:

				if len( atoms ) == 0:
					atoms = inPDBFull.select( 'segment %s' % hel ).copy()
					
				else:
					atoms += inPDBFull.select( 'segment %s' % hel ).copy()
			atoms.setTitle( '-'.join( segs ) )

			## Make comparison to all trimers in this PDB, to remove exact duplicates
			# accept the first trimer
			if len( saved ) == 0:
				saved.append( atoms.select( 'ca' ).copy() )
				writePDB( newPath, atoms )
			# remember each trimer made, and compare all new ones with those saved
			else:
				fail = 0

				for k in saved:
					if compare_Trimers( k, atoms.select( 'ca' ).copy() ) > 0:
						fail +=1
						break
				# only accept new trimer if its sequences dont match >95% with an existing trimer
				if fail ==0:
					saved.append( atoms.select( 'ca' ).copy() )
					writePDB( newPath, atoms )
				

	return

from Bio import pairwise2
from Bio.pairwise2 import format_alignment	

###  extract 'non-redundant' set of trimers, each with unique sequences

if len( os.listdir( outDir ) ) < 2:
	pairwise_contacts( storeDict )

###


def extract_windows( inPDB):

	## renumber residues by segment
	h_len  		= np.array( [0,0,0] )
	segments 	= []
	segN = 0
	for h in inPDB.iterSegments():
		if segN == 0:
			h_len[segN] = len( h.getCoords() ) - 1
		else:
			h_len[segN] = len( h.getCoords() ) + h_len[segN - 1]

		segN += 1


	# emumerate 14-residue helices in sets, determine distance matrix

	# all possible starting indices for the 3 helices: A, B, C
	A = np.arange(  0, h_len[0] - 13 )
	B = np.arange(  h_len[0] +1, h_len[1] - 13 )
	C = np.arange(  h_len[1] +1, h_len[2] - 13 )

#	print A
#	print B
#	print C
	# Score matrix to note the median inverse distance 
	#vMat 	= np.zeros( ( len(A), len(B), len(C) ) )
	#mMat 	= np.zeros( ( len(A), len(B), len(C) ) )
	minVals = []
	
	a, b, c, = 0,0,0		# counters to track the matrix of scores
	start = time.time()

	for h1 in A:
		for h2 in B:
			for h3 in C:
				arng 	= np.arange( h1, h1 + 14 ) 
				brng 	= np.arange( h2, h2 + 14 )
				crng 	= np.arange( h3, h3 + 14 )

				dMatAB 	= buildDistMatrix( inPDB[arng[0]:arng[-1] + 1], inPDB[brng[0]:brng[-1] + 1] )
				AB 		= np.mean( 1.0/ dMatAB )

				dMatAC 	= buildDistMatrix( inPDB[arng[0]:arng[-1] + 1], inPDB[crng[0]:crng[-1] + 1] )
				AC 		= np.mean( 1.0/ dMatAC )

				dMatBC 	= buildDistMatrix( inPDB[brng[0]:brng[-1] + 1], inPDB[crng[0]:crng[-1] + 1] )
				BC 		= np.mean( 1.0/ dMatBC )

				minVals.append( (np.mean( [ AB, AC, BC] ) , ( arng[0], brng[0], crng[0] ) ) )


				c += 1
			b +=1
		a += 1
	#print 'time elapsed:', time.time() - start

	validTrimers = [  ]  # sets of starting points at each position that are still okay (non-overlapping)

	for i in sorted( minVals, reverse = True ):
		
		# Skip loose bundles
		if i[0] < 0.065: break

		# take the best one, and record those positions that would consitute and overlap



		h1, h2, h3 = i[1]
		if h1 in A and h2 in B and h3 in C:

			#### ideal helices only?  
			# can do a loose phi/psi filter here if necessary, but would require full backbone 
			####
			validTrimers.append( i )
			print i

		else:
			continue


		#print 'c:', [ x for x in C if x <= crng[4] - 13 ].extendcrng[-5:]
		A2, B2, C2 = [], [], []
		for a in A:
			if a <= h1 - 9 or a >= h1 + 9: A2.append(a)
		for b in B:
			if b <= h2 - 9 or b >= h2 + 9: B2.append(b)
		for c in C:
			if c <= h3 - 9 or c >= h3 + 9: C2.append(c)
		A, B, C = A2, B2, C2


	return validTrimers

################# FUNCTIONS END  ##################


### work through the trimers, extracting highly interacting 12 residue windows 
cnt = 0
innerDir 	= os.path.join( outDir, 'bestWindows14/' )
if not os.path.exists( innerDir ):
		os.mkdir( innerDir )

for f in os.listdir( outDir ):
	if f[-3:] != 'pdb': continue

	pdbPath 	= os.path.join( outDir, f )
	print 'Entering:', pdbPath


	# renumber each segment, record segment length
	#inPDBFull 	= parsePDB( pdbPath )
	fullPDB 	= parsePDB( pdbPath )
	r = 0 
	for resi in fullPDB.iterResidues():
		resi.setResnum( r )
		r += 1 
	inPDB 		= fullPDB.select( 'ca' ).copy()
	
	# For this PDB look at all possible 12 residue windows
	windows = extract_windows( inPDB )


	for w in windows:
		a_r = ' '.join( [ str(x) for x in np.arange( w[1][0], w[1][0] + 14 )] )
		b_r = ' '.join( [ str(x) for x in np.arange( w[1][1], w[1][1] + 14 )] )
		c_r = ' '.join( [ str(x) for x in np.arange( w[1][2], w[1][2] + 14 )] )
		selStr		= 'resnum %s %s %s' % ( a_r, b_r, c_r )
		oldTitle 	= fullPDB.getTitle()
		newTitle	= oldTitle + '_%d-%d-%d' % ( w[1] )
		outPath 	= os.path.join( innerDir, newTitle + '.pdb' ) 
		atoms = fullPDB.select( selStr ).copy()
		atoms.setTitle( newTitle + ' Mean_inverse_interhelical_distance: %s' % (str( round( w[0], 4 ) )) )
		chain = ['A', 'B', 'C']
		stp = 0
		# Renumber segment by chain
		for h in atoms.iterSegments():
			h.setChids( [ chain[stp] for x in h.iterAtoms() ] )
			stp +=1 

		moveAtoms( atoms.ca, to=np.zeros(3), ag=True)
		writePDB( outPath, atoms )
		print 'wrote window', newTitle


	print
	#if cnt == 10: sys.exit()

#	cnt += 1

	








###