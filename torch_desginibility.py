# make term-like fragments of model, search database for structural matches
# 

# python ~/bin/torch_desginibility.py ~/peptideAmyloid/OsakaModels/model_1.pdb ~/peptideAmyloid/ ~/peptideAmyloid/scores_rnd1


from prody import *
import sys, os, numpy as np, random, shutil
from collections import Counter, defaultdict

def selStr( lst ):			# return string from list
	return ' '.join( [ str(x) for x in lst] )


## given a number of contacts, detemine which of the residues sets in mSet match the contact profile found
def matchSet( cnts, mSet, m=False ):

	mtch 	= {}			# count the number of contacts in cnts are residues in each continuous residue set in array mSet
	inx 	= 0				# record which index of the match set, mSet, is being evaluated
	pick 	= (-1, -1)		# to sort current best index,  pick[0], and match value at that index, pick[1]

	# quick eval of middle helix makes deciscion of middle residue set if only middle residues found, despite overlaps
	if m and 11 not in cnts and 12 not in cnts:
		mSet = mSet[1:]

	# for each array of 7 residue sets in mSet, count how many residues in cnt are in that resdiue set
	# remember which of these 7 residue sets contained the most residues in cnts --> return that residue set
	for m in mSet:
		mtch[ inx ] = 0
		for k in cnts:
			if k[0] in m:
				mtch[ inx ] += k[1]

		if mtch[ inx ] > pick[1]:
			pick = ( inx, mtch[ inx ] )

		elif mtch[ inx ] == pick[1]:	

			# Do quick eval 
			if m and 21 in cnts:
				return mid3

			# randomly select winner to proceed between two
			else:
				val = random.choice( [0, 1] )
				sel = [inx, pick[0]][val]
				pick = ( sel, mtch[ sel ] )
#				print "random selection!"
#		print pick
		inx += 1
#	print
#	print cnts
#	print 
	return mSet[ pick[0] ]


class Resi:

	def __init__( self, resnum ):
		self.dist 	= 10
		self.resnum = resnum
		self.cnts	= []

	def __repr__(self):
		return str( self.resnum )

class Frag:

	def __init__(self, obj):
		self.atoms = obj

		self.resis 		= []
		self.contacts	= {}


inPDB 	= parsePDB( sys.argv[1] )
intF 	= [10,12,14,16,18,21,23]

intFset	= inPDB.select( '(not element H) and chain A B C D E F G H I J K resnum ' + selStr( intF ) ).copy()
osaka 	= inPDB.select( '(not element H) and chain A B C D E F G H I J K' ).copy()
helices = [ inPDB.select( '(not element H) chain %s' % (c) ) for c in ['X', 'Y', 'Z'] ]

## Abort model if any backbone heavy atoms of helix clash with amyloid backbone + C-beta

osakaCsh= osaka.select( 'bb or name CB' ).copy()
torchSet= inPDB.select( '(not element H) and bb and chain X Y Z' ).copy()

clashes = osakaCsh.select( 'within 3.2 of helices', helices=torchSet )

if clashes:

	txt = 'HH SH1 SH2 SH3 # CLASH!\n0 0 0 0 0'
	scorePath = os.path.join( sys.argv[3], os.path.basename( sys.argv[1] )[:-3] + 'sc'  )

	# if clash/contact found to osaka Carbon atom, drop model 
	for i in clashes:
		if i.getElement() not in ['O', 'N']:

			print 'Helix - sheet clash found. Aborting model!'
			oFile = open( scorePath, 'w' ) 
			oFile.write( txt )
#			os.remove( sys.argv[1] )
			sys.exit()

		# if one is potential H-bonder, find what contact is in  helix, see if H-bonder
		for a in iterNeighbors( i, 3.2, torchSet ):
			if a.getElement() not in ['O', 'N']:
				print 'Helix - sheet clash found. Aborting model!'
				oFile = open( scorePath, 'w' ) 
				oFile.write( txt )
#				os.remove( sys.argv[1] )

				sys.exit()

del osakaCsh
del torchSet


# make helix-sheet fragments, search
hSels	= [ np.arange( 7, 15 ), np.arange( 12, 20 ), np.arange( 17, 25 ) ]
hFrags 	= [ Frag( helices[1].select( 'resnum ' + selStr(h) ).copy() ) for h in hSels  ] 

end1    = [10, 11, 12, 13, 14, 15]
end2	= [16, 17, 18, 19, 20, 21]

mid1	=[11, 12, 13, 14, 15, 16]
mid2	=[13, 14, 15, 16, 17, 18]
mid3	=[15, 16, 17, 18, 19, 20]

Eset 	= [ end1, end2 ]
Mset 	= [ mid1, mid2, mid3 ]
Allset 	= [ end1, mid1, mid2, mid3, end2 ]	

tmp_dir = os.path.join( sys.argv[2] ) 
if not os.path.exists( tmp_dir ):
	os.mkdir( tmp_dir )
else:
	shutil.rmtree( tmp_dir )
	os.mkdir(  tmp_dir )

# Make search fragments
print
for fr in hFrags:

	# find three CB's closest to a interface residue 
	seeds 	= []
	stp 	= 0
	ends 	= False 		# flag to match helix termini 
	index	= hFrags.index( fr )

	if index in [0,2] 	:		# determine which helix to find
		ends = True

	for a in fr.atoms.select( 'name CB' ):
		fr.resis.append( Resi( a.getResnum() ) )
		minV = 10

		for k in iterNeighbors( a, 7.3, intFset): 
			minV = min( minV, k[2] )

			r 	= ( k[1].getResnum(), k[1].getChid() )
			rz 	= ( k[1].getResnum(), k[1].getChid(), a.getResnum() )
			if r not in fr.resis[stp].cnts:
				fr.resis[stp].cnts.append( r )
			
			try:
				fr.contacts[ rz ] = min( k[2], fr.contacts[ r ] )
			except KeyError:
				fr.contacts[ rz ] = k[2]   

		if minV == 10:
			fr.resis.pop(-1)
			stp -= 1
		else:
			fr.resis[stp].dist = minV


		stp += 1

	
	if len( fr.resis ) < 3 : # abort if too few potential contacts for frags
		print 'Too few possible helix - sheet interactions; no fragment made for frag', index
		print
		continue

	## Given contacts, make sure to pick adjacent strands with most contacts

	top = Counter( [ a[1] for a in fr.contacts.keys() ] ).most_common(3)
	c2  = ('Z', 0)		# this should be replaced or operation was a failure and should quit 
	for c in top:
		if top.index( c ) == 0:
			c1 = c
			continue

		# look left and right in ASCII integers for next chain number by interacting strand 
		if ord( c[0] )  in [  ord( c1[0] ) + 1, ord( c1[0] ) - 1   ]:
			c2 = c
			break

	if c2[0] == 'Z' or c1[1] < 3 or c2[1] < 3:
		print 'Too few possible helix - sheet interactions per sheet; no fragment made for frag', index
		print
		continue

	# Find residue set and make fragment
	seeds 	= defaultdict(list)
	# Given the contact 
	for c in fr.contacts.keys():
		if c[1] in [ c1[0], c2[0] ] :
#			k = ( c[0], c[1], fr.contacts[c] )
			if c[0] == 23:
				seeds[ c[1] ].append(  ( 21 , 1.0 / (fr.contacts[c] + 3.5 ) )  )
			elif c[0] == 10 and not ends:
				seeds[ c[1] ].append(  ( 11, 1.0 / ( fr.contacts[c] + 3.5 ) )  )
			else:
				seeds[ c[1] ].append(  ( c[0], 1.0 / fr.contacts[c] ) )

	resiByCh = {}

	for k, v in seeds.items():
		# Screen out if less then three contacts per sheet

		# find which pre-determined residue set best matches contact profile found
		if ends:
			sStr 	= 'chain %s resnum %s' % ( k, selStr(matchSet( v, Allset ) ) )
		else:
			sStr 	= 'chain %s resnum %s' % ( k, selStr(matchSet( v, Mset, True ) ) )

		fr.atoms += osaka.select( sStr ).copy()

	fragPath = os.path.join( tmp_dir, 'fragSH_%d.pdb' % (index) ) 
	writePDB( fragPath, fr.atoms )

	print 'wrote', fragPath, '\n'


#### If and only if at least two fragments are made...

if len( [ x for x in os.listdir( tmp_dir ) if x[-3:] == 'pdb' ] ) < 3:

	print 'not enough helix-sheet fragments made, abort!'
	shutil.rmtree( tmp_dir )
	sys.exit()

# Extract each helical interface, write to file
helices = [ inPDB.select( '(not element H) chain %s resnum %s' % (c, selStr( np.arange( 10, 22 ) )) ).copy() for c in ['X', 'Y'] ]

pair = helices[0] + helices[1]
fragPath = os.path.join( tmp_dir, 'fragHelices.pdb') 
writePDB( fragPath, pair )
print 'wrote', fragPath, '\n'