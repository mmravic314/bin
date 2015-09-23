# python org

import sys, cPickle as pic
from collections import Counter, defaultdict
from prody import *

#   python ~/bin/orgDiZnCoordResiInfo.py biPairs_byPDB.pkl 
metalDict 	= pic.load( open( sys.argv[1], 'rb') )
validZn 	= [ x.rstrip() for x in open(sys.argv[2], 'rU').readlines() if len(x.rstrip()) ==4  ]

#Collect all sites that only ahve asp, glu, and his, with at least one his per residue
hga 	= []
ideal 	= [ 'HIS', 'GLU', 'ASP' ]

for k, v in sorted( metalDict.items() ):
	if k not in validZn: 
		continue
	#print k

	# cnt 
	

	for p in v:
		siteD 	= defaultdict(list)
		dieFlg 	= 0			# Suicide if any non ASP, GLU, HIS ligands present
		# Di Zinc only
		if p.name[:2] =='ZN' and p.name.split('+')[1][:2] =='ZN':
		#print p, '\t', ' '.join( [ r[2:5] for r in p.contacts] )
			#print p
			for x in p.contacts: 

				for meInd in x.split('=')[1:]:
					if x[2:5] in ideal:
						siteD[ meInd ].append( x[2:5] )
					else:
						dieFlg += 1 						
				#print x[2:5], x.split('=')[1:]

				#if x.split('=')[1:] 

			#print
		else:
			dieFlg += 1
		


		#Ignore site if non his, glu, asp ligand present
		if dieFlg > 0:
			continue 
		else:
			for m,lig in siteD.items():
				cnt = Counter(lig)
				# at last one 'his' per metal
				if cnt['HIS'] == 0:
					dieFlg +=1
				# at least 2 ligands per metal (since at most one water/ligand and maybe one non-natural AA allowed)
				if sum( cnt.values() ) == 0:
					dieFlg +=1 

		if dieFlg > 0:
			continue 

		#These PDB's made the cut
		if k not in hga:
			hga.append( k )
			#print k 

sharedWater = []
#import subprocess as sp 
for i in hga:
	path 	= './znDBr/%s.pdb'%(i)
	print i
	pdb 	= parsePDB( path )
	cont 	= Contacts( pdb )

	for k in metalDict[i]:
		print k
		sel =  [ j[-4:] for j in k.name.split('+') ] 


		waterShared = []
		sFlg 		= 0 
		#m = pdb.select( 'serial %s' %( ' '.join( sel ) ) )
		for a in cont.select( 2.7, pdb.select('serial %s' % (sel[0] ) ).getCoords() ):
			if a.getResname() == 'HOH':
				waterShared.append( a )

		for a in cont.select( 2.7, pdb.select('serial %s' % (sel[-1] ) ).getCoords() ):
			#if a.getResname() == 'HOH':
			if a in waterShared:
					sFlg += 1
					break
				#print a.getResname(), a

		sharedWater.append( k )

		#print sel[0]
		#for a in findNeighbors( m.select( 'serial %s' % (sel[0]) ), 2.7 ):
		#	print a

		#print sel[-1]
		#for a in findNeighbors( m.select( 'serial %s' % (sel[-1]) ), 2.7 ):
		#	print a
	print 
	sys.exit()
	 



