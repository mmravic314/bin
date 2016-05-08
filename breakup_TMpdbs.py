
## following pairIntructions.txt, with 


# parse TM definintions file and rewrite fastas

import os, sys, re, cPickle as pic
from prody import *
from collections import defaultdict, Counter
from PDButil import UnNatAA


# xray@dgl-xray:~/jan2016TMpdb$ python ~/bin/breakup_TMpdbs.py pdbDimerList.txt extendedDefinitions.txt selectedOPM/ TM_redundant/ TM_redundant_Fasta



# extendeddefinintions.txt
storeDict 		= defaultdict(list)
storeDictFull 	= defaultdict(list)
with open( sys.argv[2] ) as file:
	for i in file:
				
		#write prody selection string if segment > 8, hash by pdbID
		if int( i[18:23].strip() ) - int( i[13:18].strip() ) > 8:

			strk = i[5], i[13:18].strip(), i[18:23].strip()
			storeDict[ i[:4] ].append( 'chain %s resnum %s to %s'% ( i[5], i[13:18].strip(), i[18:23].strip() ) )

		storeDictFull[ i[:4] ].append( 'chain %s resnum %s to %s'% ( i[5], i[13:18].strip(), i[18:23].strip() ) )

motif = defaultdict(list)

# access each tertiary pdb
pdbList = []
with open( sys.argv[1] ) as file:
	for i in file:
		print "Entering", i.rstrip(), '\n'

		
		print i
		# fetch resolultion
		with open(  )

		sys.exit()
		inPDB 	= parsePDB(   os.path.join( sys.argv[3], i.rstrip() + '.pdb' )	 )

		for sel in storeDict[ i.rstrip() ]:

			# grab each
			obj = inPDB.select( sel ).copy()

			fasta = ''.join( [ UnNatAA[x.getResname()] for x in inPDB.select( sel ).copy().iterResidues() if x.getResname() in UnNatAA.keys()] )

			

			#find right-parallel GAS-x3-gas-x3
			match = re.search( r'[GAS]\w\w\w[GAS]\w\w\w[GAS]', fasta  )
			if match:
				# remove lo complexity 
				if 'GGG' in match.group() or 'AAA' in match.group() or 'G' not in match.group(): continue 

				print match.group(), obj[0].getChid(), i.rstrip()


				motif[ i.rstrip() ].append( match.group() )


		#	print fasta, obj[0].getChid()

		print 

for k, v in motif.items():

	if len(v) == 0: continue
	c = Counter( v )


	# find single pass multimers only
	if len( c.keys() ) > 1 : continue
	if c.most_common()[0][1] < 2: continue 

	print k, c.most_common()[0][0], c.most_common()[0][1]

	print 





