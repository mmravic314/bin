
## following pairIntructions.txt, with 


# parse TM definintions file and rewrite fastas

import os, sys
from prody import *
from collections import defaultdict
from PDButil import UnNatAA

# xray@dgl-xray:~/jan2016TMpdb$ python ~/bin/breakup_TMpdbs.py pdbDimerList.txt extendedDefinitions.txt selectedOPM/



# extendeddefinintions.txt
storeDict 		= defaultdict(list)
storeDictFull 	= defaultdict(list)
with open( sys.argv[2] ) as file:
	for i in file:
				
		#write prody selection string if segment > 8, hash by pdbID
		if int( i[18:23].strip() ) - int( i[13:18].strip() ) > 8:

			storeDict[ i[:4] ].append( 'chain %s resnum %s to %s'% ( i[5], i[13:18].strip(), i[18:23].strip() ) )

		storeDictFull[ i[:4] ].append( 'chain %s resnum %s to %s'% ( i[5], i[13:18].strip(), i[18:23].strip() ) )

# access each tertiary pdb
pdbList = []
with open( sys.argv[1] ) as file:
	for i in file:
		print "Entering", i.rstrip(), '\n'

		inPDB = parsePDB(   os.path.join( sys.argv[3], i.rstrip() + '.pdb' )	 )

		for sel in storeDict[ i.rstrip() ]:
			obj = inPDB.select( sel ).copy()
			fasta = ''.join( [ UnNatAA[x] for x in inPDB.select( sel ).copy().iterResidues() ] )



