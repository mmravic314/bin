#quickAlignOsk.py
## python ~/bin/quickAlignOsk.py ~/peptideAmyloid/parameterization/FullInterface.pdb ~/peptideAmyloid/parameterization/zeroAlignment.pdb 
# For alignment of the OSAKA amyloid interface  to 5 points in a new "centered" coordinate frame

import os, sys
from prody import * 

inPDB  	= parsePDB( sys.argv[1] )
outPath	= sys.argv[1][:-4] + 'CENTERED.pdb'

frame 	= parsePDB( sys.argv[2] ) 

atoms 	= []
serials = []

subset 	= inPDB.select( 'ca chain B D resnum 15' ).copy() + inPDB.select( 'ca chain C resnum 10 15 20' ).copy()

for i in subset.iterAtoms():
	print i.getResnum(), i.getChid(), i.getName()
trans = superpose( subset, frame )[1]
print calcRMSD( subset, frame )
writePDB( outPath, applyTransformation( trans, inPDB )  )