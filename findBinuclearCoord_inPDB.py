## Marco Mravic UCSF Biophysics DeGrado Lab August 2015 ##
# Look through all structures in the Protein Data Bank (Downloaded August 6, 2015) 
# skip entries of resolution < 3.2, and all non-xray non-NMR.
# write list of pdbIDs that have binuclear coordination sites for divalent metals.

from prody import *
import os, sys, gzip

#Example command line:  python bin/cleanHydrolasePDB.py ~/tertiaryBuilding/hydrolase_pdb/
"""
'BA', 'barium'
'BE', 'beryllium'
'CD', 'cadmium'
'CA', 'calcium'
'CR', 'chromium'
'CO', 'cobalt'
'CU', 'copper'
'EU', 'europium'
'GD', 'gadolinium'
'GE', 'germanium'
'FE', 'iron'
'LA', 'lanthanum'
'PD', 'lead'
'MG', 'magnesium'
'MN', 'manganese'
'HG', 'mercury'
'NI', 'nickel'
'OS', 'osmium'
'PT', 'platinum'
'RU', 'ruthenium'
'SR', 'strontium'
'SN', 'tin'
'U', 'uranium'
'V', 'vanadium'
'Y', 'yttrium'
'ZN', 'zinc'
"""

from collections import Counter
divalent_metals = [
'BA',
'BE', 
'CD', 
'CA', 
'CR', 
'CO', 
'CU', 
'EU', 
'GD', 
'GE', 
'FE', 
'LA', 
'PD', 
'MG', 
'MN', 
'HG', 
'NI', 
'OS', 
'PT', 
'RU', 
'SR', 
'SN',
'U', 
'V', 
'Y', 
'ZN' ]

found_metals = []
count = 0
inFiles = os.listdir(  os.path.abspath( sys.argv[1] ) ) 
binuclearFlg = 0
binuclearPDBs =[]
print 

for f in inFiles:
	
	if os.path.splitext(f)[-1] != '.pdb':
		continue

	#print os.path.join( os.path.abspath( sys.argv[1] ), f )
	f_stream = open( os.path.join( os.path.abspath( sys.argv[1] ), f ), 'r')

	metalFlg = 0
	for i in f_stream:

		# X-ray and NMR only
		if i[:6] == 'EXPDTA':
			if 'NMR' in i.rsplit() or 'X-RAY' in i.rsplit():
				count +=1
			else: 
				break

		if i[:6] == 'HETNAM': 

			if i.split()[1].strip() in divalent_metals:
				found_metals.append( i.split()[1].strip() )
				metalFlg +=1  #contains divalent metal

		if i[:4] == 'ATOM': break


	# If a divalent metal is present, go back in and check for any divalent metals within 4 angstroms of each other 
	if metalFlg > 0:

		metals = []		#Array containing metals found in this PDB
		pdb = parsePDB( os.path.join( os.path.abspath( sys.argv[1] ), f ))
		# Find all potential divalent metals in the 
		for a in pdb.iterAtoms():
			if a.getElement() in divalent_metals:
				metals.append( a )

		# Deep look into pdb's with multiple metals, find pairs of divalent metals < 4.3 Angstroms away from each other
		if len(metals) < 2: continue
		else:
			#print f 
			cont = Contacts( pdb )
			for a in metals:
				#print "inspecting", a 
				for c in cont.select( 4.3, a.getCoords() ):

					#Ignore self atoms and non-divalent metals
					if c.getElement() not in divalent_metals: continue
					if c.getIndex() == a.getIndex(): continue

					# If at least one binuclear site found, store pdb and break loop to next pdb
					binuclearPDBs.append( f.split('.')[0] )
					binuclearFlg = 1
					break
				if binuclearFlg > 0: break 

				

			#print binuclearPDBs

			#sys.exit()   

print len(binuclearPDBs)
# One pdb's have been retrieved, output a report of the resdiues/chains nearby 
# as well as a list of uniprot ID's to download(uniprot-pdb mapping in report too)

#Uniprot module written by Kortemme lab Shane O'Connor in tools/bio/ where tools/ is from Klab guybrush.ucsf.edu server & is located in my /home/xray/bin/ 
from uniprot import *


for p in binuclearPDBs:
	print p
	uniprot_map('PDB_ID', 'ID', p)


#os.remove( "/home/xray/tertiaryBuilding/tmp.txt" )
c = Counter( found_metals )
#print sorted( c.items() ) 	
#print count
