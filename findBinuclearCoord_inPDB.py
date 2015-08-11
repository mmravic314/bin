## Marco Mravic UCSF Biophysics DeGrado Lab August 2015 ##
# Look through all structures in the Protein Data Bank (Downloaded August 6, 2015) 
# skip entries of resolution < 3.2, and all non-xray non-NMR.
# write list of pdbIDs that have binuclear coordination sites for divalent metals.

from prody import *
import os, sys, numpy as np#, gzip

#Example command line:  xray@dgl-xray:~/tertiaryBuilding$    python     ~/bin/findBinuclearCoord_inPDB.py     ~/pdb_080614/
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

#from collections import Counter#, defaultdict
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

############### Class def
# mapped PDB chain containing which pdb residue corresponds to each uniprot entry. 
# Also holds uniprot with download option
class PdbChain():
	# does mapping at initialization given a DBREF line from a PDB file
	def __init__(self, pdbLine):
		self.name = pdbLine[7:11] + pdbLine[12].lower() 
		self.uniprot = pdbLine[33:39]
		self.pdbRes = np.arange( int( pdbLine[14:18].strip() ) , int ( pdbLine[19:24].strip() ) + 1 )
		self.uniRes  = np.arange( int( pdbLine[55:60].strip() ) , int ( pdbLine[62:68].strip() ) + 1 )
		self.pdbMetalCoords = []

	def __repr__(self):
		return '%s, %s' % (self.name, self.uniprot) 

	def pdb2uni(self, resNum):

		try:
			val = list( self.pdbRes ).index( resNum )
		except ValueError:
			print "Requested pdb residue %d not in chain %s" % (resNum, self.name)
			sys.exit()

		try:
			val = self.uniRes[val]
		except IndexError:
			print  "Requested pdb residue %d not in uniprot mapping of chain %s" % (val, self.name) 
			sys.exit()
		
		return val


	def uni2pdb(self, resNum):
		
		try:
			val = list( self.uniRes ).index( resNum )
		except ValueError:
			print "Requested uniprot residue %d not in uniprot entry %s" % (resNum, self.uniprot)
			sys.exit()

		try:
			val = self.pdbRes[val]
		except IndexError:
			print  "Requested uniprot residue %d not in pdb mapping of chain %s" % (val, self.name) 
			sys.exit()
		
		return val
####################### Class definitions end #####################


found_metals = []
count = 0
inFiles = os.listdir(  os.path.abspath( sys.argv[1] ) ) 
binuclearFlg = 0
binuclearPDBs = {}
manualMappings = ''
print 

for f in inFiles:
	
	if os.path.splitext(f)[-1] != '.pdb':
		continue

	#print os.path.join( os.path.abspath( sys.argv[1] ), f )
	f_stream = open( os.path.join( os.path.abspath( sys.argv[1] ), f ), 'r')

	metalFlg = 0
	uniprotList = [] 
	exceptions = ''
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

		if i[:6] == 'DBREF ':
			#skip non-uniprot DBREF records 
			if i[26:32].strip() != 'UNP': continue

			#print i
			uniprotList.append( PdbChain( i ) )

		#Exceptions made for DEBREF split lines, parse these manually
		elif i[:5] == 'DBREF':
			try:
				int(i[6])
				exceptions += "DBREF SPLIT LINES:\n"
				exceptions += i
				

			except ValueError:
				pass


		


	# If a divalent metal is present, go back in and check for any divalent metals within 4 angstroms of each other 
	if metalFlg > 0:

		metals = {}	#Array containing metals found in this PDB
		pdb = parsePDB( os.path.join( os.path.abspath( sys.argv[1] ), f ))
		# Find all potential divalent metals in the 
		for a in pdb.iterAtoms():
			if a.getElement() in divalent_metals:
				metals[a] = []

		# Deep look into pdb's with multiple metals, find pairs of divalent metals < 4.3 Angstroms away from each other
		if len(metals) < 2: continue
		else:
			#print f 
			cont = Contacts( pdb )
			for a in metals.keys():
				#print "inspecting", a 
				for c in cont.select( 4.3, a.getCoords() ):

					#Grab coordinating residues ( )					
					#print c.getChid(), c.getResnum(), 
					#if calcDistance(c, a) < 2.4:
					#	idStr =  "%s, %d" % (c.getChid(), c.getResnum())
					#	if idStr not in metals[a]:
					#		metals[a].append( idStr )

					#Ignore self atoms and non-divalent metals
					if c.getElement() not in divalent_metals: continue
					if c.getIndex() == a.getIndex(): continue

					# If at least one binuclear site found, store pdb and break loop to next pdb
					binuclearPDBs[ f.split('.')[0] ] = uniprotList
					binuclearFlg = 1

					if len(exceptions) > 2:
						manualMappings += exceptions

					#For second metal, grab residues coordinating (distance R = < 2.4)
					#for d in cont.select( 4, c.getCoords() ):
					#	print d, calcDistance(c,d)
					#	if d == a: continue

					break


				if binuclearFlg > 0: break 

	#if binuclearFlg > 0: break				

			#print binuclearPDBs

			#sys.exit()   

print len(binuclearPDBs.keys())
# One pdb's have been retrieved, output a report of the resdiues/chains nearby 
# as well as a list of uniprot ID's to download(uniprot-pdb mapping in report too)

#Uniprot module written by Kortemme lab Shane O'Connor in tools/bio/ where tools/ is from Klab guybrush.ucsf.edu server & is located in my /home/xray/bin/ 
#from uniprot import *


for p,l in binuclearPDBs.items():
	print  p, l 
	#, #uniprot_map('PDB_ID', 'ID', p)
	#print p
	#uniprot_map('PDB_ID', 'ID', p)
print manualMappings



#os.remove( "/home/xray/tertiaryBuilding/tmp.txt" )
#c = Counter( found_metals )
#print sorted( c.items() ) 	
#print count
