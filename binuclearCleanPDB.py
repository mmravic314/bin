## Marco Mravic UCSF Biophysics DeGrado Lab August 2015 ##
# Input list of PDB's with PDB chain mappings to UniProt entries
# Download all chains once, comment out that section
# Determine which pdb chains have zinc binding capabilities in nature (despite crystallization metals) from uniprot
# Create hash table out of chain mappings to compare coordination ligands and UniProt descriptions
# Output list of PDBs with natural zinc2+ (and Mg2+/Mn2+?) to cull for seq ID
## Details:
# rEMOVE dna/rna ENTRIES
# If no uniProt, assume denovo design or catalytic antibody. Gotta hand eval, don't disclude

from PDButil import *  ## this is a homemade module
from prody import *	 ## http://prody.csb.pitt.edu/
import sys, os, numpy as np, cPickle as pic

#Input file format 
# 4BBJ [4BBJa, Q5ZWR1]
# 4ORB [4ORBa, P63328, 4ORBb, Q63810]
# 2HZH []
# 3F79 [3F79a, Q9I045, 3F79b, Q9I045, 3F79c, Q9I045, 3F79d, Q9I045, 3F79e, Q9I045, 3F79f, Q9I045]

# Example command line
# Run from directory containing pickle files (e.g. ~/tertBuilding/)
# python ~/bin/binuclearCleanPDB.py ./binuclearDivalentPDB.txt ~/pdb_080615/ ~/binucUniProt/ ~/pdb_chain_uniprot.csv 
# python ~/bin/binuclearCleanPDB.py path2PDBchainUniprotList path2PDBfiles path2uniprot path2SIFTS_pdb2uniProt_mapping

## This function creates a hash (saved to cPickle binary file) containing PdbChain objects (see PDButil.py in my local bin)
def SIFT_mapping(pdbList):
	SIFTS_pdb_uniP_Map = {}

	with open( sys.argv[4] ) as file:
		prvCh = ''	
		segments = []
		for i in file:
			if i[:4].upper() not in pdbList:
				continue

			chId = ''.join( [ i[:4].upper(), i[5].lower() ] )

			if prvCh == '': 
				segments.append(i)
				prvCh = chId
				continue

			# If new chain is reached, initialize & store the previous. Then, begin the next
			if prvCh != chId:
				SIFTS_pdb_uniP_Map[ prvCh ] = PdbChain( 'SIFTS', segments )
				segments = []
				segments.append(i)
			#For segments of the same chain, keep adding
			else:
				segments.append(i)

			prvCh = chId

	return SIFTS_pdb_uniP_Map




####### This section commeted out but ran once to get all UniProt entries and save a hash of uniprot-PDB mapping (chain level)
#inMap = {}
#from collections import defaultdict
#import subprocess as sp
#pdbs = {}
#uniProtList = []
#with open(sys.argv[1]) as file:
#	for i in file:
#		pdbs.append(i[:4])
#
#		for ent in [ x.strip(']').strip('[').strip(',') for x in i[4:].split() ]:
#			if len(ent) < 2: continue 
#
			# remember pdb chain
#			elif ent[:4] == i[:4]:
#				storeKey = ent
#
			# dump pdb chain with uniprot in hash inMap
#			else:
#				try: 
#					inMap[storeKey] = ent
#					#if ent not in uniProtList: 
						#uniProtList.append( ent )
						#cmd = 'wget -P /home/xray/tertBuilding/binucUniProt/ http://www.uniprot.org/uniprot/%s.txt' % (ent)
						#sp.call( cmd.split() ) 
#
#				except NameError:
#					print "HEY ERRORZ! pdb chain %s doesn't match pdb entry for %s" % ( ent, i[:4] )
#					sys.exit()

#pic.dump( inMap, open('/home/xray/tertBuilding/binucUniProtPDBmap.pkl', 'wb') )
#pic.dump( pdbs, open('/home/xray/tertBuilding/binucPDB.pkl', 'wb') )
#
################

inMap 	= pic.load( open(   os.path.abspath( os.path.join('./' , 'binucUniProtPDBmap.pkl') ) ,  'rb') )
pdbs 	= pic.load( open(  os.path.abspath( os.path.join('./', 'binucPDB.pkl') ), 'rb') ) 
## Parse each PDB and UniProt entries to create lists of annotated/validated 'natural' binuclear zincs and less strict di-zinc sites
## Save each PDB chain object in hash with metal pair objects detailing contacting residues
chainObj = {}



########## Metal pair class  ##########
#Initialize with prody atom objects for metals and list of prody residue objects contacting
# Tihs class slims down info each atom carries to only store essentials (reduce filesize/RAM used)
class biMsite():
	def __init__(self, atom1, atom2, residue_matesList, pdb):
		self.name = "%s_%s+%s_%s" % ( atom1.getElement(), str( atom1.getIndex() ), atom2.getElement(), str( atom2.getIndex() ) )
		self.pdb = pdb
		self.contacts = []
		#For coordinating residues, loss of metal-res link, but bidentate should be store as 'A_ASP180-OD2+OD1'
		#self.contacts = [ '%s_%s%s-%s' % ( mate.getChid(), mate.getResname(), str( mate.getResnum() ), mate.getName() ) for mate in residue_matesList ] 
		for mate in residue_matesList:
			biFlg = 0
			fullID 	= '%s_%s%s-%s' % ( mate.getChid(), mate.getResname(), str( mate.getResnum() ), mate.getName() )	 
			if len( self.contacts ) == 0:
				self.contacts.append( fullID )
			else:
				rBaseId = '%s_%s%s' % ( mate.getChid(), mate.getResname(), str( mate.getResnum() ) ) 

				for c in self.contacts:
					if c.split('-')[0] == rBaseId:
						new = c + '+' + mate.getName()
						old = c
						biFlg +=1
				if biFlg > 0 :
					self.contacts.remove( old )
					self.contacts.append( new )
				else:
					self.contacts.append( fullID )


				

	def __repr__(self):
		return self.name

###### class end ######################


######### Find binuclear sites and store surrounding information #############	


def find_sites(pdbs):
 
 pairsByPdb = {}
 metalFlg = 0 
 for s in pdbs:

	pdbPath = os.path.join( sys.argv[2], s + '.pdb' )
	if not os.path.exists(pdbPath):
			print "ERROR Path not found %s" % (pdbPath)
			sys.exit(0)
		

	pdb = parsePDB(pdbPath)
	cont = Contacts( pdb )
	metals = []
	skip = []
	pairs = {}
	final_PairsV2 = []
	finalPairs = {}
	for atom in pdb.iterAtoms():
		if atom.getElement() in divalent_metals:
			metals.append( atom )
	
	for atom in metals:
		
		if atom in skip: continue   		#skip atoms known to be in sites with >2 metals
		metal_mates = []
		residue_mates = []
		for mate in cont.select( 4.3, atom.getCoords() ):
			
			if mate == atom: continue #ignore self


			# Store other metals
			if mate.getElement() in divalent_metals:
				#print 'FOUND pair', atom, mate, '\n'
				metal_mates.append( mate )
				continue

			#store coordinating residues with loose distance cut off (should be metal dependent really)
			if calcDistance( mate, atom ) <= 2.7:
				# Exclude carbons
				if mate.getElement() =='C': continue

				coID = '%s, %s%s' % ( mate.getChid(), mate.getResname(), str( mate.getResnum() ) ) 
				if mate.getResname() not in natAA.keys(): #only allow natural amino acids coordinating for now
					continue
				if coID not in residue_mates:
					residue_mates.append( mate )

		# Bi nuclear only (drop singles and trinuclear coord sites)
		if len( metal_mates ) == 1:
			# Go back and add second metal's coordinating residues (mates) 			
			for mate2 in cont.select( 2.7, metal_mates[0].getCoords() ):
				# Exclude carbons
				if mate2.getElement() =='C': continue

				coID = '%s, %s%s' % ( mate2.getChid(), mate2.getResname(), str( mate2.getResnum() ) ) 
				if mate2.getResname() not in natAA.keys(): #only allow natural amino acids coordinating for now
					continue
				if coID not in residue_mates:
					residue_mates.append( mate2 )


			# sort indices and compare to those already found so not double storing atom pairs (e.g. labeled 101-104 and 104-101) 
			pair = tuple( sorted( [ metal_mates[0], atom ], key = lambda x: x.getIndex() ) )
			

			if pair not in pairs.keys():
				pairs[pair] = residue_mates 

		# Ignore any atom in a >2 metal site with the blacklist 'skip'
		if len( metal_mates ) > 1 :
			for m in metal_mates:
				skip.append(m)
			skip.append(atom)
	
	# skip pairs that were found to be parts of tri/quqd metal sites
	for k, v in pairs.items():
		if k[0] in skip: continue
		if k[1] in skip: continue
		final_PairsV2.append( biMsite( k[0], k[1], v, s )  )

	if len( final_PairsV2 ) != 0:
		pairsByPdb[s] = final_PairsV2

 return pairsByPdb

#pairsByPdb = find_sites( pdbs )
#pic.dump(pairsByPdb, open( '/home/xray/tertBuilding/biPairs_byPDB.pkl', 'wb') )

### Created this hash in the 'find binuclear sites...' section, which is now commented out
pairsByPdb = pic.load( open( os.path.abspath( os.path.join('./', 'biPairs_byPDB.pkl') ), 'rb') )

#Create mapping of multinuclear sites (from inMap) OR just load the pre-made pickle
## PdbChain objects indexed in hash by chain
#SIFTS_pdb_uniP_Map = SIFT_mapping( pairsByPdb.keys() )
#pic.dump( SIFTS_pdb_uniP_Map, open( '/home/xray/tertBuilding/SIFTS_mapping.pkl' , 'wb' ) )
#sys.exit()
SIFTS_pdb_uniP_Map = pic.load(  open( '/home/xray/tertBuilding/SIFTS_mapping.pkl' , 'rb' ) )



### Functions for the following section
# PDB Mapping from pdb residue index to uniprot residue index
# try proper PDB file, then path to the hash (made above) from EMBL-EBI SIFTS (https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html)
# Return a list  containing PdbChain objects (see PDButil.py module in my local bin) represented by string '(pdbChain, uniprot)'
def find_mapping( path2pdb, chain, SIFTS_flg = False):
	tmp = ''
	with open(path2pdb) as file:
		for i in file:
			if i[:6] == 'DBREF ':
				tmp = PdbChain( 'pdb', i )
				if tmp.name == chain:
					break
				else:
					tmp = ''

			if i[:4] == 'ATOM': break

	# SIFTS_flg just forces replacement of PDB-DBREF mapping to SIFTS mapping
	if tmp != chain or not SIFTS_flg:
		try: 
			tmp = SIFTS_pdb_uniP_Map[ chain ] 
		except KeyError:
			print "No Valid uniprot - PDB mapping available for", chain

	return tmp


# Look for contact to specified elements (str type) annotated in uniprot [between input residues(takes number only)]
def parseUniprot(path, atom, resList, pChain ):
	mention = False
	reviewed = False
	metal_seen = False
	add2_DB = False
	coordinating = []
	divalents = []
	#Parse
	with open(path) as file:
		for i in file:
			if i[:2] == 'ID': 
				if 'Reviewed' in i:
					reviewed = True

			## If metal lines found and Zinc not mentioned, Turn false positive to false negative
		
			if i[0:10] == 'FT   METAL':
				coordinating.append( i )
				if 'DIVALENT METAL CATION' in i.upper():
					if i.split()[2] != i.split()[3]: 	#If string of residues bind, hand this annotation differently
						span = [ str(x) for x in np.arange( int(i.split()[2]), int( i.split()[3] ) + 1  ) ]
						divalents.extend(span)
					else:										#This is normal way to record binding sites in metal binding annotations
						divalents.append( int( i.split()[2] ) ) 

			if 'ZN' in i.upper() or 'ZINC' in i.upper():
				#print i.rstrip()
				mention = True

	#### Evaluate if this should be included in zinc data base 
				
	# If entry has not be reviewed, and there is a dizinc site  --> accept as true site
	if not reviewed:
		add2_DB = True
		print "not reviewed"
		return add2_DB


	# If entry has been reviewed but no mention of zinc, and there is a dizinc site  --> accept as true site
	if reviewed:
		if len( coordinating ) == 0:
			add2_DB = True
			print 'Reviewed, but No coordinating residue info'
			return add2_DB

		# if reviewed and ZN has residue info, check if enough info and if it matches PDB contacts found previously
		else:
			annotatedZn = []
			non_Zn		= []
			for ln in coordinating:

				if 'ZN' in ln.upper() or 'ZINC' in ln.upper():
					if ln.split()[2] != ln.split()[3]: 	#If string of residues bind, hand this annotation differently
						span = [ x for x in np.arange( int(ln.split()[2]), int( ln.split()[3] ) + 1  ) ]
						annotatedZn.extend(span)
					else:										#This is normal way to record binding sites in metal binding annotations
						annotatedZn.append( int( ln.split()[2] ) ) 
				
				else:
					non_Zn.append( ln )

			# If only Zn annotations in reviewed uniprot, and binuclear Zn present  --> accept as true
			if len(non_Zn) == 0 and len(annotatedZn) >= 0:
				add2_DB = True
				print "Reviewed, only coordinating zincs"
				return add2_DB
				#print non_Zn

			# If no Zn, but "divalent metal binding" listed, check if residues of coordinating sites match			
			if len( divalents ) > 0:
				# Get residues contacting the metal in question and convert residue numbering to uniprot numbering via mapping function

				metal_cont = [ int( filter(lambda x: x.isdigit(), res.split('-')[0] ) ) for res in resList]

				# Use mapping function to get chain object to do mapping conversion
				path2pdb = os.path.abspath( os.path.join(sys.argv[2], '%s.pdb' % chain[:4] ) )
				ChObj = find_mapping( path2pdb, pChain )


				#Convert pdb numbering to uniprot residue indices
				foundSites = []
				for res in metal_cont:
					try: 
						foundSites.append( int( ChObj.pdb2uni( res )) )   
					except ValueError:
						pass
				mtch = 0

				for i in foundSites:
					if i in divalents:
						mtch +=1 
						continue
				# If >= 80% of observed ZN coordinated residues are listed in vague metal sites, accept as valid ZN site. 
				correct = float( mtch ) /len(foundSites)
				if correct >= 0.8:
					add2_DB = True
					print 'Reviewed, lists \"divalent metal cation\" and observed ZN sites match these general metal sites'
					return add2_DB


			# If multiple metals, including those residues coordinating Zn, check if residues match
			if len(non_Zn) > 0 and len(annotatedZn) > 0:


				metal_cont = [ int( filter(lambda x: x.isdigit(), res.split('-')[0] ) ) for res in resList]

				# Use mapping function to get chain object to do mapping conversion
				path2pdb = os.path.abspath( os.path.join(sys.argv[2], '%s.pdb' % chain[:4] ) )
				ChObj = find_mapping( path2pdb, pChain )


				#Convert pdb numbering to uniprot residue indices
				foundSites = []
				for res in metal_cont:
					try: 
						foundSites.append( int( ChObj.pdb2uni( res )) )   
					except ValueError:
						pass
				mtch = 0

				for i in foundSites:
					if i in annotatedZn:
						mtch +=1 
						continue
				# If >= 80% of observed ZN coordinated residues are listed in vague metal sites, accept as valid ZN site. 
				correct = float( mtch ) /len(foundSites)
				if correct >= 0.8:
					add2_DB = True
					print 'Reviewed, lists \"divalent metal cation\" and observed ZN sites match these general metal sites'
					return add2_DB



		if not mention:
			if len(coordinating) == 0:
				add2_DB = True
				print "Reviewed, but not mentioned and no annotated metals found"
				return add2_DB

	print add2_DB
	return add2_DB



############ map coordinating residues to uniprot numbering and see if corresponds to known zinc sites  ############3

zn_cryst		=	[]
zn_annotated	=	[]
zn_dbase		=	[]
hitchains 		= 	[]
cnt = 0

from collections import defaultdict


for s, c in sorted( pairsByPdb.items() ):
#	print s, c
	# for each binuclear site in the set of pdbs
	
	for m in c:
		a1 = m.name.split('_')[0]
		a2 = m.name.split('+')[1].split('_')[0]

		#option 1: both are crystallographic zincs  -> next check if this is confirmed by uniprot/text 	 
		if a1 == 'ZN' and a2 =='ZN':
			
			if s not in zn_cryst:
				print s,m
				zn_cryst.append(s)
				
				# One a double zinc containing structure is found, look at its coordinating residues
				# make sure it's in uniprot as zinc binding
				chProt = defaultdict(list)	#chain/uniprots to check ( saving residues )
				for res in m.contacts:

					chain = s + res.split('_')[0].lower()
					chProt[chain].append( res )
				
				if len( chProt[chain] ) == 0 : 
					print 'no contacts to chain'
					continue		# no contacts, skip .


				for u in chProt.keys():
					try:
						uPath = os.path.abspath( os.path.join(sys.argv[3], '%s.txt' % inMap[u] ) )
						print "entering...", uPath 
						if parseUniprot( uPath,  'ZN' , chProt[u], u ): 
							if s not in zn_dbase:
								zn_dbase.append(s)
								print
								break
						print
					except KeyError:
						print 'no uniprot for this chain!!!!!!!!!!!!!!!!!!!!!!!', u
						if s not in zn_dbase:
							zn_dbase.append(s)
							print 
							break
			else:
				break		# quit iterating if pdb qualifies

			print
			#sys.exit()

print '\nTRUE DI-ZINC PDBs FOUND', len( zn_dbase )
print 'out of potential di-Zinc:', len( zn_cryst )

for i in zn_dbase:
	print i

	#print s, os.path.exists( uniPath)
	#if os.path.exists( uniPath)

	
