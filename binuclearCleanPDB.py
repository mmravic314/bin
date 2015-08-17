## Marco Mravic UCSF Biophysics DeGrado Lab August 2015 ##
# Input list of PDB's with PDB chain mappings to UniProt entries
# Download all chains once, comment out that section
# Determine which pdb chains have zinc binding capabilities in nature (despite crystallization metals) from uniprot
# Create hash table out of chain mappings to compare coordination ligands and UniProt descriptions
# Output list of PDBs with natural zinc2+ (and Mg2+?) to cull for seq ID
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
# python ~/bin/binuclearCleanPDB.py ./binuclearDivalentPDB.txt ~/pdb_080615 ~/tertBuilding/binucUniProt/
# python ~/bin/binuclearCleanPDB.py path2PDBchainUniprotList path2PDBfiles path2uniprot


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

# Look for contact to specified elements (str type) annotated in uniprot [between input residues(takes number only)]
def parseUniprot(path, atom, resList):
	with open(path) as file:
		for i in file:
			#print i.upper()
			if 'ZN' in i.upper():
				print i.rstrip()

	return

### Created this hash in the 'find binuclear sites...' section, which is now commented out
pairsByPdb = pic.load( open( os.path.abspath( os.path.join('./', 'biPairs_byPDB.pkl') ), 'rb') )


############ map coordinating residues to uniprot numbering and see if corresponds to known zinc sites  ############3

zn_cryst		=	[]
zn_annotated	=	[]
zn_dbase		=	[]
hitchains 		= 	[]

from collections import defaultdict

for s, c in sorted( pairsByPdb.items() ):
	print s, c
	# for each binuclear site in the set of pdbs
	for m in c:
		a1 = m.name.split('_')[0]
		a2 = m.name.split('+')[0].split('_')[0]
		#option 1: both are crystallographic zincs  -> next check if this is confirmed by uniprot/text 	 
		if a1 == 'ZN' and a2 =='ZN':
			print m

			if s not in zn_cryst:
				zn_cryst.append(s)
				
				# One a double zinc containing structure is found, look at its coordinating residues
				# make sure it's in uniprot as zinc binding

				chProt = defaultdict(list)	#chain/uniprots to check ( saving residues )
				for res in m.contacts:
					chain = s + res.split('_')[0].lower()
					chProt[chain].append( res )
				
				if len(chProt) == 0 : continue		# no contacts, skip .
				for u in chProt.keys():
					try:
						uPath = os.path.abspath( os.path.join(sys.argv[3], '%s.txt' % inMap[u] ) )
						print "entering...", uPath 
						parseUniprot( uPath,  'ZN' , chProt[u])
						print
					except KeyError:
						# What to do when there is no uniprot file!!!


			#sys.exit()




	#print s, os.path.exists( uniPath)
	#if os.path.exists( uniPath)

	
