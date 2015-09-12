## Marco Mravic UCSF Biophysics DeGrado Lab August 2015 ##
## Using frequencies calculates by potential rotamer - ligand contacts, follow cut off and build term residues
## Build term (if each ligand has at least two potentials) for each ligand in PDB and remove redundants later
# Input a list of potential cutoff's and manipulations (e.g '+' means if residue close to both ligands, add frequencies)
# input directory of 
# input hash (stored in binary pickle) mapping pdb to metal ligand objects (for more info if necessary)
# Writes eligible TERMs fragment PDB files into new directory in working dir named after paramter label (first column)

from prody import *
import sys,os, numpy as np

#Example command line
# 	python ~/bin/writeTERMpdbs.py   diZNtermParams.txt ~/tertBuilding/ZN_db/  zN_freq/ --lig
#	python  ~/bin/writeTERMpdbs.py  path2Paramsfile	 path2ZnDatabase	 path2freqDir  --flag2includeLigandsinTERMpdb


# abspath to original pdb file from first line of local file path list
pdbDir = os.path.dirname( sys.argv[2] )
if not os.path.exists( pdbDir ):
	print "ERROR: path to PDB database directory does not exist"
	sys.exit()


class Parameter():

	def __init__( self, name, fcut, comboKey, dirPath ):
		self.comboDict 		= { '0': '', '+': 'add' }
		self.name 			= name
		self.fcut 			= float( fcut )
		self.comboMethod 	= self.comboDict[ comboKey ]
		self.dirPath		= dirPath

	def __repr__(self):
		return self.name

class TERM():

	def __init__(self, rList, oPath):
		self.rList 		= rList
		self.oPath 		= oPath
		self.atomGroup 	= AtomGroup()


	def __repr__(self):
		return self.oPath

	#For residue list in term, write a selection string to query presence of prody atoms
	def resList2selectionStr(self, ligandStr = '' ):
		byCh	= defaultdict(list)
		sel 	= ''

		for r in self.rList:
			ch 		= r.split(',')[0]
			resnum 	= r.split(',')[1]
			byCh[ch].append( resnum )
		for c, resis in byCh.items():
			sel += 'chain %s resnum ' % (c) + ' '.join( resis ) + ' or '

		#Remove final 'and'
		if len( ligandStr ) > 0:
			sel += ligandStr
		else:
			sel = sel[:-3] 


		return sel




#Establish Parameter sets
params = []
with open( sys.argv[1] ) as file:
	for p in file:
		if p.split()[0] == 'label':
			continue

		init 		= p.split()
		dirPath 	= os.path.abspath( './p' + init[0] + '/' ) 
		init.append(  dirPath )
		params.append( Parameter( *init   ) )

		if not os.path.exists( dirPath ):
			os.mkdir( dirPath )
			print 'wrote', dirPath

from collections import Counter, defaultdict

# Enter freq files to find residues influenced by the ligand (diZn sites)
for f in os.listdir( sys.argv[3] ):
	# Save (residue listing) per param set for each 
	metalSite = {}

	txt = open(   os.path.join( sys.argv[3],  f ) , 'rU' ).readlines()
	for p in params:
		resSet = {}
		lig = defaultdict(list)
		ligA = []

		#Rememeber for h=which ligand residues are being added ()
		for ln in txt:
			if ln[0] != '(':
				if ln[0] == 'R':
					continue
				
				ligA.append( ln.split()[0] )
				index = len( ligA ) - 1
				continue
			

			res 	= ','.join( ln.split(',')[0:2] ).strip('(')
			freq 	= float( ln.rstrip().split(',')[-1].strip(')') )
			#print res, freq

			# Record frequency if above user-input freq cut off
			# combine frequencies of residue contacting both ligands, if applicable
			try: 
				resSet[ res ]

				#  linearly combine frequencies. If this increases freq for a residue, +1 ligand counts for each lig
				if p.comboMethod == 'add':
					
					prvLigA = ligA[ index -1 ]
						
					if freq + resSet[ res ] > p.fcut:
							if prvLigA not in lig[res]:
								lig[res].extend( ligA[index-1 :index+1] )
							else:
								lig[res].append( ligA[0] )

					resSet[ res ] = freq + resSet[ res ]


				else:
					resSet[ res ] = max( freq, resSet[ res ] )
					if freq > p.fcut:
						lig[res].append( ligA[index] )

			except KeyError:
				resSet[ res ] = freq
				if resSet[ res ] > p.fcut:
					lig[res].append( ligA[index] )

		ligCnt = []

		for r,fr in resSet.items():		
			if fr <= p.fcut :
				del resSet[r]
			else:
				ligCnt.extend( lig[r] )

		ligCnt = Counter( ligCnt )




		# Only build a term if there is two or more residues contacting each ligand
		cnt = 0
		for x in ligCnt.values():
			if x > 1: cnt += 1
		if cnt != 2:
			continue   #skip to next parameter set

		# Expand residues to +/- of interacting residues
		resList = [   ]
		for r in resSet.keys():
			
			IDup	= ','.join(  [ r.split(',')[0], str(  int( r.split(',')[-1] ) + 1  ) ]  )
			if IDup not in resList:
				resList.append( IDup )

			IDup2	= ','.join(  [ r.split(',')[0], str(  int( r.split(',')[-1] ) + 2  ) ]  )
			if IDup2 not in resList:
				resList.append( IDup2 )
			
			IDdown	= int( r.split(',')[-1] ) - 1
			if IDdown >= 0:
				IDdown2 = ','.join(  [ r.split(',')[0], str(IDdown) ]    ) 
				if IDdown2 not in resList:
					resList.append( IDdown2 )

			IDdownV2	= int( r.split(',')[-1] ) - 2
			if IDdownV2 >= 0:
				IDdownV3 = ','.join(  [ r.split(',')[0], str(IDdownV2) ]    ) 
				if IDdownV3 not in resList:
					resList.append( IDdownV3 )

			if r not in resList:
				resList.append( r )

		# Store (sorted) potential contacting residues and path to write that PDB fragment file to 
		oPath = os.path.join( p.dirPath, os.path.splitext( f )[0] + '.pdb' )

		metalSite[ p.name ] = TERM( sorted( resList, key = lambda x: ( x.split(',')[0], int( x.split(',')[1] ) ) ), oPath ) 



	# Write PDB

	pPath = os.path.join( pdbDir, os.path.splitext( f )[0][:4] + '.pdb' )
	inPdb = parsePDB( pPath )
	#print pPath
	#print inPdb.select( 'chain A resnum 158 159 160 161 or chain D resnum 1')

	

	#write a prody selection string for each, and pass selection object to be written 
	print 'Writing', f, 
	for n, term in metalSite.items() :
		print n,

		#print term.rList
		#print term.resList2selectionStr()
		#print inPdb.select( term.resList2selectionStr() )
		
		if sys.argv[4] == '--lig':
			ligStr = 'serial ' + ' '.join(  os.path.splitext( f)[0].split('-')[-1].split('_') )
			writePDB( term.oPath, inPdb.select( term.resList2selectionStr( ligStr ) ) )

		else:
			writePDB( term.oPath, inPdb.select( term.resList2selectionStr() ) )


	print 
	


