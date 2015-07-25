## Marco Mravic UCSF Biophysics DeGrado Lab July 2015 ##
#Example command line
##python ~/bin/add_ligand_contacts_Confind.py --p 3A4J/3A4J.pdb --rout 3A4J/3A4J.rout --lig 2930,2931 --rcut 7 --freq 0 --c 3A4J/3A4J.contacts 

helpStr = """
 Modification step to TERM (Zhang et al 2015, Structure) analysis after CONFIND and before generating TERMS, 
 consider that potential side chain interactions with ligand atoms
 influence that residues identity/conformation and those also contacting ligand.
	Input 1: path to original pdb file
	Input 2: path to Confind --rout rotamer output file (Confind.cpp by Gregoryian from Molecular Software Library, MSL)
	Input 3: Selection string of Atom serial numbers from ligand (e.g 1-10,35,46,60-100)
	Input 4: --rcut distance shell to consider around ligand, DEFAULT: 7 A (per ligand --rcut feature can be added later) 
	Input 5: --freqcut joint freqency cutt off of 2 residues contacting ligand DEFAULT: 0.1 (must be >0.1)
	Input 6: path to contact list file 

	Do: find residues with heavy atoms some distance "rcut" from each specified ligand atom calculate fration that interact (similar to Zhang et al 2015) rewrite contact list given an interaction frequency cut off (expandes TERM fragment)

"""

from prody import *
import os, sys, numpy as np, cPickle as pic, argparse

# Thorough help options and argument parsing 
helpOptions = ['-h', '-help', '--h', '--help']
for i in sys.argv:
	if i in helpOptions:
		print helpStr
		sys.exit()

par = argparse.ArgumentParser()
par.add_argument('--p', required = True, help = 'input PDB file')
par.add_argument('--c', required = True, help = 'path to .cmap contact list file')
par.add_argument('--lig', required = True, help = 'ligand atom(s) serial selection string (inclusive)')
par.add_argument('--rcut', default = 7.0, help = 'ligand interaction shell cutoff')
par.add_argument('--freqcut', default = 0.1, help = 'joint ligand interaction frequency cutoff')
par.add_argument('--rout', required = True, help = 'path to contact file to rewrite')
args = par.parse_args()

#####################################################################
################## Helper classes and functions  ####################
#####################################################################

# Turn valid atom selection string into a integer list
# complain & quit if invalid string
def parseLigandStr( lig ):
	l = []
	try: 
		lig = lig.split(',')
		for i in lig: 
			if '-' in i:
				if len( i.split('-') ) != 2: 
					print "\nERROR: ligand string syntax/formatting incorrect (1-2-300 failure), try --help for more information\n"
					sys.exit()
				else: 
					# Add each range to end of list
					l.extend( list( np.arange( int( i.split('-')[0]) , int( i.split('-')[-1]) + 1 ) ) ) 
			else: 
				l.append( int(i) ) 

	except (ValueError, TypeError, AttributeError, NameError ):
		print "\nERROR: ligand string syntax/formatting incorrect, try --help for more information\n"
		sys.exit()

	return l

class Ligand():
	def __init__(self, serial, coords ):
		self.name = serial
		self.coords = coords
		contacts = []

	def __repr__(self):
		return self.name

class Residue():
	def __init__(self, serial, coords ):
		self.name = serial
		self.coords = coords
		contacts = []

	def __repr__(self):
		return self.name



########################  MAIN  #####################################


ROTLIB = '~/termanal/support.default/rotlib/RR2000.rotlib'
routf	= args.rout
pdbf	= args.p
rcut	= args.rcut 
freqcu 	= args.freqcut
cmap 	= args.c
ligStr 	= args.lig
ligStr 	= parseLigandStr( ligStr )
print

### Parse Pdb, using ProDy package, and get input ligand coordinates
##  make list of residues (chain,resnum) with heavy atoms within 7 angstroms
pdbGroup = parsePDB( pdbf )
print "\n"

contactList = []
for i in ligStr:
	lig_coords = pdbGroup.getBySerial( i, i+1 ).getCoords()[0]

	if len( lig_coords ) < 3:
		print "\nWARNING: No valid coordinates of ligand of serial number", i


	cont = Contacts( pdbGroup )
	for j in cont.select( rcut , lig_coords ):

		ID = ( j.getChid() , str( j.getResnum())  )
		if ID not in contactList:
			contactList.append( ID )


########### temp section for handling rotamer library
######## Saving rotamer probability look-up table a pickle 
#RR = open( '/home/xray/bin/RR2000.rotlib', 'rU')
#flg = 0
#rotVals = {}
#for i in RR:
#
#	if i[:4] == 'RESI':
#		resi = i.split()[1]
#		flg +=1
#	if i[:4] == "PROB":
#		step = 1
#		for p in i.split()[1:]:
#			rotVals[ resi + ',' + str( step ) ] = float(p) 
#			step += 1
#
#pic.dump( rotVals, open('/home/xray/termanal/support.default/rotlib/rotVals.pkl', 'wb') )

########### Dictionary for amino acid frequency in data base taking from CONFIND.cpp July 24, 2015
## should alter these values to reflect database, maybe supply equivilant hash in pickle 
aaProp = {}
aaProp["ALA"] = 7.73; aaProp["CYS"] = 1.84; aaProp["ASP"] = 5.82; aaProp["GLU"] = 6.61; aaProp["PHE"] = 4.05;
aaProp["GLY"] = 7.11; aaProp["HIS"] = 2.35; aaProp["HSD"] = 2.35; aaProp["ILE"] = 5.66; aaProp["LYS"] = 6.27;
aaProp["LEU"] = 8.83; aaProp["MET"] = 2.08; aaProp["ASN"] = 4.50; aaProp["PRO"] = 4.52; aaProp["GLN"] = 3.94;
aaProp["ARG"] = 5.03; aaProp["SER"] = 6.13; aaProp["THR"] = 5.53; aaProp["VAL"] = 6.91; aaProp["TRP"] = 1.51; aaProp["TYR"] = 3.54;


### See fraction of potential rotamers of residues nearby ligand reach within 3 Angstroms

with open( routf ) as file:

	record_lines = []
	stream_path = './tmp'
	# Write a tmp file for each rotamer to parse w/ prody
	# Call Auxillary pickle dictionary looking up rotamer library info/probailities 
	rotVals = pic.load( open('/home/xray/termanal/support.default/rotlib/rotVals.pkl', 'rb' ))
	# And another to read amino acid propensities


	contactFlag = 0
	for i in file:

		#Reset frequency count for each chain,res  
		prvResID = ''

		# if end of rotamer model, build rotamer object and refresh line list
		if i[:3] == 'REM': 
			resID =  tuple( i.split()[1].split(',')[:2] )
			res = i.split()[1].split(',')[2]
			rotID = i.split()[1].split(',')[2] + ',' +  i.split()[-1] 

			# record rotamer model if in ligand contact list
			if resID in contactList:
				contactFlag = 1
				tempStream = open( stream_path , 'w' )
				tempStream.write( i )

			else: 
				contactFlag = 0

		else: 
			if contactFlag > 0:
				if i[:3] == 'END':
					tempStream.write(i)
					tempStream.close()

					# build rotamer model
					print resID, rotID, 'rotamer probability:', rotVals[ rotID ], 'amino acid frq:', aaProp[res]

					rotamer = parsePDB( stream_path )


					#

				
				else: 
					tempStream.write( i )




			else: continue




			#sys.exit()






