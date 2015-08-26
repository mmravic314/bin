## Marco Mravic UCSF Biophysics DeGrado Lab July 2015 ##
#Example command line
##python ~/bin/add_ligand_contacts_Confind.py --p 3A4J/3A4J.pdb --rout 3A4J/3A4J.rout --lig 2930,2931 --rcut 7 --freq 0 --c 3A4J/3A4J.contacts 
# python ~/bin/add_ligand_contacts_Confind.py --p ZN_db/1KN4.pdb --rout zN_rOUT/1KN4.rout --lig 2930,2931 --rcut 7 --freq 0 --c zN_rOUT/1KN4.cmap



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
par.add_argument('--o', required = True, help = 'pathname for .cmap contact map output file')
par.add_argument('--c', required = True, help = 'path to .cmap contact list file')
par.add_argument('--lig', required = True, help = 'ligand atoms serial selection string (inclusive)')
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

		##################### Helper classes #########################

class Ligand():
	def __init__(self, serial, coords, chain, number, identity ):
		self.name 		= 	str( serial )
		self.coords 	= 	coords
		self.chain 		= 	chain 
		self.number 	= 	number
		self.identity 	= 	identity
		self.contacts 	= 	[]

	def __repr__(self):
		return self.name

class Residue():
	def __init__(self, chain, number, identity):
		self.chain 		= 		chain 
		self.number 	= 		number
		self.identity 	= 		identity
		# Contain each rotamer as list (value) hashed within each amino acid type
		self.rotamers 	= 		{} 

		# key should be serial of ligand, value should be a len 20 array w/ number of rotamers
		# contacting ligand per possible identity
		self. frequency = 0  
	
	def __repr__(self):
		return self.chain + ',' + self.number


########### Dictionary for amino acid frequency in Di Zinc data base August 17, 2015
## should alter these values to reflect each database, maybe supply equivilant hash in pickle 
## Determined by module in my local bin 'PDButil.py': function freqAA()
aaProp = {
'CYS': 1.19, 'ILE': 5.92, 'SER': 5.34, 'GLN': 3.35, 'LYS': 5.64, 'ASN': 3.91, 'PRO': 4.77, 'THR': 5.51, 
'PHE': 3.93, 'ALA': 8.62, 'HIS': 2.99, 'HSD': 2.99, 'GLY': 7.75, 'ASP': 6.27, 'LEU': 9.24, 'ARG': 4.52,
'TRP': 1.34, 'VAL': 7.1, 'GLU': 6.86, 'TYR': 3.46, 'MET': 2.29 }


########################  MAIN  #####################################


ROTLIB = '~/termanal/support.default/rotlib/RR2000.rotlib'
routf	= args.rout
pdbf	= args.p
rcut	= args.rcut 
freqcu 	= args.freqcut
cmap 	= args.c
ligStr 	= args.lig
ligStr 	= parseLigandStr( ligStr )
outf 	= args.o
print 
print "Entering input PDB", pdbf, ligStr

### Parse Pdb, using ProDy package, and get input ligand coordinates
##  make list of residues (chain,resnum) with heavy atoms within 7 angstroms
pdbGroup = parsePDB( pdbf )


contactList = []		# joint array for any residues within 7 angstroms from any ligand heavy atom
ligands = []			# array of ligand objects (coords + serial)
unique = []				# Lazy way to recall & index residues "nearby" ligands
for i in ligStr:
	lig_coords 		= pdbGroup.getBySerial( i, i+1 ).getCoords()[0]
	lig_chain		= pdbGroup.getBySerial( i, i+1 ).getChids()[0]
	lig_number		= pdbGroup.getBySerial( i, i+1 ).getResnums()[0]
	lig_resID		= pdbGroup.getBySerial( i, i+1 ).getResnames()[0]

	ligands.append( Ligand( i, lig_coords, lig_chain, str( lig_number ), lig_resID )  )

	if len( lig_coords ) < 3:
		print "\nWARNING: No valid coordinates of ligand with serial number", i


	cont = Contacts( pdbGroup )
	# Store all residues (as objects) within rcut (default: 7 Angstroms) from each ligand 
	for j in cont.select( rcut , lig_coords ):

		ID = j.getChid() + ',' + str( j.getResnum() ) 

		#Only consider natural amino acids (optional?) 
		if j.getResname() not in aaProp.keys(): continue

		if ID not in unique: 
			unique.append(ID)
			contactList.append( Residue( j.getChid(), str( j.getResnum()), j.getResname() ) )


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


### See fraction of potential rotamers of residues nearby ligand reach within 3 Angstroms
print "\nRotamers file"
with open( routf ) as file:

	record_lines = []
	tmp_path = './tmp'
	# Write a tmp file for each rotamer to parse w/ prody

	# Call Auxillary pickle dictionary looking up rotamer library info/probailities 
	rotVals = pic.load( open('/home/xray/termanal/support.default/rotlib/rotVals.pkl', 'rb' ))

	# Amino acid propensities in hash above

	contactFlag = 0

	# store all non-clashing rotamers of Residues in potential contact with ligand atoms 
	for i in file:

		# if end of rotamer model, build rotamer object and refresh line list
		if i[:3] == 'REM': 
			
			resID =  ','.join( i.split()[1].split(',')[:2] )			
			rotID = i.split()[1].split(',')[2] + ',' +  i.split()[-1] 
			# record rotamer model if in ligand contact list
			if resID in unique:
				contactFlag = 1
				tempStream = open( tmp_path , 'w' )
				tempStream.write( i )

			else: 
				contactFlag = 0

		else: 
			if contactFlag > 0:
				if i[:3] == 'END':
					tempStream.write(i)
					tempStream.close()

					# build rotamer model
					rotamer = parsePDB( tmp_path )

					# Find index of this residue's ID (e.g "A,56"; chain,resNumber) corresponding to Residue object in contactList
					# Store this rotamer AtomGroup (ProDy  object) in Residue object's hash by its rotamer ID (e.g. "ARG,11")
					contactList[ unique.index( resID ) ].rotamers[ rotID ] = rotamer
				
				else: 
					# heavy atoms only
					if i[13] == 'H': continue
					# Store 'ATOM' lines in a tmp file containing rotamer object
					tempStream.write( i )

			else: continue



print "\nDone parsing rotamers and pdb, %s\nCalculating frequency of rotamers contacting ligands...\n" % ( os.path.splitext( os.path.basename( pdbf ) )[0] )
##############################    Frequency Calculation        #######################################
# For each residue i, start a frequency sum for each ligand atom; checking if if each rotamer is within 3 Angstroms
# r_i 				= 	potential rotamer of i
# lig_j 			= 	ligand of atom serial number j
# P( a ) 			= 	Probability (in percent) of amino acid of that identity, "a", in dataset
# p( r_i ) 			= 	Probability of rotamer in Lovell/Richardson rotamer library
# C( r_i, lig_j ) 	= 	binary indicator if any heady atom of that rotamer is within 3 Ang of ligand j
################ Frequency =   sum( sum[ C(r_i, lig_j)*p(r_i)*P(a) ]-over-all-r_i-of-a )-over-all-a  /  sum( sum[ p(r_i)*P(a) ]-over-all-r_i-of-a )-over-all-a
#####################################################################
for resi in contactList:
	for l in ligands:
		fr 		= 0.0
		maxV 	= 0.0
		for name , rotObj in resi.rotamers.items():		
			# Increase max possible count
			maxV += rotVals[ name ] * aaProp[ name.split(',')[0] ]
			
			# Check if there a contact between any heavy atom in the current rotamer, r_i and current ligand, lig_j
			# pass if not, count it if so. Break loop and stop looking at rotamer if contact (< 3 Angstrom) is found
			for a in rotObj.iterAtoms():
				if calcDistance( a.getCoords(), l.coords) <= 3:
					fr += rotVals[ name ] * aaProp[ name.split(',')[0] ]
					break
		try: 
			l.contacts.append( ( resi , round( fr/maxV, 5 ) ) )
		except ZeroDivisionError:
			l.contacts.append( ( resi , 0.0 ) )



# Write to output file
oFile = open( outf, 'w' )
for l in ligands:
	oFile.write( "Residues and frequency of potential rotamers contacting ligand atom " + l.name + '\n' ) 
	oFile.write( '%s %s %s %s\n' % ( l.name, l.identity, l.number, l.chain ) )
	for k in l.contacts:
		if k[1] != 0:
			oFile.write( str( k ) + '\n' )  
print 'Done with %s\n' % (  os.path.splitext( os.path.basename( outf ) )[0]  )  


#print r, k, 'rotamer probability:', rotVals[ k ], 'amino acid frq:', aaProp[ k.split(',')[0] ]
# Clean up temporary files
os.remove( tmp_path )



