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
		print lig
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

# Parse Pdb, using ProDy package, and get input ligand coordinates
pdbGroup = parsePDB( pdbf )
print "\n"

contactList = []
for i in ligStr:
	lig_coords = pdbGroup.getBySerial( i, i+1 ).getCoords()[0]

	if len( lig_coords ) < 3:
		print "\nWARNING: No valid coordinates of ligand of serial number", i


	#print pdbGroup.select( pdbGroup.getBySerial( i, i+1 ).getCoords[0] , 7 )

	cont = Contacts( pdbGroup )
	print cont.select( 7, lig_coords )
#	print cont.select( pdbGroup.getBySerial( i, i+1 ).getCoords[0] , 7 ) 





