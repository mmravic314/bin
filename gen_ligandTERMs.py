## Marco Mravic UCSF Biophysics DeGrado Lab August 2015 ##
# Use this to define the ligand atoms in a group and determine which residues nearby could be in potential contact with this ligand
# Calls much heavier script '~/bin/add_ligand_contacts_ConfindV2.py' to do rotamer analysis once there ligand groups per PDB is established
# Input: list of paths to PDbs to analysis, and definitions of ligand groups for each PDB 

# example command line
# python  gen_ligandTERMs.py ~/tertBuilding/localZNdbfiles.txt   /home/xray/tertBuilding/biPairs_byPDB.pkl    ~/tertBuilding/zN_rOUT/
# python  gen_ligandTERMs.py   path2pathlist    path2metal-PDBhash   path2ConfindOUTfiles


import cPickle as pic, os, sys, subprocess as sp
from PDButil import biMsite

### Created this hash in binuclearCleanPDB.py. Contains lists of pairs of ZN indexed by PDB
try:
	pairsByPdb = pic.load( open( sys.argv[2], 'rb') )
except IOError:
	print 'ERROR: cPickle binary file containing hash for pdb - diZn mapping has incorrect path\n'
	sys.exit()

scriptPath = os.path.abspath( '/home/xray/bin/add_ligand_contacts_ConfindV2.py' ) 

with open( sys.argv[1] ) as file:
	for i in file:
		print 'Entering...', i
		pID 		= os.path.splitext( os.path.basename( i.rstrip() )  )[0]
		pdbPath 	= i.rstrip()
		cmapPath 	= os.path.join(  os.path.abspath( sys.argv[3] ), '%s.cmap' % (pID) ) 
		rOutPath 	= os.path.join(  os.path.abspath( sys.argv[3] ), '%s.rout' % (pID) )


		for site in pairsByPdb[ pID ]:
			freq = '0.1'
			metals = [ filter(lambda x: x.isalpha(), metal ) for metal in site.name.split('+') ] 

			#Skip the occassional non-diZinc site
			if metals[0] != metals[1] or metals[0].upper() != 'ZN':
				continue

			#print site
			#  This ligand string aquisition is based on string splitting of strings formatted like: ZN_5458+ZN_5459 (objects with .name property) 
			# But should be switch able for other formattings like small molecule ligands. Just grabbing PDB atom index
			pathStr 	= '-' + '_'.join( [ filter(lambda x: x.isdigit(), metal ) for metal in site.name.split('+') ] )
			siteStr 	= ','.join( [ filter(lambda x: x.isdigit(), metal ) for metal in site.name.split('+') ] )
			outputPath 	= os.path.join(  os.path.abspath( sys.argv[3] ), '%s.freq' % ( pID + pathStr ) )

			cmd = [ 'python', scriptPath, '--p', pdbPath, '--rout', rOutPath, '--lig', siteStr, '--rcut', '10', '--freq', freq, '--c', cmapPath, '--o', outputPath ]
			#print cmd
			sp.call( cmd )
		
		#print



		

# python ~/bin/add_ligand_contacts_Confind.py --p ZN_db/1KN4.pdb --rout zN_rOUT/1KN4.rout --lig 2930,2931 --rcut 7 --freq 0 --c zN_rOUT/1KN4.cmap